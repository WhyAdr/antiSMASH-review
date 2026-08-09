"""Dependency-free deterministic domain-content clustering."""

from __future__ import annotations

import re
from dataclasses import dataclass

from .cohort import CohortMember


@dataclass(slots=True)
class _ClusterNode:
    members: tuple[int, ...]
    height: float
    left: _ClusterNode | None = None
    right: _ClusterNode | None = None


def domain_jaccard_distances(members: list[CohortMember]) -> list[list[float]]:
    """Return a symmetric Jaccard distance matrix for binary domain presence."""

    domain_sets = [set(member.domain_counts) for member in members]
    matrix = [[0.0 for _ in members] for _ in members]
    for left_index, left_set in enumerate(domain_sets):
        for right_index in range(left_index + 1, len(domain_sets)):
            right_set = domain_sets[right_index]
            union = left_set | right_set
            similarity = len(left_set & right_set) / len(union) if union else 1.0
            distance = 1.0 - similarity
            matrix[left_index][right_index] = distance
            matrix[right_index][left_index] = distance
    return matrix


def _member_key(node: _ClusterNode, names: list[str]) -> tuple[str, ...]:
    return tuple(sorted(names[index] for index in node.members))


def _branch_length(value: float) -> str:
    rendered = f"{max(value, 0.0):.6f}".rstrip("0").rstrip(".")
    return rendered or "0"


def _newick_label(name: str) -> str:
    if re.fullmatch(r"[A-Za-z0-9_.-]+", name):
        return name
    return "'" + name.replace("'", "''") + "'"


def _render_newick(node: _ClusterNode, names: list[str], parent_height: float | None) -> str:
    if node.left is None or node.right is None:
        rendered = _newick_label(names[node.members[0]])
    else:
        rendered = (
            "("
            + ",".join(
                (
                    _render_newick(node.left, names, node.height),
                    _render_newick(node.right, names, node.height),
                )
            )
            + ")"
        )
    if parent_height is not None:
        rendered += ":" + _branch_length(parent_height - node.height)
    return rendered


def average_linkage_domain_clustering(
    members: list[CohortMember],
) -> tuple[list[list[float]], list[str], str]:
    """Cluster members with average linkage and deterministic lexical tie-breaking."""

    names = [member.name for member in members]
    distances = domain_jaccard_distances(members)
    if not members:
        return distances, [], ";"
    if len(members) == 1:
        node = _ClusterNode((0,), 0.0)
        return distances, [names[0]], _render_newick(node, names, None) + ";"

    active = [_ClusterNode((index,), 0.0) for index in range(len(members))]
    while len(active) > 1:
        best_pair: tuple[int, int] | None = None
        best_key: tuple[float, tuple[str, ...], tuple[str, ...]] | None = None
        for left_index, left in enumerate(active):
            for right_index in range(left_index + 1, len(active)):
                right = active[right_index]
                left_key = _member_key(left, names)
                right_key = _member_key(right, names)
                pair_key = tuple(sorted((left_key, right_key)))
                distance = sum(
                    distances[left_member][right_member]
                    for left_member in left.members
                    for right_member in right.members
                ) / (len(left.members) * len(right.members))
                candidate = (distance, pair_key[0], pair_key[1])
                if best_key is None or candidate < best_key:
                    best_key = candidate
                    best_pair = (left_index, right_index)

        assert best_pair is not None and best_key is not None
        left_index, right_index = best_pair
        left = active[left_index]
        right = active[right_index]
        children = sorted((left, right), key=lambda node: _member_key(node, names))
        merged = _ClusterNode(
            members=tuple(sorted(left.members + right.members)),
            height=best_key[0] / 2.0,
            left=children[0],
            right=children[1],
        )
        active = [
            node for index, node in enumerate(active) if index not in {left_index, right_index}
        ]
        active.append(merged)
        active.sort(key=lambda node: _member_key(node, names))

    root = active[0]
    order: list[int] = []

    def visit(node: _ClusterNode) -> None:
        if node.left is None or node.right is None:
            order.append(node.members[0])
            return
        visit(node.left)
        visit(node.right)

    visit(root)
    return distances, [names[index] for index in order], _render_newick(root, names, None) + ";"
