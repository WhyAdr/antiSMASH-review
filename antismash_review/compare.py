from __future__ import annotations

import statistics
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from antismash_review.models import Record
from antismash_review.review import review_record


@dataclass(slots=True, frozen=True)
class DiagnosticFingerprint:
    code: str
    severity: str
    message: str
    feature_index: int | None


@dataclass(slots=True, frozen=True)
class IntergenicSummary:
    gap_count: int
    total_bp: int
    mean_bp: float | None
    median_bp: float | None
    max_bp: int | None
    circular_wrap_included: bool


@dataclass(slots=True, frozen=True)
class CoordinateMatchEvidence:
    overlap_bp: int
    left_span_bp: int
    right_span_bp: int
    left_overlap_fraction: float
    right_overlap_fraction: float


@dataclass(slots=True)
class RecordComparison:
    left_record_id: str
    right_record_id: str
    match_key: str
    left_region_count: int
    right_region_count: int
    left_gene_count: int
    right_gene_count: int
    left_domain_count: int
    right_domain_count: int
    left_nrps_pks_count: int
    right_nrps_pks_count: int
    gained_products: list[str]
    lost_products: list[str]
    new_diagnostics: list[DiagnosticFingerprint]
    resolved_diagnostics: list[DiagnosticFingerprint]
    left_intergenic: IntergenicSummary
    right_intergenic: IntergenicSummary
    coordinate_evidence: CoordinateMatchEvidence | None = None


@dataclass(slots=True)
class ComparisonResult:
    left_input: Path
    right_input: Path
    match_method: str
    shared_coordinate_system_assumed: bool
    min_reciprocal_overlap: float | None
    matched: list[RecordComparison]
    unmatched_left: list[str]
    unmatched_right: list[str]


def intergenic_summary(record: Record) -> IntergenicSummary:
    intervals = sorted(
        (part.start, part.end) for gene in record.genes for part in gene.location.parts
    )
    if not intervals:
        return IntergenicSummary(0, 0, None, None, None, False)

    merged: list[list[int]] = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    gaps = [
        right[0] - left[1]
        for left, right in zip(merged, merged[1:], strict=False)
        if right[0] > left[1]
    ]
    circular = record.topology is not None and record.topology.casefold() == "circular"
    if circular:
        wrap_gap = (record.length - merged[-1][1]) + merged[0][0]
        if wrap_gap > 0:
            gaps.append(wrap_gap)

    return IntergenicSummary(
        gap_count=len(gaps),
        total_bp=sum(gaps),
        mean_bp=statistics.fmean(gaps) if gaps else None,
        median_bp=statistics.median(gaps) if gaps else None,
        max_bp=max(gaps) if gaps else None,
        circular_wrap_included=circular,
    )


def _feature_span_intervals(record: Record) -> list[tuple[int, int]]:
    raw_intervals: list[tuple[int, int]] = []
    collections = (
        record.regions + record.candidate_clusters + record.protoclusters + record.proto_cores
    )
    for col in collections:
        for part in col.location.parts:
            raw_intervals.append((part.start, part.end))
    for gene in record.genes:
        for part in gene.location.parts:
            raw_intervals.append((part.start, part.end))

    if not raw_intervals:
        return []

    sorted_intervals = sorted(raw_intervals)
    merged: list[list[int]] = []
    for start, end in sorted_intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(s, e) for s, e in merged]


def _calculate_overlap_evidence(left: Record, right: Record) -> CoordinateMatchEvidence:
    left_intervals = _feature_span_intervals(left)
    right_intervals = _feature_span_intervals(right)

    left_span = sum(end - start for start, end in left_intervals)
    right_span = sum(end - start for start, end in right_intervals)

    overlap_bp = 0
    for l_start, l_end in left_intervals:
        for r_start, r_end in right_intervals:
            o_start = max(l_start, r_start)
            o_end = min(l_end, r_end)
            if o_end > o_start:
                overlap_bp += o_end - o_start

    left_frac = overlap_bp / left_span if left_span > 0 else 0.0
    right_frac = overlap_bp / right_span if right_span > 0 else 0.0

    return CoordinateMatchEvidence(
        overlap_bp=overlap_bp,
        left_span_bp=left_span,
        right_span_bp=right_span,
        left_overlap_fraction=left_frac,
        right_overlap_fraction=right_frac,
    )


def _compare_single_pair(
    left: Record,
    right: Record,
    match_key: str,
    evidence: CoordinateMatchEvidence | None = None,
) -> RecordComparison:
    left_products = [p for r in left.regions for p in r.products]
    right_products = [p for r in right.regions for p in r.products]

    left_p_counter = Counter(left_products)
    right_p_counter = Counter(right_products)

    gained = sorted((right_p_counter - left_p_counter).elements())
    lost = sorted((left_p_counter - right_p_counter).elements())

    left_diags = [
        DiagnosticFingerprint(d.code, d.severity.value, d.message, d.feature_index)
        for d in review_record(left)
    ]
    right_diags = [
        DiagnosticFingerprint(d.code, d.severity.value, d.message, d.feature_index)
        for d in review_record(right)
    ]

    left_d_counter = Counter(left_diags)
    right_d_counter = Counter(right_diags)

    new_d = sorted(
        (right_d_counter - left_d_counter).elements(),
        key=lambda d: (d.code, d.severity, d.message, d.feature_index or 0),
    )
    resolved_d = sorted(
        (left_d_counter - right_d_counter).elements(),
        key=lambda d: (d.code, d.severity, d.message, d.feature_index or 0),
    )

    return RecordComparison(
        left_record_id=left.record_id,
        right_record_id=right.record_id,
        match_key=match_key,
        left_region_count=len(left.regions),
        right_region_count=len(right.regions),
        left_gene_count=len(left.genes),
        right_gene_count=len(right.genes),
        left_domain_count=len(left.domains),
        right_domain_count=len(right.domains),
        left_nrps_pks_count=len(left.nrps_pks_domains),
        right_nrps_pks_count=len(right.nrps_pks_domains),
        gained_products=gained,
        lost_products=lost,
        new_diagnostics=new_d,
        resolved_diagnostics=resolved_d,
        left_intergenic=intergenic_summary(left),
        right_intergenic=intergenic_summary(right),
        coordinate_evidence=evidence,
    )


def compare_records(
    left_records: list[Record],
    right_records: list[Record],
    *,
    left_input: Path,
    right_input: Path,
    match_method: str = "record_id",
    assume_shared_coordinate_system: bool = False,
    min_reciprocal_overlap: float = 0.80,
) -> ComparisonResult:
    matched: list[RecordComparison] = []
    unmatched_left: list[str] = []
    unmatched_right: list[str] = []

    if match_method == "record_id":
        left_ids = [r.record_id for r in left_records]
        right_ids = [r.record_id for r in right_records]
        if len(left_ids) != len(set(left_ids)):
            raise ValueError("Duplicate record ID found on left input")
        if len(right_ids) != len(set(right_ids)):
            raise ValueError("Duplicate record ID found on right input")

        right_by_id = {r.record_id: r for r in right_records}
        matched_right_ids = set()
        for left in left_records:
            if left.record_id in right_by_id:
                right = right_by_id[left.record_id]
                matched_right_ids.add(right.record_id)
                matched.append(_compare_single_pair(left, right, match_key=left.record_id))
            else:
                unmatched_left.append(left.record_id)

        for right in right_records:
            if right.record_id not in matched_right_ids:
                unmatched_right.append(right.record_id)

    elif match_method == "record_region":
        # Every record must contain exactly one numbered region
        left_keys: dict[tuple[str, int], Record] = {}
        for r in left_records:
            nums = [reg.number for reg in r.regions if reg.number is not None]
            if len(nums) != 1:
                raise ValueError(
                    f"Record {r.record_id} does not contain exactly one numbered region "
                    f"(found {len(nums)})"
                )
            key = (r.record_id, nums[0])
            if key in left_keys:
                raise ValueError(f"Duplicate (record_id, region_number) on left: {key}")
            left_keys[key] = r

        right_keys: dict[tuple[str, int], Record] = {}
        for r in right_records:
            nums = [reg.number for reg in r.regions if reg.number is not None]
            if len(nums) != 1:
                raise ValueError(
                    f"Record {r.record_id} does not contain exactly one numbered region "
                    f"(found {len(nums)})"
                )
            key = (r.record_id, nums[0])
            if key in right_keys:
                raise ValueError(f"Duplicate (record_id, region_number) on right: {key}")
            right_keys[key] = r

        matched_right_keys = set()
        for key, left in sorted(left_keys.items()):
            if key in right_keys:
                right = right_keys[key]
                matched_right_keys.add(key)
                match_str = f"{key[0]}:region_{key[1]}"
                matched.append(_compare_single_pair(left, right, match_key=match_str))
            else:
                unmatched_left.append(f"{key[0]}:region_{key[1]}")

        for key, _right_rec in sorted(right_keys.items()):
            if key not in matched_right_keys:
                unmatched_right.append(f"{key[0]}:region_{key[1]}")

    elif match_method == "single_record":
        if len(left_records) != 1 or len(right_records) != 1:
            raise ValueError(
                f"single_record matching requires exactly one record on each side "
                f"(left has {len(left_records)}, right has {len(right_records)})"
            )
        left = left_records[0]
        right = right_records[0]
        match_str = f"{left.record_id} <-> {right.record_id}"
        matched.append(_compare_single_pair(left, right, match_key=match_str))

    elif match_method == "coordinate_overlap":
        if not assume_shared_coordinate_system:
            raise ValueError(
                "coordinate_overlap matching requires explicit "
                "--assume-shared-coordinate-system confirmation"
            )
        if not 0 < min_reciprocal_overlap <= 1.0:
            raise ValueError("min_reciprocal_overlap must be in the interval (0, 1]")

        # Find candidates for each left record
        candidates_for_left: list[
            tuple[Record, list[tuple[int, Record, CoordinateMatchEvidence]]]
        ] = []
        for left in left_records:
            candidates: list[tuple[int, Record, CoordinateMatchEvidence]] = []
            for right_index, right in enumerate(right_records):
                evidence = _calculate_overlap_evidence(left, right)
                if (
                    evidence.left_overlap_fraction >= min_reciprocal_overlap
                    and evidence.right_overlap_fraction >= min_reciprocal_overlap
                ):
                    candidates.append((right_index, right, evidence))
            candidates_for_left.append((left, candidates))

        matched_right_records: dict[int, str] = {}  # right_index -> left_id
        for left, cand_list in candidates_for_left:
            if len(cand_list) == 0:
                unmatched_left.append(left.record_id)
            elif len(cand_list) > 1:
                raise ValueError(
                    f"Ambiguous coordinate match: left record {left.record_id} matched "
                    f"{len(cand_list)} right records above threshold {min_reciprocal_overlap}"
                )
            else:
                right_index, right, evidence = cand_list[0]
                if right_index in matched_right_records:
                    prev_left = matched_right_records[right_index]
                    raise ValueError(
                        f"Non-one-to-one coordinate match: right record {right.record_id} matched "
                        f"by both {prev_left} and {left.record_id}"
                    )
                matched_right_records[right_index] = left.record_id
                match_str = (
                    f"{left.record_id} <-> {right.record_id} (overlap {evidence.overlap_bp} bp)"
                )
                matched.append(
                    _compare_single_pair(left, right, match_key=match_str, evidence=evidence)
                )

        matched_right_set = set(matched_right_records)
        for right_index, right in enumerate(right_records):
            if right_index not in matched_right_set:
                unmatched_right.append(right.record_id)

    else:
        raise ValueError(f"Unknown match method: {match_method}")

    return ComparisonResult(
        left_input=left_input.resolve(),
        right_input=right_input.resolve(),
        match_method=match_method,
        shared_coordinate_system_assumed=assume_shared_coordinate_system,
        min_reciprocal_overlap=min_reciprocal_overlap
        if match_method == "coordinate_overlap"
        else None,
        matched=matched,
        unmatched_left=unmatched_left,
        unmatched_right=unmatched_right,
    )
