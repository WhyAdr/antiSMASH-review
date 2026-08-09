"""Reusable, conservative feature fingerprints for comparison and cohort analysis."""

from __future__ import annotations

import unicodedata
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass

from .models import Record
from .review import review_record


@dataclass(slots=True, frozen=True)
class DiagnosticFingerprint:
    """Stable value object used to compare evidence-scoped diagnostics."""

    code: str
    severity: str
    message: str
    feature_index: int | None


def normalize_label(value: str) -> str:
    """Return a comparison key without changing the source/display token."""

    return unicodedata.normalize("NFKC", value.strip()).casefold()


def _product_values(records: Sequence[Record]) -> list[str]:
    return [
        product for record in records for region in record.regions for product in region.products
    ]


def product_counter(
    records: Sequence[Record],
    *,
    normalized: bool = False,
) -> Counter[str]:
    """Count region products using raw tokens or explicit normalized keys.

    Raw mode is the compatibility mode used by pairwise comparison.  Normalized mode is
    intended for new matrix identities; it does not provide a biological synonym map.
    """

    values = _product_values(records)
    if normalized:
        return Counter(normalize_label(value) for value in values)
    return Counter(values)


def _domain_label(name: str | None, domain_id: str | None) -> str:
    if name and name.strip():
        return name
    if domain_id and domain_id.strip():
        return domain_id
    return "<unnamed-domain>"


def domain_counter(
    records: Sequence[Record],
    *,
    normalized: bool = False,
) -> Counter[str]:
    """Count adapted domain names, retaining a deterministic fallback for unnamed domains."""

    values = [
        _domain_label(domain.name, domain.domain_id)
        for record in records
        for domain in record.domains
    ]
    if normalized:
        return Counter(normalize_label(value) for value in values)
    return Counter(values)


def domain_presence(
    records: Sequence[Record],
    *,
    normalized: bool = False,
) -> frozenset[str]:
    """Return the domain-name set represented by ``records``."""

    return frozenset(domain_counter(records, normalized=normalized))


def diagnostic_counter(records: Sequence[Record]) -> Counter[DiagnosticFingerprint]:
    """Count diagnostics with the same value semantics as the comparison layer."""

    return Counter(
        DiagnosticFingerprint(
            diagnostic.code,
            diagnostic.severity.value,
            diagnostic.message,
            diagnostic.feature_index,
        )
        for record in records
        for diagnostic in review_record(record)
    )
