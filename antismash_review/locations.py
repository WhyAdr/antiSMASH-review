"""Shared location-overlap helpers for scoped analysis and export."""

from __future__ import annotations

from .models import CollectionFeature, Domain, Location, Module, Record


def overlaps(left: Location, right: Location) -> bool:
    """Return whether any half-open location parts overlap."""

    return any(
        left_part.start < right_part.end and right_part.start < left_part.end
        for left_part in left.parts
        for right_part in right.parts
    )


def containing_regions(record: Record, location: Location) -> list[CollectionFeature]:
    return [region for region in record.regions if overlaps(region.location, location)]


def domains_in_collection(record: Record, collection: CollectionFeature) -> list[Domain]:
    return [domain for domain in record.domains if overlaps(domain.location, collection.location)]


def modules_in_collection(record: Record, collection: CollectionFeature) -> list[Module]:
    return [module for module in record.modules if overlaps(module.location, collection.location)]
