"""Shared feature enumeration for deterministic GFF3 and BED tracks."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from antismash_review.locations import containing_regions
from antismash_review.models import CollectionFeature, Domain, Gene, Location, Module, Record
from antismash_review.review import review_findings


@dataclass(slots=True, frozen=True)
class TrackFeature:
    record_id: str
    seqid: str
    feature_type: str
    location: Location
    ordinal: int
    preferred_id: str | None
    attributes: tuple[tuple[str, str], ...]


def sequence_id(record: Record) -> str:
    """Return a stable sequence name without merging rebased region records."""

    if ".region" in record.source_path.name.casefold():
        return f"{record.source_path.stem}:{record.record_id}"
    return record.record_id


def stable_feature_id(
    record_id: str,
    feature_type: str,
    ordinal: int,
    preferred: str | None,
) -> str:
    token = preferred.strip() if preferred and preferred.strip() else f"{ordinal:04d}"
    return f"{record_id}:{feature_type}:{token}"


def _joined(values: Sequence[object]) -> str:
    return ",".join(str(value) for value in values)


def _collection_feature(
    record: Record, seqid: str, feature: CollectionFeature, ordinal: int
) -> TrackFeature:
    return TrackFeature(
        record_id=record.record_id,
        seqid=seqid,
        feature_type=feature.feature_type,
        location=feature.location,
        ordinal=ordinal,
        preferred_id=str(feature.number) if feature.number is not None else None,
        attributes=tuple(
            (key, value)
            for key, value in (
                ("record_id", record.record_id),
                ("product", _joined(feature.products)),
                ("number", str(feature.number) if feature.number is not None else ""),
                (
                    "contig_edge",
                    str(feature.contig_edge).lower() if feature.contig_edge is not None else "",
                ),
                ("creating_tool", feature.creating_tool or ""),
            )
            if value
        ),
    )


def _gene_feature(record: Record, seqid: str, gene: Gene, ordinal: int) -> TrackFeature:
    region_numbers = [
        region.number
        for region in containing_regions(record, gene.location)
        if region.number is not None
    ]
    return TrackFeature(
        record_id=record.record_id,
        seqid=seqid,
        feature_type="CDS",
        location=gene.location,
        ordinal=ordinal,
        preferred_id=gene.locus_tag,
        attributes=tuple(
            (key, value)
            for key, value in (
                ("record_id", record.record_id),
                ("locus_tag", gene.locus_tag or ""),
                ("gene", gene.gene or ""),
                ("product", gene.product or ""),
                ("gene_kind", gene.gene_kind),
                ("region_numbers", _joined(region_numbers)),
            )
            if value
        ),
    )


def _domain_feature(record: Record, seqid: str, domain: Domain, ordinal: int) -> TrackFeature:
    return TrackFeature(
        record_id=record.record_id,
        seqid=seqid,
        feature_type="aSDomain",
        location=domain.location,
        ordinal=ordinal,
        preferred_id=domain.domain_id or domain.name,
        attributes=tuple(
            (key, value)
            for key, value in (
                ("record_id", record.record_id),
                ("Name", domain.name or ""),
                ("domain_id", domain.domain_id or ""),
                ("tool", domain.tool or ""),
                ("locus_tag", domain.locus_tag or ""),
                ("subtypes", _joined(domain.subtypes)),
                ("specificity", _joined(domain.specificity)),
            )
            if value
        ),
    )


def _module_feature(record: Record, seqid: str, module: Module, ordinal: int) -> TrackFeature:
    return TrackFeature(
        record_id=record.record_id,
        seqid=seqid,
        feature_type="aSModule",
        location=module.location,
        ordinal=ordinal,
        preferred_id=module.locus_tags[0] if module.locus_tags else None,
        attributes=tuple(
            (key, value)
            for key, value in (
                ("record_id", record.record_id),
                ("type", module.module_type or ""),
                ("locus_tags", _joined(module.locus_tags)),
                ("complete", str(module.complete).lower()),
                ("starter", str(module.starter).lower()),
                ("final", str(module.final).lower()),
                ("iterative", str(module.iterative).lower()),
                ("monomer_pairings", _joined(module.monomer_pairings)),
            )
            if value
        ),
    )


def iter_track_features(records: list[Record]) -> list[TrackFeature]:
    """Enumerate parsed entities and localized findings in stable source order."""

    features: list[TrackFeature] = []
    for record in records:
        seqid = sequence_id(record)
        collections = (
            *record.regions,
            *record.candidate_clusters,
            *record.protoclusters,
            *record.proto_cores,
        )
        ordinals: dict[str, int] = {}
        for feature in collections:
            ordinals[feature.feature_type] = ordinals.get(feature.feature_type, 0) + 1
            features.append(
                _collection_feature(record, seqid, feature, ordinals[feature.feature_type])
            )
        for ordinal, gene in enumerate(record.genes, start=1):
            features.append(_gene_feature(record, seqid, gene, ordinal))
        for ordinal, domain in enumerate(record.domains, start=1):
            features.append(_domain_feature(record, seqid, domain, ordinal))
        for ordinal, module in enumerate(record.modules, start=1):
            features.append(_module_feature(record, seqid, module, ordinal))
        finding_ordinal = 0
        for finding in review_findings(record):
            if finding.location is None:
                continue
            finding_ordinal += 1
            preferred = finding.diagnostic.code
            if finding.entity_id:
                preferred = f"{preferred}:{finding.entity_id}"
            features.append(
                TrackFeature(
                    record_id=record.record_id,
                    seqid=seqid,
                    feature_type="review_finding",
                    location=finding.location,
                    ordinal=finding_ordinal,
                    preferred_id=preferred,
                    attributes=tuple(
                        (key, value)
                        for key, value in (
                            ("record_id", record.record_id),
                            ("code", finding.diagnostic.code),
                            ("severity", finding.diagnostic.severity.value),
                            ("message", finding.diagnostic.message),
                            ("entity_type", finding.entity_type or ""),
                            ("entity_id", finding.entity_id or ""),
                        )
                        if value
                    ),
                )
            )
    return features
