from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Sequence
from pathlib import Path

from antismash_review.models import (
    ClusterBlastHit,
    ClusterBlastPairing,
    ClusterBlastResult,
    ClusterBlastSearchType,
    Record,
)

_REGION_NUM_RE = re.compile(r"_c(?P<region>\d+)\.txt$", re.IGNORECASE)
_HEADER_RE = re.compile(r"^ClusterBlast scores for\s+(?P<record_id>\S+)", re.IGNORECASE)
_SIG_HIT_RE = re.compile(r"^(?P<rank>\d+)\.\s+(?P<accession>\S+)(?:\t(?P<description>.*))?$")
_DETAILS_HEADER_RE = re.compile(r"^(?P<rank>\d+)\.\s+(?P<accession>\S+)$")


class ClusterBlastParseError(RuntimeError):
    """A ClusterBlast sidecar was recognized but could not be parsed safely."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_clusterblast_text(
    path: Path,
    *,
    search_type: ClusterBlastSearchType,
) -> ClusterBlastResult:
    """Parse one antiSMASH 8 ClusterBlast text result."""
    path = Path(path)
    match_region = _REGION_NUM_RE.search(path.name)
    if not match_region:
        raise ClusterBlastParseError(f"Could not extract region number from filename {path.name}")
    region_number = int(match_region.group("region"))

    try:
        content = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ClusterBlastParseError(f"Could not read {path}: {exc}") from exc

    lines = content.splitlines()
    if not lines:
        raise ClusterBlastParseError(f"Empty ClusterBlast file: {path}")

    header_match = _HEADER_RE.match(lines[0].strip())
    if not header_match:
        raise ClusterBlastParseError(
            f"Missing or invalid ClusterBlast header in {path}: {lines[0]!r}"
        )
    record_id = header_match.group("record_id")

    # Locate "Significant hits:" and "Details:"
    sig_index = -1
    details_index = -1
    for idx, raw_line in enumerate(lines):
        line = raw_line.strip()
        if line == "Significant hits:":
            sig_index = idx
        elif line == "Details:":
            details_index = idx

    if sig_index == -1 or details_index == -1 or details_index < sig_index:
        raise ClusterBlastParseError(
            f"Malformed structure in {path}: Significant hits or Details section missing"
        )

    # Parse Significant hits
    sig_hits: list[tuple[int, str, str]] = []
    expected_rank = 1
    for raw_line in lines[sig_index + 1 : details_index]:
        line = raw_line.strip()
        if not line:
            continue
        sig_match = _SIG_HIT_RE.match(line)
        if not sig_match:
            raise ClusterBlastParseError(f"Malformed significant hit line in {path}: {line!r}")
        rank = int(sig_match.group("rank"))
        if rank != expected_rank:
            raise ClusterBlastParseError(
                f"Non-consecutive rank in {path}: expected {expected_rank}, got {rank}"
            )
        accession = sig_match.group("accession")
        description = sig_match.group("description") or ""
        sig_hits.append((rank, accession, description))
        expected_rank += 1

    # Parse Details section
    details_lines = lines[details_index + 1 :]
    # Split details into blocks separated by ">>"
    details_blocks: list[list[str]] = []
    current_block: list[str] = []
    for raw_line in details_lines:
        line = raw_line.rstrip()
        if line.strip() == ">>":
            if current_block:
                details_blocks.append(current_block)
            current_block = []
        else:
            current_block.append(line)
    if current_block:
        details_blocks.append(current_block)

    # Filter out empty blocks
    details_blocks = [b for b in details_blocks if any(line.strip() for line in b)]

    if len(details_blocks) != len(sig_hits):
        raise ClusterBlastParseError(
            f"Mismatched hit count in {path}: {len(sig_hits)} significant hits, "
            f"{len(details_blocks)} detail blocks"
        )

    rankings: list[ClusterBlastHit] = []
    for (sig_rank, sig_accession, sig_description), block in zip(
        sig_hits, details_blocks, strict=True
    ):
        # First non-empty line of block must be "<rank>. <accession>"
        non_empty = [line.strip() for line in block if line.strip()]
        if not non_empty:
            raise ClusterBlastParseError(f"Empty detail block in {path}")
        det_header = _DETAILS_HEADER_RE.match(non_empty[0])
        if not det_header:
            raise ClusterBlastParseError(f"Malformed detail header in {path}: {non_empty[0]!r}")
        det_rank = int(det_header.group("rank"))
        det_accession = det_header.group("accession")
        if det_rank != sig_rank or det_accession != sig_accession:
            raise ClusterBlastParseError(
                f"Detail header {det_rank}. {det_accession} does not match "
                f"significant hit {sig_rank}. {sig_accession} in {path}"
            )

        cluster_type: str | None = None
        num_hits: int | None = None
        blast_score: float | None = None
        pairings: list[ClusterBlastPairing] = []
        in_blast_table = False

        for raw_line in block:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("Type:"):
                cluster_type = line.split(":", 1)[1].strip() or None
            elif line.startswith("Number of proteins with BLAST hits to this cluster:"):
                try:
                    num_hits = int(line.split(":", 1)[1].strip())
                except ValueError as exc:
                    raise ClusterBlastParseError(f"Invalid hit count in {path}: {line!r}") from exc
            elif line.startswith("Cumulative BLAST score:"):
                try:
                    blast_score = float(line.split(":", 1)[1].strip())
                except ValueError as exc:
                    raise ClusterBlastParseError(
                        f"Invalid blast score in {path}: {line!r}"
                    ) from exc
            elif line.startswith("Table of Blast hits"):
                in_blast_table = True
            elif in_blast_table:
                cols = raw_line.split("\t")
                if len(cols) == 6:
                    try:
                        pairings.append(
                            ClusterBlastPairing(
                                query_gene=cols[0].strip(),
                                subject_gene=cols[1].strip(),
                                percent_identity=float(cols[2]),
                                blast_score=float(cols[3]),
                                percent_coverage=float(cols[4]),
                                evalue=float(cols[5]),
                            )
                        )
                    except ValueError as exc:
                        raise ClusterBlastParseError(
                            f"Invalid numeric value in BLAST hit row in {path}: {raw_line!r}"
                        ) from exc
                elif cols and cols[0].strip():
                    raise ClusterBlastParseError(
                        f"Expected 6 tab-separated BLAST hit columns in {path}, "
                        f"got {len(cols)}: {raw_line!r}"
                    )

        rankings.append(
            ClusterBlastHit(
                rank=sig_rank,
                accession=sig_accession,
                description=sig_description,
                cluster_type=cluster_type,
                num_hits=num_hits,
                core_gene_hits=None,
                blast_score=blast_score,
                synteny_score=None,
                core_bonus=None,
                similarity=None,
                pairings=pairings,
            )
        )

    return ClusterBlastResult(
        record_id=record_id,
        region_number=region_number,
        search_type=search_type,
        total_hits=None,
        rankings=rankings,
        source_path=path.resolve(),
        source_sha256=_sha256(path),
        source_format="text",
    )


def parse_clusterblast_json(path: Path) -> list[ClusterBlastResult]:
    """Parse only ClusterBlast modules from a native antiSMASH JSON sidecar."""
    path = Path(path)
    try:
        content = path.read_text(encoding="utf-8")
        data = json.loads(content)
    except (OSError, json.JSONDecodeError) as exc:
        raise ClusterBlastParseError(f"Could not parse JSON {path}: {exc}") from exc

    if not isinstance(data, dict) or "records" not in data:
        raise ClusterBlastParseError(f"Invalid antiSMASH JSON root in {path}")

    source_sha256 = _sha256(path)
    results: list[ClusterBlastResult] = []

    type_mapping: tuple[tuple[str, ClusterBlastSearchType], ...] = (
        ("general", "clusterblast"),
        ("knowncluster", "knownclusterblast"),
        ("subcluster", "subclusterblast"),
    )

    for record_idx, rec in enumerate(data.get("records", [])):
        rec_id = rec.get("id")
        modules = rec.get("modules", {})
        if "antismash.modules.clusterblast" not in modules:
            continue

        cb_mod = modules["antismash.modules.clusterblast"]
        if not isinstance(cb_mod, dict):
            raise ClusterBlastParseError(
                f"ClusterBlast module in record {record_idx} is not a dict in {path}"
            )

        module_schema_version = cb_mod.get("schema_version")
        if module_schema_version != 2:
            raise ClusterBlastParseError(
                f"Unsupported ClusterBlast module schema version {module_schema_version!r} "
                f"(expected 2) in {path}"
            )

        cb_record_id = cb_mod.get("record_id")
        if not isinstance(cb_record_id, str):
            raise ClusterBlastParseError(
                f"ClusterBlast module record_id is missing or not a string in {path}"
            )
        if rec_id is not None and cb_record_id != rec_id:
            raise ClusterBlastParseError(
                f"ClusterBlast module record_id {cb_record_id!r} does not match "
                f"record id {rec_id!r} in {path}"
            )

        for sec_key, search_type in type_mapping:
            if sec_key not in cb_mod:
                continue
            section = cb_mod[sec_key]
            if not isinstance(section, dict):
                continue

            result_schema_version = section.get("schema_version")
            if result_schema_version != 5:
                raise ClusterBlastParseError(
                    f"Unsupported ClusterBlast {sec_key} result schema version "
                    f"{result_schema_version!r} (expected 5) in {path}"
                )

            for region_res in section.get("results", []):
                region_number = region_res["region_number"]
                total_hits = region_res.get("total_hits")
                raw_rankings = region_res.get("ranking", [])

                rankings: list[ClusterBlastHit] = []
                for rank_idx, rank_entry in enumerate(raw_rankings, start=1):
                    if not isinstance(rank_entry, (list, tuple)) or len(rank_entry) != 2:
                        raise ClusterBlastParseError(
                            f"Invalid ranking entry format in {path}: {rank_entry!r}"
                        )
                    hit_info, hit_details = rank_entry
                    if not isinstance(hit_info, dict) or not isinstance(hit_details, dict):
                        raise ClusterBlastParseError(
                            f"Invalid ranking elements in {path}: {rank_entry!r}"
                        )

                    pairings: list[ClusterBlastPairing] = []
                    for pairing_entry in hit_details.get("pairings", []):
                        if not isinstance(pairing_entry, (list, tuple)) or len(pairing_entry) != 3:
                            raise ClusterBlastParseError(
                                f"Invalid pairing entry in {path}: {pairing_entry!r}"
                            )
                        query_str, subject_idx, pairing_dict = pairing_entry
                        parts = query_str.split("|", 5)
                        if len(parts) < 5:
                            raise ClusterBlastParseError(
                                f"Malformed query string in {path}: {query_str!r}"
                            )
                        query_gene = parts[4]
                        pairings.append(
                            ClusterBlastPairing(
                                query_gene=query_gene,
                                subject_gene=str(pairing_dict["name"]),
                                percent_identity=float(pairing_dict["perc_ident"]),
                                blast_score=float(pairing_dict["blastscore"]),
                                percent_coverage=float(pairing_dict["perc_coverage"]),
                                evalue=float(pairing_dict["evalue"]),
                                subject_protein_id=pairing_dict.get("locus_tag"),
                                subject_index=subject_idx if isinstance(subject_idx, int) else None,
                            )
                        )

                    raw_blast_score = hit_details.get("blast_score")
                    rankings.append(
                        ClusterBlastHit(
                            rank=rank_idx,
                            accession=str(hit_info["accession"]),
                            description=str(hit_info["description"]),
                            cluster_type=hit_info.get("cluster_type"),
                            num_hits=hit_details.get("hits"),
                            core_gene_hits=hit_details.get("core_gene_hits"),
                            blast_score=float(raw_blast_score)
                            if raw_blast_score is not None
                            else None,
                            synteny_score=hit_details.get("synteny_score"),
                            core_bonus=hit_details.get("core_bonus"),
                            similarity=hit_details.get("similarity"),
                            pairings=pairings,
                        )
                    )

                results.append(
                    ClusterBlastResult(
                        record_id=cb_record_id,
                        region_number=region_number,
                        search_type=search_type,
                        total_hits=total_hits,
                        rankings=rankings,
                        source_path=path.resolve(),
                        source_sha256=source_sha256,
                        source_format="json",
                        module_schema_version=module_schema_version,
                        result_schema_version=result_schema_version,
                    )
                )

    return results


def _result_key(
    result: ClusterBlastResult,
) -> tuple[str, int, ClusterBlastSearchType]:
    return (result.record_id, result.region_number, result.search_type)


def merge_clusterblast_results(
    text_results: Sequence[ClusterBlastResult],
    json_results: Sequence[ClusterBlastResult],
) -> list[ClusterBlastResult]:
    selected: dict[
        tuple[str, int, ClusterBlastSearchType],
        ClusterBlastResult,
    ] = {}
    for result in json_results:
        key = _result_key(result)
        if key in selected:
            raise ClusterBlastParseError(f"duplicate JSON ClusterBlast result: {key}")
        selected[key] = result
    for result in text_results:
        key = _result_key(result)
        if key in selected and selected[key].source_format == "text":
            raise ClusterBlastParseError(f"duplicate text ClusterBlast result: {key}")
        selected[key] = result
    return [selected[key] for key in sorted(selected)]


def attach_clusterblast_results(
    records: list[Record],
    results: list[ClusterBlastResult],
) -> None:
    for result in results:
        candidates = [
            record
            for record in records
            if result.record_id in {record.record_id, record.name}
            and any(region.number == result.region_number for region in record.regions)
        ]
        if len(candidates) != 1:
            raise ClusterBlastParseError(
                f"expected one GenBank target for {_result_key(result)}, found {len(candidates)}"
            )
        candidates[0].clusterblast_results.append(result)
