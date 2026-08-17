from __future__ import annotations

import csv
import json
import os
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

from antismash_review.assemblyline import predict_assembly_lines
from antismash_review.discovery import discover
from antismash_review.loading import LoadedReviewInput, load_review_input
from antismash_review.review import review_record
from tests.fixtures.build_fixture import write_synthetic_cross_cds_monomer_genbank

_FINGERPRINT_SCHEMA_VERSION = "0.1.0"


@dataclass(slots=True, frozen=True)
class _IntegrationCase:
    name: str
    input_path: Path
    fingerprint_path: Path
    recursive: bool
    lenient: bool


def _parse_bool(value: str, *, field: str, row_number: int) -> bool:
    normalized = value.strip().casefold()
    if normalized in {"", "0", "false", "no"}:
        return False
    if normalized in {"1", "true", "yes"}:
        return True
    raise ValueError(f"integration manifest row {row_number}: invalid {field} value {value!r}")


def _resolve_manifest_path(value: str, manifest_dir: Path) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (manifest_dir / path).resolve()


def _read_integration_manifest(path: Path) -> list[_IntegrationCase]:
    if path.suffix.casefold() != ".tsv":
        raise ValueError("integration manifest must be a TSV file")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"name", "path", "fingerprint"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("integration manifest requires name, path, and fingerprint columns")
        cases: list[_IntegrationCase] = []
        seen_names: set[str] = set()
        for row_number, row in enumerate(reader, start=2):
            name = (row.get("name") or "").strip()
            input_value = (row.get("path") or "").strip()
            fingerprint_value = (row.get("fingerprint") or "").strip()
            if not name or not input_value or not fingerprint_value:
                raise ValueError(
                    f"integration manifest row {row_number}: name, path, and fingerprint "
                    "must be non-empty"
                )
            if name in seen_names:
                raise ValueError(f"integration manifest row {row_number}: duplicate name {name!r}")
            seen_names.add(name)
            cases.append(
                _IntegrationCase(
                    name=name,
                    input_path=_resolve_manifest_path(input_value, path.parent),
                    fingerprint_path=_resolve_manifest_path(fingerprint_value, path.parent),
                    recursive=_parse_bool(
                        row.get("recursive") or "", field="recursive", row_number=row_number
                    ),
                    lenient=_parse_bool(
                        row.get("lenient") or "", field="lenient", row_number=row_number
                    ),
                )
            )
    if not cases:
        raise ValueError("integration manifest contains no cases")
    return cases


def _counter_dict(values: list[str]) -> dict[str, int]:
    return dict(sorted(Counter(values).items()))


def _integration_fingerprint(loaded: LoadedReviewInput) -> dict[str, Any]:
    predictions = [
        prediction for record in loaded.records for prediction in predict_assembly_lines(record)
    ]
    modules = [module for prediction in predictions for module in prediction.modules]
    diagnostics = [
        diagnostic for record in loaded.records for diagnostic in review_record(record)
    ] + list(loaded.input_diagnostics)
    sidecars = sorted(
        (
            {
                "record_id": result.record_id,
                "region_number": result.region_number,
                "search_type": result.search_type,
                "source_name": result.source_path.name,
                "source_format": result.source_format,
                "module_schema_version": result.module_schema_version,
                "result_schema_version": result.result_schema_version,
                "data_version": result.data_version,
            }
            for record in loaded.records
            for result in record.clusterblast_results
        ),
        key=lambda item: (
            item["record_id"],
            item["region_number"],
            item["search_type"],
            item["source_name"],
        ),
    )
    mass_states = [
        {
            "record_id": prediction.record_id,
            "region_number": prediction.region_number,
            "linear_non_null": bool(
                prediction.mass and prediction.mass.linear_core_mass_da is not None
            ),
            "cyclic_non_null": bool(
                prediction.mass
                and prediction.mass.head_to_tail_cyclic_candidate_mass_da is not None
            ),
        }
        for prediction in predictions
    ]
    knowncluster_accessions = sorted(
        {
            hit.accession
            for record in loaded.records
            for result in record.clusterblast_results
            if result.search_type == "knownclusterblast"
            for hit in result.rankings
        }
    )
    return {
        "schema_version": _FINGERPRINT_SCHEMA_VERSION,
        "record_ids": [record.record_id for record in loaded.records],
        "antismash_versions": [record.antismash_version for record in loaded.records],
        "sidecars": sidecars,
        "region_count": sum(len(record.regions) for record in loaded.records),
        "product_counts": _counter_dict(
            [
                product
                for record in loaded.records
                for region in record.regions
                for product in region.products
            ]
        ),
        "domain_count": sum(len(record.domains) for record in loaded.records),
        "nrps_pks_domain_count": sum(len(record.nrps_pks_domains) for record in loaded.records),
        "module_count": sum(len(record.modules) for record in loaded.records),
        "assembly_line_count": len(predictions),
        "raw_monomer_pairing_count": sum(module.raw_pairing_count for module in modules),
        "interpreted_incorporation_count": sum(len(prediction.chain) for prediction in predictions),
        "pairing_status_counts": _counter_dict([module.pairing_status for module in modules]),
        "assemblyline_integrity_flag_counts": _counter_dict(
            [flag for module in modules for flag in module.integrity_flags]
        ),
        "mass_candidate_states": mass_states,
        "diagnostic_code_counts": _counter_dict([diagnostic.code for diagnostic in diagnostics]),
        "knowncluster_accessions": knowncluster_accessions,
    }


def test_integration_fingerprint_locks_cross_cds_duplicate_semantics(tmp_path: Path) -> None:
    write_synthetic_cross_cds_monomer_genbank(tmp_path / "cross_cds.gbk")
    fingerprint = _integration_fingerprint(load_review_input(discover(tmp_path)))

    assert fingerprint["raw_monomer_pairing_count"] == 2
    assert fingerprint["interpreted_incorporation_count"] == 1
    assert fingerprint["pairing_status_counts"] == {"identical_duplicate": 1}
    assert fingerprint["assemblyline_integrity_flag_counts"] == {
        "cross_cds_duplicate_monomer_pairing": 1,
        "duplicate_monomer_pairing": 1,
    }
    assert fingerprint["mass_candidate_states"] == [
        {
            "record_id": "SYNTH_CROSS.1",
            "region_number": 1,
            "linear_non_null": False,
            "cyclic_non_null": False,
        }
    ]


def test_integration_manifest_tsv_contract(tmp_path: Path) -> None:
    input_dir = tmp_path / "private-case"
    input_dir.mkdir()
    fingerprint_path = tmp_path / "fingerprints" / "case.json"
    fingerprint_path.parent.mkdir()
    fingerprint_path.write_text("{}", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "name\tpath\tfingerprint\trecursive\tlenient\n"
        "case-a\tprivate-case\tfingerprints/case.json\tyes\t0\n",
        encoding="utf-8",
    )

    cases = _read_integration_manifest(manifest)
    assert cases == [
        _IntegrationCase(
            name="case-a",
            input_path=input_dir.resolve(),
            fingerprint_path=fingerprint_path.resolve(),
            recursive=True,
            lenient=False,
        )
    ]


def test_optional_real_data_integration_manifest() -> None:
    manifest_value = os.environ.get("ANTISMASH_REVIEW_INTEGRATION_MANIFEST")
    if not manifest_value:
        pytest.skip(
            "ANTISMASH_REVIEW_INTEGRATION_MANIFEST not set; skipping private integration suite"
        )

    manifest_path = Path(manifest_value).resolve()
    if not manifest_path.is_file():
        pytest.fail(f"integration manifest does not exist: {manifest_path}")

    try:
        cases = _read_integration_manifest(manifest_path)
    except (OSError, ValueError) as exc:
        pytest.fail(str(exc))

    for case in cases:
        if not case.input_path.exists():
            pytest.fail(f"integration case {case.name!r} input does not exist: {case.input_path}")
        if not case.fingerprint_path.is_file():
            pytest.fail(
                f"integration case {case.name!r} fingerprint does not exist: "
                f"{case.fingerprint_path}"
            )
        loaded = load_review_input(
            discover(case.input_path, recursive=case.recursive),
            lenient=case.lenient,
        )
        expected = json.loads(case.fingerprint_path.read_text(encoding="utf-8"))
        assert _integration_fingerprint(loaded) == expected, case.name
