from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

GENBANK_SUFFIXES = (".gbk", ".gb", ".gbff")


@dataclass(slots=True, frozen=True)
class InputManifest:
    root: Path
    json_files: tuple[Path, ...]
    aggregate_genbanks: tuple[Path, ...]
    region_genbanks: tuple[Path, ...]
    ignored_files: tuple[Path, ...]


def discover(path: Path, *, recursive: bool = False) -> InputManifest:
    path = Path(path).resolve()
    if path.is_file():
        suffix = path.suffix.casefold()
        normalized_name = path.name.casefold()
        is_json = suffix == ".json"
        is_gbk = suffix in GENBANK_SUFFIXES
        if not is_json and not is_gbk:
            raise ValueError(f"Unsupported input: {path}")
        return InputManifest(
            root=path.parent,
            json_files=(path,) if is_json else (),
            aggregate_genbanks=() if is_json or ".region" in normalized_name else (path,),
            region_genbanks=(path,) if is_gbk and ".region" in normalized_name else (),
            ignored_files=(),
        )

    if not path.is_dir():
        raise FileNotFoundError(path)

    iterator = path.rglob("*") if recursive else path.glob("*")
    files = tuple(sorted(candidate for candidate in iterator if candidate.is_file()))
    json_files = tuple(candidate for candidate in files if candidate.suffix.casefold() == ".json")
    genbanks = tuple(
        candidate for candidate in files if candidate.suffix.casefold() in GENBANK_SUFFIXES
    )
    recognized = set(json_files) | set(genbanks)
    return InputManifest(
        root=path,
        json_files=json_files,
        aggregate_genbanks=tuple(
            item for item in genbanks if ".region" not in item.name.casefold()
        ),
        region_genbanks=tuple(item for item in genbanks if ".region" in item.name.casefold()),
        ignored_files=tuple(item for item in files if item not in recognized),
    )
