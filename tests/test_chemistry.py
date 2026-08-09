from __future__ import annotations

import pytest

from antismash_review.chemistry import (
    FREE_MONOMER_MASSES,
    canonical_monomer_name,
    free_monomer_mass,
)


@pytest.mark.parametrize("name", sorted(FREE_MONOMER_MASSES))
def test_proteinogenic_registry_has_formula_derived_mass_and_stereo_alias(name: str) -> None:
    assert canonical_monomer_name(name) == name
    assert canonical_monomer_name(f"D-{name}") == name
    assert free_monomer_mass(name) == FREE_MONOMER_MASSES[name]
    assert FREE_MONOMER_MASSES[name] > 0


@pytest.mark.parametrize("value", ["X", "NH2", "ccmal", "D-Orn", None, ""])
def test_unmodeled_antiSMASH_calls_remain_unresolved(value: str | None) -> None:
    assert canonical_monomer_name(value) is None
    assert free_monomer_mass(value) is None
