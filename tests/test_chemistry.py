from __future__ import annotations

import pytest

from antismash_review.chemistry import (
    FREE_MONOMER_MASSES,
    WATER_MASS,
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


def test_water_mass_constant() -> None:
    """Water: H2O = 2*1.00782503223 + 15.99491461957 ≈ 18.01056."""
    assert abs(WATER_MASS - 18.01056) < 0.001


def test_glycine_free_mass() -> None:
    """Free glycine: C2H5NO2 = 75.03203 Da (monoisotopic)."""
    mass = free_monomer_mass("Gly")
    assert mass is not None
    assert abs(mass - 75.03203) < 0.001


def test_serine_free_mass() -> None:
    """Free serine: C3H7NO3 = 105.04259 Da (monoisotopic)."""
    mass = free_monomer_mass("Ser")
    assert mass is not None
    assert abs(mass - 105.04259) < 0.001


def test_gly_gly_dipeptide_mass() -> None:
    """Gly-Gly linear dipeptide: 2*Gly_free - 1*H2O = 132.05350 Da."""
    gly = free_monomer_mass("Gly")
    assert gly is not None
    linear = 2 * gly - WATER_MASS
    assert abs(linear - 132.05350) < 0.001


def test_ser_leu_dipeptide_mass() -> None:
    """Ser-Leu linear dipeptide: Ser_free + Leu_free - H2O = 218.12627 Da."""
    ser = free_monomer_mass("Ser")
    leu = free_monomer_mass("Leu")
    assert ser is not None
    assert leu is not None
    linear = ser + leu - WATER_MASS
    assert abs(linear - 218.12627) < 0.001
