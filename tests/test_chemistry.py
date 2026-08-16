from __future__ import annotations

import pytest

from antismash_review.chemistry import (
    FREE_MONOMER_MASSES,
    WATER_MASS,
    canonical_monomer_name,
    free_monomer_mass,
)

# Golden neutral monoisotopic masses (Da) computed independently from IUPAC exact isotopes:
# C=12.0, H=1.00782503223, N=14.00307400443, O=15.99491461957, S=31.9720711744.
# These values serve as external chemical truth and must not be computed from production tables.
GOLDEN_H2O = 18.01056468403
GOLDEN_FREE_GLY = 75.03202840472
GOLDEN_FREE_SER = 105.04259308875
GOLDEN_FREE_LEU = 131.09462866256
GOLDEN_LINEAR_GLY_GLY = 132.05349212541
GOLDEN_LINEAR_SER_LEU = 218.12665706728
GOLDEN_CYCLIC_SER_LEU = 200.11609238325
HISTORICAL_BUGGY_SER_LEU = 182.10552769922


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
    """Water: H2O = 2*1.00782503223 + 15.99491461957 = 18.01056468403 Da."""
    assert pytest.approx(GOLDEN_H2O, abs=1e-6) == WATER_MASS


def test_free_glycine_mass() -> None:
    """Free glycine: C2H5NO2 = 75.03202840472 Da."""
    mass = free_monomer_mass("Gly")
    assert mass is not None
    assert mass == pytest.approx(GOLDEN_FREE_GLY, abs=1e-6)


def test_free_serine_mass() -> None:
    """Free serine: C3H7NO3 = 105.04259308875 Da."""
    mass = free_monomer_mass("Ser")
    assert mass is not None
    assert mass == pytest.approx(GOLDEN_FREE_SER, abs=1e-6)


def test_free_leucine_mass() -> None:
    """Free leucine: C6H13NO2 = 131.09462866256 Da."""
    mass = free_monomer_mass("Leu")
    assert mass is not None
    assert mass == pytest.approx(GOLDEN_FREE_LEU, abs=1e-6)


def test_gly_gly_dipeptide_mass() -> None:
    """Gly-Gly linear dipeptide: 2*Gly_free - 1*H2O = 132.05349212541 Da."""
    gly = free_monomer_mass("Gly")
    assert gly is not None
    linear = 2 * gly - WATER_MASS
    assert linear == pytest.approx(GOLDEN_LINEAR_GLY_GLY, abs=1e-6)


def test_ser_leu_dipeptide_mass() -> None:
    """Ser-Leu linear dipeptide: Ser_free + Leu_free - H2O = 218.12665706728 Da."""
    ser = free_monomer_mass("Ser")
    leu = free_monomer_mass("Leu")
    assert ser is not None
    assert leu is not None
    linear = ser + leu - WATER_MASS
    assert linear == pytest.approx(GOLDEN_LINEAR_SER_LEU, abs=1e-6)


def test_cyclic_ser_leu_dipeptide_mass() -> None:
    """Ser-Leu cyclic dipeptide: Ser_free + Leu_free - 2*H2O = 200.11609238325 Da."""
    ser = free_monomer_mass("Ser")
    leu = free_monomer_mass("Leu")
    assert ser is not None
    assert leu is not None
    cyclic = ser + leu - 2 * WATER_MASS
    assert cyclic == pytest.approx(GOLDEN_CYCLIC_SER_LEU, abs=1e-6)


def test_linear_dipeptide_does_not_double_apply_dehydration() -> None:
    """Regression test against the historical residue-formula bug (~182.1055 Da)."""
    ser = free_monomer_mass("Ser")
    leu = free_monomer_mass("Leu")
    assert ser is not None and leu is not None
    linear = ser + leu - WATER_MASS
    assert abs(linear - HISTORICAL_BUGGY_SER_LEU) > 30.0  # Error was exactly 2*H2O (~36.02 Da)
