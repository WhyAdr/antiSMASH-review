"""Small, explicit monomer registry for conservative core-mass candidates."""

from __future__ import annotations

import re
import unicodedata

# Monoisotopic neutral atom masses use the standard exact isotope masses for C-12,
# H-1, N-14, O-16, and S-32.  Keeping formulas here makes the fixture calculations
# independently reproducible and avoids a cheminformatics dependency.
_ATOMIC_MASSES = {
    "C": 12.0,
    "H": 1.00782503223,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "S": 31.9720711744,
}

_FORMULAS: dict[str, dict[str, int]] = {
    "Ala": {"C": 3, "H": 5, "N": 1, "O": 1},
    "Arg": {"C": 6, "H": 12, "N": 4, "O": 1},
    "Asn": {"C": 4, "H": 6, "N": 2, "O": 2},
    "Asp": {"C": 4, "H": 5, "N": 1, "O": 3},
    "Cys": {"C": 3, "H": 5, "N": 1, "O": 1, "S": 1},
    "Gln": {"C": 5, "H": 8, "N": 2, "O": 2},
    "Glu": {"C": 5, "H": 7, "N": 1, "O": 3},
    "Gly": {"C": 2, "H": 3, "N": 1, "O": 1},
    "His": {"C": 6, "H": 7, "N": 3, "O": 1},
    "Ile": {"C": 6, "H": 11, "N": 1, "O": 1},
    "Leu": {"C": 6, "H": 11, "N": 1, "O": 1},
    "Lys": {"C": 6, "H": 12, "N": 2, "O": 1},
    "Met": {"C": 5, "H": 9, "N": 1, "O": 1, "S": 1},
    "Phe": {"C": 9, "H": 9, "N": 1, "O": 1},
    "Pro": {"C": 5, "H": 7, "N": 1, "O": 1},
    "Ser": {"C": 3, "H": 5, "N": 1, "O": 2},
    "Thr": {"C": 4, "H": 7, "N": 1, "O": 2},
    "Trp": {"C": 11, "H": 10, "N": 2, "O": 1},
    "Tyr": {"C": 9, "H": 9, "N": 1, "O": 2},
    "Val": {"C": 5, "H": 9, "N": 1, "O": 1},
}

_ONE_LETTER = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
}
_STEREO_PREFIX_RE = re.compile(r"^(?:l|d)[- ](?P<base>.+)$", re.IGNORECASE)


def _formula_mass(formula: dict[str, int]) -> float:
    return sum(_ATOMIC_MASSES[element] * count for element, count in formula.items())


FREE_MONOMER_MASSES: dict[str, float] = {
    name: _formula_mass(formula) for name, formula in _FORMULAS.items()
}


def _key(value: str) -> str:
    return unicodedata.normalize("NFKC", value.strip()).casefold()


_ALIASES = {_key(name): name for name in _FORMULAS}
_ALIASES.update({_key(letter): name for letter, name in _ONE_LETTER.items()})
for _name in _FORMULAS:
    _ALIASES[_key(f"L-{_name}")] = _name
    _ALIASES[_key(f"D-{_name}")] = _name


def canonical_monomer_name(value: str | None) -> str | None:
    """Return a supported canonical amino-acid name without changing source text."""

    if value is None or not value.strip():
        return None
    key = _key(value)
    direct = _ALIASES.get(key)
    if direct is not None:
        return direct
    stereo_match = _STEREO_PREFIX_RE.match(value.strip())
    if stereo_match:
        return _ALIASES.get(_key(stereo_match.group("base")))
    return None


def free_monomer_mass(value: str | None) -> float | None:
    """Return a formula-derived free-monomer mass for an exact supported alias."""

    canonical = canonical_monomer_name(value)
    return FREE_MONOMER_MASSES.get(canonical) if canonical is not None else None
