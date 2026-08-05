"""Design of experiments — designs that state their own alias structure.

The legacy `build_128run_fractional_16factor` assigned all nine added factors to
*two-factor* interaction columns (H=AB, I=AC, ... P=BE). Enumerating the defining
relation gives twelve words of length three, i.e. **Resolution III**: every main
effect is aliased with two-factor interactions, and `T_m_C` is inseparable from
`rho_m_s*h_m` while `h_m` is inseparable from `T_m_C*rho_m_s` -- the same
interaction, so no effect-sparsity argument separates them.

Assigning the same nine factors to *odd-order* columns gives Resolution IV in
exactly the same 128 runs. This module does that, and -- more importantly --
computes and carries the defining relation so it can be printed in the paper
instead of assumed.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# Resolution IV generators for 2^(16-9), verified: shortest word length 4.
RES4_GENERATORS_16 = {
    "H": "ABC", "I": "ABD", "J": "ABE", "K": "ABF", "L": "ABG",
    "M": "ACD", "N": "ACE", "O": "ACF", "P": "ACG",
}

# The legacy assignment, kept so the old design can be reproduced for comparison.
RES3_GENERATORS_16 = {
    "H": "AB", "I": "AC", "J": "AD", "K": "AE", "L": "AF",
    "M": "AG", "N": "BC", "O": "BD", "P": "BE",
}


@dataclass(frozen=True)
class AliasStructure:
    resolution: int
    word_lengths: dict
    words: tuple
    main_effect_aliases: dict

    def summary(self) -> str:
        lines = [f"Resolution {'I' * self.resolution if self.resolution < 4 else 'IV'} "
                 f"(shortest word length {self.resolution})",
                 f"word-length distribution: {self.word_lengths}"]
        short = [w for w in self.words if len(w) == self.resolution]
        lines.append(f"defining words of length {self.resolution} "
                     f"({len(short)}): {' '.join(short)}")
        return "\n".join(lines)


def _columns(n_base: int, generators: dict, base_names: str):
    base = np.array(np.meshgrid(*[[-1, 1]] * n_base)).T.reshape(-1, n_base)
    cols = {base_names[i]: base[:, i] for i in range(n_base)}
    for letter, word in generators.items():
        c = np.ones(base.shape[0], dtype=int)
        for ch in word:
            c = c * cols[ch]
        cols[letter] = c
    return cols


def alias_structure(cols: dict, max_word: int = 8) -> AliasStructure:
    """Enumerate the defining relation and the two-factor aliases of each main effect."""
    letters = list(cols)
    n = len(next(iter(cols.values())))
    ones = np.ones(n, dtype=int)
    words = []
    for r in range(1, max_word + 1):
        for combo in itertools.combinations(letters, r):
            p = ones.copy()
            for c in combo:
                p = p * cols[c]
            if np.all(p == 1):
                words.append("".join(combo))
    lengths = {}
    for w in words:
        lengths[len(w)] = lengths.get(len(w), 0) + 1
    res = min(lengths) if lengths else max_word + 1

    aliases = {}
    for m in letters:
        al = []
        for a, b in itertools.combinations(letters, 2):
            if a == m or b == m:
                continue
            if np.array_equal(cols[m], cols[a] * cols[b]):
                al.append((a, b))
        aliases[m] = al
    return AliasStructure(res, dict(sorted(lengths.items())), tuple(words), aliases)


@dataclass
class FractionalFactorial:
    """2^(k-p) design that knows and reports its own alias structure.

    factor_levels: {name: (low, high)} in physical units, in factor order.
    """

    factor_levels: dict
    n_base: int = 7
    generators: dict = field(default_factory=lambda: dict(RES4_GENERATORS_16))
    n_centre: int = 1          # deterministic model: replicates carry zero information
    study_type: str = "FractionalFactorial"

    def __post_init__(self):
        self.factors = list(self.factor_levels)
        self.letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[: len(self.factors)]
        self._cols = _columns(self.n_base, self.generators, self.letters)
        self.alias = alias_structure(self._cols)
        if self.alias.resolution < 4:
            raise ValueError(
                f"Design is Resolution {self.alias.resolution}: main effects are "
                "aliased with two-factor interactions. Use odd-order generators "
                "(RES4_GENERATORS_16) or pass generators explicitly and accept it "
                "by constructing with FractionalFactorial.allow_low_resolution()."
            )

    @classmethod
    def allow_low_resolution(cls, *a, **kw):
        obj = cls.__new__(cls)
        object.__setattr__ if False else None
        # bypass the resolution guard, for reproducing the legacy design only
        dataclass_fields = {"factor_levels", "n_base", "generators", "n_centre", "study_type"}
        for k in dataclass_fields:
            setattr(obj, k, kw.get(k, getattr(cls, k, None)))
        if a:
            obj.factor_levels = a[0]
        obj.factors = list(obj.factor_levels)
        obj.letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[: len(obj.factors)]
        obj._cols = _columns(obj.n_base, obj.generators, obj.letters)
        obj.alias = alias_structure(obj._cols)
        return obj

    @property
    def coded(self) -> np.ndarray:
        return np.column_stack([self._cols[l] for l in self.letters])

    def design(self) -> pd.DataFrame:
        coded = self.coded
        rows = []
        for r in coded:
            row = {}
            for j, f in enumerate(self.factors):
                lo, hi = self.factor_levels[f]
                row[f] = 0.5 * (lo + hi) + 0.5 * (hi - lo) * r[j]
                row["coded_" + f] = float(r[j])
            row["is_centre"] = False
            rows.append(row)
        for _ in range(self.n_centre):
            row = {}
            for f in self.factors:
                lo, hi = self.factor_levels[f]
                row[f] = 0.5 * (lo + hi)
                row["coded_" + f] = 0.0
            row["is_centre"] = True
            rows.append(row)
        return pd.DataFrame(rows)

    def metadata(self) -> dict:
        return {
            "study_type": self.study_type,
            "n_factors": len(self.factors),
            "n_runs": int(self.coded.shape[0]) + self.n_centre,
            "n_centre": self.n_centre,
            "resolution": self.alias.resolution,
            "generators": {k: v for k, v in self.generators.items()},
            "factor_letters": dict(zip(self.factors, self.letters)),
            "defining_words_len_min": [w for w in self.alias.words
                                       if len(w) == self.alias.resolution],
            "word_length_distribution": self.alias.word_lengths,
            "main_effect_2fi_aliases": {
                self.factors[self.letters.index(m)]: [
                    (self.factors[self.letters.index(a)],
                     self.factors[self.letters.index(b)]) for a, b in v
                ] for m, v in self.alias.main_effect_aliases.items()
            },
            "factor_levels": {k: list(v) for k, v in self.factor_levels.items()},
            "note": (
                "Centre points are NOT replicates: the simulator is deterministic, so "
                "repeated centre runs give exactly zero pure error and cannot support "
                "an F test for curvature. Curvature must be assessed against a "
                "lack-of-fit criterion, not replicate error."
            ),
        }
