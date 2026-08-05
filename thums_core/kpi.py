"""The KPI registry — every metric declared once, computed in one place.

Why a registry rather than ad-hoc columns:

* `direction` makes a genuine Pareto front computable. The legacy
  `Pareto_for_THUMS_DoE` notebook contains no dominance test at all.
* `depends_on` records algebraic dependence between metrics, so a
  "consensus across KPIs" cannot silently count the same physical mechanism
  two or three times (see `independent_subset`).
* `label` is the single LaTeX string used by every figure and table, so the
  manuscript cannot drift from the code.

Deliberate change from the legacy code: no KPI is baked into another. The legacy
`weighted_RTE = RTE_real * m_m_eff` is split back into `eta_rte` and `eps_pcm`.
That product is not an efficiency of any control volume -- it falls when PCM
utilization falls even though electricity in and out are unchanged. The
efficiency/utilization trade-off belongs on a Pareto front, not in a
scalarisation with an unstated weighting.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Iterable

MAX = "max"
MIN = "min"


@dataclass(frozen=True)
class KPI:
    name: str
    label: str
    unit: str
    direction: str
    group: str
    doc: str = ""
    depends_on: tuple[str, ...] = ()
    fn: Callable | None = None

    def __post_init__(self):
        if self.direction not in (MAX, MIN):
            raise ValueError(f"{self.name}: direction must be 'max' or 'min'")


def _kpi(*args, **kw):
    k = KPI(*args, **kw)
    REGISTRY[k.name] = k
    return k


REGISTRY: dict[str, KPI] = {}

# ----------------------------------------------------------------- thermodynamic
_kpi("eta_rte", r"$\eta_{\mathrm{RTE}}$", "-", MAX, "thermodynamic",
     doc="Round-trip efficiency: (W_out - W_pump,dc)*t_dc / ((W_in + W_pump,ch)*t_ch). "
         "The actual electricity-in to electricity-out ratio, nothing else folded in.")
_kpi("eta_rte_nopump", r"$\eta_{\mathrm{RTE}}^{0}$", "-", MAX, "thermodynamic",
     doc="Round-trip efficiency excluding parasitic pumping. Reported for comparison "
         "with the conference paper, which used this quantity.")
_kpi("cop_hp", r"$\mathrm{COP}_{\mathrm{HP}}$", "-", MAX, "thermodynamic",
     doc="HTHP coefficient of performance.")
_kpi("eta_orc", r"$\eta_{\mathrm{ORC}}$", "-", MAX, "thermodynamic",
     doc="ORC thermal efficiency.")
_kpi("eta_ex", r"$\eta_{\mathrm{ex}}$", "-", MAX, "thermodynamic",
     doc="Exergy efficiency of the complete power-to-heat-to-power chain.")

# ---------------------------------------------------------------------- storage
_kpi("eps_pcm", r"$\varepsilon_{\mathrm{PCM}}$", "-", MAX, "storage",
     doc="PCM utilization: melted mass / total mass. Legacy name: m_m_eff. "
         "Reported on its own -- never multiplied into an efficiency.")
_kpi("E_well", r"$E_{\mathrm{well}}$", "MWh", MAX, "storage",
     doc="Energy stored per well.")
_kpi("rho_E", r"$\rho_E$", "kWh/m$^3$", MAX, "storage",
     depends_on=("E_well",),
     doc="Energy density per unit well volume.")
_kpi("self_discharge", r"$\dot{\sigma}$", "%/day", MIN, "storage",
     doc="Self-discharge rate. Meaningful only once formation heat loss is modelled; "
         "with the legacy flat 5% surplus this is structurally zero.")

# -------------------------------------------------------------- sizing / parasitics
_kpi("N_wells", r"$N_{\mathrm{wells}}$", "-", MIN, "sizing",
     doc="Wells required for the stated duty (charging side).")
_kpi("wells_per_MW", r"$N/\dot{W}$", "wells/MW", MIN, "sizing",
     depends_on=("N_wells",),
     doc="Specific well count.")
_kpi("f_pump", r"$f_{\mathrm{pump}}$", "-", MIN, "sizing",
     doc="Pumping work as a fraction of gross work output.")

# --------------------------------------------------------------------- economic
_kpi("c_pcm", r"$c_{\mathrm{PCM}}$", "USD/kWh", MIN, "economic",
     doc="PCM cost per kWh of stored energy.")


def names(group: str | None = None) -> list[str]:
    return [k.name for k in REGISTRY.values() if group is None or k.group == group]


def labels(kpis: Iterable[str]) -> dict[str, str]:
    return {n: REGISTRY[n].label for n in kpis}


def independent_subset(kpis: Iterable[str]) -> list[str]:
    """Drop KPIs that are algebraically derived from others in the same set.

    Guards the "consensus ranking across KPIs" analysis. Averaging normalised
    effects over {eta_rte, eps_pcm, eta_rte*eps_pcm} weights the same mechanism
    twice; this returns a set where that cannot happen.
    """
    kpis = list(kpis)
    keep = []
    for n in kpis:
        deps = set(REGISTRY[n].depends_on)
        if deps & set(kpis):
            continue
        keep.append(n)
    return keep


def check_independent(kpis: Iterable[str]) -> None:
    kpis = list(kpis)
    ok = independent_subset(kpis)
    dropped = [n for n in kpis if n not in ok]
    if dropped:
        raise ValueError(
            "KPIs are algebraically dependent and cannot be averaged into a "
            f"consensus ranking: {dropped}. Use kpi.independent_subset() or "
            "report them separately."
        )
