"""
iRBC Digital Twin - immortalized RBC manufacturing simulation.

Stages from KCABP ontology (graphify-out/graph.json):
  1. Expansion      - bioreactor scale-up (1L->7L, ArtRed/HUDEP-2)
  2. Differentiation - erythroid induction, GPA up, cKit down
  3. Maturation     - hemoglobinization via ferric citrate + apo-transferrin
  4. Enucleation    - HDACi (Fluoro-SAHA/M344/TSA) + micropore extrusion
  5. QC harvest     - P50, 2,3-DPG, ATP, enucleation rate gate

Parameter sources (ontology community labels):
  - ArtBlood Bioreactor Scale-Up: 1e7/mL @1L, 3-5e6/mL @7L
  - Iron Metabolism & Hemoglobinization: deferiprone, ferric citrate
  - HDACi-Mediated Enucleation: ~35% current baseline (DRAQ-/CD235a+)
  - Manufacturing Process Ontology: stage gate QC criteria
"""

from __future__ import annotations
import math, random
from dataclasses import dataclass, field
from typing import Literal

Stage = Literal["expansion", "differentiation", "maturation", "enucleation", "harvest"]

@dataclass
class SimState:
    time_h: float = 0.0
    volume_L: float = 1.0          # working volume (L)
    density: float = 1e5           # cells/mL
    gpa: float = 0.05              # CD235a (0-1)
    ckit: float = 0.90             # CD117  (0-1)
    hemoglobin_pg: float = 0.5     # pg/cell, adult target ~27 pg
    enuc_rate: float = 0.0         # fraction DRAQ-/CD235a+
    viability: float = 0.98
    p50_mmhg: float = 0.0
    dpg: float = 0.0               # 2,3-DPG umol/gHb
    atp: float = 0.0               # umol/gHb
    stage: Stage = "expansion"
    history: list = field(default_factory=list, repr=False)

    def snapshot(self) -> dict:
        return {
            "time_h": round(self.time_h, 1),
            "stage": self.stage,
            "volume_L": round(self.volume_L, 1),
            "density_M": round(self.density / 1e6, 3),
            "viability": round(self.viability, 3),
            "GPA": round(self.gpa, 3),
            "cKit": round(self.ckit, 3),
            "Hb_pg": round(self.hemoglobin_pg, 2),
            "enuc_rate": round(self.enuc_rate, 3),
        }


PARAMS = {
    "expansion": {
        "doubling_time_h": 18.0,     # ArtRed/HUDEP-2 typical
        "max_density": 1.2e7,        # 1e7/mL cap (ArtBlood 1L data)
        "scale_trigger": 5e6,        # seed 7L when >5e6/mL
        "scale_volume_L": 7.0,
        "duration_h": 7 * 24,
    },
    "differentiation": {
        "duration_h": 6 * 24,
        "gpa_rate": 0.003,           # per hour, sigmoid
        "ckit_decay": 0.003,
        "hb_rate": 0.05,             # pg/cell/h
        "viability_loss": 0.001,
    },
    "maturation": {
        # ferric citrate + apo-transferrin iron loading (Lee Ji-yoon paper)
        "duration_h": 4 * 24,
        "hb_rate": 0.20,             # pg/cell/h, boosted by iron
        "hb_target": 27.0,           # adult RBC ~27 pg
        "viability_loss": 0.0005,
    },
    "enucleation": {
        # HDACi baseline: Fluoro-SAHA or M344 (Community 9)
        # KCABP current ~35%, target >80% clinical
        "duration_h": 3 * 24,
        "hdaci_k": 0.018,            # rate constant /h (Fluoro-SAHA)
        "max_enuc": 0.38,            # ponytail: current ceiling 35-38%
        "viability_loss": 0.0008,
    },
    "qc": {
        # QC panel: Lee Dong-jun community (P50, 2,3-DPG, ATP, blood type)
        "p50_mean": 26.5, "p50_sd": 0.8,
        "dpg_mean": 5.2,  "dpg_sd": 0.4,
        "atp_mean": 4.8,  "atp_sd": 0.3,
        "pass_enuc": 0.15,
        "pass_viability": 0.70,
        "pass_hb": 10.0,
    },
}

STAGE_ORDER = ["expansion", "differentiation", "maturation", "enucleation", "harvest"]


def _step(s: SimState, dt: float, p: dict) -> None:
    if s.stage == "expansion":
        gr = math.log(2) / p["expansion"]["doubling_time_h"]
        K  = p["expansion"]["max_density"]
        s.density += gr * s.density * (1 - s.density / K) * dt
        s.density  = max(s.density, 1e3)
        if s.volume_L < p["expansion"]["scale_volume_L"] and s.density >= p["expansion"]["scale_trigger"]:
            s.density *= s.volume_L / p["expansion"]["scale_volume_L"]
            s.volume_L = p["expansion"]["scale_volume_L"]

    elif s.stage == "differentiation":
        ep = p["differentiation"]
        s.gpa = min(1.0, s.gpa + ep["gpa_rate"] * s.gpa * (1 - s.gpa) * dt)
        s.ckit = max(0.0, s.ckit - ep["ckit_decay"] * s.ckit * dt)
        s.hemoglobin_pg = min(10.0, s.hemoglobin_pg + ep["hb_rate"] * dt)
        s.viability = max(0.0, s.viability - ep["viability_loss"] * dt)

    elif s.stage == "maturation":
        mp = p["maturation"]
        rem = mp["hb_target"] - s.hemoglobin_pg
        s.hemoglobin_pg += mp["hb_rate"] * (rem / mp["hb_target"]) * dt
        s.viability = max(0.0, s.viability - mp["viability_loss"] * dt)

    elif s.stage == "enucleation":
        ep = p["enucleation"]
        s.enuc_rate += ep["hdaci_k"] * (ep["max_enuc"] - s.enuc_rate) * dt
        s.viability = max(0.0, s.viability - ep["viability_loss"] * dt)


def _qc(s: SimState, p: dict) -> dict:
    rng = random.Random(42)
    qp = p["qc"]
    s.p50_mmhg = rng.gauss(qp["p50_mean"], qp["p50_sd"])
    s.dpg = rng.gauss(qp["dpg_mean"], qp["dpg_sd"])
    s.atp = rng.gauss(qp["atp_mean"], qp["atp_sd"])
    gates = {
        "enucleation": s.enuc_rate >= qp["pass_enuc"],
        "viability":   s.viability >= qp["pass_viability"],
        "hemoglobin":  s.hemoglobin_pg >= qp["pass_hb"],
        "p50":         24.0 <= s.p50_mmhg <= 29.0,
        "2,3-DPG":     s.dpg >= 4.0,
    }
    return {"pass": all(gates.values()), "gates": gates}


class IRBCTwin:
    """
    Discrete-time digital twin for immortalized RBC manufacturing.

        twin = IRBCTwin()
        result = twin.run()
        twin.plot(result)

    Override any stage parameter:
        twin = IRBCTwin({"enucleation": {"max_enuc": 0.85}})  # next-gen HDACi
    """
    def __init__(self, overrides: dict | None = None, dt_h: float = 1.0):
        self.p = {k: dict(v) for k, v in PARAMS.items()}
        if overrides:
            for stage, vals in overrides.items():
                self.p.setdefault(stage, {}).update(vals)
        self.dt = dt_h

    def run(self, seed_density: float = 1e5, seed_volume_L: float = 1.0) -> dict:
        s = SimState(density=seed_density, volume_L=seed_volume_L)
        durations = {st: self.p[st]["duration_h"] for st in STAGE_ORDER if st != "harvest"}
        idx, elapsed = 0, 0.0
        s.history.append(s.snapshot())

        while s.stage != "harvest":
            _step(s, self.dt, self.p)
            s.time_h += self.dt
            elapsed  += self.dt
            if s.time_h % 24 < self.dt:
                s.history.append(s.snapshot())
            if elapsed >= durations.get(STAGE_ORDER[idx], 0):
                idx += 1
                elapsed = 0.0
                s.stage = STAGE_ORDER[idx]

        qc = _qc(s, self.p)
        s.history.append(s.snapshot())
        total = s.density * s.volume_L * 1000
        return {
            "qc": qc,
            "final": {
                "days":          round(s.time_h / 24, 1),
                "volume_L":      s.volume_L,
                "total_cells":   f"{total:.2e}",
                "enuc_cells":    f"{total * s.enuc_rate:.2e}",
                "enuc_pct":      round(s.enuc_rate * 100, 1),
                "hb_pg":         round(s.hemoglobin_pg, 1),
                "viability_pct": round(s.viability * 100, 1),
                "p50_mmhg":      round(s.p50_mmhg, 1),
                "dpg":           round(s.dpg, 2),
                "atp":           round(s.atp, 2),
            },
            "history": s.history,
        }

    def plot(self, result: dict) -> None:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("pip install matplotlib to plot"); return

        h = result["history"]
        days  = [r["time_h"] / 24 for r in h]
        dens  = [r["density_M"]   for r in h]
        hb    = [r["Hb_pg"]       for r in h]
        enuc  = [r["enuc_rate"] * 100 for r in h]
        gpa   = [r["GPA"] * 100   for r in h]
        viab  = [r["viability"] * 100 for r in h]

        bounds = [0,
                  self.p["expansion"]["duration_h"] / 24,
                  (self.p["expansion"]["duration_h"] + self.p["differentiation"]["duration_h"]) / 24,
                  (self.p["expansion"]["duration_h"] + self.p["differentiation"]["duration_h"] + self.p["maturation"]["duration_h"]) / 24,
                  days[-1]]
        names  = ["Expansion", "Differentiation", "Maturation", "Enucleation"]
        colors = ["#d4edda", "#cce5ff", "#fff3cd", "#f8d7da"]

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle("iRBC Digital Twin - ArtRed/HUDEP-2 Manufacturing Simulation (KCABP ontology)", fontsize=10)

        def shade(ax):
            for i, (a, b) in enumerate(zip(bounds, bounds[1:])):
                ax.axvspan(a, b, alpha=0.25, color=colors[i])
            for d in bounds[1:-1]:
                ax.axvline(d, color="gray", lw=0.7, ls="--")

        ax = axes[0, 0]
        ax.plot(days, dens, "b-", lw=2); shade(ax)
        ax.set(ylabel="Density (x10^6/mL)", title="Cell Density", xlabel="Day")

        ax = axes[0, 1]
        ax.plot(days, hb, "r-", lw=2); shade(ax)
        ax.axhline(27, color="red", lw=1, ls=":", label="Adult target")
        ax.axhline(self.p["qc"]["pass_hb"], color="orange", lw=1, ls=":", label="QC min")
        ax.legend(fontsize=7)
        ax.set(ylabel="Hemoglobin (pg/cell)", title="Hemoglobinization", xlabel="Day")

        ax = axes[1, 0]
        ax.plot(days, gpa,  "g-",  lw=2, label="GPA (CD235a %)")
        ax.plot(days, viab, "k--", lw=1.5, label="Viability %")
        shade(ax); ax.legend(fontsize=8)
        ax.set(ylabel="%", title="Differentiation Markers & Viability", xlabel="Day")

        ax = axes[1, 1]
        ax.plot(days, enuc, "m-", lw=2); shade(ax)
        ax.axhline(self.p["enucleation"]["max_enuc"] * 100, color="purple", lw=1, ls=":", label=f"Ceiling {self.p['enucleation']['max_enuc']*100:.0f}%")
        ax.axhline(self.p["qc"]["pass_enuc"] * 100, color="orange", lw=1, ls=":", label=f"QC min {self.p['qc']['pass_enuc']*100:.0f}%")
        ax.legend(fontsize=7)
        ax.set(ylabel="%", title="Enucleation Rate (HDACi)", xlabel="Day")

        handles = [plt.Rectangle((0,0),1,1,fc=c,alpha=0.5) for c in colors]
        axes[0,0].legend(handles, names, fontsize=7, loc="upper left")

        plt.tight_layout()
        out = "graphify-out/irbc_twin_plot.png"
        plt.savefig(out, dpi=150, bbox_inches="tight")
        plt.show()
        print("Saved:", out)


# -- scenarios -----------------------------------------------------------------
# Each returns a fresh IRBCTwin with a description attribute.

def _twin(overrides: dict, label: str) -> IRBCTwin:
    t = IRBCTwin(overrides)
    t.label = label
    return t


def scenario_baseline() -> IRBCTwin:
    # Current KCABP state: ArtRed/HUDEP-2, 7L bioreactor, Fluoro-SAHA HDACi ~35% enuc
    return _twin({}, "Baseline (current KCABP)")


def scenario_next_gen_hdaci(_twin_unused=None) -> IRBCTwin:
    # Fluoro-SAHA optimized + micropore filter extrusion (Kim Hyun-kyung BONE-V work)
    # Enucleation ceiling raised toward clinical target >80%
    return _twin({
        "enucleation": {**PARAMS["enucleation"], "max_enuc": 0.85, "hdaci_k": 0.025},
    }, "Next-gen HDACi + micropore extrusion")


def scenario_gmp_200l(_twin_unused=None) -> IRBCTwin:
    # 200L GMP scale-up as planned in KCABP Phase 2 (2028-2034)
    # Longer expansion required to seed the larger vessel
    return _twin({
        "expansion": {**PARAMS["expansion"], "scale_volume_L": 200.0, "duration_h": 10 * 24},
    }, "GMP 200L scale-up (Phase 2)")


def scenario_ferric_citrate() -> IRBCTwin:
    # Enhanced iron loading via ferric citrate + apo-transferrin (Lee Ji-yoon, 이지윤 박사 논문)
    # Faster hemoglobinization: hb_rate boosted, maturation shortened by 1 day
    # Deferiprone (Def3.Fe3+) iron-chelate pre-loading accelerates heme synthesis
    return _twin({
        "maturation": {**PARAMS["maturation"],
                       "hb_rate": 0.38,           # ~2x boost from iron chelate loading
                       "duration_h": 3 * 24,      # 3d instead of 4d
                       "viability_loss": 0.0003}, # better medium tolerance
    }, "Ferric Citrate / Apo-Transferrin Enhanced")


def scenario_crispr_enucleation() -> IRBCTwin:
    # CRISPR gene editing: Orai1-KO + SH2B3-KO hiPSC (Park Chan-young, Community 12)
    # SOCE/Orai1 pathway disruption promotes nuclear expulsion
    # Combined with chemical enucleation ceiling rises to ~65%
    return _twin({
        "enucleation": {**PARAMS["enucleation"],
                        "max_enuc": 0.65,          # Orai1/Orai2 KO raises ceiling
                        "hdaci_k":  0.022,
                        "duration_h": 4 * 24},     # extra day for gene-edited line
        "differentiation": {**PARAMS["differentiation"],
                             "viability_loss": 0.0007},  # KO lines more robust
    }, "CRISPR Orai1/SH2B3-KO Enhanced Enucleation")


def scenario_phase2_clinical() -> IRBCTwin:
    # Full Phase 2 target: 200L + next-gen HDACi + ferric citrate + CRISPR
    # Represents the 2028-2034 clinical-grade production target
    return _twin({
        "expansion":   {**PARAMS["expansion"],
                        "scale_volume_L": 200.0, "duration_h": 10 * 24},
        "maturation":  {**PARAMS["maturation"],
                        "hb_rate": 0.38, "duration_h": 3 * 24,
                        "viability_loss": 0.0003},
        "enucleation": {**PARAMS["enucleation"],
                        "max_enuc": 0.88, "hdaci_k": 0.028,
                        "viability_loss": 0.0005},
        "differentiation": {**PARAMS["differentiation"],
                             "viability_loss": 0.0006},
    }, "Phase 2 Clinical-Grade (2028 target)")


def scenario_rapid_protocol() -> IRBCTwin:
    # Speed-optimized protocol: faster-growing line, shorter stages
    # Trade-off: lower enucleation ceiling, acceptable for research-grade product
    # Inspired by BEL-A (Bristol) rapid ex vivo protocol
    return _twin({
        "expansion":       {**PARAMS["expansion"],
                            "doubling_time_h": 14.0,  # BEL-A grows faster
                            "duration_h": 5 * 24},
        "differentiation": {**PARAMS["differentiation"],
                            "duration_h": 4 * 24,
                            "viability_loss": 0.0015},
        "maturation":      {**PARAMS["maturation"],
                            "duration_h": 2 * 24, "hb_rate": 0.15},
        "enucleation":     {**PARAMS["enucleation"],
                            "duration_h": 2 * 24, "max_enuc": 0.28},
    }, "Rapid Protocol (BEL-A inspired, 13d)")


# Registry: name -> factory function (no args needed)
SCENARIOS = {
    s().label: s for s in [
        scenario_baseline,
        scenario_next_gen_hdaci,
        scenario_ferric_citrate,
        scenario_crispr_enucleation,
        scenario_gmp_200l,
        scenario_phase2_clinical,
        scenario_rapid_protocol,
    ]
}


# -- self-check ----------------------------------------------------------------

def demo():
    twin = IRBCTwin()
    r = twin.run()
    f = r["final"]
    print("=== iRBC Digital Twin - baseline run ===")
    for k, v in f.items():
        print(f"  {k:20s}: {v}")
    print(f"  QC result           : {'PASS' if r['qc']['pass'] else 'FAIL'}")
    print(f"  QC gates            : {r['qc']['gates']}")

    ng = scenario_next_gen_hdaci(twin)
    ng_r = ng.run()
    print(f"\n=== Next-gen HDACi scenario ===")
    print(f"  enuc_pct            : {ng_r['final']['enuc_pct']}%  (target >80%)")
    print(f"  enuc_cells          : {ng_r['final']['enuc_cells']}")

    assert r["final"]["days"] > 15
    assert r["final"]["enuc_pct"] > 20
    assert ng_r["final"]["enuc_pct"] > r["final"]["enuc_pct"]
    assert 24 <= r["final"]["p50_mmhg"] <= 29
    print("\nAll assertions passed.")


if __name__ == "__main__":
    demo()
    twin = IRBCTwin()
    twin.plot(twin.run())