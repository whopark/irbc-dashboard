# iRBC Manufacturing Dashboard
# Run: streamlit run irbc_dashboard.py  (from 김교수님2/ or repo root)
import sys, pathlib, json
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from irbc_digital_twin import IRBCTwin, SCENARIOS
from irbc_charts import process_fig, comparison_fig, sensitivity_fig

st.set_page_config(page_title="iRBC Manufacturing Monitor", page_icon="🩸", layout="wide")
st.title("iRBC Manufacturing Process Monitor")
st.caption("ArtRed / HUDEP-2 immortalized RBC — KCABP ontology-calibrated digital twin")

# ── sidebar ──────────────────────────────────────────────────────────────────

_DESC = {
    "Baseline (current KCABP)":
        "ArtRed/HUDEP-2, 7L bioreactor, Fluoro-SAHA HDACi. Current program state (~35% enucleation).",
    "Next-gen HDACi + micropore extrusion":
        "Optimized Fluoro-SAHA + micropore filter (Kim Hyun-kyung). Targets >80% enucleation ceiling.",
    "Ferric Citrate / Apo-Transferrin Enhanced":
        "Iron-chelate loading (Lee Ji-yoon paper): ferric citrate + Def3.Fe3+. Maturation shortened to 3d.",
    "CRISPR Orai1/SH2B3-KO Enhanced Enucleation":
        "CRISPR gene editing of SOCE/Orai1 pathway (Park Chan-young). Enucleation ceiling ~65%.",
    "GMP 200L scale-up (Phase 2)":
        "200L bioreactor — KCABP Phase 2 (2028-2034). Longer expansion to seed vessel.",
    "Phase 2 Clinical-Grade (2028 target)":
        "Full stack: 200L + ferric citrate + next-gen HDACi + KO cell line.",
    "Rapid Protocol (BEL-A inspired, 13d)":
        "Speed-optimized based on Bristol BEL-A. 13-day total, research-grade enucleation.",
}

with st.sidebar:
    st.header("Scenario")
    scenario_name = st.selectbox("Select scenario", list(SCENARIOS.keys()))
    st.info(_DESC.get(scenario_name, ""))
    st.divider()
    st.header("Parameter Overrides")
    st.caption("Adjust on top of selected scenario")
    doubling_h  = st.slider("Doubling time (h)",      12.0, 30.0, 18.0, 0.5)
    max_enuc    = st.slider("HDACi enuc ceiling (%)", 20.0, 95.0, 38.0, 1.0)
    hb_target   = st.slider("Hb target (pg/cell)",    20.0, 35.0, 27.0, 0.5)
    st.divider()
    st.header("QC Gate Thresholds")
    qc_enuc_min = st.slider("Min enucleation (%)", 10.0, 80.0, 15.0, 1.0)
    qc_hb_min   = st.slider("Min hemoglobin (pg)",  5.0, 27.0, 10.0, 0.5)
    qc_viab_min = st.slider("Min viability (%)",   50.0, 99.0, 70.0, 1.0)

# ── run simulation ────────────────────────────────────────────────────────────

base_twin = SCENARIOS[scenario_name]()
twin = IRBCTwin({
    **base_twin.p,
    "expansion":   {**base_twin.p["expansion"],   "doubling_time_h": doubling_h},
    "enucleation": {**base_twin.p["enucleation"], "max_enuc": max_enuc / 100},
    "maturation":  {**base_twin.p["maturation"],  "hb_target": hb_target},
    "qc":          {**base_twin.p["qc"], "pass_enuc": qc_enuc_min / 100,
                    "pass_hb": qc_hb_min, "pass_viability": qc_viab_min / 100},
})
result = twin.run()
f, qc  = result["final"], result["qc"]

# ── KPI row ──────────────────────────────────────────────────────────────────

c1, c2, c3, c4, c5, c6 = st.columns(6)
c1.metric("Duration",    f"{f['days']} d")
c2.metric("Volume",      f"{f['volume_L']} L")
c3.metric("Total Cells", f["total_cells"])
c4.metric("Enucleation", f"{f['enuc_pct']}%",      "PASS" if qc["gates"]["enucleation"] else "FAIL",
          delta_color="normal" if qc["gates"]["enucleation"] else "inverse")
c5.metric("Hemoglobin",  f"{f['hb_pg']} pg",       "PASS" if qc["gates"]["hemoglobin"]  else "FAIL",
          delta_color="normal" if qc["gates"]["hemoglobin"]  else "inverse")
c6.metric("Viability",   f"{f['viability_pct']}%", "PASS" if qc["gates"]["viability"]   else "FAIL",
          delta_color="normal" if qc["gates"]["viability"]   else "inverse")

st.divider()

# ── stage timeline ────────────────────────────────────────────────────────────

_STAGE_KEYS   = ("expansion", "differentiation", "maturation", "enucleation")
_STAGE_NAMES  = ["Expansion", "Differentiation", "Maturation", "Enucleation"]
_STAGE_COLORS = ["#28a745", "#007bff", "#fd7e14", "#dc3545"]
stage_days    = [twin.p[k]["duration_h"] / 24 for k in _STAGE_KEYS]
total_days    = sum(stage_days)

st.subheader("Stage Timeline")
for col, name, days, color in zip(st.columns(4), _STAGE_NAMES, stage_days, _STAGE_COLORS):
    col.markdown(
        f"<div style='background:{color};color:white;border-radius:8px;"
        f"padding:12px;text-align:center;font-weight:bold'>"
        f"{name}<br><span style='font-size:1.4em'>{days:.0f}d</span>"
        f"<br><small>{days / total_days * 100:.0f}% of total</small></div>",
        unsafe_allow_html=True,
    )
st.divider()

# ── tabs ──────────────────────────────────────────────────────────────────────

tab_proc, tab_qc, tab_cmp, tab_sens = st.tabs(
    ["📈 Process Charts", "🔬 QC Report", "📊 Scenario Compare", "🎛 Sensitivity"]
)

with tab_proc:
    fig = process_fig(twin, result, hb_target, qc_hb_min, qc_enuc_min, qc_viab_min)
    st.pyplot(fig); plt.close(fig)

with tab_qc:
    left, right = st.columns([1, 2])
    with left:
        st.markdown(f"### {'🟢 PASS' if qc['pass'] else '🔴 FAIL'}")
        for gate, passed in qc["gates"].items():
            st.markdown(f"{'✅' if passed else '❌'} **{gate}**")
    with right:
        st.dataframe(pd.DataFrame({
            "P50 (mmHg)":       [f["p50_mmhg"],       "24–29",                "✅" if qc["gates"]["p50"]         else "❌"],
            "2,3-DPG":          [f["dpg"],             ">5.0 umol/gHb",        "✅" if qc["gates"]["2,3-DPG"]     else "❌"],
            "ATP":              [f["atp"],             ">4.0 umol/gHb",        "—"],
            "Enucleation":      [f"{f['enuc_pct']}%",  f">{qc_enuc_min:.0f}%", "✅" if qc["gates"]["enucleation"] else "❌"],
            "Hemoglobin":       [f"{f['hb_pg']} pg",   f">{qc_hb_min:.0f} pg", "✅" if qc["gates"]["hemoglobin"]  else "❌"],
            "Viability":        [f"{f['viability_pct']}%", f">{qc_viab_min:.0f}%", "✅" if qc["gates"]["viability"] else "❌"],
            "Enucleated Cells": [f["enuc_cells"],      "—",                    "—"],
        }, index=["Value", "Target", "Status"]).T, use_container_width=True)
    st.download_button("⬇ Download QC JSON",
                       json.dumps({"scenario": scenario_name, "final": f, "qc": qc}, indent=2),
                       "irbc_qc_report.json", "application/json")

# ponytail: cache all-scenario runs — 7 simulations only on first load, not every slider move
@st.cache_data
def _compare_all():
    rows = []
    for name, factory in SCENARIOS.items():
        r  = factory().run()
        f2 = r["final"]
        rows.append({"Scenario": name, "Days": f2["days"], "Volume (L)": f2["volume_L"],
                     "Enuc (%)": f2["enuc_pct"], "Hb (pg)": f2["hb_pg"],
                     "Viability (%)": f2["viability_pct"], "Total Cells": f2["total_cells"],
                     "Enuc Cells": f2["enuc_cells"], "QC": "PASS" if r["qc"]["pass"] else "FAIL"})
    return pd.DataFrame(rows).set_index("Scenario")

with tab_cmp:
    df = _compare_all()
    def _hl(row):
        bg = "#d4edda" if row["QC"] == "PASS" else "#f8d7da"
        tx = "#155724" if row["QC"] == "PASS" else "#721c24"
        s  = [f"background-color:{bg};color:{tx}"] * len(row)
        qc_style = ("background-color:#28a745;color:white;font-weight:bold" if row["QC"] == "PASS"
                    else "background-color:#dc3545;color:white;font-weight:bold")
        s[row.index.get_loc("QC")] = qc_style
        return s
    st.dataframe(df.style.apply(_hl, axis=1).format(
        {"Days": "{:.0f}", "Volume (L)": "{:.0f}", "Enuc (%)": "{:.1f}",
         "Hb (pg)": "{:.1f}", "Viability (%)": "{:.1f}"}), use_container_width=True)
    fig2 = comparison_fig(df); st.pyplot(fig2); plt.close(fig2)

with tab_sens:
    st.caption("Sweep HDACi ceiling to find the enucleation break-even point for the selected scenario and QC thresholds.")
    sweep    = [v / 100 for v in range(20, 96, 5)]
    enuc_sw, qc_sw = [], []
    for v in sweep:
        t = IRBCTwin({**base_twin.p,
                      "expansion":   {**base_twin.p["expansion"],   "doubling_time_h": doubling_h},
                      "enucleation": {**base_twin.p["enucleation"], "max_enuc": v},
                      "maturation":  {**base_twin.p["maturation"],  "hb_target": hb_target},
                      "qc":          {**base_twin.p["qc"], "pass_enuc": qc_enuc_min / 100,
                                      "pass_hb": qc_hb_min, "pass_viability": qc_viab_min / 100}})
        r = t.run()
        enuc_sw.append(r["final"]["enuc_pct"])
        qc_sw.append(r["qc"]["pass"])
    fig3 = sensitivity_fig([v * 100 for v in sweep], enuc_sw, qc_sw, max_enuc, qc_enuc_min)
    st.pyplot(fig3); plt.close(fig3)
    st.dataframe(pd.DataFrame({"HDACi Ceiling (%)": [v * 100 for v in sweep],
                                "Final Enuc (%)": enuc_sw,
                                "QC": ["PASS" if p else "FAIL" for p in qc_sw]}),
                 use_container_width=True, hide_index=True)
