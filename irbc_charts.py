"""Matplotlib figure builders for iRBC dashboard — no Streamlit calls."""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

_SHADE  = ["#28a74522", "#007bff22", "#fd7e1422", "#dc354522"]
_NAMES  = ["Expansion", "Differentiation", "Maturation", "Enucleation"]
_COLORS = ["#28a745",   "#007bff",   "#fd7e14",   "#dc3545"]


def _bounds(twin):
    b = [0]
    for k in ("expansion", "differentiation", "maturation", "enucleation"):
        b.append(b[-1] + twin.p[k]["duration_h"] / 24)
    return b


def _dark(ax):
    ax.set_facecolor("#1a1d27")
    ax.tick_params(colors="white")
    for item in [ax.xaxis.label, ax.yaxis.label, ax.title]:
        item.set_color("white")
    for spine in ax.spines.values():
        spine.set_edgecolor("#444")


def _shade(ax, bounds):
    for i, (a, b) in enumerate(zip(bounds, bounds[1:])):
        ax.axvspan(a, b, alpha=0.35, color=_SHADE[i])
    for d in bounds[1:-1]:
        ax.axvline(d, color="#555", lw=0.8, ls="--")


def process_fig(twin, result, hb_target, qc_hb_min, qc_enuc_min, qc_viab_min):
    h  = result["history"]
    ts = [r["time_h"] / 24 for r in h]
    b  = _bounds(twin)

    fig, axes = plt.subplots(2, 2, figsize=(13, 7))
    fig.patch.set_facecolor("#0e1117")
    for ax in axes.flat:
        _dark(ax)

    ax = axes[0, 0]
    ax.plot(ts, [r["density_M"] for r in h], color="#4fc3f7", lw=2)
    _shade(ax, b)
    ax.legend(handles=[mpatches.Patch(fc=c, alpha=0.5, label=n)
                       for c, n in zip(_COLORS, _NAMES)],
              fontsize=7, facecolor="#1a1d27", labelcolor="white")
    ax.set(ylabel="x10^6 / mL", title="Cell Density", xlabel="Day")

    ax = axes[0, 1]
    ax.plot(ts, [r["Hb_pg"] for r in h], color="#ef9a9a", lw=2)
    _shade(ax, b)
    ax.axhline(hb_target, color="#ff5252", lw=1, ls=":", label=f"Target {hb_target:.0f} pg")
    ax.axhline(qc_hb_min, color="#ffa726", lw=1, ls=":", label=f"QC min {qc_hb_min:.0f} pg")
    ax.legend(fontsize=7, facecolor="#1a1d27", labelcolor="white")
    ax.set(ylabel="pg / cell", title="Hemoglobinization", xlabel="Day")

    ax = axes[1, 0]
    ax.plot(ts, [r["GPA"] * 100 for r in h],       color="#a5d6a7", lw=2, label="GPA / CD235a (%)")
    ax.plot(ts, [r["viability"] * 100 for r in h], color="#eeeeee", lw=1.5, ls="--", label="Viability (%)")
    ax.axhline(qc_viab_min, color="#ffa726", lw=1, ls=":", label=f"QC min {qc_viab_min:.0f}%")
    _shade(ax, b)
    ax.legend(fontsize=7, facecolor="#1a1d27", labelcolor="white")
    ax.set(ylabel="%", title="Differentiation Markers & Viability", xlabel="Day")

    ax = axes[1, 1]
    ax.plot(ts, [r["enuc_rate"] * 100 for r in h], color="#ce93d8", lw=2)
    _shade(ax, b)
    ax.axhline(twin.p["enucleation"]["max_enuc"] * 100, color="#9c27b0", lw=1, ls=":",
               label=f"HDACi ceiling {twin.p['enucleation']['max_enuc'] * 100:.0f}%")
    ax.axhline(qc_enuc_min, color="#ffa726", lw=1, ls=":", label=f"QC min {qc_enuc_min:.0f}%")
    ax.legend(fontsize=7, facecolor="#1a1d27", labelcolor="white")
    ax.set(ylabel="%", title="Enucleation Rate (HDACi)", xlabel="Day")

    plt.tight_layout()
    return fig


def comparison_fig(df):
    metrics    = ["Days", "Enuc (%)", "Viability (%)", "Hb (pg)"]
    bar_colors = ["#28a745" if q == "PASS" else "#dc3545" for q in df["QC"]]
    labels     = [s[:28] + "…" if len(s) > 28 else s for s in df.index]

    fig, axes = plt.subplots(1, 4, figsize=(18, max(3, len(df) * 0.55)))
    fig.patch.set_facecolor("#f8f9fa")
    for ax, metric in zip(axes, metrics):
        ax.set_facecolor("#ffffff")
        vals = df[metric].tolist()
        bars = ax.barh(labels, vals, color=bar_colors, edgecolor="none", height=0.6)
        ax.bar_label(bars, fmt="%.0f" if metric == "Days" else "%.1f",
                     padding=3, fontsize=8, color="#333")
        ax.set_xlabel(metric, fontsize=9)
        ax.tick_params(labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(0, max(vals) * 1.18)
        ax.invert_yaxis()
    fig.tight_layout(pad=1.5)
    return fig


def sensitivity_fig(ceilings, enuc_out, qc_pass, current, qc_min):
    fig, ax = plt.subplots(figsize=(10, 4))
    fig.patch.set_facecolor("#0e1117")
    _dark(ax)
    ax.plot(ceilings, enuc_out, color="#ce93d8", lw=1.5, zorder=2)
    ax.scatter(ceilings, enuc_out,
               c=["#28a745" if p else "#dc3545" for p in qc_pass], zorder=3, s=60)
    ax.axhline(qc_min,  color="#ffa726", lw=1, ls=":", label=f"QC min {qc_min:.0f}%")
    ax.axvline(current, color="#4fc3f7", lw=1, ls="--", label=f"Current ceiling {current:.0f}%")
    ax.legend(fontsize=8, facecolor="#1a1d27", labelcolor="white")
    ax.set(xlabel="HDACi Ceiling (%)", ylabel="Final Enucleation (%)",
           title="Sensitivity: HDACi Ceiling → Enucleation Rate")
    plt.tight_layout()
    return fig
