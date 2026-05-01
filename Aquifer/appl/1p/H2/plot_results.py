#!/usr/bin/env python3
"""
Plot DuMuX UHS simulation results from the JSON output file.

Usage:
    python3 plot_results.py                    # looks for H2.json in cwd
    python3 plot_results.py path/to/file.json  # explicit path

Produces 3 figures:
  1. Reservoir state    — pressure (all + gas-bearing) and mass balance error
  2. Well performance   — recovery factor, H2 purity, water cut
  3. H2 budget          — injection/production rates and cumulative H2 balance
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def load_json(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def safe_get(data, *keys, default=None):
    """Nested dict lookup that returns default if any key is missing."""
    node = data
    for k in keys:
        if isinstance(node, dict) and k in node:
            node = node[k]
        else:
            return default
    return node


def to_days(seconds):
    return np.asarray(seconds) / 86400.0


def add_cycle_shading(ax, t_days, cycle_numbers):
    """Shade background by cycle number for visual reference."""
    if cycle_numbers is None or len(cycle_numbers) == 0:
        return
    cn = np.asarray(cycle_numbers)
    unique_cycles = sorted(set(cn))
    colors = ["#f0f0f0", "#e0e8f0"]
    for i, c in enumerate(unique_cycles):
        mask = cn == c
        idx = np.where(mask)[0]
        if len(idx) < 2:
            continue
        ax.axvspan(t_days[idx[0]], t_days[idx[-1]],
                   alpha=0.3, color=colors[i % 2], zorder=0)


# ── Figure 1: Reservoir State ──────────────────────────────────────────
def plot_reservoir_state(data, t_days, cycle_numbers, outdir):
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

    # (a) Pressure
    ax = axes[0]
    add_cycle_shading(ax, t_days, cycle_numbers)
    p_all = safe_get(data, "averageReservoirPressure")
    p_gas = safe_get(data, "averageReservoirPressureGasBearing")
    if p_all:
        ax.plot(t_days, np.array(p_all) / 1e6, label="All cells", lw=1.5)
    if p_gas:
        ax.plot(t_days, np.array(p_gas) / 1e6, label="Gas-bearing cells", lw=1.5, ls="--")
    ax.set_ylabel("Pressure [MPa]")
    ax.set_title("Average Reservoir Pressure")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (b) Mass balance error
    ax = axes[1]
    add_cycle_shading(ax, t_days, cycle_numbers)
    mbe_abs = safe_get(data, "massBalanceError", "absolute")
    mbe_rel = safe_get(data, "massBalanceError", "relative")
    if mbe_abs:
        ax.plot(t_days, mbe_abs, label="Absolute", lw=1.2, color="C3")
    ax.set_ylabel("Mass balance error [mol]")
    ax.set_xlabel("Time [days]")
    ax.set_title("Mass Balance Error")
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.3)
    if mbe_rel:
        ax2 = ax.twinx()
        ax2.plot(t_days, np.array(mbe_rel) * 100, label="Relative", lw=1.2,
                 color="C4", ls="--")
        ax2.set_ylabel("Relative error [%]")
        ax2.legend(fontsize=9, loc="upper right")

    fig.suptitle("Reservoir State", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(outdir / "fig1_reservoir_state.png", dpi=150)
    plt.close(fig)
    print(f"  Saved fig1_reservoir_state.png")


# ── Figure 2: Well Performance ─────────────────────────────────────────
def plot_well_performance(data, t_days, cycle_numbers, outdir):
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)

    # (a) Recovery Factor
    ax = axes[0]
    add_cycle_shading(ax, t_days, cycle_numbers)
    rf_cycle = safe_get(data, "recoveryFactor", "perCycle")
    rf_cum = safe_get(data, "recoveryFactor", "cumulative")
    if rf_cycle:
        ax.plot(t_days, rf_cycle, label="Per cycle", lw=1.5, color="C0")
    if rf_cum:
        ax.plot(t_days, rf_cum, label="Cumulative", lw=1.5, color="C1", ls="--")
    ax.set_ylabel("Recovery Factor [-]")
    ax.set_title("H₂ Recovery Factor")
    ax.set_ylim(bottom=-0.05)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (b) H2 Purity
    ax = axes[1]
    add_cycle_shading(ax, t_days, cycle_numbers)
    purity = safe_get(data, "wellPurity", "H2MoleFraction")
    if purity:
        ax.plot(t_days, np.array(purity) * 100, lw=1.5, color="C2")
    ax.set_ylabel("H₂ mole fraction [%]")
    ax.set_title("H₂ Purity at Well")
    ax.set_ylim(-5, 105)
    ax.grid(True, alpha=0.3)

    # (c) Water Cut
    ax = axes[2]
    add_cycle_shading(ax, t_days, cycle_numbers)
    wcut = safe_get(data, "wellWaterCut")
    if wcut:
        ax.plot(t_days, np.array(wcut) * 100, lw=1.5, color="C5")
    ax.set_ylabel("Water cut [%]")
    ax.set_xlabel("Time [days]")
    ax.set_title("Water Cut at Well")
    ax.set_ylim(-5, 105)
    ax.grid(True, alpha=0.3)

    fig.suptitle("Well Performance", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(outdir / "fig2_well_performance.png", dpi=150)
    plt.close(fig)
    print(f"  Saved fig2_well_performance.png")


# ── Figure 3: H2 Budget ────────────────────────────────────────────────
def plot_h2_budget(data, t_days, cycle_numbers, outdir):
    fig, axes = plt.subplots(2, 2, figsize=(14, 8))

    # (a) Instantaneous injection/production rates (H2 only)
    ax = axes[0, 0]
    add_cycle_shading(ax, t_days, cycle_numbers)
    inj_h2 = safe_get(data, "InjectionValues", "H2")
    prod_h2 = safe_get(data, "ProductionValues", "H2")
    if inj_h2:
        ax.plot(t_days, inj_h2, label="H₂ injection", lw=1.2, color="C0")
    if prod_h2:
        ax.plot(t_days, prod_h2, label="H₂ production", lw=1.2, color="C3")
    ax.set_ylabel("Rate [mol/s]")
    ax.set_title("H₂ Injection & Production Rate")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (b) Cumulative H2 injected vs produced
    ax = axes[0, 1]
    add_cycle_shading(ax, t_days, cycle_numbers)
    cum_inj = safe_get(data, "recoveryFactor", "cumulativeH2Injected")
    cum_prod = safe_get(data, "recoveryFactor", "cumulativeH2Produced")
    if cum_inj:
        ax.plot(t_days, cum_inj, label="Cumulative injected", lw=1.5, color="C0")
    if cum_prod:
        ax.plot(t_days, cum_prod, label="Cumulative produced", lw=1.5, color="C3")
    ax.set_ylabel("Cumulative H₂ [mol]")
    ax.set_title("Cumulative H₂ Balance")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (c) Component inventory
    ax = axes[1, 0]
    add_cycle_shading(ax, t_days, cycle_numbers)
    inv = safe_get(data, "inventory")
    if inv:
        for comp, vals in inv.items():
            if comp == "SimpleH2O":
                continue  # skip water — dominates the scale
            ax.plot(t_days, vals, label=comp, lw=1.2)
    ax.set_ylabel("Inventory [mol]")
    ax.set_xlabel("Time [days]")
    ax.set_title("Gas Component Inventory (excl. H₂O)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (d) Per-component mass balance error
    ax = axes[1, 1]
    add_cycle_shading(ax, t_days, cycle_numbers)
    mbe = safe_get(data, "materialBalanceError")
    if mbe:
        for comp, vals in mbe.items():
            if comp == "SimpleH2O":
                continue
            ax.plot(t_days, vals, label=comp, lw=1.0)
    ax.set_ylabel("Error [mol]")
    ax.set_xlabel("Time [days]")
    ax.set_title("Per-Component Mass Balance Error")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.suptitle("H₂ Budget & Inventory", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(outdir / "fig3_h2_budget.png", dpi=150)
    plt.close(fig)
    print(f"  Saved fig3_h2_budget.png")


# ── Main ───────────────────────────────────────────────────────────────
def main():
    # Find JSON file
    if len(sys.argv) > 1:
        json_path = Path(sys.argv[1])
    else:
        # Default: look for H2.json in the matching build-cmake H2 folder.
        json_path = (
            Path(__file__).resolve().parents[3]
            / "build-cmake"
            / "appl"
            / "1p"
            / "H2"
            / "H2.json"
        )

    print(f"Loading: {json_path}")
    data = load_json(json_path)

    # Time vector
    t_sec = data.get("time", [])
    if not t_sec:
        print("ERROR: No 'time' array found in JSON.")
        sys.exit(1)

    t_days = to_days(t_sec)
    n = len(t_days)
    print(f"  Time steps: {n}")
    print(f"  Time range: {t_days[0]:.1f} – {t_days[-1]:.1f} days")

    cycle_numbers = safe_get(data, "recoveryFactor", "cycleNumber")

    # Output directory
    outdir = json_path.parent
    print(f"  Output dir: {outdir}")
    print()

    # Quick summary
    p_all = safe_get(data, "averageReservoirPressure")
    if p_all:
        print(f"  Pressure range: {min(p_all)/1e6:.2f} – {max(p_all)/1e6:.2f} MPa")

    rf_cum = safe_get(data, "recoveryFactor", "cumulative")
    if rf_cum:
        print(f"  Final cumulative RF: {rf_cum[-1]:.4f}")

    purity = safe_get(data, "wellPurity", "H2MoleFraction")
    if purity:
        nonzero = [p for p in purity if p > 0]
        if nonzero:
            print(f"  H2 purity range: {min(nonzero)*100:.1f}% – {max(nonzero)*100:.1f}%")

    wcut = safe_get(data, "wellWaterCut")
    if wcut:
        nonzero = [w for w in wcut if w > 0]
        if nonzero:
            print(f"  Water cut range: {min(nonzero)*100:.1f}% – {max(nonzero)*100:.1f}%")

    print()

    # Generate figures
    plot_reservoir_state(data, t_days, cycle_numbers, outdir)
    plot_well_performance(data, t_days, cycle_numbers, outdir)
    plot_h2_budget(data, t_days, cycle_numbers, outdir)

    print(f"\nDone. 3 figures saved to {outdir}/")


if __name__ == "__main__":
    main()
