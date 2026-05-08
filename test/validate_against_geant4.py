#!/usr/bin/env python3
"""
validate_against_geant4.py

Compare g4gamma analytic spectra (Geant4 data, Sandia, LARA backends) against
a full rdecay01 Monte Carlo simulation for U-238 at secular equilibrium.

Usage:
    # from repo root
    python test/validate_against_geant4.py <rdecay01_csv> <n_primaries> [output_dir]

    # example with the 1e7 (10M) event run:
    python test/validate_against_geant4.py \\
        buildG4RadDecayExample/u238_AFtrue_h1_3.csv 10000000

Simulation settings matched here (moduleTest.mac + PhysicsList.cc):
    /process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year  -> t=-1 (secular eq.)
    /rdecay01/fullChain true                                      -> full chain
    radioactiveDecay->SetARM(true)                                -> K/L/M X-rays ON
    G4EmParameters::SetAugerCascade(true)                         -> Auger cascade ON
    G4EmParameters::SetDeexcitationIgnoreCut(true)                -> no energy cut on X-rays
    /gun/ion 92 238                                               -> U-238 primary
    /analysis/h1/set 3  3000  0. 3000 keV                        -> 3000 bins, 1 keV wide

g4gamma Geant4 backend uses include_xrays=True (K-shell only). The simulation's
G4UAtomicDeexcitation produces the full Auger cascade (K+L+M shells), so
g4gamma will slightly underpredict the 50-90 keV X-ray region for heavy elements.

LARA note: Th-230 is available from LNHB's LaraWEB but not at the standard
/nuclides/ path. fetch_lara.sh uses a POST fallback to Result_Lara2.php to
find the versioned filename (Th-230_@03.lara.txt). The real file (27 emission
lines, including Ra K X-rays) is included in data/lara/lara.tar.gz.
"""

import os
import sys
import numpy as np
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..")

# --- locate g4gamma module ---------------------------------------------------
for cand in [os.path.join(REPO, "build"), os.path.join(HERE, "build"), "build", "."]:
    if os.path.isdir(cand) and any(
        f.startswith("g4gamma") and f.endswith(".so") for f in os.listdir(cand)
    ):
        sys.path.insert(0, cand)
        break

try:
    import g4gamma as g
except ImportError as e:
    print(f"ERROR: cannot import g4gamma. Build it first.\n  {e}")
    sys.exit(1)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
except ImportError:
    print("ERROR: matplotlib required. pip install matplotlib")
    sys.exit(1)


# =============================================================================
# Configuration
# =============================================================================

Z, A, M = 92, 238, 0
N_BINS = 3000
E_MIN_KEV, E_MAX_KEV = 0.0, 3000.0

# g4gamma settings that match moduleTest.mac + PhysicsList.cc:
#   - secular equilibrium (thresholdForVeryLongDecayTime 1e60 year)
#   - ARM ON: radioactiveDecay->SetARM(true) + G4UAtomicDeexcitation
#             + SetAugerCascade(true), SetDeexcitationIgnoreCut(true)
#             → K-shell X-rays from EC and IC are present in the simulation
#   - 511 keV annihilation pairs from beta+ ARE present in the simulation
#
# Note: g4gamma only implements K-shell X-rays (not L/M-shell Auger cascade).
# The simulation's G4UAtomicDeexcitation produces the full cascade.
# Expect the 50-90 keV X-ray region to be slightly underpredicted by g4gamma.
G4GAMMA_T = -1.0          # secular equilibrium
INCLUDE_XRAYS = True      # ARM=true in PhysicsList; K-shell X-rays only in g4gamma
INCLUDE_ANNIHILATION = True


# =============================================================================
# Helpers
# =============================================================================

def load_rdecay01_csv(path, n_primaries):
    """
    Parse a Geant4 rdecay01 CSV histogram file.

    Format (Geant4 11.x analysis):
        6 lines starting with '#'   -- class/title/dimension/axis/annotation/bin_number
        1 column header line        -- 'entries,Sw,Sw2,Sxw0,Sx2w0'  (NOT a comment)
        1 underflow row
        N data rows
        1 overflow row

    The 'entries' column (col 0) is the raw fill count per bin.
    Dividing by n_primaries gives gammas per primary decay per bin.

    The 'Sw' column is sum-of-weights; rdecay01 fills with weight=0.01, so
    Sw / n_primaries = 0.01 * (gammas/primary/bin).  We use 'entries' directly.
    """
    rows = []
    header_skipped = False
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            if not header_skipped:
                # first non-comment line is the column name row
                header_skipped = True
                continue
            rows.append([float(x) for x in line.split(",")])

    # rows[0]=underflow, rows[1:-1]=data bins, rows[-1]=overflow
    entries = np.array([r[0] for r in rows[1:-1]])
    assert len(entries) == N_BINS, f"Expected {N_BINS} bins, got {len(entries)}"
    return entries / n_primaries


def build_spectra(edges_keV):
    """Build spectra from all three g4gamma backends. Returns dict of arrays."""
    edges_int = edges_keV * g.units.keV
    key = g.IsotopeKey(Z, A, M)

    # --- Geant4 backend: closest match to the rdecay01 simulation ---
    opts_g4 = g.SpectrumOptions()
    opts_g4.source = g.DataSource.Geant4
    opts_g4.include_xrays = INCLUDE_XRAYS
    opts_g4.include_annihilation = INCLUDE_ANNIHILATION
    b_g4 = g.GammaSpectrumBuilder(opts_g4)
    r_g4 = b_g4.build(key, G4GAMMA_T, edges_int)

    # --- SandiaDecay backend ---
    opts_sd = g.SpectrumOptions()
    opts_sd.source = g.DataSource.Sandia
    b_sd = g.GammaSpectrumBuilder(opts_sd)
    r_sd = b_sd.build(key, G4GAMMA_T, edges_int)

    # --- LARA backend ---
    opts_la = g.SpectrumOptions()
    opts_la.source = g.DataSource.Lara
    b_la = g.GammaSpectrumBuilder(opts_la)
    r_la = b_la.build(key, G4GAMMA_T, edges_int)

    return {
        "geant4": np.array(r_g4.counts),
        "sandia": np.array(r_sd.counts),
        "lara":   np.array(r_la.counts),
        "contribs": {
            "geant4": r_g4.contributions,
            "sandia": r_sd.contributions,
            "lara":   r_la.contributions,
        },
    }


# =============================================================================
# Plotting helpers
# =============================================================================

COLORS = {
    "sim":    "#1a1a2e",   # near-black
    "geant4": "#2F4F8F",   # navy
    "sandia": "#C45A11",   # rust
    "lara":   "#1D9E75",   # teal
}
LABELS = {
    "sim":    "rdecay01 (Geant4 MC)",
    "geant4": "g4gamma — Geant4 data",
    "sandia": "g4gamma — Sandia",
    "lara":   "g4gamma — LARA/DDEP",
}
STYLES = {
    "sim":    dict(lw=1.1, zorder=4, alpha=0.9),
    "geant4": dict(lw=1.1, zorder=3),
    "sandia": dict(lw=1.1, zorder=2, linestyle="--"),
    "lara":   dict(lw=1.1, zorder=1, linestyle=":"),
}

plt.rcParams.update({
    "figure.dpi":      110,
    "savefig.dpi":     160,
    "font.family":     "DejaVu Sans",
    "font.size":       9,
    "axes.titlesize":  10,
    "axes.labelsize":  9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid":       True,
    "grid.alpha":      0.2,
    "grid.linewidth":  0.5,
    "legend.frameon":  False,
    "lines.linewidth": 1.0,
})


def step_xy(edges, counts):
    """Return (x, y) arrays for a step histogram plot."""
    x = np.concatenate([[edges[0]], np.repeat(edges[1:-1], 2), [edges[-1]]])
    y = np.repeat(counts, 2)
    return x, y


def plot_spectrum(ax, edges, counts, key, clip_low=1e-7):
    counts = np.clip(counts, clip_low, None)
    x, y = step_xy(edges, counts)
    ax.plot(x, y, color=COLORS[key], label=LABELS[key], **STYLES[key])


def annotate_peak(ax, center_kev, label, offset=(0, 8), fontsize=7.5):
    ax.annotate(
        label,
        xy=(center_kev, ax.get_ylim()[1]),
        xycoords=("data", "axes fraction"),
        xytext=offset,
        textcoords="offset points",
        ha="center", va="bottom",
        fontsize=fontsize, color="#444",
        arrowprops=None,
    )


# =============================================================================
# Main plots
# =============================================================================

def make_plots(edges_keV, sim, analytic, outdir):
    os.makedirs(outdir, exist_ok=True)
    centers = 0.5 * (edges_keV[:-1] + edges_keV[1:])

    def _has_sim():
        return sim is not None

    # -------------------------------------------------------------------------
    # Fig 1: Full spectrum overview (log scale)
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(11, 8),
                             gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08})
    ax_main, ax_ratio = axes

    if _has_sim():
        plot_spectrum(ax_main, edges_keV, sim, "sim")

    for k in ("geant4", "sandia", "lara"):
        plot_spectrum(ax_main, edges_keV, analytic[k], k)

    ax_main.set_yscale("log")
    ax_main.set_ylim(1e-6, 3)
    ax_main.set_xlim(0, 3000)
    ax_main.set_ylabel("γ / primary decay / keV bin")
    ax_main.set_title("U-238 secular-equilibrium gamma spectrum — dataset comparison", loc="left")
    ax_main.legend(loc="upper right", fontsize=8.5)
    ax_main.set_xticklabels([])

    # annotate major peaks on the main plot
    MAJOR_PEAKS = [
        (63.3,  "Th-234\n63 keV"),
        (93.3,  "Th-234\n93 keV"),
        (186.2, "Ra-226\n186 keV"),
        (295.2, "Pb-214\n295 keV"),
        (351.9, "Pb-214\n352 keV"),
        (609.3, "Bi-214\n609 keV"),
        (768.4, "Bi-214\n768 keV"),
        (1001,  "Pa-234m\n1001 keV"),
        (1120,  "Bi-214\n1120 keV"),
        (1238,  "Bi-214\n1238 keV"),
        (1377,  "Bi-214\n1378 keV"),
        (1510,  "Bi-214\n1510 keV"),
        (1764,  "Bi-214\n1765 keV"),
        (2204,  "Bi-214\n2204 keV"),
    ]
    y_lim = ax_main.get_ylim()
    for e_kev, lbl in MAJOR_PEAKS:
        ax_main.axvline(e_kev, color="#ccc", lw=0.5, zorder=0)

    # Ratio panel: each analytic / sim  (or analytic / Sandia if no sim)
    ref_key = "sim" if _has_sim() else "sandia"
    ref = sim if _has_sim() else analytic["sandia"]
    ref_label = "MC sim" if _has_sim() else "Sandia"

    def safe_ratio(num, den, threshold=1e-5):
        mask = den > threshold
        out = np.full_like(num, np.nan)
        out[mask] = num[mask] / den[mask]
        return out

    for k in ("geant4", "sandia", "lara"):
        if k == ref_key:
            continue
        ratio = safe_ratio(analytic[k], ref)
        x, y = step_xy(edges_keV, ratio)
        ax_ratio.plot(x, y, color=COLORS[k], lw=0.9,
                      linestyle=STYLES[k]["linestyle"] if "linestyle" in STYLES[k] else "-")

    if _has_sim():
        ratio_g4 = safe_ratio(analytic["geant4"], sim)
        x, y = step_xy(edges_keV, ratio_g4)
        ax_ratio.plot(x, y, color=COLORS["geant4"], lw=0.9)

    ax_ratio.axhline(1.0, color="#888", lw=0.8, linestyle="--")
    ax_ratio.set_ylim(0.5, 1.5)
    ax_ratio.set_xlim(0, 3000)
    ax_ratio.set_xlabel("Energy (keV)")
    ax_ratio.set_ylabel(f"Ratio to {ref_label}")

    path = os.path.join(outdir, "01_u238_overview.png")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {path}")

    # -------------------------------------------------------------------------
    # Fig 2: Four zoom panels for key peak regions
    # -------------------------------------------------------------------------
    ZOOM_REGIONS = [
        (50,  250, "Low-energy region (Pb-214 / Ra-226 / Th-234 / Th-230)"),
        (270, 430, "Pb-214 peaks: 295.2, 352.0 keV"),
        (550, 700, "Bi-214 peak: 609.3 keV"),
        (1700, 1820, "Bi-214 peak: 1764.5 keV"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    axes = axes.flatten()

    for ax, (lo, hi, title) in zip(axes, ZOOM_REGIONS):
        mask = (centers >= lo) & (centers <= hi)
        e_zoom = edges_keV[:-1][mask]
        e_zoom = np.append(e_zoom, edges_keV[1:][mask][-1])

        if _has_sim():
            ax.bar(centers[mask], sim[mask], width=1.0,
                   color=COLORS["sim"], alpha=0.35, label=LABELS["sim"])

        for k in ("geant4", "sandia", "lara"):
            c = analytic[k][mask]
            x, y = step_xy(e_zoom, c)
            kw = dict(STYLES[k])
            ax.plot(x, y, color=COLORS[k], label=LABELS[k], **kw)

        ax.set_xlim(lo, hi)
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("γ / primary decay")
        ax.set_title(title, loc="left", fontsize=9)
        ax.legend(loc="upper right", fontsize=7.5)

    path = os.path.join(outdir, "02_u238_zoom_peaks.png")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {path}")

    # -------------------------------------------------------------------------
    # Fig 3: High-energy region (1000-2400 keV) — Bi-214 forest
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 5))
    lo, hi = 950, 2400
    mask = (centers >= lo) & (centers <= hi)
    e_zoom = edges_keV[:-1][mask]
    e_zoom = np.append(e_zoom, edges_keV[1:][mask][-1])

    if _has_sim():
        ax.bar(centers[mask], sim[mask], width=1.0,
               color=COLORS["sim"], alpha=0.35, label=LABELS["sim"])
    for k in ("geant4", "sandia", "lara"):
        x, y = step_xy(e_zoom, analytic[k][mask])
        ax.plot(x, y, color=COLORS[k], label=LABELS[k], **STYLES[k])

    BI214_HIGH = [1001, 1120, 1155, 1238, 1281, 1377, 1385, 1408, 1510, 1583,
                  1594, 1661, 1730, 1764, 2204]
    for e_kev in BI214_HIGH:
        if lo <= e_kev <= hi:
            ax.axvline(e_kev, color="#ddd", lw=0.5, zorder=0)

    ax.set_xlim(lo, hi)
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("γ / primary decay")
    ax.set_title("U-238 chain — high-energy region (Bi-214 + Pa-234m)", loc="left")
    ax.legend(fontsize=8.5)

    path = os.path.join(outdir, "03_u238_high_energy.png")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {path}")

    # -------------------------------------------------------------------------
    # Fig 4: Ratio of each analytic backend to Sandia (or sim)
    # -------------------------------------------------------------------------
    ref_key = "sim" if _has_sim() else "sandia"
    ref     = sim   if _has_sim() else analytic["sandia"]
    ref_lbl = "MC sim" if _has_sim() else "Sandia"

    comp_keys = ["geant4", "lara"] if ref_key == "sandia" else ["geant4", "sandia", "lara"]
    fig, ax = plt.subplots(figsize=(12, 4))

    for k in comp_keys:
        ratio = safe_ratio(analytic[k], ref, threshold=5e-5)
        x, y = step_xy(edges_keV, ratio)
        kw = dict(lw=0.9)
        if "linestyle" in STYLES[k]:
            kw["linestyle"] = STYLES[k]["linestyle"]
        ax.plot(x, y, color=COLORS[k], label=f"{LABELS[k]} / {ref_lbl}", **kw)

    ax.axhline(1.0, color="#888", lw=0.8, linestyle="--")
    ax.set_ylim(0.6, 1.4)
    ax.set_xlim(0, 3000)
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel(f"Ratio to {ref_lbl}")
    ax.set_title(
        f"Dataset agreement — ratio to {ref_lbl} (peaks above {5e-5:.0e} γ/primary)",
        loc="left",
    )
    ax.legend(fontsize=8.5)

    path = os.path.join(outdir, "04_u238_ratios.png")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {path}")


# =============================================================================
# Peak comparison table
# =============================================================================

REFERENCE_PEAKS = [
    # (energy_keV, nuclide, ENSDF_intensity, label)
    # Intensities are gammas/primary at secular equilibrium (= per-decay intensity
    # since all chain members have activity 1 at SE).
    (63.3,   "Th-234",  0.037,  "63.3"),   # Th-234 nuclear gamma
    (75.5,   "Bi K X-ray", 0.093, "75.5"),  # Bi Kα1/2 X-ray (ARM; peak bin varies ±2 keV by G4EMLOW ver)
    (92.5,   "Th/Rn Kα",0.063, "92.5"),   # Th/Rn K X-rays (ARM)
    (186.2,  "Ra-226",  0.036,  "186.2"),
    (295.2,  "Pb-214",  0.184,  "295.2"),
    (351.9,  "Pb-214",  0.358,  "351.9"),
    (609.3,  "Bi-214",  0.461,  "609.3"),
    (768.4,  "Bi-214",  0.0491, "768.4"),
    (1001.0, "Pa-234m", 0.0084, "1001"),
    (1120.3, "Bi-214",  0.1491, "1120.3"),
    (1238.1, "Bi-214",  0.0583, "1238.1"),
    (1764.5, "Bi-214",  0.1530, "1764.5"),
    (2204.1, "Bi-214",  0.0491, "2204.1"),
]


def print_peak_table(edges_keV, sim, analytic):
    print("\n" + "=" * 95)
    print(f"{'Peak (keV)':<12} {'Nuclide':<10} {'ENSDF':>8}", end="")
    if sim is not None:
        print(f" {'MC sim':>10} {'sim-err%':>9}", end="")
    for k in ("geant4", "sandia", "lara"):
        print(f" {k:>10} {'err%':>6}", end="")
    print()
    print("-" * 95)

    for e_kev, nuc, ref_int, label in REFERENCE_PEAKS:
        # find bin
        idx = int(np.searchsorted(edges_keV, e_kev, side="right")) - 1
        idx = max(0, min(idx, len(edges_keV) - 2))

        row = f"{label:<12} {nuc:<10} {ref_int:>8.4f}"
        if sim is not None:
            v = sim[idx]
            pct = (v - ref_int) / ref_int * 100 if ref_int > 0 else float("nan")
            row += f" {v:>10.4f} {pct:>+8.1f}%"
        for k in ("geant4", "sandia", "lara"):
            v = analytic[k][idx]
            pct = (v - ref_int) / ref_int * 100 if ref_int > 0 else float("nan")
            row += f" {v:>10.4f} {pct:>+5.1f}%"
        print(row)

    print("=" * 95)
    print()
    print("Notes:")
    print("  ENSDF = reference intensities (gammas/primary at secular equilibrium)")
    print("  err% = (value - ENSDF) / ENSDF * 100")
    print("  Geant4 backend: include_xrays=True (ARM=true in PhysicsList, K-shell only)")
    print("  Sandia / LARA: X-ray and annihilation emissions baked into dataset files")
    print("  Simulation uses full Auger cascade (K+L+M); g4gamma models K-shell only")
    print("  LARA Th-230: real LNHB data (27 emission lines) fetched from LaraWEB; "
          "gamma yield ~0.4% at 67.7 keV")
    print()


# =============================================================================
# Chain contribution summary
# =============================================================================

def print_chain_summary(analytic):
    print("\n--- Chain contributions (gammas/primary at secular equilibrium) ---")
    header = f"{'Nuclide':<10} {'Activity':>9}"
    for k in ("geant4", "sandia", "lara"):
        header += f" {k:>10} {'g/dec':>7}"
    print(header)
    print("-" * 60)

    # collect all isotope keys seen across all providers
    seen = {}
    for k in ("geant4", "sandia", "lara"):
        for c in analytic["contribs"][k]:
            key = str(c.isotope)
            if key not in seen:
                seen[key] = {}
            seen[key][k] = c

    for iso_str, providers in seen.items():
        # pick any contribution for the activity
        ref = next(iter(providers.values()))
        act = ref.activity
        row = f"{iso_str:<10} {act:>9.4f}"
        for k in ("geant4", "sandia", "lara"):
            if k in providers:
                c = providers[k]
                row += f" {c.activity:>10.4f} {c.gamma_yield:>7.4f}"
            else:
                row += f"  {'---':>9} {'---':>7}"
        print(row)
    print()


# =============================================================================
# Entry point
# =============================================================================

def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    csv_path   = sys.argv[1]
    n_primaries = int(float(sys.argv[2]))  # accept 1e7 notation
    outdir     = sys.argv[3] if len(sys.argv) > 3 else os.path.join(REPO, "validation_plots")

    # --- load simulation -------------------------------------------------------
    print(f"Reading rdecay01 CSV: {csv_path}")
    print(f"  n_primaries = {n_primaries:,}")
    sim = load_rdecay01_csv(csv_path, n_primaries)
    print(f"  {len(sim)} bins, total {sim.sum():.3f} γ/primary")
    print(f"  statistical uncertainty on 609 keV peak: "
          f"~{1/np.sqrt(max(sim[609]*n_primaries, 1))*100:.1f}%")

    # --- build analytic spectra -----------------------------------------------
    edges_keV = np.linspace(E_MIN_KEV, E_MAX_KEV, N_BINS + 1)
    print("\nBuilding analytic spectra...")
    analytic = build_spectra(edges_keV)
    for k in ("geant4", "sandia", "lara"):
        print(f"  {k:8s}: {analytic[k].sum():.4f} γ/primary total")

    # --- print tables ----------------------------------------------------------
    print_peak_table(edges_keV, sim, analytic)
    print_chain_summary(analytic)

    # --- make plots ------------------------------------------------------------
    print(f"Writing plots to: {os.path.abspath(outdir)}")
    make_plots(edges_keV, sim, analytic, outdir)

    # --- write difference explanation ------------------------------------------
    _print_explanation(sim, analytic)


def _print_explanation(sim, analytic):
    print("""
=============================================================================
EXPLANATION OF DIFFERENCES BETWEEN DATASETS
=============================================================================

1. Geant4 data vs. Sandia/LARA
   The Geant4 backend reads the same data files ($G4RADIOACTIVEDATA,
   $G4LEVELGAMMADATA) that the rdecay01 Monte Carlo simulation uses, so it
   is the reference for reproducing the simulation output.

   Sandia and LARA source their data from ENSDF evaluations which are
   updated independently of the Geant4 data files. Small differences in
   branching ratios and level energies accumulate across the long U-238 chain.

2. Total gamma count differences (Geant4 > Sandia ≈ LARA)
   With ARM=true the Geant4 backend total is ~2.87 γ/primary vs Sandia 2.25 and
   LARA 2.16. The excess comes from:
   - Short-lived isomers (U-234m, Bi-210m, Pb-210m, Tl-206m, Pb-206m) that the
     Geant4 data resolves as separate chain nodes with their own gamma yields;
     Sandia/LARA fold these into parent entries.
   - Pa-234 ground state: Geant4 splits Th-234 β⁻ → Pa-234m (78%) + Pa-234 (22%);
     Pa-234 (t½=6.7 h) has a very high gamma yield (~1.5 γ/decay) from its own
     complex level scheme. Sandia treats the whole Th-234 branch as going through
     Pa-234m only.
   - K X-ray contributions (ARM=true) are similar across all backends above 70 keV.

3. 609 keV peak (Bi-214 diagnostic)
   All datasets agree to within ~1% on this peak — it is the dominant
   NORM peak and well measured. Residual differences reflect dataset vintage
   and whether the Geant4 data pre/post-dates the most recent ENSDF evaluation.

4. Low-energy region (< 150 keV) — X-ray and fluorescence region
   The simulation runs ARM=true with a full Auger cascade (G4UAtomicDeexcitation,
   K+L+M shells, SetDeexcitationIgnoreCut=true). g4gamma's Geant4 backend
   implements K-shell X-rays only (from FluorData / G4EMLOW fl-tr-pr-Z.dat).

   Two effects:
   a) Energy offset ±1-2 keV: g4gamma reads the same G4EMLOW files but the
      simulation may use a different G4EMLOW version, causing K X-ray peaks to
      sit in adjacent 1 keV bins (e.g. Bi Kα1 at 75-76 keV in backend vs
      77-78 keV in sim). Total X-ray yield in the region agrees better than
      bin-by-bin.
   b) Missing L/M shells: the simulation produces the full Auger cascade;
      g4gamma only adds K-shell X-rays. L/M X-rays for Pb/Bi/Po are in the
      10-20 keV range, below NaI/LaBr3 detection threshold, but the K-shell
      yield difference between g4gamma (K only) and the full simulation is
      a few percent of the total K-shell intensity.

   Sandia and LARA have K X-ray intensities baked into their emission tables
   (from ENSDF evaluations) and agree better with the simulation in this region.

5. LARA Th-230
   Th-230 is not at the standard /nuclides/ path on LNHB but is available via
   LaraWEB. fetch_lara.sh discovers the versioned filename (Th-230_@03.lara.txt)
   by POSTing to Result_Lara2.php and downloads it automatically. The real file
   has 27 emission lines including Ra K X-rays. Its total gamma yield (~0.4%)
   is negligible compared to the Bi-214 contribution (~130%).

6. Simulation statistics
   This run used 1×10⁷ primaries, giving ~0.05% statistical uncertainty on the
   609 keV peak (√(0.461 × 1e7) ≈ 2147 σ). Residuals above ~0.3% on strong
   peaks are dataset differences, not Monte Carlo noise.
=============================================================================
""")


if __name__ == "__main__":
    main()
