#!/usr/bin/env python3
"""
validate_against_geant4.py

Compare g4gamma analytic spectra (Geant4 data, Sandia, LARA backends) against
a full rdecay01 Monte Carlo simulation at secular equilibrium.

Usage:
    # from repo root
    python test/validate_against_geant4.py <rdecay01_csv> <n_primaries> [output_dir]

    # U-238 example (10k-event run):
    python test/validate_against_geant4.py \\
        buildG4RadDecayExample/u238_AFtrue_h1_3.csv 10000

    # Th-232 example (1M-event run):
    python test/validate_against_geant4.py \\
        buildG4RadDecayExample/th232_AFtrue_h1_3.csv 1000000

The isotope (Z, A) is auto-detected from the CSV filename (e.g. th232_, u238_).
For U-238 and Th-232 the script uses isotope-specific peak annotations and
reference intensities. For other isotopes it produces unlabelled spectra.

Simulation settings matched here (moduleTest.mac + PhysicsList.cc):
    /process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year  -> t=-1 (secular eq.)
    /rdecay01/fullChain true                                      -> full chain
    radioactiveDecay->SetARM(true)                                -> K/L/M X-rays ON
    G4EmParameters::SetAugerCascade(true)                         -> Auger cascade ON
    G4EmParameters::SetDeexcitationIgnoreCut(true)                -> no energy cut on X-rays

g4gamma Geant4 backend uses full_xray_cascade=True (K→L→M fluorescence cascade,
no Auger data). The simulation's G4UAtomicDeexcitation produces the full Auger
cascade (K+L+M shells including Auger electrons), so g4gamma will slightly
undercount secondary vacancies from Auger transitions (~3% for Pb K-shell
where fluorescence yield ω_K≈0.97, larger for L-shell).
"""

import os
import re
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

N_BINS = 3000
E_MIN_KEV, E_MAX_KEV = 0.0, 3000.0

# g4gamma settings that match moduleTest.mac + PhysicsList.cc:
#   - secular equilibrium (thresholdForVeryLongDecayTime 1e60 year)
#   - ARM ON: radioactiveDecay->SetARM(true) + G4UAtomicDeexcitation
#             + SetAugerCascade(true), SetDeexcitationIgnoreCut(true)
#             → K/L/M X-rays + Auger cascade present in simulation
#   - 511 keV annihilation pairs from beta+ ARE present in the simulation
G4GAMMA_T = -1.0               # secular equilibrium
INCLUDE_XRAYS = True
FULL_XRAY_CASCADE = True       # Geant4 backend: K→L→M fluorescence cascade
INCLUDE_ANNIHILATION = True


# =============================================================================
# Isotope detection and per-isotope config
# =============================================================================

# Each config entry has keys:
#   name           : str          display name e.g. "U-238"
#   file_tag       : str          used in output filenames e.g. "u238"
#   major_peaks    : [(keV, lbl)] vertical lines on overview plot
#   zoom_regions   : [(lo, hi, title)] four zoom panels
#   high_energy    : (lo, hi, title, [keV...])  high-energy figure
#   reference_peaks: [(keV, nuclide, ENSDF_intensity, label)]
#                    gammas/primary at secular equilibrium from ENSDF

ISOTOPE_CONFIGS = {
    (92, 238, 0): {
        "name": "U-238",
        "file_tag": "u238",
        "major_peaks": [
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
        ],
        "zoom_regions": [
            (50,  250, "Low-energy region (Pb-214 / Ra-226 / Th-234 / Th-230)"),
            (270, 430, "Pb-214 peaks: 295.2, 352.0 keV"),
            (550, 700, "Bi-214 peak: 609.3 keV"),
            (1700, 1820, "Bi-214 peak: 1764.5 keV"),
        ],
        "high_energy": (950, 2400,
                        "U-238 chain — high-energy region (Bi-214 + Pa-234m)",
                        [1001, 1120, 1155, 1238, 1281, 1377, 1385,
                         1408, 1510, 1583, 1594, 1661, 1730, 1764, 2204]),
        "reference_peaks": [
            # Intensities = gammas/primary at secular equilibrium
            (63.3,   "Th-234",      0.037,  "63.3"),
            (75.5,   "Bi K X-ray",  0.093,  "75.5"),   # ARM; bin varies ±2 keV by G4EMLOW
            (92.5,   "Th/Rn Kα",   0.063,  "92.5"),   # ARM
            (186.2,  "Ra-226",      0.036,  "186.2"),
            (295.2,  "Pb-214",      0.184,  "295.2"),
            (351.9,  "Pb-214",      0.358,  "351.9"),
            (609.3,  "Bi-214",      0.461,  "609.3"),
            (768.4,  "Bi-214",      0.0491, "768.4"),
            (1001.0, "Pa-234m",     0.0084, "1001"),
            (1120.3, "Bi-214",      0.1491, "1120.3"),
            (1238.1, "Bi-214",      0.0583, "1238.1"),
            (1764.5, "Bi-214",      0.1530, "1764.5"),
            (2204.1, "Bi-214",      0.0491, "2204.1"),
        ],
    },
    (90, 232, 0): {
        "name": "Th-232",
        "file_tag": "th232",
        # At secular equilibrium each chain member has activity 1 per Th-232 primary,
        # EXCEPT Tl-208 which is fed only by the Bi-212 α-branch (35.94%), so
        # Tl-208 gammas/primary = intensity_per_Tl208_decay × 0.3594.
        "major_peaks": [
            (84.4,   "Th-228\n84 keV"),
            (238.6,  "Pb-212\n239 keV"),
            (338.3,  "Ac-228\n338 keV"),
            (463.0,  "Ac-228\n463 keV"),
            (583.1,  "Tl-208\n583 keV"),
            (727.3,  "Bi-212\n727 keV"),
            (794.9,  "Ac-228\n795 keV"),
            (911.2,  "Ac-228\n911 keV"),
            (968.9,  "Ac-228\n969 keV"),
            (2614.5, "Tl-208\n2615 keV"),
        ],
        "zoom_regions": [
            (50,  300, "Low-energy region (Th-228 / Ra-224 / Pb-212 / Ac-228)"),
            (220, 420, "Pb-212 239 keV, Ac-228 338 keV"),
            (550, 760, "Tl-208 583 keV, Bi-212 727 keV"),
            (880, 1010, "Ac-228 911 / 969 keV"),
        ],
        "high_energy": (550, 2700,
                        "Th-232 chain — high-energy region (Ac-228 / Tl-208)",
                        [583, 727, 860, 911, 969, 1588, 2614]),
        "reference_peaks": [
            # ENSDF intensities (gammas/primary at secular equilibrium).
            # At SE every member has activity 1.0 except Tl-208 (activity=0.3594,
            # fed only by Bi-212 α-branch 35.94%) and Po-212 (activity=0.6406).
            # Tl-208 intensities = per-Tl208-decay × 0.3594.
            # Bi-212 intensities = per-Bi212-decay × 1.0 (already include both branches).
            (84.4,   "Th-228",  0.0122, "84.4"),    # Geant4/sim give +25% (ICC diff)
            (238.6,  "Pb-212",  0.436,  "238.6"),
            (300.1,  "Pb-212",  0.033,  "300.1"),
            (338.3,  "Ac-228",  0.113,  "338.3"),   # Geant4 data gives +17% (data diff)
            (463.0,  "Ac-228",  0.044,  "463.0"),
            (583.1,  "Tl-208",  0.304,  "583.1"),   # 84.5% × 0.3594
            (727.3,  "Bi-212",  0.0665, "727.3"),   # 6.65% per Bi-212 decay
            (794.9,  "Ac-228",  0.042,  "794.9"),
            (911.2,  "Ac-228",  0.258,  "911.2"),   # Geant4 data gives +3.8% (data diff)
            (968.9,  "Ac-228",  0.158,  "968.9"),
            (2614.5, "Tl-208",  0.358,  "2614.5"),  # 99.75% × 0.3594
        ],
    },
}


def detect_isotope(csv_path):
    """
    Auto-detect (Z, A, M) from the CSV filename.
    Looks for patterns like 'th232', 'u238', 'cs137' at the start of the basename.
    Returns (Z, A, 0) if recognised, or None if not.
    """
    KNOWN = {
        ('u',  238): (92, 238, 0),
        ('th', 232): (90, 232, 0),
        ('ra', 226): (88, 226, 0),
        ('cs', 137): (55, 137, 0),
        ('co',  60): (27,  60, 0),
        ('k',   40): (19,  40, 0),
    }
    base = os.path.basename(csv_path).lower()
    for (sym, mass), zam in KNOWN.items():
        # \b treats '_' as a word char, so use explicit char-class boundaries instead
        if re.search(rf'(?<![a-z0-9]){sym}{mass}(?![a-z0-9])', base):
            return zam
    return None


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


def build_spectra(edges_keV, Z, A, M):
    """Build spectra from all three g4gamma backends. Returns dict of arrays."""
    edges_int = edges_keV * g.units.keV
    key = g.IsotopeKey(Z, A, M)

    # --- Geant4 backend: closest match to the rdecay01 simulation ---
    opts_g4 = g.SpectrumOptions()
    opts_g4.source = g.DataSource.Geant4
    opts_g4.full_xray_cascade = FULL_XRAY_CASCADE
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
    "sim":    "#1a1a2e",
    "geant4": "#2F4F8F",
    "sandia": "#C45A11",
    "lara":   "#1D9E75",
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


def safe_ratio(num, den, threshold=1e-5):
    mask = den > threshold
    out = np.full_like(num, np.nan)
    out[mask] = num[mask] / den[mask]
    return out


# =============================================================================
# Main plots
# =============================================================================

def make_plots(edges_keV, sim, analytic, outdir, cfg):
    os.makedirs(outdir, exist_ok=True)
    centers = 0.5 * (edges_keV[:-1] + edges_keV[1:])
    name = cfg["name"]
    tag  = cfg["file_tag"]

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
    ax_main.set_title(
        f"{name} secular-equilibrium gamma spectrum — dataset comparison", loc="left"
    )
    ax_main.legend(loc="upper right", fontsize=8.5)
    ax_main.set_xticklabels([])

    for e_kev, _lbl in cfg["major_peaks"]:
        ax_main.axvline(e_kev, color="#ccc", lw=0.5, zorder=0)

    # Ratio panel
    ref_key = "sim" if _has_sim() else "sandia"
    ref = sim if _has_sim() else analytic["sandia"]
    ref_label = "MC sim" if _has_sim() else "Sandia"

    for k in ("geant4", "sandia", "lara"):
        if k == ref_key:
            continue
        ratio = safe_ratio(analytic[k], ref)
        x, y = step_xy(edges_keV, ratio)
        ax_ratio.plot(x, y, color=COLORS[k], lw=0.9,
                      linestyle=STYLES[k].get("linestyle", "-"))

    if _has_sim():
        ratio_g4 = safe_ratio(analytic["geant4"], sim)
        x, y = step_xy(edges_keV, ratio_g4)
        ax_ratio.plot(x, y, color=COLORS["geant4"], lw=0.9)

    ax_ratio.axhline(1.0, color="#888", lw=0.8, linestyle="--")
    ax_ratio.set_ylim(0.5, 1.5)
    ax_ratio.set_xlim(0, 3000)
    ax_ratio.set_xlabel("Energy (keV)")
    ax_ratio.set_ylabel(f"Ratio to {ref_label}")

    path = os.path.join(outdir, f"01_{tag}_overview.png")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {path}")

    # -------------------------------------------------------------------------
    # Fig 2: Four zoom panels for key peak regions
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    axes = axes.flatten()

    for ax, (lo, hi, title) in zip(axes, cfg["zoom_regions"]):
        mask = (centers >= lo) & (centers <= hi)
        e_zoom = edges_keV[:-1][mask]
        e_zoom = np.append(e_zoom, edges_keV[1:][mask][-1])

        if _has_sim():
            ax.bar(centers[mask], sim[mask], width=1.0,
                   color=COLORS["sim"], alpha=0.35, label=LABELS["sim"])

        for k in ("geant4", "sandia", "lara"):
            c = analytic[k][mask]
            x, y = step_xy(e_zoom, c)
            ax.plot(x, y, color=COLORS[k], label=LABELS[k], **dict(STYLES[k]))

        ax.set_xlim(lo, hi)
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("γ / primary decay")
        ax.set_title(title, loc="left", fontsize=9)
        ax.legend(loc="upper right", fontsize=7.5)

    path = os.path.join(outdir, f"02_{tag}_zoom_peaks.png")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {path}")

    # -------------------------------------------------------------------------
    # Fig 3: High-energy region
    # -------------------------------------------------------------------------
    lo, hi, he_title, he_peaks = cfg["high_energy"]
    fig, ax = plt.subplots(figsize=(12, 5))
    mask = (centers >= lo) & (centers <= hi)
    e_zoom = edges_keV[:-1][mask]
    e_zoom = np.append(e_zoom, edges_keV[1:][mask][-1])

    if _has_sim():
        ax.bar(centers[mask], sim[mask], width=1.0,
               color=COLORS["sim"], alpha=0.35, label=LABELS["sim"])
    for k in ("geant4", "sandia", "lara"):
        x, y = step_xy(e_zoom, analytic[k][mask])
        ax.plot(x, y, color=COLORS[k], label=LABELS[k], **STYLES[k])

    for e_kev in he_peaks:
        if lo <= e_kev <= hi:
            ax.axvline(e_kev, color="#ddd", lw=0.5, zorder=0)

    ax.set_xlim(lo, hi)
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("γ / primary decay")
    ax.set_title(he_title, loc="left")
    ax.legend(fontsize=8.5)

    path = os.path.join(outdir, f"03_{tag}_high_energy.png")
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

    path = os.path.join(outdir, f"04_{tag}_ratios.png")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {path}")


# =============================================================================
# Peak comparison table
# =============================================================================

def print_peak_table(edges_keV, sim, analytic, cfg):
    reference_peaks = cfg["reference_peaks"]
    print("\n" + "=" * 95)
    print(f"{'Peak (keV)':<12} {'Nuclide':<12} {'ENSDF':>8}", end="")
    if sim is not None:
        print(f" {'MC sim':>10} {'sim-err%':>9}", end="")
    for k in ("geant4", "sandia", "lara"):
        print(f" {k:>10} {'err%':>6}", end="")
    print()
    print("-" * 95)

    for e_kev, nuc, ref_int, label in reference_peaks:
        idx = int(np.searchsorted(edges_keV, e_kev, side="right")) - 1
        idx = max(0, min(idx, len(edges_keV) - 2))

        row = f"{label:<12} {nuc:<12} {ref_int:>8.4f}"
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
    print("  Geant4 backend: full_xray_cascade=True (K→L→M fluorescence cascade, no Auger data)")
    print("  Sandia / LARA: X-ray and annihilation emissions baked into dataset files")
    print("  Simulation uses full Auger cascade; g4gamma drops Auger secondary vacancies (~3% for Pb K)")
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

    seen = {}
    for k in ("geant4", "sandia", "lara"):
        for c in analytic["contribs"][k]:
            key = str(c.isotope)
            if key not in seen:
                seen[key] = {}
            seen[key][k] = c

    for iso_str, providers in seen.items():
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

    csv_path    = sys.argv[1]
    n_primaries = int(float(sys.argv[2]))   # accept 1e7 notation
    outdir      = sys.argv[3] if len(sys.argv) > 3 else os.path.join(REPO, "validation_plots")

    # --- detect isotope -------------------------------------------------------
    zam = detect_isotope(csv_path)
    if zam is None:
        print("WARNING: could not detect isotope from filename; defaulting to U-238.")
        zam = (92, 238, 0)
    Z, A, M = zam
    cfg = ISOTOPE_CONFIGS.get((Z, A, M))
    if cfg is None:
        print(f"WARNING: no config for ({Z},{A},{M}); using generic (no peak annotations).")
        cfg = {
            "name": f"Z={Z} A={A}",
            "file_tag": f"z{Z}a{A}",
            "major_peaks": [],
            "zoom_regions": [
                (0, 500, "0–500 keV"),
                (500, 1000, "500–1000 keV"),
                (1000, 2000, "1000–2000 keV"),
                (2000, 3000, "2000–3000 keV"),
            ],
            "high_energy": (1000, 3000, f"Z={Z} A={A} — high-energy region", []),
            "reference_peaks": [],
        }

    print(f"Isotope: {cfg['name']} (Z={Z}, A={A}, M={M})")

    # --- load simulation -------------------------------------------------------
    print(f"Reading rdecay01 CSV: {csv_path}")
    print(f"  n_primaries = {n_primaries:,}")
    sim = load_rdecay01_csv(csv_path, n_primaries)
    print(f"  {len(sim)} bins, total {sim.sum():.3f} γ/primary")

    # Reference bin for statistical uncertainty: largest peak in REFERENCE_PEAKS,
    # or fall back to 609 keV for unknown isotopes.
    ref_peak_keV = 609.3
    if cfg["reference_peaks"]:
        ref_peak_keV = max(cfg["reference_peaks"], key=lambda r: r[2])[0]
    edges_keV_tmp = np.linspace(E_MIN_KEV, E_MAX_KEV, N_BINS + 1)
    ref_idx = int(np.searchsorted(edges_keV_tmp, ref_peak_keV, side="right")) - 1
    print(f"  statistical uncertainty on {ref_peak_keV:.1f} keV peak: "
          f"~{1/np.sqrt(max(sim[ref_idx]*n_primaries, 1))*100:.1f}%")

    # --- build analytic spectra -----------------------------------------------
    edges_keV = edges_keV_tmp
    print(f"\nBuilding analytic spectra for {cfg['name']}...")
    analytic = build_spectra(edges_keV, Z, A, M)
    for k in ("geant4", "sandia", "lara"):
        print(f"  {k:8s}: {analytic[k].sum():.4f} γ/primary total")

    # --- print tables ----------------------------------------------------------
    if cfg["reference_peaks"]:
        print_peak_table(edges_keV, sim, analytic, cfg)
    print_chain_summary(analytic)

    # --- make plots ------------------------------------------------------------
    print(f"Writing plots to: {os.path.abspath(outdir)}")
    make_plots(edges_keV, sim, analytic, outdir, cfg)

    # --- write difference explanation ------------------------------------------
    _print_explanation(sim, analytic, cfg)


def _print_explanation(sim, analytic, cfg):
    name = cfg["name"]
    print(f"""
=============================================================================
EXPLANATION OF DIFFERENCES BETWEEN DATASETS  [{name}]
=============================================================================

1. Geant4 data vs. Sandia/LARA
   The Geant4 backend reads the same data files ($G4RADIOACTIVEDATA,
   $G4LEVELGAMMADATA) that the rdecay01 Monte Carlo simulation uses, so it
   is the reference for reproducing the simulation output.

   Sandia and LARA source their data from ENSDF evaluations which are
   updated independently of the Geant4 data files. Small differences in
   branching ratios and level energies accumulate across long decay chains.

2. Total gamma count differences (Geant4 > Sandia ≈ LARA)
   The Geant4 backend total with full_xray_cascade=True is generally higher
   than Sandia/LARA because it resolves short-lived isomers as separate
   chain nodes with their own gamma yields; Sandia/LARA fold these into
   parent entries. The Geant4 total typically agrees with simulation to
   within ~1%.

3. Low-energy region (< 150 keV) — X-ray and fluorescence region
   The simulation runs ARM=true with a full Auger cascade (G4UAtomicDeexcitation,
   K+L+M shells, SetDeexcitationIgnoreCut=true). g4gamma's Geant4 backend
   computes a K→L→M fluorescence cascade (full_xray_cascade=True), propagating
   secondary vacancies from each fluorescence transition.

   Remaining differences:
   a) No Auger data: g4gamma reads fl-tr-pr-Z.dat (fluorescence only) and
      does not read au-tr-pr-Z.dat (Auger transitions). When an Auger
      transition occurs instead of X-ray emission it creates two secondary
      vacancies that g4gamma drops. For Pb K-shell ω_K≈0.97 so ~3% error.
   b) Energy offset ±1-2 keV: different G4EMLOW versions can shift X-ray
      peaks into adjacent 1 keV bins. Total X-ray yield in the region
      agrees better than bin-by-bin.

   Sandia and LARA have X-ray intensities baked into their emission tables
   from ENSDF evaluations and agree well with the simulation in this region.

4. Strong-peak agreement
   For well-measured ENSDF peaks (> 5% intensity), all datasets typically
   agree with the simulation to within ±2% at n=10M events. Residual
   differences above ~3σ are dataset differences, not Monte Carlo noise.

5. Simulation statistics
   Statistical uncertainty on the largest peak scales as 1/√(I × n_primaries):
   n=10k → ~1.5%, n=100k → ~0.5%, n=1M → ~0.15%, n=10M → ~0.05%.
=============================================================================
""")


if __name__ == "__main__":
    main()
