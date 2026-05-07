#!/usr/bin/env python3
"""
showcase.py — comprehensive feature demo for g4gamma against real Geant4 data.

Produces a series of PNG plots demonstrating each feature:
  Fig 1.  Cs-137 spectrum at SE (the hello-world)
  Fig 2.  K-40 spectrum: with vs without B+ annihilation gammas
  Fig 3.  K-40 spectrum: with vs without K X-rays
  Fig 4.  Co-60 spectrum (two-photon cascade)
  Fig 5.  Cs-137 build-up: A(Ba-137m) vs time
  Fig 6.  U-238 chain: full secular-equilibrium spectrum + chain dump
  Fig 7.  Th-232 chain: full secular-equilibrium spectrum + chain dump
  Fig 8.  Spectrum format effects: linear vs log binning
  Fig 9.  Na-22 spectrum: 511 + 1274 keV with annihilation
  Fig 10. Multi-isotope comparison panel (NORM-relevant)

Run from anywhere; auto-locates the build directory:
    python3 test/showcase.py [output_dir]

Requires matplotlib + numpy.
"""

import os
import sys
import numpy as np

# ----------- locate the g4gamma module --------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
for cand in [
    os.path.join(HERE, "..", "build"),
    os.path.join(HERE, "build"),
    "build",
    ".",
]:
    if not os.path.isdir(cand):
        continue
    if any(f.startswith("g4gamma") and f.endswith(".so") for f in os.listdir(cand)):
        sys.path.insert(0, cand)
        break

try:
    import g4gamma as g
except ImportError as e:
    print(f"ERROR: cannot import g4gamma. Build it first.\n  {e}")
    sys.exit(1)

try:
    import matplotlib
    matplotlib.use("Agg")  # headless
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, NullFormatter
except ImportError:
    print("ERROR: matplotlib required. pip install matplotlib")
    sys.exit(1)


# ----------- output dir -----------------------------------------------------
OUTDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "showcase_plots")
os.makedirs(OUTDIR, exist_ok=True)
print(f"Writing plots to: {os.path.abspath(OUTDIR)}\n")


# ----------- detect available data sources --------------------------------
SOURCES = {"geant4": True}
# SandiaDecay XML location (auto-detect, allow override via env)
SANDIA_XML = os.environ.get("SANDIA_DECAY_XML", "")
if not SANDIA_XML:
    for cand in [
        os.path.join(HERE, "..", "data", "sandia", "sandia.decay.nocoinc.min.xml.gz"),
        os.path.join(HERE, "..", "data", "sandia", "sandia.decay.nocoinc.min.xml"),
        os.path.join(HERE, "..", "data", "sandia", "sandia.decay.xml.gz"),
        os.path.join(HERE, "..", "data", "sandia", "sandia.decay.xml"),
        os.path.join(HERE, "..", "data", "sandia", "sandia.decay.min.xml"),
    ]:
        if os.path.isfile(cand):
            SANDIA_XML = cand
            break
SOURCES["sandia"] = bool(SANDIA_XML)

# LARA data directory
LARA_DIR = os.environ.get("LARA_DATA_DIR", "")
if not LARA_DIR:
    cand = os.path.join(HERE, "..", "data", "lara")
    if os.path.isdir(cand):
        has_loose   = any(f.endswith(".lara.txt") for f in os.listdir(cand))
        has_tarball = (os.path.isfile(os.path.join(cand, "lara.tar.gz")) or
                       os.path.isfile(os.path.join(cand, "lara.tar")))
        if has_loose or has_tarball:
            LARA_DIR = cand
SOURCES["lara"] = bool(LARA_DIR)

print(f"Available sources: {[s for s, ok in SOURCES.items() if ok]}")
if SANDIA_XML: print(f"  Sandia XML: {SANDIA_XML}")
if LARA_DIR:   print(f"  LARA dir:   {LARA_DIR}")
print()


# ----------- plot styling ---------------------------------------------------
plt.rcParams.update({
    "figure.dpi":         110,
    "savefig.dpi":        160,
    "font.family":        "DejaVu Sans",
    "font.size":          10,
    "axes.titlesize":     11,
    "axes.labelsize":     10,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          True,
    "grid.alpha":         0.25,
    "grid.linewidth":     0.5,
    "legend.frameon":     False,
    "lines.linewidth":    1.0,
})

# Consistent colour palette (Tim, your thesis is going to want a consistent style)
C = {
    "primary":   "#2F4F8F",   # navy
    "secondary": "#C45A11",   # rust
    "accent":    "#1D9E75",   # teal
    "warning":   "#BA7517",   # amber
    "muted":     "#666666",
    "fill":      "#B5D4F4",   # pale blue
}


# ----------- helpers --------------------------------------------------------
def safe_build(primary, t, edges, source="geant4", **kwargs):
    """Build a spectrum, returning (counts, contributions) or (None, msg) on failure."""
    if source == "sandia":
        src = g.DataSource.Sandia
    elif source == "lara":
        src = g.DataSource.Lara
    else:
        src = g.DataSource.Geant4
    extra = {}
    if source == "sandia" and SANDIA_XML:
        extra["sandia_xml"] = SANDIA_XML
    if source == "lara" and LARA_DIR:
        extra["lara_dir"] = LARA_DIR
    try:
        res = g.build_spectrum(primary, t, edges, source=src, **extra, **kwargs)
        return np.array(res.counts), list(res.contributions)
    except Exception as e:
        return None, str(e)


def plot_spectrum(ax, edges, counts, label=None, color=C["primary"],
                   fill_alpha=0.25, line_kwargs=None):
    """Plot a histogram-style spectrum onto ax. Edges in keV, counts per keV bin."""
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths  = np.diff(edges)
    line_kwargs = line_kwargs or {}
    # pre-zero-pad so step looks correct at the start
    e = np.concatenate([[edges[0]], np.repeat(edges[1:-1], 2), [edges[-1]]])
    c = np.repeat(counts, 2)
    ax.plot(e, c, color=color, label=label, **line_kwargs)
    if fill_alpha > 0:
        ax.fill_between(e, 0, c, color=color, alpha=fill_alpha, linewidth=0)
    return centers, widths


def annotate_peaks(ax, edges, counts, max_peaks=5, min_count=1e-3,
                   y_factor=1.05, fontsize=8, color=None):
    """Label the largest peaks with their bin-centre energies."""
    centers = 0.5 * (edges[:-1] + edges[1:])
    nz = np.where(counts > min_count)[0]
    if len(nz) == 0:
        return
    # take the top N by intensity
    idx = nz[np.argsort(counts[nz])[::-1][:max_peaks]]
    idx.sort()
    for i in idx:
        ax.annotate(f"{centers[i]:.0f} keV",
                    xy=(centers[i], counts[i]),
                    xytext=(0, 6), textcoords="offset points",
                    ha="center", va="bottom",
                    fontsize=fontsize,
                    color=color or C["muted"])


def style_ax(ax, title=None, xlabel="Energy (keV)", ylabel="γ / primary decay",
              ylog=False):
    if title:  ax.set_title(title, loc="left")
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    if ylog:
        ax.set_yscale("log")
    ax.tick_params(direction="out", length=3)


def save(fig, name):
    path = os.path.join(OUTDIR, name)
    # tight_layout() warns on gridspec'd figures; suppress that one warning
    # rather than restructuring all the figure-building code.
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {name}")


# Default bin grid: 1 keV bins from 0-3000 keV
EDGES_KEV = np.linspace(0, 3000, 3001)
EDGES = EDGES_KEV * g.units.keV


# ===========================================================================
# Fig 1 — Cs-137 secular equilibrium
# ===========================================================================
print("Fig 1: Cs-137 secular equilibrium")

counts, contrib = safe_build(g.IsotopeKey(55, 137, 0), -1.0, EDGES)
if counts is not None:
    fig, ax = plt.subplots(figsize=(9, 4))
    plot_spectrum(ax, EDGES_KEV, counts, label="Cs-137 (SE)", color=C["primary"])
    annotate_peaks(ax, EDGES_KEV, counts, max_peaks=3)
    ax.set_xlim(0, 800)
    ax.set_ylim(0, max(counts.max() * 1.15, 1e-6))
    style_ax(ax, "Cs-137 secular-equilibrium γ spectrum")

    # chain caption — bottom right, away from the peaks
    text = "chain: " + " → ".join(
        f"{c.isotope}(A={c.activity:.3f})" for c in contrib if c.activity > 0 or not c.isotope.M)
    ax.text(0.99, 0.05, text, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8, color=C["muted"], family="monospace")
    save(fig, "01_cs137_se.png")
else:
    print(f"  SKIPPED: {contrib}")


# ===========================================================================
# Fig 2 — Annihilation toggle (Na-22 fallback to K-40)
# ===========================================================================
print("Fig 2: annihilation toggle")

# Try Na-22 first (heavy B+ emitter -- best demo). Fall back to K-40 (which
# has a tiny 0.001% B+ branch, making the toggle a near-no-op).
ann_primary, ann_label = None, None
for primary, label in [(g.IsotopeKey(11, 22, 0), "Na-22"),
                       (g.IsotopeKey(19, 40, 0), "K-40")]:
    test, _ = safe_build(primary, -1.0, EDGES)
    if test is not None and test.sum() > 0:
        ann_primary, ann_label = primary, label
        break

if ann_primary is not None:
    c_with, _    = safe_build(ann_primary, -1.0, EDGES, include_annihilation=True)
    c_without, _ = safe_build(ann_primary, -1.0, EDGES, include_annihilation=False)

    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    plot_spectrum(axes[0], EDGES_KEV, c_with,    label="includes 511 keV pairs",
                  color=C["primary"])
    plot_spectrum(axes[1], EDGES_KEV, c_without, label="cascade γ only",
                  color=C["secondary"])
    xmax = 1700 if ann_label == "K-40" else 1500
    for ax in axes:
        ax.set_xlim(0, xmax)
        ax.set_ylim(bottom=0)
        ax.legend(loc="upper right")
    style_ax(axes[0], f"{ann_label} — annihilation pairs included", xlabel="")
    style_ax(axes[1], f"{ann_label} — annihilation pairs excluded")
    annotate_peaks(axes[0], EDGES_KEV, c_with,    max_peaks=2)
    annotate_peaks(axes[1], EDGES_KEV, c_without, max_peaks=2)

    diff = c_with.sum() - c_without.sum()
    fig.suptitle(f"{ann_label}: B+ annihilation contributes {diff:.4g} γ per primary decay",
                 fontsize=10, y=1.00)
    save(fig, "02_annihilation_toggle.png")
else:
    print("  SKIPPED")


# ===========================================================================
# Fig 3 — K-40, with vs without K X-rays
# ===========================================================================
print("Fig 3: K-40, K X-ray toggle")

c_no, _   = safe_build(g.IsotopeKey(19, 40, 0), -1.0, EDGES, include_xrays=False)
c_yes, _  = safe_build(g.IsotopeKey(19, 40, 0), -1.0, EDGES, include_xrays=True)

if c_no is not None and c_yes is not None:
    fig, ax = plt.subplots(figsize=(9, 4.5))

    # Diff highlight
    diff = c_yes - c_no
    plot_spectrum(ax, EDGES_KEV, c_yes, label="with X-rays (full)",
                  color=C["primary"], fill_alpha=0.15)
    plot_spectrum(ax, EDGES_KEV, c_no,  label="cascade only (no X-rays)",
                  color=C["secondary"], fill_alpha=0)
    ax.set_xlim(0, 1700)
    style_ax(ax, "K-40: optional K-shell X-rays from EC vacancies", ylog=True)
    ax.set_ylim(1e-5, 1)
    ax.legend(loc="upper right")

    # Inset for the X-ray region
    inset = ax.inset_axes([0.08, 0.50, 0.30, 0.40])
    plot_spectrum(inset, EDGES_KEV, c_yes, color=C["primary"])
    plot_spectrum(inset, EDGES_KEV, c_no,  color=C["secondary"], fill_alpha=0)
    inset.set_xlim(0, 10)
    ymax = max(c_yes[:10].max(), 1e-4) * 1.3
    inset.set_ylim(0, ymax)
    inset.tick_params(labelsize=8)
    inset.set_title("Ar K X-rays (zoom)", fontsize=8, loc="left")
    inset.set_xlabel("keV", fontsize=8)
    inset.grid(alpha=0.2)
    save(fig, "03_k40_xrays_toggle.png")
else:
    print(f"  SKIPPED")


# ===========================================================================
# Fig 4 — Co-60 (two-line cascade, the classic)
# ===========================================================================
print("Fig 4: Co-60 two-line cascade")

counts, contrib = safe_build(g.IsotopeKey(27, 60, 0), -1.0, EDGES)
if counts is not None:
    fig, ax = plt.subplots(figsize=(9, 4))
    plot_spectrum(ax, EDGES_KEV, counts, label="Co-60 (SE)", color=C["primary"])
    annotate_peaks(ax, EDGES_KEV, counts, max_peaks=4, min_count=1e-3)
    ax.set_xlim(0, 1500)
    ax.set_ylim(bottom=0)
    style_ax(ax, "Co-60 secular-equilibrium γ spectrum (1173 + 1332 keV cascade)")

    total = counts.sum()
    ax.text(0.99, 0.95, f"total: {total:.4f} γ / primary decay",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=9, color=C["muted"])
    save(fig, "04_co60.png")
else:
    print(f"  SKIPPED: {contrib}")


# ===========================================================================
# Fig 5 — Cs-137 build-up: A(Ba-137m) vs time
# ===========================================================================
print("Fig 5: Cs-137 / Ba-137m build-up vs time")

# Sweep across times. Ba-137m T_half = 153.12s; full SE reached well before 1 hr.
times_s = np.concatenate([[0.0],
                          np.geomspace(1, 3600, 100)])
peak_661 = np.empty_like(times_s)
A_cs     = np.empty_like(times_s)
A_ba     = np.empty_like(times_s)

ok = True
for i, t_s in enumerate(times_s):
    t = t_s * g.units.s if t_s > 0 else 0.0
    counts, contrib = safe_build(g.IsotopeKey(55, 137, 0), t, EDGES)
    if counts is None:
        print(f"  SKIPPED: {contrib}"); ok = False; break
    peak_661[i] = counts[661] if 661 < len(counts) else 0
    cs = next((c.activity for c in contrib if c.isotope.Z == 55), 0.0)
    ba = next((c.activity for c in contrib if c.isotope.Z == 56 and c.isotope.M == 1), 0.0)
    A_cs[i] = cs; A_ba[i] = ba

if ok:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    # Left: activities
    axes[0].plot(times_s, A_cs, color=C["primary"], label="A(Cs-137)")
    axes[0].plot(times_s, A_ba, color=C["secondary"], label="A(Ba-137m)")
    axes[0].axvline(153.12, color=C["muted"], ls=":", lw=0.8)
    axes[0].text(153.12, 0.5, "T½(Ba-137m)", rotation=90,
                 fontsize=8, color=C["muted"], va="center", ha="right")
    axes[0].set_xscale("log")
    axes[0].set_xlim(1, 3600)
    axes[0].set_ylim(0, 1.05)
    axes[0].legend(loc="center right")
    style_ax(axes[0], "Cs-137 → Ba-137m build-up to transient equilibrium",
             xlabel="time since A(Cs-137) = 1 Bq  (s)", ylabel="activity (Bq)")

    # Right: 661 keV peak intensity
    axes[1].plot(times_s, peak_661, color=C["accent"])
    axes[1].axhline(0.852385, color=C["muted"], ls=":", lw=0.8)
    axes[1].text(2, 0.852, "secular eqm", fontsize=8, color=C["muted"], va="bottom")
    axes[1].set_xscale("log")
    axes[1].set_xlim(1, 3600)
    axes[1].set_ylim(0, 1.0)
    style_ax(axes[1], "Time-resolved 661.7 keV peak intensity",
             xlabel="time (s)", ylabel="γ / primary decay")
    save(fig, "05_cs137_buildup.png")


# ===========================================================================
# Fig 6 — U-238 chain
# ===========================================================================
print("Fig 6: U-238 full chain at SE")

# Try Geant4 first; if no data, fall back to Sandia (we ship the data).
counts, contrib = safe_build(g.IsotopeKey(92, 238, 0), -1.0, EDGES, source="geant4")
src_used = "Geant4"
if (counts is None or counts.sum() == 0) and SOURCES["sandia"]:
    counts, contrib = safe_build(g.IsotopeKey(92, 238, 0), -1.0, EDGES, source="sandia")
    src_used = "SandiaDecay"

if counts is not None and counts.sum() > 0:
    fig = plt.figure(figsize=(11, 6))
    gs = fig.add_gridspec(2, 2, height_ratios=[2, 1.5], hspace=0.45, wspace=0.25)
    ax_spec = fig.add_subplot(gs[0, :])
    ax_chain_lin = fig.add_subplot(gs[1, 0])
    ax_chain_log = fig.add_subplot(gs[1, 1])

    plot_spectrum(ax_spec, EDGES_KEV, counts, color=C["primary"])
    ax_spec.set_xlim(0, 3000)
    style_ax(ax_spec, f"U-238 chain at secular equilibrium  ({src_used})", ylog=True)
    ax_spec.set_ylim(1e-5, max(counts.max()*2, 1))

    # Annotate top peaks with isotope guesses
    annotate_peaks(ax_spec, EDGES_KEV, counts, max_peaks=8, min_count=1e-3)

    # Chain bar charts: members by activity (linear + log)
    members = [c for c in (contrib or []) if c.activity > 0]
    members.sort(key=lambda c: c.activity, reverse=True)
    if len(members) > 0:
        names = [str(c.isotope) for c in members][:25]
        acts  = [c.activity for c in members][:25]
        y = np.arange(len(names))[::-1]
        for ax_c, scale, title in [(ax_chain_lin, "linear", "Chain activity (linear)"),
                                    (ax_chain_log, "log",    "Chain activity (log)")]:
            ax_c.barh(y, acts, color=C["fill"], edgecolor=C["primary"], linewidth=0.5)
            ax_c.set_yticks(y); ax_c.set_yticklabels(names, fontsize=8, family="monospace")
            ax_c.set_xlabel("relative activity")
            ax_c.set_title(title, loc="left", fontsize=10)
            if scale == "log":
                ax_c.set_xscale("log")
                ax_c.set_xlim(max(min(acts)/3, 1e-12), 2)
    else:
        for ax_c in (ax_chain_lin, ax_chain_log):
            ax_c.text(0.5, 0.5, "no active chain members", transform=ax_c.transAxes,
                      ha="center", va="center", color=C["muted"])
            ax_c.set_xticks([]); ax_c.set_yticks([])

    save(fig, "06_u238_chain.png")
else:
    msg = contrib if isinstance(contrib, str) else "no decay data for U-238 (check $G4RADIOACTIVEDATA)"
    print(f"  SKIPPED: {msg}")


# ===========================================================================
# Fig 7 — Th-232 chain
# ===========================================================================
print("Fig 7: Th-232 full chain at SE")

counts, contrib = safe_build(g.IsotopeKey(90, 232, 0), -1.0, EDGES, source="geant4")
src_used = "Geant4"
if (counts is None or counts.sum() == 0) and SOURCES["sandia"]:
    counts, contrib = safe_build(g.IsotopeKey(90, 232, 0), -1.0, EDGES, source="sandia")
    src_used = "SandiaDecay"

if counts is not None and counts.sum() > 0:
    fig, ax = plt.subplots(figsize=(10, 4.5))
    plot_spectrum(ax, EDGES_KEV, counts, color=C["primary"])
    ax.set_xlim(0, 3000)
    style_ax(ax, f"Th-232 chain at secular equilibrium  ({src_used})", ylog=True)
    ax.set_ylim(1e-5, max(counts.max()*2, 1))
    annotate_peaks(ax, EDGES_KEV, counts, max_peaks=10, min_count=1e-3)

    nm = sum(1 for c in (contrib or []) if c.activity > 0)
    ax.text(0.99, 0.95, f"{nm} active chain members",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=9, color=C["muted"])
    save(fig, "07_th232_chain.png")
else:
    msg = contrib if isinstance(contrib, str) else "no decay data for Th-232"
    print(f"  SKIPPED: {msg}")


# ===========================================================================
# Fig 8 — Binning options
# ===========================================================================
print("Fig 8: linear vs log binning")

# Pick a primary that's available in the dataset. U-238 is ideal (broad
# spectrum), Cs-137 is the fallback.
binning_primary = None
binning_label = None
for primary, label in [(g.IsotopeKey(92, 238, 0), "U-238 chain"),
                       (g.IsotopeKey(55, 137, 0), "Cs-137 SE"),
                       (g.IsotopeKey(19,  40, 0), "K-40 SE")]:
    test_counts, _ = safe_build(primary, -1.0, EDGES)
    if test_counts is not None and test_counts.sum() > 0:
        binning_primary = primary
        binning_label = label
        break

if binning_primary is not None:
    grids = {
        "1 keV linear":         np.linspace(0, 3000, 3001),
        "10 keV linear":        np.linspace(0, 3000,  301),
        "log (50 bins/decade)": np.logspace(np.log10(10), np.log10(3000), 250),
    }
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    for ax, (name, edges_kev) in zip(axes, grids.items()):
        edges = edges_kev * g.units.keV
        counts, _ = safe_build(binning_primary, -1.0, edges)
        if counts is None or counts.sum() == 0: continue
        centers = 0.5 * (edges_kev[:-1] + edges_kev[1:])
        ax.bar(centers, counts, width=np.diff(edges_kev),
               align="center", color=C["fill"], edgecolor=C["primary"], linewidth=0.4)
        ax.set_xlim(10, 3000)
        if "log" in name: ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(1e-5, max(counts.max()*2, 1))
        style_ax(ax, name, ylog=True)
    fig.suptitle(f"{binning_label} — same data, different bin grids", y=1.02)
    save(fig, "08_binning_grids.png")
else:
    print("  SKIPPED — no isotope available")


# ===========================================================================
# Fig 9 — Na-22 (B+ emitter)
# ===========================================================================
print("Fig 9: Na-22 with annihilation")

# Try with annihilation, fallback to without
counts, contrib = safe_build(g.IsotopeKey(11, 22, 0), -1.0, EDGES, include_annihilation=True)
if counts is not None and counts.sum() > 0:
    fig, ax = plt.subplots(figsize=(9, 4))
    plot_spectrum(ax, EDGES_KEV, counts, color=C["primary"])
    annotate_peaks(ax, EDGES_KEV, counts, max_peaks=2, min_count=1e-2)
    ax.set_xlim(0, 1500)
    ax.set_ylim(bottom=0)
    style_ax(ax, "Na-22 — B+ emitter with 511 + 1274 keV signature")
    ax.text(0.99, 0.95,
            f"total: {counts.sum():.4f} γ / decay\n"
            f"(2× 511 keV pairs from B+ + 1274.5 keV from Ne-22 cascade)",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color=C["muted"])
    save(fig, "09_na22.png")
else:
    print("  SKIPPED (Na-22 not in dataset?)")


# ===========================================================================
# Fig 10 — Multi-isotope NORM panel
# ===========================================================================
print("Fig 10: NORM-relevant isotope panel")

panel_isos = [
    (g.IsotopeKey(19,  40,  0), "K-40",    C["primary"]),
    (g.IsotopeKey(55, 137,  0), "Cs-137",  C["secondary"]),
    (g.IsotopeKey(27,  60,  0), "Co-60",   C["accent"]),
    (g.IsotopeKey(92, 238,  0), "U-238",   C["warning"]),
    (g.IsotopeKey(90, 232,  0), "Th-232",  "#933"),
    (g.IsotopeKey(88, 226,  0), "Ra-226",  "#525"),
]

n_avail = []
for primary, label, color in panel_isos:
    counts, _ = safe_build(primary, -1.0, EDGES)
    if counts is not None and counts.sum() > 0:
        n_avail.append((primary, label, color, counts))

if len(n_avail) > 0:
    n = len(n_avail)
    ncols = 2
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 2.8*nrows), sharex=True)
    axes = np.atleast_1d(axes).ravel()

    for ax, (primary, label, color, counts) in zip(axes, n_avail):
        plot_spectrum(ax, EDGES_KEV, counts, color=color)
        ax.set_xlim(0, 3000)
        ax.set_yscale("log")
        ax.set_ylim(1e-5, max(counts.max()*2, 1))
        ax.set_title(f"{label} (Σγ = {counts.sum():.3f}/decay)",
                     loc="left", fontsize=10)
        ax.grid(alpha=0.25)

    # Hide unused subplots
    for ax in axes[len(n_avail):]:
        ax.axis("off")

    for ax in axes[-ncols:]:
        ax.set_xlabel("Energy (keV)")
    for i in range(0, len(axes), ncols):
        if i < len(n_avail):
            axes[i].set_ylabel("γ / decay")

    fig.suptitle("Comparison: NORM-relevant isotope spectra at secular equilibrium",
                 y=1.00)
    save(fig, "10_norm_panel.png")
else:
    print("  SKIPPED — no isotopes available")


# ===========================================================================
# Fig 11 — Tri-source cross-validation (Geant4 vs Sandia vs LARA)
# ===========================================================================
print("Fig 11: Geant4 vs Sandia vs LARA cross-validation")

active_sources = [s for s in ("geant4", "sandia", "lara") if SOURCES[s]]
if len(active_sources) >= 2:
    iso_panel = [
        (g.IsotopeKey(55, 137, 0), "Cs-137",  (0, 800)),
        (g.IsotopeKey(19,  40, 0), "K-40",    (0, 1700)),
        (g.IsotopeKey(27,  60, 0), "Co-60",   (0, 1500)),
        (g.IsotopeKey(92, 238, 0), "U-238",   (0, 3000)),
    ]

    color_map = {
        "geant4": C["primary"],
        "sandia": C["secondary"],
        "lara":   C["accent"],
    }
    label_map = {"geant4": "Geant4", "sandia": "SandiaDecay", "lara": "LARA/DDEP"}
    style_map = {"geant4": "-", "sandia": "--", "lara": ":"}

    panel_data = []
    for primary, label, xrange in iso_panel:
        rows = []
        for src in active_sources:
            counts, _ = safe_build(primary, -1.0, EDGES, source=src)
            if counts is not None and counts.sum() > 0:
                rows.append((src, counts))
        if len(rows) >= 2:
            panel_data.append((primary, label, xrange, rows))

    if panel_data:
        n = len(panel_data)
        fig, axes = plt.subplots(n, 1, figsize=(10, 2.6 * n + 0.5), sharex=False)
        if n == 1: axes = [axes]
        for ax, (primary, label, xrange, rows) in zip(axes, panel_data):
            for src, counts in rows:
                color = color_map[src]
                fill_alpha = 0.20 if src == "geant4" else 0
                lk = {"linestyle": style_map[src], "linewidth": 1.2}
                plot_spectrum(ax, EDGES_KEV, counts,
                              label=f"{label_map[src]} (Σ={counts.sum():.4f})",
                              color=color, fill_alpha=fill_alpha, line_kwargs=lk)
            ax.set_xlim(*xrange)
            ax.set_yscale("log")
            ymax = max(c.max() for _, c in rows) * 2
            ax.set_ylim(1e-5, max(ymax, 1))
            style_ax(ax, label, ylog=True)
            ax.legend(loc="upper right", fontsize=8)
        fig.suptitle(f"Cross-validation: {' vs '.join(label_map[s] for s in active_sources)} — same nuclide, independent ENSDF evaluations",
                     y=1.00)
        save(fig, "11_tri_source_cross_validation.png")

        # Residual scatter plot vs Geant4 baseline
        if "geant4" in active_sources and len(active_sources) >= 2:
            fig, ax = plt.subplots(figsize=(10, 4))
            for primary, label, xrange, rows in panel_data:
                gd = next((c for s, c in rows if s == "geant4"), None)
                if gd is None or gd.sum() == 0: continue
                centers = 0.5 * (EDGES_KEV[:-1] + EDGES_KEV[1:])
                mask = gd > 0.005
                if mask.sum() == 0: continue
                for src, counts in rows:
                    if src == "geant4": continue
                    frac = (counts - gd) / np.maximum(gd, 1e-12)
                    marker = "s" if src == "sandia" else "^"
                    ax.scatter(centers[mask], frac[mask] * 100, s=20,
                               marker=marker, alpha=0.7,
                               color=color_map[src],
                               label=f"{label} vs Geant4 ({label_map[src]})" if label == panel_data[0][1] else None)
            ax.axhline(0, color=C["muted"], lw=0.5)
            ax.set_xlim(0, 3000)
            ax.set_ylim(-5, 5)
            style_ax(ax, "Fractional differences vs Geant4 baseline (peaks > 0.5%)",
                     xlabel="Energy (keV)", ylabel="(other − Geant4) / Geant4  (%)")
            ax.legend(loc="upper right", fontsize=8)
            save(fig, "11b_cross_validation_residuals.png")
    else:
        print("  SKIPPED -- no panel data with multiple sources")
else:
    print("  SKIPPED -- need at least 2 sources, only have:", active_sources)
print(f"\nDone. Plots written to: {os.path.abspath(OUTDIR)}")
print(f"Files:")
for f in sorted(os.listdir(OUTDIR)):
    full = os.path.join(OUTDIR, f)
    size = os.path.getsize(full)
    print(f"  {f}  ({size//1024} KB)")
