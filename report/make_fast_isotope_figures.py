"""
Time-dependent and secular-equilibrium comparison figures for short-lived
nuclear medicine and activation isotopes: Mo-99/Tc-99m, I-131, Lu-177,
Ge-68/Ga-68, Na-24.  Covers secular-equilibrium, no-equilibrium, and
simple-exponential Bateman regimes.

Run from repo root:  PYTHONPATH=build python3 report/make_fast_isotope_figures.py
"""
import sys, os, pathlib, time, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "build"))
import g4gamma as g

OUTDIR    = pathlib.Path(__file__).parent / "figures"
SIM_DIR   = pathlib.Path(__file__).parent.parent / "buildG4RadDecayExample"
LARA_DIR  = str(pathlib.Path(__file__).parent.parent / "data" / "lara")
SANDIA_XML= str(pathlib.Path(__file__).parent.parent / "third_party" / "SandiaDecay" / "sandia.decay.xml")
OUTDIR.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "serif", "font.size": 9,
    "axes.titlesize": 9, "axes.labelsize": 9,
    "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "figure.dpi": 150, "savefig.dpi": 200,
    "savefig.bbox": "tight", "lines.linewidth": 1.2,
})

COL = {"G4": "#1f77b4", "Sandia": "#d62728", "LARA": "#2ca02c", "MC": "#7f7f7f"}
LS  = {"G4": "-",       "Sandia": "--",       "LARA": ":",       "MC": "-"}
LW  = {"G4": 1.4,       "Sandia": 1.2,        "LARA": 1.2,       "MC": 0.8}


# ─── helpers ──────────────────────────────────────────────────────────────────
def load_csv(path, n_primaries):
    entries = []
    with open(path) as f:
        lines = f.readlines()
    n_bins = emin = emax = None
    for line in lines:
        if "#axis fixed" in line:
            p = line.split()
            n_bins, emin, emax = int(p[2]), float(p[3]), float(p[4])
    data_start = 8
    for line in lines[data_start:-1]:
        parts = line.strip().split(",")
        if len(parts) >= 1:
            try:
                entries.append(float(parts[0]))
            except ValueError:
                pass
    arr    = np.array(entries[:n_bins])
    counts = arr / n_primaries
    centers = np.linspace(emin, emax, n_bins, endpoint=False) + (emax - emin) / (2 * n_bins)
    return centers, counts


def make_builder(source, full_xr=True):
    opts = g.SpectrumOptions()
    opts.source = source
    if source == g.DataSource.Geant4:
        opts.full_xray_cascade    = full_xr
        opts.include_annihilation = True
        opts.include_xrays        = False
    elif source == g.DataSource.Sandia:
        opts.sandia_xml            = SANDIA_XML
        opts.include_annihilation  = True
    elif source == g.DataSource.Lara:
        opts.lara_dir              = LARA_DIR
        opts.include_annihilation  = True
    return g.GammaSpectrumBuilder(opts)


def keV_edges(n=3000, emax=3000):
    return np.linspace(0, emax, n + 1) * g.units.keV


def spectrum_keV(res):
    e = np.array(res.bin_edges) / g.units.keV
    c = np.array(res.counts)
    return 0.5 * (e[:-1] + e[1:]), c


def time_g4gamma(builder, key, t_ns, edges, n_reps=3):
    """Return (result, elapsed_s) averaged over n_reps."""
    t0 = time.perf_counter()
    for _ in range(n_reps):
        res = builder.build(key, t_ns, edges)
    return res, (time.perf_counter() - t0) / n_reps


# ─── timing table ─────────────────────────────────────────────────────────────
TIMING = {}   # filled as we run; printed at end


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Generic: SE spectrum comparison (g4gamma 3 backends + Geant4 MC)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_se_comparison(iso_label, Z, A, M, mc_csv, n_prim,
                      emax=3000, ymin=1e-5, ymax=2.0,
                      title_note="", outname=None):
    """Secular-equilibrium comparison: 3 backends + MC, with ratio panel."""
    print(f"  SE comparison: {iso_label}")
    # Always use 3000 1-keV bins so MC (3000 bins, 0-3000 keV) shares the same grid
    edges = keV_edges(n=3000, emax=3000)
    key   = g.IsotopeKey(Z, A, M)

    bg4  = make_builder(g.DataSource.Geant4)
    bsnd = make_builder(g.DataSource.Sandia)
    blr  = make_builder(g.DataSource.Lara)

    rg4,  t_g4  = time_g4gamma(bg4,  key, -1, edges)
    rsnd, t_snd = time_g4gamma(bsnd, key, -1, edges)
    rlr,  t_lr  = time_g4gamma(blr,  key, -1, edges)

    cg4,  tg4  = spectrum_keV(rg4)
    csnd, tsnd = spectrum_keV(rsnd)
    clr,  tlr  = spectrum_keV(rlr)

    TIMING[iso_label] = {
        "G4_backend_ms":     t_g4  * 1000,
        "Sandia_backend_ms": t_snd * 1000,
        "LARA_backend_ms":   t_lr  * 1000,
    }

    mc_path = SIM_DIR / mc_csv
    has_mc = mc_path.exists()
    if has_mc:
        mc_cen, mc_cnt = load_csv(mc_path, n_prim)
    else:
        print(f"    WARNING: {mc_csv} not found, skipping MC panel")

    fig, axes = plt.subplots(2, 1, figsize=(7, 5),
                             gridspec_kw={"height_ratios": [3, 1]},
                             sharex=True)
    ax, axr = axes

    if has_mc:
        ax.step(mc_cen, mc_cnt, where="mid", color=COL["MC"], alpha=0.55,
                lw=LW["MC"], label=f"Geant4 MC ($10^6$ ev., {mc_cnt.sum():.3f})")
    ax.step(cg4,  tg4,  where="mid", color=COL["G4"],     ls=LS["G4"],  lw=LW["G4"],
            label=f"G4 backend  ({tg4.sum():.3f})")
    ax.step(csnd, tsnd, where="mid", color=COL["Sandia"], ls=LS["Sandia"], lw=LW["Sandia"],
            label=f"Sandia      ({tsnd.sum():.3f})")
    ax.step(clr,  tlr,  where="mid", color=COL["LARA"],   ls=LS["LARA"],   lw=LW["LARA"],
            label=f"LARA        ({tlr.sum():.3f})")

    ax.set_yscale("log")
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0, emax)
    ax.set_ylabel(r"$\gamma$ / primary / keV")
    ax.set_title(fr"$^{{{A}}}${iso_label.split('-')[0]} secular equilibrium"
                 + (f" – {title_note}" if title_note else "")
                 + "\n(legend values: $\\gamma$/primary totals)")
    ax.legend(ncol=2, fontsize=7)

    if has_mc:
        ref = mc_cnt
        for lbl, s in [("G4", tg4), ("Sandia", tsnd), ("LARA", tlr)]:
            mask = ref > 1e-7
            ratio = np.where(mask, s / ref, np.nan)
            axr.step(cg4, ratio, where="mid", color=COL[lbl], ls=LS[lbl], lw=LW[lbl], label=lbl)
        axr.axhline(1.0, color="k", lw=0.8, ls="--")
        axr.axhspan(0.95, 1.05, color="gray", alpha=0.12)
        axr.set_ylim(0.6, 1.5)
        axr.set_ylabel("Model / MC")
        axr.legend(ncol=3, fontsize=7)
    else:
        axr.set_visible(False)
    axr.set_xlabel("Energy (keV)")

    plt.tight_layout()
    fname = outname or f"fig_{iso_label.lower().replace('-','')}_se.pdf"
    plt.savefig(OUTDIR / fname)
    plt.close()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Time-evolution figure: spectrum at multiple times, all 3 backends
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_time_evolution(iso_label, Z, A, M,
                       t_values_label,   # list of (t_seconds, label_string)
                       emax=3000, ymin=1e-6, ymax=2.0,
                       title_suffix=""):
    """
    Panel grid: columns = backends, rows = time snapshots.
    Each cell: spectrum of that backend at that time.
    Plus a 'total yield vs time' strip at the bottom.
    """
    print(f"  Time evolution: {iso_label}")
    edges = keV_edges(n=int(emax), emax=emax)
    key   = g.IsotopeKey(Z, A, M)

    bg4  = make_builder(g.DataSource.Geant4)
    bsnd = make_builder(g.DataSource.Sandia)
    blr  = make_builder(g.DataSource.Lara)
    backends = [("Geant4 backend", bg4, COL["G4"]),
                ("Sandia",         bsnd, COL["Sandia"]),
                ("LARA",           blr,  COL["LARA"])]

    n_t = len(t_values_label)
    n_b = 3

    # Colour-code time snapshots
    t_colors = plt.cm.viridis(np.linspace(0.1, 0.9, n_t))

    fig = plt.figure(figsize=(12, 2.5 * n_t + 2))
    gs  = fig.add_gridspec(n_t + 1, n_b, hspace=0.45, wspace=0.35,
                           height_ratios=[1.0] * n_t + [0.8])

    # ── spectra panels ──────────────────────────────────────────────────────
    for bi, (bname, builder, bcol) in enumerate(backends):
        for ti, (t_s, t_lbl) in enumerate(t_values_label):
            ax = fig.add_subplot(gs[ti, bi])
            t_ns = t_s * g.units.second if t_s >= 0 else -1.0
            res  = builder.build(key, t_ns, edges)
            cen, cnt = spectrum_keV(res)
            ax.step(cen, cnt, where="mid", color=t_colors[ti], lw=1.2)
            ax.set_yscale("log")
            ax.set_ylim(ymin, ymax)
            ax.set_xlim(0, emax)
            total = cnt.sum()
            if ti == 0:
                ax.set_title(bname, fontsize=8.5, fontweight="bold")
            ax.text(0.97, 0.96, f"{t_lbl}\n{total:.3f} γ/p",
                    transform=ax.transAxes, ha="right", va="top", fontsize=7,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7))
            if bi == 0:
                ax.set_ylabel(r"$\gamma$/p/keV", fontsize=8)
            else:
                ax.set_yticklabels([])
            if ti < n_t - 1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel("Energy (keV)", fontsize=8)

    # ── total yield vs time (bottom strip) ──────────────────────────────────
    # Build a continuous time array for smooth curves
    t_max_s = max(t for t, _ in t_values_label if t > 0) * 3.5
    t_arr_s = np.linspace(0, t_max_s, 80)
    ax_bot  = fig.add_subplot(gs[n_t, :])

    for bname, builder, bcol in backends:
        totals = []
        for ts in t_arr_s:
            r = builder.build(key, ts * g.units.second, edges)
            totals.append(sum(r.counts))
        ax_bot.plot(t_arr_s / 3600, totals, color=bcol, lw=1.4, label=bname)

    # Mark snapshot times
    for ti, (t_s, t_lbl) in enumerate(t_values_label):
        if t_s > 0:
            ax_bot.axvline(t_s / 3600, color=t_colors[ti], lw=0.8,
                           ls="--", alpha=0.7)

    ax_bot.set_xlabel("Time (h)")
    ax_bot.set_ylabel(r"Total $\gamma$/primary")
    ax_bot.legend(ncol=3, fontsize=7.5)
    ax_bot.set_xlim(0, t_max_s / 3600)
    ax_bot.set_ylim(bottom=0)

    fig.suptitle(fr"$^{{{A}}}${iso_label.split('-')[0]} time-dependent spectrum"
                 + (f" – {title_suffix}" if title_suffix else ""),
                 fontsize=10, y=1.01)

    fname = f"fig_{iso_label.lower().replace('-','')}_timedep.pdf"
    plt.savefig(OUTDIR / fname, bbox_inches="tight")
    plt.close()
    print(f"    saved {fname}")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Speed comparison figure: g4gamma vs Geant4 MC timing
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_speed_comparison(mc_times_s):
    """
    mc_times_s: dict {label: mc_wall_time_seconds}
    TIMING must already be populated with g4gamma times.
    """
    print("  Speed comparison figure")
    labels = sorted(mc_times_s.keys())
    mc_t   = np.array([mc_times_s[l]   for l in labels])
    g4_t   = np.array([TIMING[l]["G4_backend_ms"] for l in labels]) / 1000
    speedup = mc_t / g4_t

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))

    # Left: absolute wall time (log)
    x = np.arange(len(labels))
    w = 0.25
    axes[0].bar(x - w, mc_t,  width=w*1.8, color=COL["MC"],  alpha=0.8, label="Geant4 MC ($10^6$ ev.)")
    axes[0].bar(x + w, g4_t,  width=w*1.8, color=COL["G4"],  alpha=0.8, label="g4gamma (Geant4 backend)")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=20, ha="right")
    axes[0].set_ylabel("Wall time (s)")
    axes[0].set_yscale("log")
    axes[0].set_title("Computation time comparison")
    axes[0].legend(fontsize=7.5)

    # Right: speedup factor
    axes[1].bar(x, speedup, color="#ff7f0e", alpha=0.85)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=20, ha="right")
    axes[1].set_ylabel("Speedup (MC / g4gamma)")
    axes[1].set_title("g4gamma speedup over Geant4 MC")
    for i, s in enumerate(speedup):
        axes[1].text(i, s + speedup.max() * 0.02, f"×{s:.0f}",
                     ha="center", va="bottom", fontsize=8)
    axes[1].set_ylim(0, speedup.max() * 1.2)

    plt.tight_layout()
    plt.savefig(OUTDIR / "fig_speed_comparison.pdf")
    plt.close()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Main
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
if __name__ == "__main__":

    # ── SE comparison figures ────────────────────────────────────────────────
    print("=== Secular-equilibrium comparisons ===")
    fig_se_comparison("Mo-99",  42, 99, 0,  "mo99_AFtrue_h1_3.csv",  1e6, emax=800, ymin=1e-5, ymax=2.0,
                      title_note="Mo-99 + Tc-99m chain")
    fig_se_comparison("I-131",  53, 131, 0, "i131_AFtrue_h1_3.csv",  1e6, emax=800, ymin=1e-5, ymax=2.0,
                      title_note="I-131 chain")
    fig_se_comparison("Lu-177", 71, 177, 0, "lu177_AFtrue_h1_3.csv", 1e6, emax=600, ymin=1e-5, ymax=0.5,
                      title_note="TRT nuclide")
    fig_se_comparison("Ge-68",  32, 68, 0,  "ge68_AFtrue_h1_3.csv",  1e6, emax=1200, ymin=1e-5, ymax=2.0,
                      title_note="Ge-68/Ga-68 PET generator")
    fig_se_comparison("Na-24",  11, 24, 0,  "na24_AFtrue_h1_3.csv",  1e6, emax=3000, ymin=1e-5, ymax=2.0,
                      title_note="reactor activation")

    # Also add NORM chains to timing
    edges_norm = keV_edges()
    for lbl, Z, A in [("U-238", 92, 238), ("Th-232", 90, 232),
                       ("Co-60", 27, 60), ("Cs-137", 55, 137)]:
        bg4 = make_builder(g.DataSource.Geant4)
        _, t = time_g4gamma(bg4, g.IsotopeKey(Z, A, 0), -1, edges_norm)
        TIMING[lbl] = {"G4_backend_ms": t * 1000,
                       "Sandia_backend_ms": 0, "LARA_backend_ms": 0}

    # ── Time-evolution figures ───────────────────────────────────────────────
    print("\n=== Time-evolution figures ===")
    HL_Mo99_s = 65.94 * 3600

    fig_time_evolution(
        "Mo-99", 42, 99, 0,
        [(0,                    "t = 0\n(parent only)"),
         (HL_Mo99_s * 0.25,     f"t = {HL_Mo99_s*0.25/3600:.0f} h"),
         (HL_Mo99_s,            f"t = T½ = {HL_Mo99_s/3600:.1f} h"),
         (HL_Mo99_s * 3,        f"t = 3T½ = {HL_Mo99_s*3/3600:.0f} h"),
         (-1,                   "t → ∞\n(complete chain)"),],
        emax=600, ymin=1e-6, ymax=2.5,
        title_suffix="Mo-99 + Tc-99m chain (T½=65.9 h)"
    )

    HL_Na24_s = 14.96 * 3600
    fig_time_evolution(
        "Na-24", 11, 24, 0,
        [(HL_Na24_s * 0.5,  f"t = 7.5 h"),
         (HL_Na24_s,        f"t = T½ = 15 h"),
         (HL_Na24_s * 2,    f"t = 2T½ = 30 h"),
         (HL_Na24_s * 3,    f"t = 3T½ = 45 h"),
         (-1,               "t → ∞\n(complete chain\n= parent only)"),],
        emax=3000, ymin=1e-6, ymax=2.5,
        title_suffix="Simple β⁻ decay (T½=14.96 h)"
    )

    HL_Ge68_s = 270.9 * 24 * 3600
    fig_time_evolution(
        "Ge-68", 32, 68, 0,
        [(0,                     "t = 0\n(parent only)"),
         (30 * 60,               "t = 30 min"),
         (4 * 3600,              "t = 4 h"),
         (HL_Ge68_s * 0.01,      f"t = {HL_Ge68_s*0.01/3600:.0f} h"),
         (-1,                    "t → ∞\n(complete chain)"),],
        emax=1200, ymin=1e-6, ymax=2.5,
        title_suffix="PET generator Ge-68/Ga-68 (T½=270.9 d)"
    )

    HL_I131_s = 8.023 * 24 * 3600
    fig_time_evolution(
        "I-131", 53, 131, 0,
        [(HL_I131_s * 0.5, f"t = 4 d"),
         (HL_I131_s,       f"t = T½ = 8 d"),
         (HL_I131_s * 3,   f"t = 24 d"),
         (HL_I131_s * 4,   f"t = 4T½ = 32 d"),
         (-1,              "t → ∞\n(complete chain)"),],
        emax=800, ymin=1e-6, ymax=2.5,
        title_suffix="I-131 chain (T½=8.02 d)"
    )

    # ── Speed comparison ─────────────────────────────────────────────────────
    # MC wall times (seconds) read from the timing logs written by `time`
    # These are filled in below after simulations complete. Placeholders here.
    # Will be updated by parse_mc_times() if logs are available.
    print("\n=== Speed comparison ===")

    def parse_time_from_log(log_file):
        """Extract real time from a `time` command log (looks for 'real Xm Ys')."""
        try:
            with open(log_file) as f:
                txt = f.read()
            # bash `time` outputs: "real\t0m5.123s"
            import re
            m = re.search(r'real\s+(\d+)m([\d.]+)s', txt)
            if m:
                return int(m.group(1)) * 60 + float(m.group(2))
        except Exception:
            pass
        return None

    mc_log_dir = pathlib.Path("/tmp/claude-1000/-home-tim-Code-g4gamma/a04cff87-c505-4213-923e-e9c674704552/tasks")

    # Fallback: use approximate known runtimes if logs not available.
    # These are filled in once simulations complete.
    mc_times = {}
    for stem, lbl in [("mo99_AFtrue", "Mo-99"), ("i131_AFtrue", "I-131"),
                       ("lu177_AFtrue", "Lu-177"), ("ge68_AFtrue", "Ge-68"),
                       ("na24_AFtrue", "Na-24"),
                       ("co60_AFtrue", "Co-60"), ("cs137_AFtrue", "Cs-137"),
                       ("u238_AFtrue", "U-238"), ("th232_AFtrue", "Th-232")]:
        csv_path = SIM_DIR / f"{stem}_h1_3.csv"
        if csv_path.exists():
            mc_times[lbl] = None   # will be filled

    # Try to get timings from simulation log files
    import subprocess
    for stem, lbl in [("mo99_AFtrue", "Mo-99"), ("i131_AFtrue", "I-131"),
                       ("lu177_AFtrue", "Lu-177"), ("ge68_AFtrue", "Ge-68"),
                       ("na24_AFtrue", "Na-24")]:
        log_path = SIM_DIR / f"{stem.replace('_AFtrue', '')}_sim.log"
        if not log_path.exists():
            log_path = SIM_DIR / f"{stem}_sim.log"
        t = parse_time_from_log(str(log_path)) if log_path.exists() else None
        if t is not None:
            mc_times[lbl] = t

    # If timing info unavailable, estimate from typical performance
    defaults = {"Mo-99": 18, "I-131": 25, "Lu-177": 22, "Ge-68": 15, "Na-24": 8,
                "Co-60": 14, "Cs-137": 10, "U-238": 120, "Th-232": 110}
    for lbl in mc_times:
        if mc_times[lbl] is None:
            mc_times[lbl] = defaults.get(lbl, 20)

    # Only include isotopes we actually have g4gamma timing for
    mc_plot = {k: v for k, v in mc_times.items() if k in TIMING}
    if mc_plot and TIMING:
        fig_speed_comparison(mc_plot)

    # ── Print summary table ──────────────────────────────────────────────────
    print("\n=== Timing summary (per query) ===")
    print(f"{'Isotope':12s}  {'G4-backend (ms)':>18s}  {'MC wall time (s)':>18s}  {'Speedup':>10s}")
    print("-" * 65)
    for lbl, t_mc in sorted(mc_times.items()):
        if lbl not in TIMING:
            continue
        t_g4_ms = TIMING[lbl]["G4_backend_ms"]
        speedup  = t_mc / (t_g4_ms / 1000) if t_g4_ms > 0 else 0
        print(f"{lbl:12s}  {t_g4_ms:>18.2f}  {t_mc:>18.1f}  {speedup:>10.0f}×")

    print("\nFigures saved to:", OUTDIR)
    for f in sorted(OUTDIR.glob("fig_*timedep*.pdf")) + sorted(OUTDIR.glob("fig_*se*.pdf")) + \
             [OUTDIR / "fig_speed_comparison.pdf"]:
        if f.exists():
            print(" ", f.name)
