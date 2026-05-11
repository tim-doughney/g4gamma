"""
Generate all production figures for the g4gamma technical report.
Run from repo root:  PYTHONPATH=build python3 report/make_figures.py
"""
import sys, os, csv, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "build"))
import g4gamma as g

OUTDIR = pathlib.Path(__file__).parent / "figures"
OUTDIR.mkdir(exist_ok=True)
SIM_DIR  = pathlib.Path(__file__).parent.parent / "buildG4RadDecayExample"
LARA_DIR = str(pathlib.Path(__file__).parent.parent / "data" / "lara")
SANDIA_XML = str(pathlib.Path(__file__).parent.parent / "third_party" / "SandiaDecay" / "sandia.decay.xml")
SANDIA_CSV = str(pathlib.Path(__file__).parent.parent / "build" / "sandia_compare_out.csv")

# ─── matplotlib style ─────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":      "serif",
    "font.size":        10,
    "axes.titlesize":   10,
    "axes.labelsize":   10,
    "legend.fontsize":  8,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "figure.dpi":       150,
    "savefig.dpi":      200,
    "savefig.bbox":     "tight",
    "lines.linewidth":  1.2,
})

COLORS = {
    "g4_af":    "#1f77b4",   # blue
    "sandia":   "#d62728",   # red
    "lara":     "#2ca02c",   # green
    "g4_mc":    "#7f7f7f",   # gray
    "cutoff":   "#ff7f0e",   # orange
    "g4_no":    "#9467bd",   # purple
}

# ─── helpers ──────────────────────────────────────────────────────────────────
def load_rdecay_csv(path, n_primaries):
    """Return (bin_centers_keV, counts_per_primary_per_bin) from rdecay01 CSV."""
    entries_col = []
    with open(path) as f:
        lines = f.readlines()
    # 6 comment lines + 1 header line + 1 underflow row  => data at index 8
    data_start = 8
    data_end   = len(lines) - 1  # skip overflow
    for line in lines[data_start:data_end]:
        parts = line.strip().split(",")
        entries_col.append(float(parts[0]))
    arr = np.array(entries_col)
    counts = arr / n_primaries
    # reconstruct bin centers from "axis fixed N Emin Emax" in header
    for line in lines:
        if "#axis fixed" in line:
            parts = line.split()
            n_bins = int(parts[2])
            e_min  = float(parts[3])
            e_max  = float(parts[4])
            break
    centers = np.linspace(e_min, e_max, n_bins, endpoint=False) + (e_max - e_min) / (2 * n_bins)
    return centers, counts


def make_builder(source, **kw):
    opts = g.SpectrumOptions()
    opts.source = source
    if source == g.DataSource.Geant4:
        opts.full_xray_cascade  = kw.get("full_xray_cascade", True)
        opts.include_xrays      = kw.get("include_xrays", False)
        opts.include_annihilation = kw.get("include_annihilation", True)
    elif source == g.DataSource.Sandia:
        opts.sandia_xml = SANDIA_XML
        opts.include_xrays = kw.get("include_xrays", False)
        opts.include_annihilation = kw.get("include_annihilation", True)
    elif source == g.DataSource.Lara:
        opts.lara_dir = LARA_DIR
        opts.include_annihilation = kw.get("include_annihilation", True)
    if "chain_cutoffs" in kw:
        opts.chain_cutoffs = kw["chain_cutoffs"]
    return g.GammaSpectrumBuilder(opts)


def keV_edges(n=3000, emax=3000):
    return np.linspace(0, emax, n + 1) * g.units.keV


def spectrum_to_keV(res):
    """Return (centers_keV, counts_per_primary_per_bin)."""
    edges  = np.array(res.bin_edges) / g.units.keV
    counts = np.array(res.counts)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, counts


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Figure 1 – U-238 SE: g4gamma (Geant4 backend, AF=true) vs Geant4 MC
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_u238_g4_vs_mc():
    print("Fig 1: U-238 SE g4gamma vs MC")
    mc_path = SIM_DIR / "u238_AFtrue_1e7_h1_3.csv"
    mc_cen, mc_cnt = load_rdecay_csv(mc_path, 1e7)

    edges = keV_edges()
    b = make_builder(g.DataSource.Geant4)
    res = b.build(g.IsotopeKey(92, 238, 0), -1, edges)
    cen, cnt = spectrum_to_keV(res)

    fig, axes = plt.subplots(2, 1, figsize=(7, 5),
                             gridspec_kw={"height_ratios": [3, 1]},
                             sharex=True)
    ax, axr = axes

    ax.step(mc_cen, mc_cnt, where="mid", color=COLORS["g4_mc"],
            alpha=0.6, linewidth=0.8, label="Geant4 MC ($10^7$ events)")
    ax.step(cen, cnt, where="mid", color=COLORS["g4_af"],
            label=r"g4gamma (Geant4 backend, full $X$-ray cascade)")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 0.5)
    ax.set_ylabel(r"$\gamma$ / primary / keV")
    ax.set_title(r"$^{238}$U secular equilibrium chain – Geant4 data backend")
    ax.legend(loc="upper right")
    ax.text(0.02, 0.97, fr"Model total: {cnt.sum():.3f} $\gamma$/primary"
            + "\n" + fr"MC total: {mc_cnt.sum():.3f} $\gamma$/primary",
            transform=ax.transAxes, va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

    # ratio on bins where MC > small threshold
    mask = mc_cnt > 1e-6
    ratio = np.where(mask, cnt / mc_cnt, np.nan)
    axr.step(cen, ratio, where="mid", color=COLORS["g4_af"])
    axr.axhline(1.0, color="k", linewidth=0.8, linestyle="--")
    axr.axhspan(0.95, 1.05, color="gray", alpha=0.15)
    axr.set_ylim(0.5, 1.8)
    axr.set_ylabel("Model / MC")
    axr.set_xlabel("Energy (keV)")

    plt.tight_layout()
    plt.savefig(OUTDIR / "fig_u238_g4_vs_mc.pdf")
    plt.close()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Figure 2 – Th-232 SE: g4gamma vs Geant4 MC
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_th232_g4_vs_mc():
    print("Fig 2: Th-232 SE g4gamma vs MC")
    mc_path = SIM_DIR / "th232_AFtrue_1e7_h1_3.csv"
    mc_cen, mc_cnt = load_rdecay_csv(mc_path, 1e7)

    edges = keV_edges()
    b = make_builder(g.DataSource.Geant4)
    res = b.build(g.IsotopeKey(90, 232, 0), -1, edges)
    cen, cnt = spectrum_to_keV(res)

    fig, axes = plt.subplots(2, 1, figsize=(7, 5),
                             gridspec_kw={"height_ratios": [3, 1]},
                             sharex=True)
    ax, axr = axes

    ax.step(mc_cen, mc_cnt, where="mid", color=COLORS["g4_mc"],
            alpha=0.6, linewidth=0.8, label="Geant4 MC ($10^7$ events)")
    ax.step(cen, cnt, where="mid", color=COLORS["g4_af"],
            label=r"g4gamma (Geant4 backend, full $X$-ray cascade)")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 0.5)
    ax.set_ylabel(r"$\gamma$ / primary / keV")
    ax.set_title(r"$^{232}$Th secular equilibrium chain – Geant4 data backend")
    ax.legend(loc="upper right")
    ax.text(0.02, 0.97, fr"Model total: {cnt.sum():.3f} $\gamma$/primary"
            + "\n" + fr"MC total: {mc_cnt.sum():.3f} $\gamma$/primary",
            transform=ax.transAxes, va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

    mask = mc_cnt > 1e-6
    ratio = np.where(mask, cnt / mc_cnt, np.nan)
    axr.step(cen, ratio, where="mid", color=COLORS["g4_af"])
    axr.axhline(1.0, color="k", linewidth=0.8, linestyle="--")
    axr.axhspan(0.95, 1.05, color="gray", alpha=0.15)
    axr.set_ylim(0.5, 1.8)
    axr.set_ylabel("Model / MC")
    axr.set_xlabel("Energy (keV)")

    plt.tight_layout()
    plt.savefig(OUTDIR / "fig_th232_g4_vs_mc.pdf")
    plt.close()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Figure 3 – U-238 SE: all three backends compared
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_u238_backends():
    print("Fig 3: U-238 all backends")
    edges = keV_edges()
    key = g.IsotopeKey(92, 238, 0)

    bg4  = make_builder(g.DataSource.Geant4)
    bsnd = make_builder(g.DataSource.Sandia)
    blr  = make_builder(g.DataSource.Lara)

    rg4  = bg4.build(key, -1, edges)
    rsnd = bsnd.build(key, -1, edges)
    rlr  = blr.build(key, -1, edges)

    c4,  cnt4  = spectrum_to_keV(rg4)
    cs,  cnts  = spectrum_to_keV(rsnd)
    cl,  cntl  = spectrum_to_keV(rlr)

    mc_path = SIM_DIR / "u238_AFtrue_1e7_h1_3.csv"
    mc_cen, mc_cnt = load_rdecay_csv(mc_path, 1e7)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.step(mc_cen, mc_cnt, where="mid", color=COLORS["g4_mc"],
            alpha=0.5, linewidth=0.8, label="Geant4 MC ($10^7$ events)")
    ax.step(c4,  cnt4,  where="mid", color=COLORS["g4_af"],  label=f"Geant4 backend  ({cnt4.sum():.3f})")
    ax.step(cs,  cnts,  where="mid", color=COLORS["sandia"],  label=f"Sandia backend  ({cnts.sum():.3f})")
    ax.step(cl,  cntl,  where="mid", color=COLORS["lara"],    label=f"LARA backend    ({cntl.sum():.3f})")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 0.5)
    ax.set_xlim(0, 3000)
    ax.set_ylabel(r"$\gamma$ / primary / keV")
    ax.set_xlabel("Energy (keV)")
    ax.set_title(r"$^{238}$U secular equilibrium – database comparison (legend shows $\gamma$/primary total)")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig_u238_backends.pdf")
    plt.close()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Figure 4 – Th-232 SE: all three backends compared
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_th232_backends():
    print("Fig 4: Th-232 all backends")
    edges = keV_edges()
    key = g.IsotopeKey(90, 232, 0)

    bg4  = make_builder(g.DataSource.Geant4)
    bsnd = make_builder(g.DataSource.Sandia)
    blr  = make_builder(g.DataSource.Lara)

    rg4  = bg4.build(key, -1, edges)
    rsnd = bsnd.build(key, -1, edges)
    rlr  = blr.build(key, -1, edges)

    c4,  cnt4  = spectrum_to_keV(rg4)
    cs,  cnts  = spectrum_to_keV(rsnd)
    cl,  cntl  = spectrum_to_keV(rlr)

    mc_path = SIM_DIR / "th232_AFtrue_1e7_h1_3.csv"
    mc_cen, mc_cnt = load_rdecay_csv(mc_path, 1e7)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.step(mc_cen, mc_cnt, where="mid", color=COLORS["g4_mc"],
            alpha=0.5, linewidth=0.8, label="Geant4 MC ($10^7$ events)")
    ax.step(c4,  cnt4,  where="mid", color=COLORS["g4_af"],  label=f"Geant4 backend  ({cnt4.sum():.3f})")
    ax.step(cs,  cnts,  where="mid", color=COLORS["sandia"],  label=f"Sandia backend  ({cnts.sum():.3f})")
    ax.step(cl,  cntl,  where="mid", color=COLORS["lara"],    label=f"LARA backend    ({cntl.sum():.3f})")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 0.5)
    ax.set_xlim(0, 3000)
    ax.set_ylabel(r"$\gamma$ / primary / keV")
    ax.set_xlabel("Energy (keV)")
    ax.set_title(r"$^{232}$Th secular equilibrium – database comparison (legend shows $\gamma$/primary total)")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig_th232_backends.pdf")
    plt.close()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Figure 5 – Co-60 and Cs-137 single isotope: all backends vs MC
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_single_isotopes():
    print("Fig 5: Co-60 and Cs-137 all backends vs MC")
    edges = keV_edges(n=3000, emax=3000)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    for ax, (sym, Z, A, mc_file, nprim) in zip(axes, [
        ("Co-60",  27, 60, "co60_AFtrue_h1_3.csv",  1e6),
        ("Cs-137", 55, 137, "cs137_AFtrue_h1_3.csv", 1e6),
    ]):
        key  = g.IsotopeKey(Z, A, 0)

        bg4  = make_builder(g.DataSource.Geant4)
        bsnd = make_builder(g.DataSource.Sandia)
        blr  = make_builder(g.DataSource.Lara) if (
            Z == 27 and os.path.exists(LARA_DIR + "/Co-60.lara.txt")) or (
            Z == 55 and os.path.exists(LARA_DIR + "/Cs-137.lara.txt")) else None

        rg4  = bg4.build(key, -1, edges)
        rsnd = bsnd.build(key, -1, edges)

        c4, cnt4   = spectrum_to_keV(rg4)
        cs, cnts   = spectrum_to_keV(rsnd)

        mc_path = SIM_DIR / mc_file
        mc_cen, mc_cnt = load_rdecay_csv(mc_path, nprim)

        ax.step(mc_cen, mc_cnt, where="mid", color=COLORS["g4_mc"],
                alpha=0.5, linewidth=0.8, label=f"Geant4 MC ($10^6$ ev., {mc_cnt.sum():.3f})")
        ax.step(c4,  cnt4, where="mid", color=COLORS["g4_af"],
                label=f"Geant4 ({cnt4.sum():.3f})")
        ax.step(cs,  cnts, where="mid", color=COLORS["sandia"],
                label=f"Sandia ({cnts.sum():.3f})")

        if blr is not None:
            rlr = blr.build(key, -1, edges)
            cl, cntl = spectrum_to_keV(rlr)
            ax.step(cl, cntl, where="mid", color=COLORS["lara"],
                    label=f"LARA ({cntl.sum():.3f})")

        ax.set_yscale("log")
        ax.set_ylim(1e-5, 2.0)
        ax.set_xlim(0, 3000)
        ax.set_ylabel(r"$\gamma$ / primary / keV")
        ax.set_xlabel("Energy (keV)")
        ax.set_title(fr"$^{{{A}}}${sym.split('-')[0]} secular equilibrium"
                     "\n(legend shows γ/primary total)")
        ax.legend(fontsize=7.5)

    plt.tight_layout()
    plt.savefig(OUTDIR / "fig_single_isotopes.pdf")
    plt.close()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Figure 6 – Time-dependent: Co-60 total gamma yield vs time, all backends
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_timedep():
    print("Fig 6: time-dependent Co-60")
    edges = keV_edges(n=3000, emax=3000)
    key   = g.IsotopeKey(27, 60, 0)

    bg4  = make_builder(g.DataSource.Geant4)
    bsnd = make_builder(g.DataSource.Sandia)
    blr  = make_builder(g.DataSource.Lara)

    half_life_s = 5.2713 * 365.25 * 24 * 3600   # Co-60 T1/2 in seconds
    t_ns_values = np.linspace(0, 20 * 365.25 * 24 * 3600,
                               80) * g.units.second   # 0..20 years in ns

    def total_counts(builder, t_ns):
        res = builder.build(key, t_ns, edges)
        return np.array(res.counts).sum()

    totals_g4  = [total_counts(bg4,  t) for t in t_ns_values]
    totals_snd = [total_counts(bsnd, t) for t in t_ns_values]
    totals_lr  = [total_counts(blr,  t) for t in t_ns_values]

    t_years = t_ns_values / g.units.year

    # Normalise to t=0 value (Geant4 backend)
    norm = totals_g4[0] if totals_g4[0] > 0 else 1.0

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t_years, np.array(totals_g4)  / norm, color=COLORS["g4_af"],
            label="Geant4 backend")
    ax.plot(t_years, np.array(totals_snd) / norm, color=COLORS["sandia"],
            linestyle="--", label="Sandia backend")
    ax.plot(t_years, np.array(totals_lr)  / norm, color=COLORS["lara"],
            linestyle=":", label="LARA backend")

    # Expected exponential decay
    t_half_y = 5.2713
    ax.plot(t_years, np.exp(-np.log(2) * t_years / t_half_y),
            color="k", linestyle="-.", linewidth=0.8,
            label=fr"Analytic $e^{{-t \ln 2 / T_{{1/2}}}}$")

    ax.set_xlabel("Time (years)")
    ax.set_ylabel(r"Normalised $\gamma$ yield")
    ax.set_title(r"$^{60}$Co total $\gamma$ yield vs time – database comparison")
    ax.legend()
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig_timedep_co60.pdf")
    plt.close()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Figure 7 – Time-dependent spectra: Co-60 at 0, 5, 10, 20 years
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_timedep_spectra():
    print("Fig 7: time-dependent Co-60 spectra")
    edges = keV_edges(n=3000, emax=3000)
    key   = g.IsotopeKey(27, 60, 0)
    bg4   = make_builder(g.DataSource.Geant4)

    times_y = [0, 5, 10, 20]
    cols    = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

    fig, ax = plt.subplots(figsize=(7, 4))
    for yr, col in zip(times_y, cols):
        t_ns = yr * g.units.year
        res  = bg4.build(key, t_ns, edges)
        cen, cnt = spectrum_to_keV(res)
        ax.step(cen, cnt, where="mid", color=col,
                label=fr"$t = {yr}$ yr  (total {cnt.sum():.3f})")

    ax.set_yscale("log")
    ax.set_ylim(1e-6, 2.0)
    ax.set_xlim(0, 2000)
    ax.set_ylabel(r"$\gamma$ / primary / keV")
    ax.set_xlabel("Energy (keV)")
    ax.set_title(r"$^{60}$Co gamma spectrum evolution over time (Geant4 backend)")
    ax.legend()
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig_timedep_spectra.pdf")
    plt.close()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Figure 8 – Rn-222 chain cutoff effect on U-238 spectrum
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_rn222_cutoff():
    print("Fig 8: Rn-222 cutoff")
    edges = keV_edges()
    key   = g.IsotopeKey(92, 238, 0)

    b_full = make_builder(g.DataSource.Geant4)
    b_cut  = make_builder(g.DataSource.Geant4,
                          chain_cutoffs=[g.IsotopeKey(86, 222, 0)])

    r_full = b_full.build(key, -1, edges)
    r_cut  = b_cut.build(key,  -1, edges)

    cf, cntf = spectrum_to_keV(r_full)
    cc, cntc = spectrum_to_keV(r_cut)

    # find which isotopes were excluded
    cut_set = {str(c.isotope) for c in r_cut.contributions if c.cutoff}
    full_set = {str(c.isotope) for c in r_full.contributions}
    excluded = full_set - {str(c.isotope) for c in r_cut.contributions}

    fig, axes = plt.subplots(2, 1, figsize=(7, 5),
                             gridspec_kw={"height_ratios": [3, 1]},
                             sharex=True)
    ax, axd = axes

    ax.step(cf, cntf, where="mid", color=COLORS["g4_af"],
            label=f"Full chain ({cntf.sum():.3f} γ/primary)")
    ax.step(cc, cntc, where="mid", color=COLORS["cutoff"], linestyle="--",
            label=f"Rn-222 cutoff ({cntc.sum():.3f} γ/primary)")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 0.5)
    ax.set_ylabel(r"$\gamma$ / primary / keV")
    ax.set_title(r"$^{238}$U secular equilibrium – effect of $^{222}$Rn chain cutoff (radon escape)")
    ax.legend()

    diff = cntf - cntc
    axd.step(cf, diff, where="mid", color="purple")
    axd.axhline(0, color="k", linewidth=0.8)
    axd.set_ylabel(r"$\Delta\gamma$ (removed)")
    axd.set_xlabel("Energy (keV)")
    axd.set_ylim(-0.001, max(diff.max() * 1.2, 0.001))

    plt.tight_layout()
    plt.savefig(OUTDIR / "fig_rn222_cutoff.pdf")
    plt.close()

    # Print chain members removed
    print(f"  Full chain members: {len(r_full.contributions)}, "
          f"Cut chain members: {len(r_cut.contributions)}")
    print(f"  Removed isotopes: {excluded}")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Figure 9 – Auger/X-ray cascade settings comparison (U-238 low-energy region)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_auger_xray():
    print("Fig 9: Auger/X-ray settings")
    edges = keV_edges()
    key   = g.IsotopeKey(92, 238, 0)

    b_none = make_builder(g.DataSource.Geant4,
                          full_xray_cascade=False, include_xrays=False,
                          include_annihilation=True)
    b_ksh  = make_builder(g.DataSource.Geant4,
                          full_xray_cascade=False, include_xrays=True,
                          include_annihilation=True)
    b_full = make_builder(g.DataSource.Geant4,
                          full_xray_cascade=True, include_xrays=False,
                          include_annihilation=True)

    r_none = b_none.build(key, -1, edges)
    r_ksh  = b_ksh.build(key,  -1, edges)
    r_full = b_full.build(key, -1, edges)

    cn, tn = spectrum_to_keV(r_none)
    ck, tk = spectrum_to_keV(r_ksh)
    cf, tf = spectrum_to_keV(r_full)

    mc_path = SIM_DIR / "u238_AFtrue_1e7_h1_3.csv"
    mc_cen, mc_cnt = load_rdecay_csv(mc_path, 1e7)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Left: full range log
    ax = axes[0]
    ax.step(mc_cen, mc_cnt, where="mid", color=COLORS["g4_mc"],
            alpha=0.5, linewidth=0.8, label=f"Geant4 MC ({mc_cnt.sum():.3f})")
    ax.step(cn, tn, where="mid", color="#9467bd", linestyle=":",
            label=f"No X-rays ({tn.sum():.3f})")
    ax.step(ck, tk, where="mid", color="#e377c2", linestyle="--",
            label=f"K-shell only ({tk.sum():.3f})")
    ax.step(cf, tf, where="mid", color=COLORS["g4_af"],
            label=f"Full cascade ({tf.sum():.3f})")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 0.5)
    ax.set_xlim(0, 3000)
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel(r"$\gamma$ / primary / keV")
    ax.set_title(r"$^{238}$U SE – X-ray/Auger settings (full range)")
    ax.legend(fontsize=7.5)

    # Right: low energy zoom
    ax2 = axes[1]
    mask_lo = mc_cen < 150
    ax2.step(mc_cen[mask_lo], mc_cnt[mask_lo], where="mid",
             color=COLORS["g4_mc"], alpha=0.5, linewidth=0.8, label="Geant4 MC")
    mask_lo2 = cn < 150
    ax2.step(cn[mask_lo2], tn[mask_lo2], where="mid", color="#9467bd",
             linestyle=":", label="No X-rays")
    ax2.step(ck[mask_lo2], tk[mask_lo2], where="mid", color="#e377c2",
             linestyle="--", label="K-shell only")
    ax2.step(cf[mask_lo2], tf[mask_lo2], where="mid",
             color=COLORS["g4_af"], label="Full cascade")
    ax2.set_yscale("log")
    ax2.set_ylim(1e-5, 0.5)
    ax2.set_xlim(0, 150)
    ax2.set_xlabel("Energy (keV)")
    ax2.set_ylabel(r"$\gamma$ / primary / keV")
    ax2.set_title(r"$^{238}$U SE – low-energy region (X-ray / Auger lines)")
    ax2.legend(fontsize=7.5)

    plt.tight_layout()
    plt.savefig(OUTDIR / "fig_auger_xray.pdf")
    plt.close()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Figure 10 – Bateman solver validation vs SandiaDecay
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_bateman_validation():
    print("Fig 10: Bateman validation")
    if not os.path.exists(SANDIA_CSV):
        print("  sandia_compare_out.csv not found – skipping")
        return

    # Parse sandia_compare CSV
    # Columns: isotope, chain_member, time_s, g4gamma_activity, sandia_activity, rel_err_pct, pass
    data = {}
    with open(SANDIA_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            iso = row["isotope"]
            dau = row["chain_member"]
            key = (iso, dau)
            if key not in data:
                data[key] = {"t": [], "g4": [], "sd": []}
            data[key]["t"].append(float(row["time_s"]))
            data[key]["g4"].append(float(row["g4gamma_activity"]))
            data[key]["sd"].append(float(row["sandia_activity"]))

    # Pick interesting chains: Mo-99/Tc-99m and Co-60/Ni-60
    chains_of_interest = {}
    for (iso, dau), d in data.items():
        if iso in ("Mo-99",) and dau in ("Tc-99m", "Tc-99"):
            chains_of_interest[f"{iso}→{dau}"] = d
        if iso == "Co-60":
            chains_of_interest[f"{iso}→{dau}"] = d

    if not chains_of_interest:
        # Fallback: pick first 4 entries
        for k in list(data.keys())[:4]:
            chains_of_interest[f"{k[0]}→{k[1]}"] = data[k]

    fig, axes = plt.subplots(1, min(len(chains_of_interest), 4),
                             figsize=(12, 4))
    if len(chains_of_interest) == 1:
        axes = [axes]

    for ax, (label, d) in zip(axes, list(chains_of_interest.items())[:4]):
        t_h = np.array(d["t"]) / 3600
        g4  = np.array(d["g4"])
        sd  = np.array(d["sd"])
        norm = sd[0] if sd[0] > 0 else 1.0
        ax.semilogy(t_h, g4 / norm, "o-", color=COLORS["g4_af"], markersize=4,
                    label="g4gamma")
        ax.semilogy(t_h, sd / norm, "s--", color=COLORS["sandia"], markersize=4,
                    label="SandiaDecay")
        ax.set_xlabel("Time (h)")
        ax.set_ylabel("Normalised activity")
        ax.set_title(label)
        ax.legend(fontsize=7.5)

    plt.tight_layout()
    plt.savefig(OUTDIR / "fig_bateman_validation.pdf")
    plt.close()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Figure 11 – Chain contributions pie / bar chart for U-238 SE
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def fig_chain_contributions():
    print("Fig 11: chain contributions")
    edges = keV_edges()
    key   = g.IsotopeKey(92, 238, 0)
    b     = make_builder(g.DataSource.Geant4)
    res   = b.build(key, -1, edges)

    names  = []
    yields = []
    for c in res.contributions:
        if c.gamma_yield > 0.005:
            names.append(str(c.isotope))
            yields.append(c.gamma_yield * c.activity)

    idx  = np.argsort(yields)[::-1][:12]
    nms  = [names[i] for i in idx]
    ylds = [yields[i] for i in idx]

    fig, ax = plt.subplots(figsize=(8, 4))
    bars = ax.barh(range(len(nms)), ylds, color=COLORS["g4_af"], alpha=0.8)
    ax.set_yticks(range(len(nms)))
    ax.set_yticklabels(nms)
    ax.invert_yaxis()
    ax.set_xlabel(r"Contribution to total $\gamma$ yield (activity × yield)")
    ax.set_title(r"$^{238}$U SE – top contributing chain members (Geant4 backend)")
    for i, (b2, v) in enumerate(zip(bars, ylds)):
        ax.text(v + 0.002, i, f"{v:.3f}", va="center", fontsize=8)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig_chain_contributions.pdf")
    plt.close()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Run all
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
if __name__ == "__main__":
    fig_u238_g4_vs_mc()
    fig_th232_g4_vs_mc()
    fig_u238_backends()
    fig_th232_backends()
    fig_single_isotopes()
    fig_timedep()
    fig_timedep_spectra()
    fig_rn222_cutoff()
    fig_auger_xray()
    fig_bateman_validation()
    fig_chain_contributions()
    print("\nAll figures written to", OUTDIR)
    for f in sorted(OUTDIR.glob("*.pdf")):
        print(" ", f.name)
