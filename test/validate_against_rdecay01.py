#!/usr/bin/env python3
"""
validate_against_rdecay01.py

Compare g4gamma output against a CSV file produced by an rdecay01 simulation
(the user's reference / ground truth). For instance, the K40.csv produced by
the macro the user posted, with /analysis/h1/set 3 3000 0. 3000 keV.

Usage:
    python validate_against_rdecay01.py <isotope> <Z> <A> <M> <rdecay01.csv> <num_primaries>

The rdecay01 ROOT file gets converted to CSV by Geant4 (via /analysis/setDefaultFileType csv).
Histograms are stored as a 1D table; we read column 3 (gamma) which corresponds to /analysis/h1/set 3.

Example:
    # for K40.csv produced from the macro you posted, with 1e7 primaries:
    python validate_against_rdecay01.py K40 19 40 0 K40_h1_3.csv 10000000
"""
import sys
import os
import numpy as np
try:
    import g4gamma as g
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "build"))
    import g4gamma as g


def main():
    if len(sys.argv) < 7:
        print(__doc__)
        sys.exit(1)
    label   = sys.argv[1]
    Z       = int(sys.argv[2])
    A       = int(sys.argv[3])
    M       = int(sys.argv[4])
    csv_path = sys.argv[5]
    n_prim  = int(sys.argv[6])

    # rdecay01 h1.set 3 3000 0. 3000 keV gives 3000 bins from 0 to 3000 keV.
    # The CSV file has 1 row per bin: column 0 is bin lower edge, column 1 is count.
    # (Geant4 csv format depends on options; adjust as needed.)
    print(f"Reading rdecay01 reference: {csv_path}")
    data = np.loadtxt(csv_path, comments="#", delimiter=",")[1:-1,0]
    if data.ndim == 1:
        # single column of bin counts
        rd_counts = data
        nbins = len(rd_counts)
        rd_edges_keV = np.linspace(0, 3000, nbins + 1)
    else:
        rd_edges_keV = np.append(data[:, 0], data[-1, 0] + (data[1, 0] - data[0, 0]))
        rd_counts = data[:, 1]

    rd_pdf = rd_counts / n_prim  # gammas/primary/bin
    print(f"  {len(rd_pdf)} bins, total {rd_pdf.sum():.4f} gammas/primary "
          f"in [{rd_edges_keV[0]:.0f}, {rd_edges_keV[-1]:.0f}] keV")

    # g4gamma forward
    print(f"\nComputing g4gamma analytic spectrum for {label} (Z={Z} A={A} M={M})...")
    edges = rd_edges_keV * g.units.keV
    opts = g.SpectrumOptions()
    opts.include_annihilation = True
    opts.include_xrays = True  # rdecay01 with ARM=false doesn't emit X-rays
    builder = g.GammaSpectrumBuilder(opts)
    res = builder.build(g.IsotopeKey(Z, A, M), -1.0, edges)
    g4g_counts = np.array(res.counts)
    print(f"  total {g4g_counts.sum():.4f} gammas/primary")

    # Print non-zero analytic peaks
    print("\nAnalytic peaks (g4gamma):")
    for i in np.where(g4g_counts > 1e-6)[0]:
        print(f"  {rd_edges_keV[i]:.0f}-{rd_edges_keV[i+1]:.0f} keV: {g4g_counts[i]:.6f}")

    # Compare bin-by-bin where rdecay01 has signal
    print("\nBin-by-bin comparison (top 30 bins by rdecay01 count):")
    top = np.argsort(rd_pdf)[::-1][:30]
    print(f"{'bin (keV)':>15s}  {'rdecay01':>14s}  {'g4gamma':>14s}  {'rel diff':>10s}")
    for i in sorted(top):
        rd = rd_pdf[i]
        g4 = g4g_counts[i]
        if rd > 1e-8:
            rel = (g4 - rd) / rd
            print(f"  {rd_edges_keV[i]:>5.0f}-{rd_edges_keV[i+1]:>5.0f}  {rd:>14.6e}  {g4:>14.6e}  {rel:>+10.2%}")

    print("\nChain breakdown:")
    for c in res.contributions:
        if c.activity > 1e-9 or c.gamma_yield > 0:
            print(f"  {c.isotope}  A={c.activity:.4f}  gammas/decay={c.gamma_yield:.4f}")


if __name__ == "__main__":
    main()
