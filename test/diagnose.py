#!/usr/bin/env python3
"""
diagnose.py -- run this if g4gamma returns all-zero histograms.

It will tell you which step is failing:
  1. Are env vars / GEANT4_DATA_DIR set?
  2. Did g4gamma resolve them to real directories?
  3. Do the per-isotope decay files actually exist?
  4. What does the chain look like?
  5. What's the output for a known-good isotope (Cs-137)?
"""
import sys
import os

# Adjust this if your build dir is elsewhere
if not any('g4gamma' in os.path.basename(p) for p in sys.path):
    for cand in ('build', '../build', '.'):
        if os.path.isdir(cand) and any(f.startswith('g4gamma') and f.endswith('.so')
                                       for f in os.listdir(cand)):
            sys.path.insert(0, cand)
            break
import g4gamma as g
import numpy as np


def main():
    print("=" * 60)
    print("Step 1: Environment variables")
    print("=" * 60)
    for v in ('G4RADIOACTIVEDATA', 'G4LEVELGAMMADATA', 'G4LEDATA',
              'GEANT4_DATA_DIR'):
        val = os.environ.get(v, '<unset>')
        print(f"  {v} = {val}")

    print()
    print("=" * 60)
    print("Step 2: Path resolution by g4gamma")
    print("=" * 60)
    rad = evap = ledata = None
    try:
        rad = g.locate_radioactive_data()
        print(f"  RadioactiveDecay -> {rad}")
    except Exception as e:
        print(f"  ERROR resolving RadioactiveDecay: {e}")
        print(f"  Tip: try export G4RADIOACTIVEDATA=/path/to/RadioactiveDecay6.1.2")
        sys.exit(1)
    try:
        evap = g.locate_photon_evap_data()
        print(f"  PhotonEvaporation -> {evap}")
    except Exception as e:
        print(f"  ERROR resolving PhotonEvaporation: {e}")
        sys.exit(1)
    try:
        ledata = g.locate_le_data()
        print(f"  G4LEDATA -> {ledata}")
    except Exception as e:
        print(f"  WARNING: G4LEDATA not found ({e}). X-rays will be disabled.")

    print()
    print("=" * 60)
    print("Step 3: Check that data files exist")
    print("=" * 60)
    test_isotopes = [
        (55, 137, "Cs-137"),
        (19, 40, "K-40"),
        (92, 238, "U-238"),
        (56, 137, "Ba-137"),  # daughter of Cs-137
    ]
    for Z, A, label in test_isotopes:
        rad_path = f"{rad}/z{Z}.a{A}"
        evap_path = f"{evap}/z{Z}.a{A}"
        rad_ok = os.path.isfile(rad_path)
        evap_ok = os.path.isfile(evap_path)
        rad_size = os.path.getsize(rad_path) if rad_ok else 0
        evap_size = os.path.getsize(evap_path) if evap_ok else 0
        print(f"  {label}: rad={rad_ok} ({rad_size}B)  evap={evap_ok} ({evap_size}B)")
        if rad_ok and rad_size > 0:
            with open(rad_path) as f:
                first_lines = []
                for i, line in enumerate(f):
                    if i >= 5: break
                    first_lines.append(line.rstrip())
            print(f"    first 5 lines of {rad_path}:")
            for ln in first_lines:
                print(f"      [{len(ln):3d}] {ln!r}")

    print()
    print("=" * 60)
    print("Step 4: Cs-137 chain dump (verbose)")
    print("=" * 60)
    edges = np.linspace(0, 1500, 1501) * g.units.keV
    res = g.build_spectrum(g.IsotopeKey(55, 137, 0), -1.0, edges, verbose=1)
    counts = np.array(res.counts)
    print()
    print(f"  Counts sum: {counts.sum():.6f}")
    print(f"  Non-zero bins: {np.where(counts > 0)[0].tolist()}")

    if counts.sum() == 0:
        print()
        print("  *** ALL ZERO. Check above output -- usually means ***")
        print("  *** the file format wasn't parsed correctly. Inspect the ***")
        print("  *** first 5 lines of the rad file printed above. The format ***")
        print("  *** should look like:")
        print("       P  <excitation_keV>  <floatFlag>  <halfLife_s>")
        print("       <mode>  <0>  <total_BR>            (3 cols, line < 72 chars)")
        print("       <mode>  <daughter_keV>  <flag>  <BR%>  <Q_keV>   (5 cols, line >= 72 chars)")

    print()
    print("=" * 60)
    print("Step 5: K-40 (single peak at 1460.8 keV)")
    print("=" * 60)
    res = g.build_spectrum(g.IsotopeKey(19, 40, 0), -1.0, edges, verbose=1)
    counts = np.array(res.counts)
    print(f"  Counts sum: {counts.sum():.6f}  (expect ~0.107)")
    nz = np.where(counts > 0)[0]
    for i in nz[:10]:
        print(f"    bin {i} ({edges[i]/g.units.keV:.0f}-{edges[i+1]/g.units.keV:.0f} keV): {counts[i]:.6f}")


if __name__ == "__main__":
    main()
