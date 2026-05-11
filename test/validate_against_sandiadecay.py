#!/usr/bin/env python3
"""
test/validate_against_sandiadecay.py
=====================================
Drives the sandia_compare C++ binary and reports pass/fail.

The binary compares g4gamma's analytic Bateman solver (using SandiaProvider)
against SandiaDecay's independent Bateman implementation.  Both sides read the
same sandia.decay.xml so any disagreement is a solver bug, not a data mismatch.

Usage
-----
  # From repo root (after building):
  PYTHONPATH=build python test/validate_against_sandiadecay.py

  # Explicit paths:
  python test/validate_against_sandiadecay.py \\
      --binary build/sandia_compare \\
      --xml data/sandia/sandia.decay.nocoinc.min.xml \\
      --out comparison.csv

  # Tighter tolerance (default is 1%):
  python test/validate_against_sandiadecay.py --tol 0.001

Exit code
---------
  0  all comparisons pass
  1  one or more comparisons fail (or binary missing / build error)

Test cases (ANSTO Opal-reactor relevant)
-----------------------------------------
Isotope      Regime                           ANSTO relevance
------------ -------------------------------- --------------------------------
Mo-99        secular-eq (Mo-99 >> Tc-99m)    Primary medical product (SPECT)
I-131        minor isomeric branch Xe-131m   Fission product / thyroid therapy
Lu-177       simple chain                     Cancer therapy (TRT)
Ge-68        rapid secular-eq (>> Ga-68)      PET generator production
Na-24        simple exponential               Reactor activation product
Co-60        simple chain                     Structural activation product
Sr-90        secular-eq (Sr-90 >> Y-90)      Fission product / Y-90 therapy
Cs-137       rapid secular-eq (>> Ba-137m)   Fission product / isotropy calibr.

Literature references
---------------------
Bateman (1910) Proc. Cambridge Phil. Soc. 15:423 — original chain equations.
Cetnar (2006)  Ann. Nucl. Energy 33:640 — general branched-chain solution.
"""
import argparse
import csv
import os
import subprocess
import sys


def find_binary(candidates):
    for p in candidates:
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    return None


def find_xml(candidates):
    for p in candidates:
        if os.path.isfile(p):
            return p
    return None


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(here)

    parser = argparse.ArgumentParser(
        description="Validate g4gamma Bateman solver against SandiaDecay")
    parser.add_argument("--binary", default=None,
                        help="Path to sandia_compare binary")
    parser.add_argument("--xml", default=None,
                        help="Path to sandia.decay.xml (uncompressed)")
    parser.add_argument("--out", default=None,
                        help="Output CSV path for detailed results")
    parser.add_argument("--tol", type=float, default=0.01,
                        help="Relative tolerance (default 0.01 = 1%%)")
    args = parser.parse_args()

    # --- Locate binary ------------------------------------------------------
    binary = args.binary or find_binary([
        os.path.join(repo, "build", "sandia_compare"),
        os.path.join(repo, "sandia_compare"),
        "./sandia_compare",
        "./build/sandia_compare",
    ])
    if not binary:
        print("ERROR: sandia_compare binary not found.")
        print("  Build it with:")
        print("    cd build && cmake .. && make sandia_compare")
        print("  (requires git submodule update --init)")
        return 1

    # --- Locate XML ---------------------------------------------------------
    xml = args.xml or os.environ.get("SANDIA_DECAY_XML") or find_xml([
        os.path.join(repo, "data", "sandia", "sandia.decay.nocoinc.min.xml"),
        os.path.join(repo, "third_party", "SandiaDecay",
                     "sandia.decay.nocoinc.min.xml"),
    ])
    if not xml:
        print("ERROR: sandia.decay.xml not found. Pass --xml <path> or set "
              "SANDIA_DECAY_XML.")
        return 1

    # --- Run comparison binary ----------------------------------------------
    out_csv = args.out or os.path.join(repo, "build", "sandia_compare_out.csv")
    cmd = [binary, "--xml", xml, "--out", out_csv, "--tol", str(args.tol)]
    print("Running:", " ".join(cmd))
    print()

    result = subprocess.run(cmd, text=True)
    print()  # blank line after binary output

    if result.returncode == 2:
        print("ERROR: binary could not find XML or other setup error.")
        return 1

    # --- Parse CSV for summary ----------------------------------------------
    if os.path.isfile(out_csv):
        rows = []
        with open(out_csv) as f:
            for row in csv.DictReader(f):
                rows.append(row)

        if rows:
            # Per-isotope summary
            isotopes = {}
            for r in rows:
                iso = r["isotope"]
                isotopes.setdefault(iso, {"pass": 0, "fail": 0, "max_err": 0.0})
                if r["pass"] == "1":
                    isotopes[iso]["pass"] += 1
                else:
                    isotopes[iso]["fail"] += 1
                try:
                    isotopes[iso]["max_err"] = max(
                        isotopes[iso]["max_err"], float(r["rel_err_pct"]))
                except ValueError:
                    pass

            print("\nPer-isotope summary:")
            print(f"  {'Isotope':<10} {'Pass':>5} {'Fail':>5} {'MaxErr%':>8}")
            print(f"  {'-'*10} {'-'*5} {'-'*5} {'-'*8}")
            all_ok = True
            for iso, s in sorted(isotopes.items()):
                ok = s["fail"] == 0
                if not ok:
                    all_ok = False
                mark = "   " if ok else "***"
                print(f"  {mark} {iso:<10} {s['pass']:>5} {s['fail']:>5} "
                      f"{s['max_err']:>8.4f}")
            print()

    # Final verdict
    if result.returncode == 0:
        print("=== ALL PASS ===")
        return 0
    else:
        print(f"=== FAIL (exit code {result.returncode}) ===")
        return 1


if __name__ == "__main__":
    sys.exit(main())
