#!/usr/bin/env python3
"""
test/run_tests.py -- self-contained validation suite for g4gamma.

Runs against synthetic data files derived from real ENSDF parameters for
Cs-137, K-40, Na-22, and a synthetic U-238 mini-chain. Each test checks
either an exact analytical prediction or a known reference value.

Run from inside build/:   python ../test/run_tests.py
Or specify the module dir: PYTHONPATH=build python test/run_tests.py
"""
import sys
import os
import numpy as np

# Test data must be at /tmp/g4testdata relative to wherever the script lives.
# This script writes the data files itself so it's portable.
TESTDATA = "/tmp/g4testdata"
RAD = f"{TESTDATA}/RadioactiveDecay"
EVAP = f"{TESTDATA}/PhotonEvaporation"
LE = f"{TESTDATA}/G4EMLOW"
FLUOR = f"{LE}/fluor"


def write_test_data():
    for d in (RAD, EVAP, FLUOR):
        os.makedirs(d, exist_ok=True)

    # ---- Cs-137 in REAL Geant4 format (string mode names, no daughter
    # ---- ground-state P block in Ba-137 file) -----------------------------
    # Uses subdir /real to avoid clobbering the integer-mode test data.
    real_rad = f"{TESTDATA}/RealRad"
    real_evap = f"{TESTDATA}/RealEvap"
    real_ensdf = f"{TESTDATA}/RealEnsdf"
    for d in (real_rad, real_evap, real_ensdf):
        os.makedirs(d, exist_ok=True)
    open(f"{real_rad}/z55.a137", "w").write("""\
# 137CS ( 30.08 Y   )
#  Excitation  flag   Halflife  Mode    Daughter Ex flag   Intensity          Q
P            0  - 9.492526e+08
                              BetaMinus            0               1
                              BetaMinus            0  -          5.3      1175.63
                              BetaMinus      661.659  -         94.7       513.97
""")
    open(f"{real_rad}/z56.a137", "w").write("""\
# 137BA ( 2.552 M   )
P      661.659  -       153.12
                                     IT            0               1
""")
    open(f"{real_evap}/z55.a137", "w").write("0  -\n0.0  0.0  3.5  0\n")
    open(f"{real_evap}/z56.a137", "w").write("""\
0  -
0.0  0.0  1.5  0
1  -
661.659  2.209e+11  5.5  1
0  661.659  1.0  4  0.0  0.111  0.0913  0.0151  0.000756  0.000746  0.00154  0.000292  0.000295  0.0  0.0  0.0
""")
    open(f"{real_ensdf}/ENSDFSTATE.dat", "w").write("""\
55 137 0.0       -    1.369e+18  7  0.0
56 137 0.0       -    -1.0       3  0.0
56 137 661.659   -    2.209e+11  11  0.0
""")

    # ---- Cs-137 (legacy integer-mode synthetic data) ---------------------
    open(f"{RAD}/z55.a137", "w").write("""\
P                  0.0       -      9.491e+08
                   1               0.000               1.000
                   1             661.659          -        94.700                1175.630
                   1               0.000          -         5.300                1175.630
""")
    open(f"{RAD}/z56.a137", "w").write("""\
P                  0.0       -      0.0
P             661.659       -      153.12
                   0               0.000               1.000
""")
    open(f"{EVAP}/z55.a137", "w").write("0  -\n0.0  0.0  3.5  0\n")
    open(f"{EVAP}/z56.a137", "w").write("""\
0  -
0.0  0.0  1.5  0
1  -
661.659  2.209e+11  5.5  1
0  661.659  1.0  4  0.0  0.111  0.0913  0.0151  0.000756  0.000746  0.00154  0.000292  0.000295  0.0  0.0  0.0
""")

    # ---- K-40 -------------------------------------------------------------
    open(f"{RAD}/z19.a40", "w").write("""\
P                  0.0       -      3.937e+16
                   1               0.000              0.8928
                   3               0.000              0.0918
                   4               0.000              0.0127
                   5               0.000              0.0010
                   1             0.0     -        100.000          1311.073
                   3             1460.820  -      100.000          1504.400
                   4             1460.820  -      100.000          1504.400
                   5             1460.820  -      100.000          1504.400
""")
    open(f"{EVAP}/z18.a40", "w").write("""\
0  -
0.0  0.0  0.0  0
1  -
1460.820  0.00162  2.0  1
0  1460.820  1.0  2  0.0  0.0
""")
    open(f"{EVAP}/z20.a40", "w").write("0  -\n0.0  0.0  0.0  0\n")

    # ---- Na-22 ------------------------------------------------------------
    open(f"{RAD}/z11.a22", "w").write("""\
P                  0.0       -      8.21e+07
                   2               0.000               0.9036
                   3               0.000               0.0964
                   2            1274.500          -        99.940              1820.700
                   2               0.000          -         0.060              1820.700
                   3            1274.500          -       100.000              1820.700
""")
    open(f"{EVAP}/z10.a22", "w").write("""\
0  -
0.0  0.0  0.0  0
1  -
1274.500  5.24e-3  2.0  1
0  1274.500  1.0  2  0.0  0.0
""")
    open(f"{EVAP}/z11.a22", "w").write("0  -\n0.0  0.0  3.0  0\n")

    # ---- Argon fluor data (for K-40 X-ray test) --------------------------
    open(f"{FLUOR}/fl-tr-pr-18.dat", "w").write("""\
1 1 1
3 0.115 0.002957
4 0.005 0.003190
5 0.002 0.003190
-1 -1 -1
3 3 3
-1 -1 -1
-2 -2 -2
""")


# ----- Test helpers --------------------------------------------------------
class Test:
    passed = 0
    failed = 0

    @classmethod
    def check(cls, name, value, expected, tol=1e-3):
        ok = abs(value - expected) <= tol * max(1.0, abs(expected))
        print(f"  {'OK' if ok else 'FAIL'}  {name}: got {value:.6f}, expected {expected:.6f}")
        if ok: cls.passed += 1
        else: cls.failed += 1


def main():
    # Make module importable
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(here, "..", "build"),
        "build",
        os.getcwd(),
    ]
    for c in candidates:
        if os.path.isfile(os.path.join(c, [f for f in os.listdir(c) if f.startswith("g4gamma") and f.endswith(".so")][0] if os.path.isdir(c) and any(f.startswith("g4gamma") and f.endswith(".so") for f in os.listdir(c)) else "_NOPE")):
            sys.path.insert(0, c)
            break
    try:
        import g4gamma as g
    except ImportError:
        for c in candidates:
            if os.path.isdir(c):
                sys.path.insert(0, c)
        import g4gamma as g

    write_test_data()
    os.environ["G4RADIOACTIVEDATA"] = RAD
    os.environ["G4LEVELGAMMADATA"] = EVAP
    os.environ["G4LEDATA"] = LE

    edges_keV = np.linspace(0, 3000, 3001)
    edges = edges_keV * g.units.keV

    builder = g.GammaSpectrumBuilder(g.SpectrumOptions())

    # ----- Cs-137 SE -------------------------------------------------------
    print("\n[1] Cs-137 secular equilibrium")
    res = builder.build(g.IsotopeKey(55, 137, 0), -1.0, edges)
    counts = np.array(res.counts)
    Test.check("661 keV peak", counts[661], 0.947 * (1.0 / 1.111))
    Test.check("total integrated", counts.sum(), 0.947 * (1.0 / 1.111))

    # ----- Cs-137 finite t = T_half(Ba137m) --------------------------------
    print("\n[2] Cs-137 at t = T(1/2) of Ba-137m")
    res = builder.build(g.IsotopeKey(55, 137, 0), 153.12 * g.units.s, edges)
    ba_act = next(c.activity for c in res.contributions if c.isotope.M == 1)
    Test.check("Ba-137m activity", ba_act, 0.947 * 0.5, tol=1e-3)

    # ----- K-40 ------------------------------------------------------------
    print("\n[3] K-40 secular equilibrium")
    res = builder.build(g.IsotopeKey(19, 40, 0), -1.0, edges)
    counts = np.array(res.counts)
    Test.check("1460 keV peak", counts[1460], 0.0918 + 0.0127 + 0.0010, tol=1e-4)

    # ----- K-40 with X-rays ------------------------------------------------
    print("\n[4] K-40 with X-rays enabled")
    opts = g.SpectrumOptions(); opts.include_xrays = True
    res = g.GammaSpectrumBuilder(opts).build(g.IsotopeKey(19, 40, 0), -1.0, edges)
    counts = np.array(res.counts)
    Test.check("K-alpha at 2.957 keV (in [2,3])", counts[2], 0.0918 * 0.115)
    Test.check("K-beta at 3.190 keV (in [3,4])", counts[3], 0.0918 * 0.007)

    # ----- Na-22 ----------------------------------------------------------
    print("\n[5] Na-22 with annihilation")
    res = builder.build(g.IsotopeKey(11, 22, 0), -1.0, edges)
    counts = np.array(res.counts)
    Test.check("511 keV (annihilation)", counts[510], 2.0 * 0.9036)
    Test.check("1274 keV", counts[1274], 0.9036 * 0.9994 + 0.0964)

    # ----- Na-22 without annihilation -------------------------------------
    print("\n[6] Na-22 without annihilation")
    opts = g.SpectrumOptions(); opts.include_annihilation = False
    res = g.GammaSpectrumBuilder(opts).build(g.IsotopeKey(11, 22, 0), -1.0, edges)
    counts = np.array(res.counts)
    Test.check("no 511 keV peak", counts[510], 0.0, tol=1e-9)
    Test.check("1274 keV unchanged", counts[1274], 0.9036 * 0.9994 + 0.0964)

    # ----- Bin edge cases -------------------------------------------------
    print("\n[7] Out-of-range peak gives zero counts")
    edges2 = np.array([200, 300, 400, 500.0]) * g.units.keV
    res = builder.build(g.IsotopeKey(55, 137, 0), -1.0, edges2)
    Test.check("counts sum (peak at 661 keV outside [200,500])",
               np.array(res.counts).sum(), 0.0, tol=1e-9)

    # ----- Real Geant4 format (string modes, no daughter ground P block) --
    print("\n[8] Cs-137 in real Geant4 file format (with ENSDFSTATE)")
    os.environ["G4RADIOACTIVEDATA"]  = f"{TESTDATA}/RealRad"
    os.environ["G4LEVELGAMMADATA"]   = f"{TESTDATA}/RealEvap"
    os.environ["G4ENSDFSTATEDATA"]   = f"{TESTDATA}/RealEnsdf"
    res = g.GammaSpectrumBuilder(g.SpectrumOptions()).build(
        g.IsotopeKey(55, 137, 0), -1.0, edges_keV * g.units.keV)
    counts = np.array(res.counts)
    Test.check("661 keV peak (real format)", counts[661], 0.947 * (1.0 / 1.111))
    # Check the chain found Ba-137m (M=1)
    has_ba137m = any(c.isotope.M == 1 and c.isotope.Z == 56 for c in res.contributions)
    Test.check("Ba-137m found in chain (real format)", float(has_ba137m), 1.0)

    print(f"\n--- {Test.passed} passed, {Test.failed} failed ---")
    sys.exit(0 if Test.failed == 0 else 1)


if __name__ == "__main__":
    main()
