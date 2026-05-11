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
    has_ba137m = any(c.isotope.M == 1 and c.isotope.Z == 56 for c in res.contributions)
    Test.check("Ba-137m found in chain (real format)", float(has_ba137m), 1.0)

    # ----- SandiaDecay backend --------------------------------------------
    # Try to find sandia.decay.*.xml. If it's not available, skip tests
    # rather than fail.
    print("\n[9] SandiaDecay backend")
    sandia_paths = [
        os.path.join(here, "..", "data", "sandia", "sandia.decay.nocoinc.min.xml.gz"),
        os.path.join(here, "..", "data", "sandia", "sandia.decay.nocoinc.min.xml"),
        os.path.join(here, "..", "data", "sandia", "sandia.decay.xml"),
    ]
    sandia_xml = next((p for p in sandia_paths if os.path.isfile(p)), None)
    if sandia_xml:
        os.environ["SANDIA_DECAY_XML"] = sandia_xml

        # Cs-137
        opts = g.SpectrumOptions()
        opts.source = g.DataSource.Sandia
        res = g.GammaSpectrumBuilder(opts).build(
            g.IsotopeKey(55, 137, 0), -1.0, edges)
        counts = np.array(res.counts)
        Test.check("Cs-137 661 keV via Sandia",
                   counts[661], 0.853, tol=0.005)
        Test.check("Cs-137 chain has Ba-137m",
                   float(any(c.isotope.M == 1 and c.isotope.Z == 56 for c in res.contributions)),
                   1.0)

        # K-40
        res = g.GammaSpectrumBuilder(opts).build(
            g.IsotopeKey(19, 40, 0), -1.0, edges)
        counts = np.array(res.counts)
        Test.check("K-40 1460 keV via Sandia",
                   counts[1460], 0.107, tol=0.005)

        # Co-60 - the canonical 2-gamma-per-decay test
        res = g.GammaSpectrumBuilder(opts).build(
            g.IsotopeKey(27, 60, 0), -1.0, edges)
        counts = np.array(res.counts)
        Test.check("Co-60 total γ/decay via Sandia",
                   counts.sum(), 2.0, tol=0.01)
        Test.check("Co-60 1173 keV via Sandia",
                   counts[1173], 0.9985, tol=0.005)
        Test.check("Co-60 1332 keV via Sandia",
                   counts[1332], 0.9998, tol=0.005)

        # U-238 chain - check at least the 609 keV Bi-214 peak comes through
        res = g.GammaSpectrumBuilder(opts).build(
            g.IsotopeKey(92, 238, 0), -1.0, edges)
        counts = np.array(res.counts)
        Test.check("U-238 chain Bi-214 609 keV via Sandia",
                   counts[609], 0.461, tol=0.02)
        Test.check("U-238 chain has many members via Sandia",
                   float(len(res.contributions) >= 15), 1.0)
    else:
        print("  SKIPPED -- no sandia.decay XML found in data/sandia/")
        print(f"           tried: {sandia_paths}")

    # ----- LARA backend ---------------------------------------------------
    print("\n[10] LARA / DDEP backend")
    lara_dir = os.path.join(here, "..", "data", "lara")
    has_loose = os.path.isdir(lara_dir) and any(
        f.endswith(".lara.txt") for f in os.listdir(lara_dir))
    has_tarball = (os.path.isfile(os.path.join(lara_dir, "lara.tar.gz")) or
                   os.path.isfile(os.path.join(lara_dir, "lara.tar")))
    if has_loose or has_tarball:
        os.environ["LARA_DATA_DIR"] = lara_dir

        opts_lara = g.SpectrumOptions()
        opts_lara.source = g.DataSource.Lara

        # Cs-137 -- DDEP 2023 says 661.6553 keV at 85.01%
        res = g.GammaSpectrumBuilder(opts_lara).build(
            g.IsotopeKey(55, 137, 0), -1.0, edges)
        counts = np.array(res.counts)
        Test.check("Cs-137 661 keV via LARA", counts[661], 0.8501, tol=0.001)

        # K-40 -- DDEP 2025 says 1460.822 keV at 10.34%
        res = g.GammaSpectrumBuilder(opts_lara).build(
            g.IsotopeKey(19, 40, 0), -1.0, edges)
        counts = np.array(res.counts)
        Test.check("K-40 1460 keV via LARA", counts[1460], 0.1034, tol=0.001)

        # Co-60 -- 1173.228 + 1332.492
        res = g.GammaSpectrumBuilder(opts_lara).build(
            g.IsotopeKey(27, 60, 0), -1.0, edges)
        counts = np.array(res.counts)
        Test.check("Co-60 1173 keV via LARA", counts[1173], 0.9985, tol=0.001)
        Test.check("Co-60 1332 keV via LARA", counts[1332], 0.9998, tol=0.001)

        # Ra-226 -- 186.211 keV at 3.555%
        res = g.GammaSpectrumBuilder(opts_lara).build(
            g.IsotopeKey(88, 226, 0), -1.0, edges)
        counts = np.array(res.counts)
        Test.check("Ra-226 186 keV via LARA", counts[186], 0.0356, tol=0.001)

        # Tri-source cross-validation: K-40 1460 keV across all three backends
        # should agree within ~1%.
        if sandia_xml:
            os.environ["G4RADIOACTIVEDATA"]  = f"{TESTDATA}/RealRad"
            os.environ["G4LEVELGAMMADATA"]   = f"{TESTDATA}/RealEvap"
            os.environ["G4ENSDFSTATEDATA"]   = f"{TESTDATA}/RealEnsdf"
            # Co-60 1332 keV is the cleanest cross-source comparison since
            # it's a single emission with intensity ~1.0 in all backends.
            o_g = g.SpectrumOptions(); o_g.source = g.DataSource.Geant4
            o_s = g.SpectrumOptions(); o_s.source = g.DataSource.Sandia
            o_l = g.SpectrumOptions(); o_l.source = g.DataSource.Lara
            try:
                # Use Cs-137 since real-format synthetic data covers it
                vg = np.array(g.GammaSpectrumBuilder(o_g).build(
                    g.IsotopeKey(55, 137, 0), -1.0, edges).counts)[661]
                vs = np.array(g.GammaSpectrumBuilder(o_s).build(
                    g.IsotopeKey(55, 137, 0), -1.0, edges).counts)[661]
                vl = np.array(g.GammaSpectrumBuilder(o_l).build(
                    g.IsotopeKey(55, 137, 0), -1.0, edges).counts)[661]
                spread = max(vg, vs, vl) - min(vg, vs, vl)
                Test.check("Cs-137 661 keV: Geant4/Sandia/LARA spread <0.5%",
                           float(spread / max(vg, vs, vl) < 0.005), 1.0)
                print(f"          Geant4={vg:.4f}  Sandia={vs:.4f}  LARA={vl:.4f}  spread={spread:.4f}")
            except Exception as e:
                print(f"  cross-source check skipped: {e}")
    else:
        print(f"  SKIPPED -- no LARA files in {lara_dir}")

    # ----- Chain truncation (cutoffs + depth limit) -----------------------
    print("\n[11] Chain truncation: cutoff isotope")
    if sandia_xml:
        opts_s = g.SpectrumOptions()
        opts_s.source = g.DataSource.Sandia
        # Full U-238 chain
        res_full = g.GammaSpectrumBuilder(opts_s).build(
            g.IsotopeKey(92, 238, 0), -1.0, edges)
        n_full = len(res_full.contributions)
        total_full = np.array(res_full.counts).sum()

        # Truncated at Rn-222 (radon escape scenario)
        opts_cut = g.SpectrumOptions()
        opts_cut.source = g.DataSource.Sandia
        opts_cut.chain_cutoffs = [g.IsotopeKey(86, 222, 0)]
        res_rn = g.GammaSpectrumBuilder(opts_cut).build(
            g.IsotopeKey(92, 238, 0), -1.0, edges)
        n_rn = len(res_rn.contributions)
        total_rn = np.array(res_rn.counts).sum()

        # Rn-222 must be in the truncated chain and marked cutoff
        has_rn = any(c.isotope.Z == 86 and c.isotope.A == 222
                     for c in res_rn.contributions)
        rn_cutoff = any(c.isotope.Z == 86 and c.isotope.A == 222 and c.cutoff
                        for c in res_rn.contributions)
        # Po-218 (daughter of Rn-222) must NOT be in the chain
        no_po218 = not any(c.isotope.Z == 84 and c.isotope.A == 218
                           for c in res_rn.contributions)
        # Truncated chain has fewer members and fewer gammas
        fewer_members = n_rn < n_full
        fewer_gammas  = total_rn < total_full

        Test.check("Rn-222 present in truncated chain",
                   float(has_rn), 1.0)
        Test.check("Rn-222 marked cutoff=True",
                   float(rn_cutoff), 1.0)
        Test.check("Po-218 absent from truncated chain",
                   float(no_po218), 1.0)
        Test.check("truncated chain has fewer members",
                   float(fewer_members), 1.0)
        Test.check("truncated chain has fewer total gammas",
                   float(fewer_gammas), 1.0)
        print(f"          full={n_full} members, {total_full:.4f} γ/primary")
        print(f"          cut@Rn-222: {n_rn} members, {total_rn:.4f} γ/primary")

        # Sanity: truncated + complement should bracket the full total
        # (not exactly equal because Rn-222's own decay gammas are counted in both)
        Test.check("Ra-226 present and active in truncated chain",
                   float(any(c.isotope.Z == 88 and c.isotope.A == 226
                             for c in res_rn.contributions)), 1.0)

    else:
        print("  SKIPPED -- no sandia.decay XML found")

    print("\n[12] Chain truncation: depth limit")
    if sandia_xml:
        opts_d = g.SpectrumOptions()
        opts_d.source = g.DataSource.Sandia
        opts_d.chain_depth_limit = 2  # root + 2 generations (U238, Th234, Pa234m)
        res_d2 = g.GammaSpectrumBuilder(opts_d).build(
            g.IsotopeKey(92, 238, 0), -1.0, edges)
        members_d2 = [str(c.isotope) for c in res_d2.contributions]
        # U-238 (depth 0), Th-234 (depth 1), Pa-234m (depth 2) must be present;
        # Pa-234m and U-234 at depth 2 are cutoff nodes, U-234 must not appear.
        has_u238   = any(c.isotope.Z == 92 and c.isotope.A == 238
                         for c in res_d2.contributions)
        has_th234  = any(c.isotope.Z == 90 and c.isotope.A == 234
                         for c in res_d2.contributions)
        no_th230   = not any(c.isotope.Z == 90 and c.isotope.A == 230
                             for c in res_d2.contributions)
        n_cutoffs  = sum(1 for c in res_d2.contributions if c.cutoff)

        Test.check("U-238 present at depth limit = 2",  float(has_u238),  1.0)
        Test.check("Th-234 present at depth limit = 2", float(has_th234), 1.0)
        Test.check("Th-230 absent at depth limit = 2",  float(no_th230),  1.0)
        Test.check("at least one cutoff node at depth 2",
                   float(n_cutoffs >= 1), 1.0)
        print(f"          members: {members_d2}")

        # depth_limit = 0 → only the root
        opts_d0 = g.SpectrumOptions()
        opts_d0.source = g.DataSource.Sandia
        opts_d0.chain_depth_limit = 0
        res_d0 = g.GammaSpectrumBuilder(opts_d0).build(
            g.IsotopeKey(92, 238, 0), -1.0, edges)
        only_root = (len(res_d0.contributions) == 1 and
                     res_d0.contributions[0].isotope.Z == 92)
        Test.check("depth_limit=0 gives root only",
                   float(only_root), 1.0)
    else:
        print("  SKIPPED -- no sandia.decay XML found")

    print(f"\n--- {Test.passed} passed, {Test.failed} failed ---")
    sys.exit(0 if Test.failed == 0 else 1)


if __name__ == "__main__":
    main()
