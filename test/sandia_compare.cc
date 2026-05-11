// test/sandia_compare.cc
//
// Validates g4gamma's Bateman solver against SandiaDecay's independent
// implementation by comparing per-nuclide activities at multiple time points
// for ANSTO Opal-reactor-relevant isotopes.
//
// Both sides use the same sandia.decay.xml data file, so the only differences
// can arise from the Bateman solver logic itself (branching, eigenvalue
// handling, chain construction).
//
// Usage:
//   ./build/sandia_compare [--xml path/to/sandia.decay.xml] [--out file.csv]
//                          [--tol 0.01]
//
// Defaults: xml = data/sandia/sandia.decay.nocoinc.min.xml (relative to repo
//           root), tol = 0.01 (1 %).
//
// Literature context:
//   The generalised Bateman equations for branched chains were established by
//   Bateman (1910) and extended by Cetnar (2006).  The test cases here cover
//   the three classical regimes:
//     1. Secular equilibrium    — parent T½ >> daughter T½
//        (Mo-99/Tc-99m, Sr-90/Y-90, Ge-68/Ga-68, Cs-137/Ba-137m, Co-60)
//     2. Transient equilibrium  — parent T½ slightly > daughter T½
//        (I-131/Xe-131m, where Xe-131m branch fraction ≈ 1.17%)
//     3. No equilibrium         — parent T½ < daughter T½
//        (dominated by a parent that decays faster than its offspring builds)
//   Na-24 and Lu-177 are single-step chains that test the trivial case and
//   confirm that simple exponential decay is reproduced exactly.

#include "g4gamma/SandiaProvider.hh"
#include "g4gamma/ChainBuilder.hh"
#include "g4gamma/Bateman.hh"
#include "g4gamma/IsotopeKey.hh"
#include "g4gamma/Units.hh"
#include "SandiaDecay.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace g4gamma;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static const char* kSymbols[] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U"
};

static std::string symbol(int Z) {
    if (Z < 1 || Z > 92) return "??";
    return kSymbols[Z - 1];
}

static std::string nuclideLabel(int Z, int A, int M) {
    std::string s = symbol(Z) + std::to_string(A);
    if (M == 1) s += "m";
    else if (M > 1) s += "m" + std::to_string(M);
    return s;
}

// ---------------------------------------------------------------------------
// Test-case definitions
// ---------------------------------------------------------------------------

struct TestCase {
    std::string label;       // human-readable name
    int    Z, A, M;          // root isotope
    std::string regime;      // Bateman regime (for documentation)
    std::string context;     // ANSTO / reactor relevance
    std::vector<double> t_s; // evaluation times in seconds
};

// Time helpers (seconds)
static constexpr double kMin = 60.0;
static constexpr double kHr  = 3600.0;
static constexpr double kDay = 86400.0;
static constexpr double kYr  = 365.2425 * kDay;

static std::vector<TestCase> buildTestCases() {
    return {
        // --- ANSTO Opal primary medical products ----------------------------
        {
            "Mo99",  42, 99, 0,
            "secular-eq (T½_parent=65.94h >> T½_Tc99m=6.007h)",
            "ANSTO primary product; Mo-99/Tc-99m generator for SPECT imaging",
            {0.1 * 65.94*kHr,   // early: Tc-99m barely built up
             0.5 * 65.94*kHr,   // half T½(Mo-99)
             1.0 * 65.94*kHr,   // 1 T½(Mo-99): Tc-99m near secular eq
             3.0 * 65.94*kHr,   // 3 T½(Mo-99)
             10.0* 65.94*kHr}   // secular eq well-established
        },
        {
            "I131",  53, 131, 0,
            "minor-branch isomer: Xe-131m (1.17% branch, T½=11.93d) and direct Xe-131",
            "ANSTO fission-product/medical; tests minor isomeric branch with longer-lived daughter",
            {0.5 * 8.025*kDay,
             1.0 * 8.025*kDay,  // 1 T½(I-131)
             2.0 * 8.025*kDay,
             5.0 * 8.025*kDay,
             10.0* 8.025*kDay}  // Xe-131m decaying after I-131 mostly gone
        },
        {
            "Lu177", 71, 177, 0,
            "simple-chain (Hf-177 stable daughter)",
            "ANSTO Lu-177 production for targeted radionuclide therapy (TRT)",
            {0.5 * 6.647*kDay,
             1.0 * 6.647*kDay,
             3.0 * 6.647*kDay,
             10.0* 6.647*kDay}
        },
        // --- PET generator (Opal produces Ge-68 for Ga-68 generators) ------
        {
            "Ge68",  32, 68, 0,
            "rapid-secular-eq (T½_parent=270.95d >> T½_Ga68=67.71min)",
            "ANSTO PET generator production; Ga-68 daughter in secular eq within hours",
            {1.0 * 67.71*kMin,   // 1 T½(Ga-68): 50% eq
             5.0 * 67.71*kMin,   // 5 T½(Ga-68): >97% eq
             24.0* kHr,          // 1 day
             30.0* kDay,         // 1 month
             1.0 * 270.95*kDay}  // 1 T½(Ge-68)
        },
        // --- Reactor activation products -----------------------------------
        {
            "Na24",  11, 24, 0,
            "simple-chain (Mg-24 stable)",
            "Reactor activation: Na-23(n,γ)Na-24; monitor for primary-loop leaks",
            {0.5 * 14.96*kHr,
             1.0 * 14.96*kHr,
             3.0 * 14.96*kHr,
             10.0* 14.96*kHr}
        },
        {
            "Co60",  27, 60, 0,
            "simple-chain (Ni-60 stable)",
            "Structural activation (Co-59(n,γ)Co-60); dominant long-lived activation product",
            {0.5 * 5.271*kYr,
             1.0 * 5.271*kYr,
             3.0 * 5.271*kYr}
        },
        // --- Fission products ----------------------------------------------
        {
            "Sr90",  38, 90, 0,
            "secular-eq (T½_Sr90=28.79yr >> T½_Y90=64.0h)",
            "Fission product; Y-90 daughter used in TRT (ANSTO produces Y-90 microspheres)",
            {1.0 * 64.0*kHr,    // 1 T½(Y-90): ~50% secular-eq reached
             5.0 * 64.0*kHr,    // >97% secular-eq
             1.0 * 365*kDay,    // 1 year
             1.0 * 28.79*kYr}   // 1 T½(Sr-90)
        },
        {
            "Cs137", 55, 137, 0,
            "rapid-secular-eq (T½_Cs137=30.08yr >> T½_Ba137m=2.552min)",
            "Fission product; Ba-137m isomer reaches eq in <15 min — tests fast-equilibrium branch",
            {0.5 * 153.12,      // 0.5 T½(Ba-137m) in seconds
             1.0 * 153.12,      // 1 T½(Ba-137m): 50% eq (94.7% branch)
             5.0 * 153.12,      // >97% eq
             1.0 * kHr,         // long after Ba-137m eq
             10.0* kYr}         // quasi-infinite
        },
    };
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    // --- Parse args ---------------------------------------------------------
    std::string xml_path;
    std::string out_path;
    double tol = 0.01;  // 1% default tolerance

    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--xml") && i + 1 < argc) {
            xml_path = argv[++i];
        } else if (!std::strcmp(argv[i], "--out") && i + 1 < argc) {
            out_path = argv[++i];
        } else if (!std::strcmp(argv[i], "--tol") && i + 1 < argc) {
            tol = std::atof(argv[++i]);
        }
    }

    // Default XML path: look for uncompressed file relative to binary
    if (xml_path.empty()) {
        const char* env = std::getenv("SANDIA_DECAY_XML");
        if (env) {
            xml_path = env;
        } else {
            // Try common relative paths (from build/ or repo root)
            for (const char* p : {
                    "../data/sandia/sandia.decay.nocoinc.min.xml",
                    "data/sandia/sandia.decay.nocoinc.min.xml",
                    "../third_party/SandiaDecay/sandia.decay.nocoinc.min.xml"}) {
                std::ifstream f(p);
                if (f.good()) { xml_path = p; break; }
            }
        }
    }

    if (xml_path.empty()) {
        std::cerr << "ERROR: no sandia.decay.xml found. Pass --xml <path> or set SANDIA_DECAY_XML.\n";
        return 2;
    }

    std::cout << "XML: " << xml_path << "\n";
    std::cout << "Tolerance: " << tol * 100 << " %\n\n";

    // --- Initialise providers ----------------------------------------------
    SandiaProvider g4g_provider(xml_path);

    SandiaDecay::SandiaDecayDataBase sd_db;
    sd_db.initialize(xml_path);

    // --- Run tests ----------------------------------------------------------
    struct Row {
        std::string isotope, member;
        double t_s, g4g_act, sd_act, rel_err;
        bool pass;
    };
    std::vector<Row> rows;
    int n_pass = 0, n_fail = 0, n_skip = 0;

    for (const auto& tc : buildTestCases()) {
        std::cout << "=== " << tc.label << " ===\n";
        std::cout << "  Regime : " << tc.regime << "\n";
        std::cout << "  Context: " << tc.context << "\n";

        // g4gamma chain
        IsotopeKey root{tc.Z, tc.A, tc.M};
        ChainBuilder cb(g4g_provider);
        cb.build(root);
        const auto& chain = cb.nodes();
        Bateman bateman(chain);

        // SandiaDecay root nuclide
        const SandiaDecay::Nuclide* sd_root =
            sd_db.nuclide(tc.Z, tc.A, tc.M);
        if (!sd_root) {
            std::cout << "  SKIP: SandiaDecay does not have "
                      << nuclideLabel(tc.Z, tc.A, tc.M) << "\n\n";
            ++n_skip;
            continue;
        }

        for (double t_s : tc.t_s) {
            // g4gamma activities (time in internal units = ns)
            auto g4g_acts = bateman.solve(t_s * units::s);

            // SandiaDecay: 1 Bq initial, query at t_s seconds
            SandiaDecay::NuclideMixture mix;
            mix.addNuclideByActivity(sd_root, 1.0);  // 1 Bq initial

            for (int i = 0; i < static_cast<int>(chain.size()); ++i) {
                // Skip genuinely stable nodes (no further gammas or decays)
                if (chain[i].stable) continue;

                double g4_a = g4g_acts[i];

                // Query SandiaDecay for this specific nuclide's activity
                double sd_a = mix.activity(
                    t_s,
                    chain[i].isotope.Z,
                    chain[i].isotope.A,
                    chain[i].isotope.M);

                // Skip if both are essentially zero
                double max_a = std::max(g4_a, sd_a);
                if (max_a < 1e-9) continue;

                double rel_err = std::abs(g4_a - sd_a) / max_a;
                bool pass = (rel_err <= tol);

                if (pass) ++n_pass; else ++n_fail;

                std::string mem_label =
                    nuclideLabel(chain[i].isotope.Z,
                                 chain[i].isotope.A,
                                 chain[i].isotope.M);

                // Human-readable time label
                char t_label[32];
                if (t_s < kMin)        std::snprintf(t_label, sizeof(t_label), "%6.1f s  ", t_s);
                else if (t_s < kHr)    std::snprintf(t_label, sizeof(t_label), "%6.2f min", t_s/kMin);
                else if (t_s < kDay)   std::snprintf(t_label, sizeof(t_label), "%6.2f hr ", t_s/kHr);
                else if (t_s < kYr)    std::snprintf(t_label, sizeof(t_label), "%6.2f d  ", t_s/kDay);
                else                   std::snprintf(t_label, sizeof(t_label), "%6.2f yr ", t_s/kYr);

                std::printf("  %s  %-9s  g4g=%8.5f  sd=%8.5f  err=%6.3f%%  %s\n",
                            t_label, mem_label.c_str(),
                            g4_a, sd_a, rel_err * 100.0,
                            pass ? "OK" : "FAIL");

                rows.push_back({tc.label, mem_label, t_s, g4_a, sd_a, rel_err, pass});
            }
        }
        std::cout << "\n";
    }

    // --- Write CSV ----------------------------------------------------------
    if (!out_path.empty()) {
        std::ofstream ofs(out_path);
        ofs << "isotope,chain_member,time_s,g4gamma_activity,sandia_activity,"
               "rel_err_pct,pass\n";
        for (const auto& r : rows) {
            ofs << r.isotope << "," << r.member << "," << r.t_s << ","
                << r.g4g_act << "," << r.sd_act << ","
                << r.rel_err * 100.0 << "," << (r.pass ? "1" : "0") << "\n";
        }
        std::cout << "CSV written to " << out_path << "\n";
    }

    // --- Summary ------------------------------------------------------------
    std::cout << "=== SUMMARY ===\n";
    std::cout << "  Pass: " << n_pass << "\n";
    std::cout << "  Fail: " << n_fail << "\n";
    std::cout << "  Skip: " << n_skip << "\n";
    if (n_fail == 0)
        std::cout << "  RESULT: ALL PASS\n";
    else
        std::cout << "  RESULT: " << n_fail << " FAILURE(S)\n";

    return (n_fail == 0) ? 0 : 1;
}
