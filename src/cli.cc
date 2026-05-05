// Simple CLI for testing without Python.
// Usage: g4gamma_cli <Z> <A> <M> <time_seconds_or_-1> <e_min_keV> <e_max_keV> <nbins>
#include "g4gamma/GammaSpectrum.hh"
#include <iostream>
#include <cstdlib>

using namespace g4gamma;

int main(int argc, char** argv) {
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0]
                  << " Z A M time_s(-1=SE) Emin_keV Emax_keV nbins\n";
        return 1;
    }
    int Z = std::atoi(argv[1]);
    int A = std::atoi(argv[2]);
    int M = std::atoi(argv[3]);
    double t_s = std::atof(argv[4]);
    double e_min_keV = std::atof(argv[5]);
    double e_max_keV = std::atof(argv[6]);
    int N = std::atoi(argv[7]);

    using namespace units;
    std::vector<double> edges(N + 1);
    double w = (e_max_keV - e_min_keV) / N;
    for (int i = 0; i <= N; ++i) edges[i] = (e_min_keV + i * w) * keV;

    SpectrumOptions opts;
    opts.includeAnnihilation = true;
    opts.includeXrays        = false;

    GammaSpectrumBuilder builder(opts);
    double t = (t_s < 0) ? -1.0 : (t_s * second);
    auto res = builder.build(IsotopeKey{Z, A, M}, t, edges);

    std::cout << "# Chain (" << res.contributions.size() << " members):\n";
    for (const auto& c : res.contributions) {
        std::cout << "#   " << c.isotope.str()
                  << "  meanLife=" << c.meanLife / second << " s"
                  << "  A=" << c.activity
                  << "  gammasPerDecay=" << c.gammaYield << "\n";
    }
    std::cout << "# bin_low_keV  bin_high_keV  counts_per_primary_decay\n";
    double total = 0.0;
    for (size_t i = 0; i < res.counts.size(); ++i) {
        std::cout << edges[i] / keV << "  " << edges[i+1] / keV
                  << "  " << res.counts[i] << "\n";
        total += res.counts[i];
    }
    std::cout << "# Total in bins: " << total << "\n";
    return 0;
}
