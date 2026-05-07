// python/bindings.cc -- pybind11 module exposing g4gamma to Python.
//
// Mirrors the C++ API while accepting numpy arrays for bin edges and
// returning numpy arrays for the histogram counts. Also exposes a `units`
// submodule with the same constants as g4gamma::units so users can write:
//
//   import g4gamma as g
//   edges = np.linspace(0, 3000, 3001) * g.units.keV
//   res = g.build_spectrum(g.IsotopeKey(55, 137, 0), -1, edges)
//
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "g4gamma/Units.hh"
#include "g4gamma/IsotopeKey.hh"
#include "g4gamma/DataPath.hh"
#include "g4gamma/GammaSpectrum.hh"

namespace py = pybind11;
using namespace g4gamma;

static SpectrumResult build_spectrum(GammaSpectrumBuilder& self,
                                      const IsotopeKey& primary,
                                      double t,
                                      py::array_t<double> edges) {
    auto buf = edges.request();
    if (buf.ndim != 1) throw std::runtime_error("edges must be 1D");
    const double* p = static_cast<const double*>(buf.ptr);
    std::vector<double> v(p, p + buf.size);
    return self.build(primary, t, v);
}

PYBIND11_MODULE(g4gamma, m) {
    m.doc() = "Direct gamma-spectrum extraction from Geant4 decay datasets, "
              "matching rdecay01-style output without running a simulation.";

    // ---- Units submodule -------------------------------------------------
    auto u = m.def_submodule("units", "Geant4/CLHEP system of units (energy in MeV, time in ns)");
    u.attr("MeV") = units::MeV;
    u.attr("eV")  = units::eV;
    u.attr("keV") = units::keV;
    u.attr("GeV") = units::GeV;
    u.attr("ns")  = units::ns;
    u.attr("ps")  = units::ps;
    u.attr("us")  = units::us;
    u.attr("ms")  = units::ms;
    u.attr("s")   = units::s;
    u.attr("second") = units::second;
    u.attr("minute") = units::minute;
    u.attr("hour")   = units::hour;
    u.attr("day")    = units::day;
    u.attr("year")   = units::year;
    u.attr("h") = units::h;
    u.attr("d") = units::d;
    u.attr("y") = units::y;
    u.attr("electron_mass_c2") = units::electron_mass_c2;

    // ---- IsotopeKey ------------------------------------------------------
    py::class_<IsotopeKey>(m, "IsotopeKey")
        .def(py::init<int,int,int>(), py::arg("Z"), py::arg("A"), py::arg("M") = 0)
        .def_readwrite("Z", &IsotopeKey::Z)
        .def_readwrite("A", &IsotopeKey::A)
        .def_readwrite("M", &IsotopeKey::M)
        .def("__str__", &IsotopeKey::str)
        .def("__repr__", [](const IsotopeKey& k){
            return "IsotopeKey(Z=" + std::to_string(k.Z) +
                   ", A=" + std::to_string(k.A) +
                   ", M=" + std::to_string(k.M) + ", name='" + k.str() + "')";
        })
        .def("__eq__", [](const IsotopeKey& a, const IsotopeKey& b){ return a == b; })
        .def("__hash__", [](const IsotopeKey& k){
            return std::hash<IsotopeKey>{}(k);
        });

    // ---- ChainContribution -----------------------------------------------
    py::class_<ChainContribution>(m, "ChainContribution")
        .def_readonly("isotope",     &ChainContribution::isotope)
        .def_readonly("activity",    &ChainContribution::activity)
        .def_readonly("gamma_yield", &ChainContribution::gammaYield)
        .def_readonly("mean_life",   &ChainContribution::meanLife)
        .def("__repr__", [](const ChainContribution& c){
            return "ChainContribution(" + c.isotope.str() +
                   ", A=" + std::to_string(c.activity) +
                   ", gammasPerDecay=" + std::to_string(c.gammaYield) + ")";
        });

    // ---- SpectrumResult --------------------------------------------------
    py::class_<SpectrumResult>(m, "SpectrumResult")
        .def_property_readonly("bin_edges", [](const SpectrumResult& r){
            return py::array_t<double>(r.binEdges.size(), r.binEdges.data());
        })
        .def_property_readonly("counts", [](const SpectrumResult& r){
            return py::array_t<double>(r.counts.size(), r.counts.data());
        })
        .def_readonly("contributions", &SpectrumResult::contributions)
        .def_readonly("source_name", &SpectrumResult::sourceName);

    // ---- DataSource enum -------------------------------------------------
    py::enum_<DataSource>(m, "DataSource")
        .value("Geant4", DataSource::Geant4)
        .value("Sandia", DataSource::Sandia)
        .value("Lara",   DataSource::Lara)
        .export_values();

    // ---- SpectrumOptions -------------------------------------------------
    py::class_<SpectrumOptions>(m, "SpectrumOptions")
        .def(py::init<>())
        .def_readwrite("source",                 &SpectrumOptions::source)
        .def_readwrite("include_annihilation",   &SpectrumOptions::includeAnnihilation)
        .def_readwrite("include_xrays",          &SpectrumOptions::includeXrays)
        .def_readwrite("max_chain_depth",        &SpectrumOptions::maxChainDepth)
        .def_readwrite("isomer_lifetime_thresh", &SpectrumOptions::isomerLifetimeThresh)
        .def_readwrite("geant4_sh",              &SpectrumOptions::geant4Sh)
        .def_readwrite("sandia_xml",             &SpectrumOptions::sandiaXml)
        .def_readwrite("lara_dir",               &SpectrumOptions::laraDir)
        .def_readwrite("level_match_tolerance",  &SpectrumOptions::levelMatchTolerance)
        .def_readwrite("verbose",                &SpectrumOptions::verbose);

    // ---- GammaSpectrumBuilder --------------------------------------------
    py::class_<GammaSpectrumBuilder>(m, "GammaSpectrumBuilder")
        .def(py::init<const SpectrumOptions&>(), py::arg("options") = SpectrumOptions{})
        .def("build", &build_spectrum,
             py::arg("primary"), py::arg("t"), py::arg("bin_edges"),
             "Build the spectrum. t<0 means secular equilibrium. Energies in internal units (MeV).")
        .def("source_name", &GammaSpectrumBuilder::sourceName);

    // ---- Top-level convenience function ----------------------------------
    m.def("build_spectrum",
          [](const IsotopeKey& iso, double t, py::array_t<double> edges,
             bool include_annihilation, bool include_xrays,
             const std::string& geant4_sh, int verbose,
             DataSource source, const std::string& sandia_xml,
             const std::string& lara_dir) {
              SpectrumOptions o;
              o.source              = source;
              o.includeAnnihilation = include_annihilation;
              o.includeXrays        = include_xrays;
              o.geant4Sh            = geant4_sh;
              o.sandiaXml           = sandia_xml;
              o.laraDir             = lara_dir;
              o.verbose             = verbose;
              GammaSpectrumBuilder b(o);
              return build_spectrum(b, iso, t, edges);
          },
          py::arg("primary"), py::arg("t"), py::arg("bin_edges"),
          py::arg("include_annihilation") = true,
          py::arg("include_xrays") = false,
          py::arg("geant4_sh") = "",
          py::arg("verbose") = 0,
          py::arg("source") = DataSource::Geant4,
          py::arg("sandia_xml") = "",
          py::arg("lara_dir") = "",
          "One-shot convenience wrapper. Set source=g.DataSource.Sandia for "
          "SandiaDecay or g.DataSource.Lara for LARA/DDEP.");

    // ---- DataPath helpers (useful for debugging path resolution) ---------
    m.def("locate_radioactive_data", &DataPath::radioactiveDecayDir,
          py::arg("geant4_sh") = "");
    m.def("locate_photon_evap_data", &DataPath::photonEvaporationDir,
          py::arg("geant4_sh") = "");
    m.def("locate_le_data",          &DataPath::lowEnergyDir,
          py::arg("geant4_sh") = "");
}
