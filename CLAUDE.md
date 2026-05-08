# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project does

g4gamma computes analytic gamma spectra from nuclear decay datasets without running Geant4 Monte Carlo. Given an isotope `(Z, A, M)`, time `t`, and bin edges, it returns a histogram of gammas/primary decay summed over the full decay chain. Used for Tim's NORM thesis work as a forward model and sanity-check against `rdecay01` output.

## Build

```bash
# With pybind11 installed (e.g. via pip):
mkdir build && cd build
cmake -Dpybind11_DIR=$(python3 -c "import pybind11; print(pybind11.get_cmake_dir())") ..
make -j

# Without pybind11 (CMake fetches it):
mkdir build && cd build
cmake ..
make -j
```

CMake produces three targets: `g4gamma_core` (static C++ lib), `g4gamma` (pybind11 `.so` module), and `g4gamma_cli` (C++ CLI binary).

## Running tests

```bash
# Unit tests (synthetic data, 21 cases)
python ../test/run_tests.py          # from inside build/
PYTHONPATH=build python test/run_tests.py  # from repo root

# Validate against rdecay01 Geant4 Monte Carlo (requires a CSV run):
#   1. cd buildG4RadDecayExample && ./rdecay01 moduleTest.mac  (1e7 events ~10 min; adjust /run/beamOn)
#   2. python test/validate_against_geant4.py \
#          buildG4RadDecayExample/u238_AFtrue_h1_3.csv 10000000
```

`test/diagnose.py` is a runtime diagnostic that prints env vars, resolved data directories, chain construction, and spectrum output — run it first when debugging data-path or zero-output issues.

### rdecay01 CSV format

Geant4 writes `<name>_h1_<id>.csv` (e.g. `u238_AFtrue_h1_3.csv`). Format:
- 6 `#...` comment lines, then one column-header line (`entries,Sw,Sw2,Sxw0,Sx2w0`) — **not a comment**
- 1 underflow row, N data rows, 1 overflow row
- **gammas/primary/bin = `entries` column (col 0) ÷ n_primaries**
- `Sw` = entries × 0.01 (rdecay01 fills with weight 0.01); do not use `Sw/n_primaries` directly

## Architecture

The codebase has a clean three-layer design:

**Layer 1 — Data providers** (`IDecayProvider` interface in `include/g4gamma/IDecayProvider.hh`)

Each provider wraps a nuclear-data library and exposes a uniform API:
- `get(IsotopeKey)` → `ParentDecayInfo*` (decay branches + pre-aggregated photon emissions)
- `name()`, `emissionsIncludeXrays()`, `emissionsIncludeAnnihilation()`, `emissionsArePerDecay()`

Concrete providers:
- `Geant4Provider` — reads `$G4RADIOACTIVEDATA` + `$G4LEVELGAMMADATA` + `$G4ENSDFSTATEDATA`; performs cascade computation internally since Geant4 data is per-level
- `SandiaProvider` — reads `sandia.decay.xml`; emissions are already aggregated per-decay
- `LaraProvider` — reads LARA/DDEP data; uses `emissionsArePerDecay()=true` convention

To add a new data source: implement `IDecayProvider` (3–4 methods) and add the `.cc` to `CMakeLists.txt`.

**Layer 2 — Chain solver** (`ChainBuilder` + `Bateman`)

`ChainBuilder` walks the provider graph from the root isotope into a topologically-ordered `std::vector<ChainNode>`. `Bateman` then solves the generalised Bateman equations to get activities:
- `solve(t)` — analytic finite-time solution
- `solveSecularEq()` — secular equilibrium (equivalent to `thresholdForVeryLongDecayTime 1e60 year` in rdecay01); pass `t=-1` in the Python API

**Layer 3 — Binner** (`GammaSpectrum.cc`)

`GammaSpectrumBuilder::build()` iterates chain nodes, multiplies emissions by activity (from Bateman), and accumulates into the user-supplied bin grid. Provider-specific flags (`emissionsArePerDecay`, `emissionsIncludeXrays`, `emissionsIncludeAnnihilation`) control whether the binner adds 511 keV pairs or X-rays itself.

## Internal units

All internal quantities use Geant4/CLHEP conventions: **energy in MeV**, **time in ns**. The `units` submodule (`include/g4gamma/Units.hh`, exposed as `g4gamma.units` in Python) provides conversion factors (`keV`, `s`, `year`, etc.).

## Data path resolution

For the Geant4 backend, paths are resolved in order:
1. Per-dataset env vars (`$G4RADIOACTIVEDATA`, `$G4LEVELGAMMADATA`, `$G4LEDATA`, `$G4ENSDFSTATEDATA`)
2. Parsing `geant4.sh` (default `/opt/geant4/bin/geant4.sh`, or `SpectrumOptions.geant4Sh`)
3. `$GEANT4_DATA_DIR` with auto-detection of versioned subdirectories

`ENSDFSTATE.dat` is required for finite-time queries (`t > 0`). Without it, all queries fall back to secular equilibrium with a warning.

## Key types

- `IsotopeKey(Z, A, M)` — identifies a nuclide; M=0 is ground state
- `SpectrumOptions` — controls source, xray/annihilation flags, verbosity, data paths
- `SpectrumResult` — `binEdges`, `counts` (gammas/primary/bin), `contributions` (per chain member), `sourceName`
- `GammaSpectrumBuilder` — persistent builder that caches provider data; preferred for scanning many isotopes

## Python API quick reference

```python
import g4gamma as g
import numpy as np

edges = np.linspace(0, 3000, 3001) * g.units.keV
res = g.build_spectrum(g.IsotopeKey(55, 137, 0), t=-1, bin_edges=edges)
# res.counts, res.bin_edges, res.source_name, res.contributions

# Persistent builder (caches provider data):
opts = g.SpectrumOptions()
opts.source = g.DataSource.Sandia
builder = g.GammaSpectrumBuilder(opts)
res = builder.build(g.IsotopeKey(92, 238, 0), -1, edges)
```

## C++ CLI

```bash
./build/g4gamma_cli <Z> <A> <M> <time_s_or_-1> <Emin_keV> <Emax_keV> <nbins> [geant4|sandia]
```

## LARA data notes

`data/lara/lara.tar.gz` contains ~44 nuclides including the full U-238 and Th-232 chains. Run `data/lara/fetch_lara.sh [--pack]` to re-download.

**Th-230** is not at the standard `/nuclides/` path on LNHB but is available via LaraWEB. `fetch_lara.sh` discovers the versioned filename (`Th-230_@03.lara.txt`) by POSTing to `Result_Lara2.php` and downloads it automatically. The real file (27 emission lines, including Ra K X-rays) is included in `lara.tar.gz`.

## rdecay01 settings → g4gamma options mapping

| rdecay01 / PhysicsList setting | g4gamma equivalent |
|---|---|
| `thresholdForVeryLongDecayTime 1e60 year` | `t = -1` (secular equilibrium) |
| `radioactiveDecay->SetARM(true)` + `G4UAtomicDeexcitation` | `include_xrays = True` |
| `SetAugerCascade(true)`, `SetDeexcitationIgnoreCut(true)` | K-shell only in g4gamma (L/M not modelled) |
| Beta+ produces 511 keV annihilation pairs | `include_annihilation = True` |

The simulation produces the **full Auger cascade** (K+L+M shells). g4gamma only models K-shell X-rays. Expect g4gamma to underpredict the 50–90 keV region for heavy chain members (Pb, Bi, Po) where L-shell X-rays at ~10–15 keV and M-shell at ~3–5 keV are below detection anyway, but K-shell Pb Kα (~75 keV) and Kβ (~85 keV) are well within NaI/LaBr3 range.
