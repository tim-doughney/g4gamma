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
# Unit tests (synthetic data, 27 cases)
python ../test/run_tests.py          # from inside build/
PYTHONPATH=build python test/run_tests.py  # from repo root

# Validate against rdecay01 Geant4 Monte Carlo (bundled 10k-event CSV):
python test/validate_against_geant4.py buildG4RadDecayExample/u238_AFtrue_h1_3.csv 10000
# Th-232 (1e7-event run):
python test/validate_against_geant4.py buildG4RadDecayExample/th232_AFtrue_1e7_h1_3.csv 10000000
# The script auto-detects the isotope from the CSV filename (e.g. u238, th232).
# For your own run: adjust n_primaries to match /run/beamOn in your macro
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

`DecayBranch` has an optional `terminals` field (vector of `(IsotopeKey, fraction)` pairs). When set by `Geant4Provider`, `ChainBuilder` uses it *instead of* `daughter` for activity routing — splitting the daughter activity by terminal fraction while the full `branchingRatio` still governs emission weight. This correctly handles the case where a β⁻ channel to a non-isomeric excitation cascades to an isomeric intermediate level (e.g., Th-234 → Pa at 166 keV → 100% Pa-234m via photon evaporation).

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

`data/lara/lara.tar.gz` contains 48 nuclides covering the full U-238 and Th-232 chains. Run `data/lara/fetch_lara.sh [--pack]` to re-download.

**Th-230** is not at the standard `/nuclides/` path on LNHB but is available via LaraWEB. `fetch_lara.sh` discovers the versioned filename (`Th-230_@03.lara.txt`) by POSTing to `Result_Lara2.php` and downloads it automatically. The real file (27 emission lines, including Ra K X-rays) is included in `lara.tar.gz`.

**Th-232 chain intermediates** (Ra-224, Rn-220, Po-216, Po-212) are now included. Ra-224 and Po-212 have LNHB data files; Rn-220 and Po-216 are pure-alpha emitters with chain-structure-only entries. These nuclides are required to connect Pb-212/Bi-212/Tl-208 to the rest of the Th-232 chain.

**Newer LARA file format (NIST 2025+):** Some files (e.g. Pb-212) use a cascade format with a 9th "Parent" column and embed gammas from all daughter nuclides in one file. `LaraProvider` detects this (by checking the Parent column) and skips daughter-nuclide emissions — they are counted separately when each daughter's own file is loaded. Without this guard, cascade-format files would double-count all daughter gammas.

## rdecay01 settings → g4gamma options mapping

| rdecay01 / PhysicsList setting | g4gamma equivalent |
|---|---|
| `thresholdForVeryLongDecayTime 1e60 year` | `t = -1` (secular equilibrium) |
| `radioactiveDecay->SetARM(true)` + `G4UAtomicDeexcitation` | `full_xray_cascade = True` (K→L→M cascade) |
| `SetAugerCascade(true)`, `SetDeexcitationIgnoreCut(true)` | `full_xray_cascade = True` (Auger/CK secondary vacancies propagated) |
| Beta+ produces 511 keV annihilation pairs | `include_annihilation = True` |

Setting `full_xray_cascade=True` propagates both fluorescence (fl-tr-pr-Z.dat) and Auger/Coster-Kronig (au-tr-pr-Z.dat) secondary vacancies through K→L→M→N shells, reproducing `SetARM(true)` + `SetAugerCascade(true)`. IC electrons are attributed to the **daughter** atom's electron cloud (`Zfluor = ch.daughter.Z`) in all decay modes. `include_xrays=True` gives K-shell fluorescence only, no Auger cascade. L X-rays (10–16 keV for Pb/Bi) are below NaI/LaBr3 threshold but included for HPGe comparisons.

## Geant4 backend cascade terminal tracking

`getCascade()` tracks where photon-evaporation probability stops at long-lived intermediate levels (τ > 1 ns). The `terminalProb` field in `Cascade` records the distribution; `ChainBuilder` uses `DecayBranch::terminals` to split daughter activity correctly. This reproduces Geant4's behavior: e.g., Th-234 β- channels to Pa-234 at 166–187 keV all cascade 100% to Pa-234m (73.92 keV), giving Pa-234m its correct SE activity of ~1 Bq while Pa-234 ground (6.7 h) receives only ~0.16% from the Pa-234m IT branch.

**Ghost isomer guard:** Before stopping a cascade at a long-lived level, `getCascade()` verifies that G4RADIOACTIVEDATA contains a P-block for that level (i.e., a `z{Z}.a{A}` file with a parent-excitation entry matching within `levelTolerance`). If no P-block exists, the level is a *ghost isomer* — ENSDFSTATE records its lifetime but Geant4 has no explicit decay data for it, so Geant4 de-excites it immediately via photon evaporation. In this case `dstIsIsomer` is forced false and the cascade continues. The corresponding guard in `compute()` checks `parentFor(br.daughter) != nullptr` before setting `deferToIsomer` — if false, the cascade fires directly rather than deferring.

The canonical example is Ra-224 level 84.372 keV (T½=7.48×10⁻¹⁰ s → τ=1.08 ns, marginally above the 1 ns threshold). ENSDFSTATE records it; no `z88.a224.m1` exists in G4RADIOACTIVEDATA. Without the guard, 26% of Th-228 decay activity was silently dropped (Ra-224 activity=0.736 instead of 1.000). With the guard: Ra-224 activity=1.000 and Th-232 SE total 4.257 γ/primary vs rdecay01 simulation 4.261 (−0.1%).

Validation results:
- **U-238 SE** (`full_xray_cascade=True`): model total 3.221 γ/primary vs rdecay01 10,000-event simulation 3.247 (−0.8%). All major gamma peaks within ±2% of ENSDF reference.
- **Th-232 SE** (`full_xray_cascade=True`): model total 4.257 γ/primary vs rdecay01 1,000,000-event simulation 4.261 (−0.1%). Key peaks: Pb-212 239 keV −0.2%, Tl-208 583 keV +1.4%, Bi-212 727 keV +0.7%, Ac-228 911 keV +3.8%, Tl-208 2615 keV +0.1%.
