# g4gamma

Direct gamma-spectrum extraction from Geant4 decay datasets, matching
`rdecay01`-style output without running a Monte Carlo simulation.

Built for Tim's NORM thesis work — computes per-bin gammas/primary at any
chain time using the same data files Geant4 11.3.0 reads at runtime.

## What it does

Given an isotope `(Z, A, M)`, a time `t`, and a bin grid, returns a histogram
of expected gammas per primary decay summed over the full decay chain. The
input matches the convention used in your rdecay01 macros: the primary nuclide
has activity 1 Bq at `t = 0`. At `t = -1` the system is in secular equilibrium
(equivalent to the `thresholdForVeryLongDecayTime 1.0e+60 year` mode).

The output approximates what you'd get from rdecay01 with infinite statistics,
without the runtime cost.

## Sources of gammas

Matches the user's `PhysicsList.cc` / macro configuration:
- Discrete nuclear gammas from PhotonEvaporation cascades following
  beta/EC/alpha decays into excited daughter levels
- Discrete gammas from IT decays of metastable parents
- 511 keV annihilation pairs from beta+ branches (controllable)
- Optional K-shell X-rays from EC decays and from internal-conversion
  vacancies in level transitions (`include_xrays=True`)
- **Not** included: atomic Auger electrons, L/M-shell X-rays (matches
  `SetARM(false)` plus a pragmatic restriction to detectable energies in
  NaI/LaBr3 spectroscopy)

## Building

The `CMakeLists.txt` uses `find_package(pybind11)` first and falls back to
the FetchContent pattern you specified. So either:

```bash
# Option A: pybind11 already installed (e.g. via pip)
mkdir build && cd build
cmake -Dpybind11_DIR=$(python3 -c "import pybind11; print(pybind11.get_cmake_dir())") ..
make -j

# Option B: let CMake fetch pybind11 from GitHub
mkdir build && cd build
cmake ..
make -j
```

On Phoenix HPC with your existing setup, option A should work straight off
since pybind11 is available in your Python env. Option B is the fallback for
machines without pybind11 installed.

## Data path resolution

In order:
1. `$G4RADIOACTIVEDATA`, `$G4LEVELGAMMADATA`, `$G4LEDATA`, `$G4ENSDFSTATEDATA`
   env vars
2. Parse a `geant4.sh` (default `/opt/geant4/bin/geant4.sh`, or the path you pass
   via `SpectrumOptions.geant4_sh`)
3. `$GEANT4_DATA_DIR` with auto-detection of versioned subdirectories
   (`RadioactiveDecay6.1.2`, `PhotonEvaporation6.1`, `G4ENSDFSTATE2.3`, etc.)

This matches your install pattern where `geant4.sh` sets `GEANT4_DATA_DIR`
but leaves the per-dataset env vars commented out.

**ENSDFSTATE.dat is required for finite-time evolution.** The half-life
column in the per-isotope `RadioactiveDecay/zZ.aA` files is a placeholder;
the actual mean lives live in `$G4ENSDFSTATEDATA/ENSDFSTATE.dat`. If
ENSDFSTATE is not found, secular-equilibrium results (`t = -1`) still work
correctly, but finite-time queries fall back to SE with a warning. So
make sure `$G4ENSDFSTATEDATA` resolves on your machine — either via the env
var, via your geant4.sh, or via `GEANT4_DATA_DIR/G4ENSDFSTATE*`.

## Diagnosing problems

If `g4gamma` returns all zeros, run the included diagnostic:

```bash
python3 test/diagnose.py
```

It prints what it found at every step (env vars, resolved directories,
per-isotope file existence, the chain it built, and the spectrum it
produced) so the failing step is obvious.

You can also enable verbose output on any builder call:

```python
res = g.build_spectrum(g.IsotopeKey(55, 137, 0), -1, edges, verbose=1)
# Prints to stderr: data dirs found, chain members and their mean lives
# (with ENSDFSTATE patches applied), warnings about missing files, etc.
```

## Python API

```python
import numpy as np
import g4gamma as g

# Bin edges in user units, multiplied by the g.units constants
edges = np.linspace(0, 3000, 3001) * g.units.keV

# Cs-137 at secular equilibrium
res = g.build_spectrum(
    primary=g.IsotopeKey(55, 137, 0),
    t=-1,                           # SE; or e.g. 60 * g.units.s for finite t
    bin_edges=edges,
    include_annihilation=True,
    include_xrays=False,
)

counts = np.array(res.counts)        # gammas per primary decay, per bin
edges  = np.array(res.bin_edges)     # internal units (MeV)

# Chain breakdown (diagnostic)
for c in res.contributions:
    print(c.isotope, c.activity, c.gamma_yield)
```

For repeated calls (e.g. scanning many isotopes), use the persistent builder
which caches level data:

```python
opts = g.SpectrumOptions()
opts.include_xrays = True
opts.geant4_sh = "/opt/geant4/bin/geant4.sh"  # only needed if env vars not set
builder = g.GammaSpectrumBuilder(opts)

for (Z, A, M) in isotopes:
    res = builder.build(g.IsotopeKey(Z, A, M), -1, edges)
    ...
```

## C++ CLI

```bash
./build/g4gamma_cli <Z> <A> <M> <time_s_or_-1> <Emin_keV> <Emax_keV> <nbins>

# K-40, secular equilibrium, 0-3000 keV in 1-keV bins
./build/g4gamma_cli 19 40 0 -1 0 3000 3000
```

## Validation

Reference tests against ENSDF data:

| isotope | peak (keV) | expected | g4gamma |
|---------|-----------|----------|---------|
| Cs-137  | 661.659   | 0.853    | 0.8524  |
| K-40    | 1460.820  | 0.1055\* | 0.1055  |

\* total EC branch ratio in the synthetic test data; the literature value is
0.1066, depending on which evaluation you use.

To validate against your own rdecay01 outputs, see
`test/validate_against_rdecay01.py`.

## Known limitations

1. **K-shell X-rays only**. L/M X-rays are not produced. For NORM gamma
   spectroscopy at NaI/LaBr3 resolution this is fine — high-Z parents have
   K X-rays at 70-100 keV which is what you actually see. For low-Z parents
   K X-rays are below 10 keV which is below detection threshold anyway.
2. **Auger electrons not produced.** They're not gammas anyway. Means the
   K-shell yield emitted here equals ω_K × N_K_vacancies (correct, because
   the fluor file probabilities already include the fluorescence yield).
3. **Isomer M assignment** is by closest excitation match (1 keV tolerance).
   Edge cases involving floating-level conventions in ENSDFSTATE are not
   read directly. Validate against rdecay01 for any new isotope you depend on.
4. **The 72-character line-length convention** for distinguishing 3-column vs
   5-column records in the radioactive decay files is honoured, matching the
   Geant4 parser exactly.
5. **Dead-end levels** (no transitions) drop their probability silently.
   Geant4's behaviour at top-of-cascade orphan levels is similar but not
   identical; cross-validation against rdecay01 recommended.

## Files

```
include/g4gamma/
    Units.hh          Geant4/CLHEP system of units (energy=MeV, time=ns)
    IsotopeKey.hh     (Z,A,M) tuple
    DataPath.hh       env-var / geant4.sh / GEANT4_DATA_DIR resolution
    DecayData.hh      RadioactiveDecay zZ.aA file parser
    PhotonEvap.hh     PhotonEvaporation zZ.aA file parser
    FluorData.hh      G4LEDATA/fluor/fl-tr-pr-Z.dat parser
    EnsdfState.hh     ENSDFSTATE.dat parser (canonical mean lives)
    ChainBuilder.hh   recursive daughter walker, isomer M assignment
    Bateman.hh        analytic generalised Bateman with branching
    GammaSpectrum.hh  main API
src/                  matching .cc files
python/bindings.cc    pybind11 module
test/run_tests.py     self-contained synthetic-data test suite (11 cases)
test/diagnose.py      runtime diagnostic for path/file/format issues
test/validate_against_rdecay01.py   compare against existing rdecay01 CSVs
TRANSCRIPT.md         development notes & design transcript
```
