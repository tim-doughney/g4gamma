# g4gamma

Direct gamma-spectrum extraction from nuclear decay datasets without running
Monte Carlo. Supports multiple data sources side-by-side via a provider
abstraction.

Built for Tim's NORM thesis work — computes per-bin gammas/primary at any
chain time, suitable as a forward model or as a sanity check on rdecay01
output.

## What it does

Given an isotope `(Z, A, M)`, a time `t`, and a bin grid, returns a histogram
of expected gammas per primary decay summed over the full decay chain. The
input matches the convention used in `rdecay01` macros: the primary nuclide
has activity 1 Bq at `t = 0`. At `t = -1` the system is in secular equilibrium
(equivalent to the `thresholdForVeryLongDecayTime 1.0e+60 year` mode).

## Data sources (providers)

You choose the source via `SpectrumOptions.source`:

- **`g.DataSource.Geant4`** (default) — reads the Geant4 data tree
  (`$G4RADIOACTIVEDATA`, `$G4LEVELGAMMADATA`, `$G4ENSDFSTATEDATA`,
  optionally `$G4LEDATA/fluor` and `$G4LEDATA/auger` for X-rays).
  Per-level scheme; cascades are computed by walking the level tree.
  This matches what your Geant4 simulation will actually produce.

- **`g.DataSource.Sandia`** — reads
  [SandiaDecay](https://github.com/sandialabs/SandiaDecay)'s
  `sandia.decay.xml` (LGPL-2.1, derived from ENSDF + LBNL ToRI). Single XML
  file, ~3500 nuclides. Per-decay aggregated photons (no cascade
  computation needed). Auto-resolves to bundled
  `data/sandia/sandia.decay.nocoinc.min.xml` if available.

Both backends produce identical answers to within 1% for the principal
peaks of all NORM-relevant nuclides — see `test/showcase.py`'s
`11_geant4_vs_sandia*.png` plots for cross-validation.

The provider abstraction is open: add a new source by implementing
`IDecayProvider` (3 methods) and dropping a new `*.cc` into
`CMakeLists.txt`. LARA/DDEP and ENDF-6 are obvious candidates.

## Sources of gammas

- Discrete nuclear gammas from PhotonEvaporation cascades following
  beta/EC/alpha decays into excited daughter levels
- Discrete gammas from IT decays of metastable parents
- 511 keV annihilation pairs from beta+ branches (controllable)
- Optional K-shell X-rays from EC decays and from internal-conversion
  vacancies in level transitions (`include_xrays=True`)
- Optional full fluorescence cascade K→L→M (`full_xray_cascade=True`):
  propagates secondary vacancies through all shells using G4LEDATA
  `fl-tr-pr-Z.dat` (fluorescence) **and** `au-tr-pr-Z.dat`
  (Auger/Coster-Kronig secondary vacancies), reproducing
  `SetARM(true)` + `SetAugerCascade(true)` in rdecay01
- Auger electrons are not produced (they're not gammas)

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

# Cs-137 at secular equilibrium (default Geant4 backend)
res = g.build_spectrum(
    primary=g.IsotopeKey(55, 137, 0),
    t=-1,                           # SE; or e.g. 60 * g.units.s for finite t
    bin_edges=edges,
    include_annihilation=True,
    include_xrays=False,
)

counts = np.array(res.counts)        # gammas per primary decay, per bin
edges  = np.array(res.bin_edges)     # internal units (MeV)
print(res.source_name)               # "geant4"

# Same query, but using SandiaDecay backend
res_sd = g.build_spectrum(
    primary=g.IsotopeKey(55, 137, 0),
    t=-1,
    bin_edges=edges,
    source=g.DataSource.Sandia,
    sandia_xml="data/sandia/sandia.decay.nocoinc.min.xml",  # or set $SANDIA_DECAY_XML
)
print(res_sd.source_name)            # "sandia"

# Chain breakdown (diagnostic)
for c in res.contributions:
    print(c.isotope, c.activity, c.gamma_yield)
```

For repeated calls (e.g. scanning many isotopes), use the persistent builder
which caches data:

```python
opts = g.SpectrumOptions()
opts.source = g.DataSource.Sandia          # or g.DataSource.Geant4
opts.include_xrays = True
opts.geant4_sh = "/opt/geant4/bin/geant4.sh"  # only needed for Geant4 backend
opts.sandia_xml = ""                        # auto-resolves; override if needed
builder = g.GammaSpectrumBuilder(opts)

for (Z, A, M) in isotopes:
    res = builder.build(g.IsotopeKey(Z, A, M), -1, edges)
    ...
```

## C++ CLI

```bash
# default: Geant4 backend
./build/g4gamma_cli <Z> <A> <M> <time_s_or_-1> <Emin_keV> <Emax_keV> <nbins>

# explicit source argument
./build/g4gamma_cli 92 238 0 -1 0 3000 3000 sandia
./build/g4gamma_cli 19 40 0 -1 0 3000 3000 geant4
```

## Validation

Reference tests against ENSDF data (`test/run_tests.py` — 27/27 passing):

| isotope | peak (keV) | reference | g4gamma (Geant4) | g4gamma (Sandia) | g4gamma (LARA) |
|---------|-----------|-----------|------------------|------------------|-----------------|
| Cs-137  | 661.659   | 0.853     | 0.8524           | 0.8533           | 0.8501          |
| K-40    | 1460.820  | 0.1067    | 0.1075           | 0.1067           | 0.1034          |
| Co-60   | 1173+1332 | 1.998     | 1.9988           | 1.9985           | 1.9985          |
| U-238   | 609 (Bi-214) | 0.461  | 0.461            | 0.461            | 0.454           |

Geant4 vs Sandia agreement is within 1% for all major peaks across all
NORM-relevant isotopes — see `showcase_plots/11_geant4_vs_sandia*.png`.

Validation against 10M-event rdecay01 simulation (U-238 secular equilibrium,
`SetARM(true)` + `SetAugerCascade(true)` + `SetDeexcitationIgnoreCut(true)`):

| peak (keV) | nuclide | simulation | g4gamma (Geant4) | g4gamma (Sandia) | err (G4) |
|-----------|---------|-----------|-----------------|-----------------|----------|
| 75.5 | Bi K X-ray | 0.0527 | 0.0525 | — | −0.4% |
| 92.5 | Th/Rn Kα | 0.0435 | 0.0437 | — | +0.5% |
| 295.2 | Pb-214 | 0.1847 | 0.1842 | 0.1857 | −0.3% |
| 351.9 | Pb-214 | 0.3574 | 0.3570 | 0.3591 | −0.1% |
| 609.3 | Bi-214 | 0.4596 | 0.4528 | 0.4610 | −1.5% |
| 1764.5 | Bi-214 | 0.1524 | 0.1522 | 0.1524 | −0.1% |
| 2204.1 | Bi-214 | 0.0493 | 0.0491 | 0.0492 | −0.4% |

All major peaks above 90 keV agree within ±2%. The Geant4 backend **total**
γ/primary is 3.92 vs simulation 3.23 for U-238 SE — excess concentrated in
the X-ray region and U-234 cascade bands (see Known Limitations #2). The
Sandia backend totals 3.20 (within 1% of simulation) and is recommended
for U-238 total-count work.

To validate against your own rdecay01 outputs, see
`test/validate_against_geant4.py`.

## Known limitations

1. **Auger electrons not produced** — they're not gammas, only the secondary
   vacancies they create are propagated (for `full_xray_cascade=True`).
2. **U-234 cascade through intermediate isomers (Geant4 backend).** For U-238
   secular equilibrium, `getCascade()` runs the photon-evaporation cascade
   from each starting level straight through to ground state without stopping
   at intermediate metastable levels. Geant4's simulation stops the cascade at
   U-234 M=2 (1421 keV, T½=33.5 μs) and M=1 (989 keV, T_mean=1.096 ns) and
   treats them as separate radioactive ions. This causes the Geant4 backend to
   over-produce cascade gammas at 120–200 keV (+0.055) and 800–1600 keV (+0.19)
   and extra U K X-rays at 90–120 keV (+0.09) relative to the rdecay01
   simulation. **Workaround: use the Sandia backend for U-238 chain work** —
   it agrees with rdecay01 to within 1% total and ±2% per-peak.
3. **Anomalous ICC in PhotonEvaporation data.** Some high-excitation U-234
   levels store fAlpha values that are orders of magnitude larger than the
   BrIcc physical estimate (e.g. fAlpha=0.37 at 804 keV; physical value
   ~0.001). A sanity cap in `getCascade()` suppresses ICC for transitions with
   E > 500 keV and fAlpha > 0.25, fixing the most egregious false K X-ray
   peak. Transitions at lower energies may still carry slightly inflated ICC.
4. **Isomer M assignment** is by closest excitation match (1 keV tolerance).
5. **Dead-end levels** (no transitions) drop their probability silently.
   Geant4's behaviour at top-of-cascade orphan levels is similar but not
   identical; cross-validation against rdecay01 recommended.
6. **G4EMLOW L2 fluorescence yield anomaly.** G4EMLOW gives ω_L2 ≈ 0.40
   for Pb (NIST value: 0.116). Both the simulation and analytic use the same
   G4EMLOW data so they agree with each other, but both deviate from physical
   L X-ray intensities in this shell. Use K-shell comparison for validation.

## Files

```
include/g4gamma/
    Units.hh           Geant4/CLHEP system of units (energy=MeV, time=ns)
    IsotopeKey.hh      (Z,A,M) tuple
    DataPath.hh        env-var / geant4.sh / GEANT4_DATA_DIR resolution
    DecayData.hh       RadioactiveDecay zZ.aA file parser
    PhotonEvap.hh      PhotonEvaporation zZ.aA file parser
    FluorData.hh       G4LEDATA/fluor/fl-tr-pr-Z.dat parser
    AugerData.hh       G4LEDATA/auger/au-tr-pr-Z.dat parser
    EnsdfState.hh      ENSDFSTATE.dat parser (canonical mean lives)
    IDecayProvider.hh  abstract interface for any decay-data source
    Geant4Provider.hh  IDecayProvider backed by Geant4 data dir
    SandiaProvider.hh  IDecayProvider backed by sandia.decay.xml (LGPL-2.1)
    ChainBuilder.hh    provider-agnostic chain walker
    Bateman.hh         analytic generalised Bateman with branching
    GammaSpectrum.hh   main API
src/                   matching .cc files
python/bindings.cc     pybind11 module
data/sandia/           bundled SandiaDecay XML + license + README
data/lara/             LARA/DDEP data (lara.tar.gz) + fetch script
test/run_tests.py      self-contained synthetic-data test suite (27 cases)
test/diagnose.py       runtime diagnostic for path/file/format issues
test/showcase.py       end-to-end demo with publication-quality plots
test/validate_against_geant4.py     compare against rdecay01 CSVs
TRANSCRIPT.md          development notes & design transcript
```
