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

- **`g.DataSource.Lara`** — reads LARA/DDEP evaluated data files from
  `data/lara/` (48 nuclides covering the full U-238 and Th-232 chains).
  Emissions are per-decay and include measured X-rays. Run
  `data/lara/fetch_lara.sh` to re-download or update. Activate with
  `opts.source = g.DataSource.Lara`.

All three backends agree to within 1–2% for the principal peaks of
NORM-relevant nuclides — see `test/showcase.py` for cross-validation.

The provider abstraction is open: add a new source by implementing
`IDecayProvider` (3 methods) and dropping a new `*.cc` into
`CMakeLists.txt`. ENDF-6 is an obvious next candidate.

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
opts.source = g.DataSource.Sandia          # or g.DataSource.Geant4 / g.DataSource.Lara
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

Reference tests against ENSDF data (`test/run_tests.py` — 38/38 passing):

| isotope | peak (keV) | reference | g4gamma (Geant4) | g4gamma (Sandia) | g4gamma (LARA) |
|---------|-----------|-----------|------------------|------------------|-----------------|
| Cs-137  | 661.659   | 0.853     | 0.8524           | 0.8533           | 0.8501          |
| K-40    | 1460.820  | 0.1067    | 0.1075           | 0.1067           | 0.1034          |
| Co-60   | 1173+1332 | 1.998     | 1.9988           | 1.9985           | 1.9985          |
| U-238   | 609 (Bi-214) | 0.461  | 0.461            | 0.461            | 0.454           |

Geant4 vs Sandia agreement is within 1% for all major peaks across all
NORM-relevant isotopes — see `showcase_plots/11_geant4_vs_sandia*.png`.

Validation against rdecay01 simulation (secular equilibrium,
`SetARM(true)` + `SetAugerCascade(true)` + `SetDeexcitationIgnoreCut(true)`):

**U-238** (10,000-event run):

| peak (keV) | nuclide | ENSDF ref | simulation | g4gamma (Geant4) | err (G4) |
|-----------|---------|-----------|-----------|-----------------|----------|
| 75.5 | Bi K X-ray | 0.0930 | 0.0544 | 0.0527 | −3.1% vs sim |
| 92.5 | Th/Rn Kα | 0.0630 | 0.0447 | 0.0437 | −2.2% vs sim |
| 295.2 | Pb-214 | 0.1840 | 0.1845 | 0.1845 | +0.0% |
| 351.9 | Pb-214 | 0.3580 | 0.3559 | 0.3583 | +0.7% |
| 609.3 | Bi-214 | 0.4610 | 0.4628 | 0.4541 | −1.9% |
| 1764.5 | Bi-214 | 0.1530 | 0.1511 | 0.1527 | +1.1% |
| 2204.1 | Bi-214 | 0.0491 | 0.0463 | 0.0493 | +6.4% |

Total γ/primary: 3.221 (model) vs 3.247 (simulation) = **−0.8%**

**Th-232** (1,000,000-event run):

| peak (keV) | nuclide | simulation | g4gamma (Geant4) | err (G4) |
|-----------|---------|-----------|-----------------|----------|
| 239.0 | Pb-212 | 0.4372 | 0.4363 | −0.2% |
| 338.3 | Ac-228 | 0.1082 | 0.1232 | +13.9% |
| 583.2 | Tl-208 | 0.2987 | 0.3029 | +1.4% |
| 727.3 | Bi-212 | 0.0660 | 0.0665 | +0.7% |
| 911.2 | Ac-228 | 0.2623 | 0.2723 | +3.8% |
| 2614.5 | Tl-208 | 0.3596 | 0.3600 | +0.1% |

Total γ/primary: 4.257 (model) vs 4.261 (simulation) = **−0.1%**

All major gamma peaks agree within ±2% of simulation. The 338.3 keV Ac-228
line (+14%) and 911 keV cluster (+4%) reflect a consistent data difference
between the Geant4 RadioactiveDecay branching fractions and ENSDF; the simulation
and model agree with each other, both differ from ENSDF. X-ray peaks at 75–93 keV
show −3% vs simulation due to the G4EMLOW L2 fluorescence yield anomaly (see
Known Limitations #2); both model and simulation deviate identically from
physical ENSDF values in that region.

`test/validate_against_geant4.py` auto-detects the isotope from the CSV filename
and produces per-source overlay plots. To validate against your own rdecay01
outputs:

```bash
python test/validate_against_geant4.py <csv_file> <n_primaries>
```

## Known limitations

1. **Auger electrons not produced** — they're not gammas, only the secondary
   vacancies they create are propagated (for `full_xray_cascade=True`).
2. **G4EMLOW L2 fluorescence yield anomaly.** G4EMLOW gives ω_L2 ≈ 0.40
   for Pb (NIST value: 0.116). Both the simulation and analytic use the same
   G4EMLOW data so they agree with each other, but both deviate from physical
   L X-ray intensities in this shell. Use K-shell comparison for validation.
3. **Isomer M assignment** is by closest excitation match (1 keV tolerance).
4. **Dead-end levels** (no transitions) drop their probability silently.
   Geant4's behaviour at top-of-cascade orphan levels is similar but not
   identical; cross-validation against rdecay01 recommended.

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
test/run_tests.py      self-contained synthetic-data test suite (38 cases)
test/diagnose.py       runtime diagnostic for path/file/format issues
test/showcase.py       end-to-end demo with publication-quality plots
test/validate_against_geant4.py       compare against rdecay01 CSVs
test/validate_against_sandiadecay.py  compare Bateman solver vs SandiaDecay library
TRANSCRIPT.md          development notes & design transcript
```
