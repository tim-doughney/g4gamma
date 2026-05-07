# Development transcript — g4gamma

This file captures the design conversation and reasoning that produced
`g4gamma`. Original conversation between Tim Doughney (PhD student, NORM
gamma spectroscopy, U Adelaide / DSTG) and Claude (Anthropic), May 2026.

The repository was built specifically to support Tim's PhD thesis work on
computational modelling of gamma radiation from naturally occurring
radioactive materials, where the existing approach was running large
rdecay01 simulations to extract per-isotope gamma PDFs. This tool gets the
same information directly from the underlying Geant4 datasets, with no Monte
Carlo statistical noise and effectively zero runtime.

---

## Initial problem statement

> I use Geant4 11.3.0. Traditionally I have been using the rdecay example to
> get the full chain gamma spectrum of isotopes, like 137Cs or 238U. I used
> 1e9 counts and then normalised to number of events to get the pdf of the
> gamma spectra per bq of primary. I want to grab the underlying probability
> distribution without the time impost of simulating large numbers of events
> so the statistics are good.

User wants to short-circuit the rdecay01 simulation by reading directly from
the Geant4 RadioactiveDecay and PhotonEvaporation datasets, producing a
binned histogram of expected gammas per primary decay. Inputs: isotope of
interest, time past which Bq was specified (with -1 meaning secular
equilibrium), energy bin boundaries. Output: a histogram. Must be CMake-built
with a pybind11 wrapper.

User provided their `PhysicsList.cc` and rdecay01 macro for Geant4 11.3.0.
Key configuration choices to mirror:
- `SetARM(false)` — no atomic rearrangement, so no X-rays from EC/IC by
  default
- `SetCorrelatedGamma(false)` — gammas sampled independently from level
  intensities (matches the per-level data in PhotonEvaporation files)
- `SetIsomerProduction(true)` with `MaxLifeTime(1*ns)` — metastable states
  with mean life > 1 ns treated as separate isotopes (Tc99m, Ba137m, etc.)
- `thresholdForVeryLongDecayTime 1.0e+60 year` — all decays in the chain
  happen regardless of half-life (this is what `t = -1` corresponds to in
  this tool)

User's geant4.sh has `GEANT4_DATA_DIR` set but the per-dataset env vars
(`G4RADIOACTIVEDATA`, `G4LEVELGAMMADATA`, `G4LEDATA`) commented out — so
path resolution has to discover the versioned subdirectories.

User confirmed:
- Include 511 keV annihilation pairs from B+ branches
- Optional X-ray inclusion (decided: K-shell only, see below)
- Hand-rolled generalised Bateman without Eigen dependency
- Inputs in Geant4 internal units, with unit constants exposed

---

## Architecture decisions

### Source of gamma lines

Two sources were identified:

1. **Discrete nuclear gammas** — born when a daughter nucleus is created in
   an excited state (after beta/EC/alpha decay), or when the parent nuclide
   is itself a metastable state that IT-decays. The transition data
   (gamma energy, internal-conversion coefficient α, branching among
   transitions) lives in `PhotonEvaporation/zZ.aA` files. The probability of
   gamma emission per de-excitation cascade was computed via dynamic
   programming: maintain `p[i] = probability of being in level i during the
   cascade`, distribute downwards, accumulate weighted gamma yields.

2. **511 keV annihilation pairs** — 2 photons per B+ decay. Cross-checked
   with Na-22 reference (90.4% B+ → 1.808 photons per decay).

### Sources of X-rays (when enabled)

Two pathways:

1. **EC decays** create a vacancy in the K, L, M, or N shell of the
   *daughter* atom (since EC removes an inner electron). Each shell-EC mode
   appears separately in the decay-data files (KshellEC, LshellEC, etc.) so
   we know directly which shell got the vacancy.

2. **Internal conversion** in level transitions creates a vacancy in the
   *parent* atom. The PhotonEvaporation file gives per-transition ICC values
   in 10 columns (K, L1, L2, L3, M1, M2, M3, M4, M5, N+) — these tell us the
   shell distribution of the IC vacancy.

X-ray fluorescence from a vacancy in shell s is computed using the
G4LEDATA/fluor/fl-tr-pr-Z.dat file (per-element, per-vacancy-shell
fluorescence transitions with energies and probabilities).

**Restriction: K-shell only.** L and M X-rays are below the meaningful
detection threshold for NaI and LaBr3 detectors at NORM-relevant energies.
For high-Z parents (Pb, Th, U) K X-rays are 70–100 keV — visible, important
peaks. Auger electrons are not produced at all (they're not gammas anyway).
This matches both the practical detector physics and the user's
`SetARM(false)` configuration.

### Decay chain handling

`ChainBuilder` walks daughters from the input nuclide via BFS, with the
following M (isomer index) assignment logic:

- For each decay channel of a parent (Z,A,M_parent), the channel specifies
  a daughter excitation E* in keV (from the 5-column line in the decay
  file). The daughter nuclide (Z',A') has its own decay-data file with one
  or more `P` blocks at different excitation energies. We match E* to the
  closest `P`-block parentExcitation; if it matches a non-ground block
  within 1 keV tolerance, the daughter is assigned that block's M index
  (1 for first metastable, 2 for second, etc.). Otherwise M=0.

This catches Ba-137m (Cs-137 daughter), Pa-234m (Th-234 daughter), Tc-99m,
etc. — long-lived metastable states that are tracked separately in the
Geant4 nuclide table.

### Bateman with branching

The user explicitly wanted no Eigen dependency. Solved with a hand-rolled
generalised Bateman:

For each chain member i (in topological order from root), the time-dependent
atom count is

    N_i(t) = sum_j n_{i,j} exp(-lambda_j t)

with the recursion

    n_{i,j} (j != i) = sum_{p in parents(i)} BR(p->i) * lambda_p * n_{p,j}
                    / (lambda_i - lambda_j)
    n_{i,i} = -sum_{j != i} n_{i,j}        (so N_i(0) = 0)

Activity is A_i(t) = lambda_i * N_i(t). Repeated eigenvalues (which never
arise in practice for real decay chains) are guarded against with a tiny
relative perturbation.

For secular equilibrium (t = -1), we walk the chain topologically and
propagate A_i = sum_p A_p * BR(p->i), with A_root = 1.

**Critical bug found during user testing**: The half-life column in the
RadioactiveDecay zZ.aA files turned out to be a placeholder — Geant4 reads
the actual values from `ENSDFSTATE.dat`. Without that, every chain member
got `meanLife = 0` and the SE solver was zeroing out anything with
`lambda == 0`, producing all-zero output. Fixed by:

1. Distinguishing "stable" (no decay channels) from "unknown lambda" in the
   SE solver — stability comes from channels, not lambda.
2. Adding an `EnsdfStateLoader` that reads the canonical mean lives from
   `$G4ENSDFSTATEDATA/ENSDFSTATE.dat`. This is the same file Geant4 itself
   reads via `G4NuclideTable`. Mean lives are patched into chain nodes
   before the Bateman solver runs.
3. For finite-time evolution when lambda is still unknown after the patch
   (e.g. ENSDFSTATE not available), falling back to secular equilibrium
   rather than returning zero.

### Path resolution

User's geant4.sh has the per-dataset env vars commented out, so resolution
falls through to:

1. Try direct env var (G4RADIOACTIVEDATA, etc.)
2. Parse `geant4.sh` (default `/opt/geant4/bin/geant4.sh`, or user-supplied)
3. Read `GEANT4_DATA_DIR` from the parsed sh, find versioned subdirectory

The sh parser handles backquote and `$()` constructs of the form
`cd ... > /dev/null ; pwd`, since that's what Geant4's installer scripts
generate. It also resolves `dirname $BASH_SOURCE` to the actual script
directory.

---

## Validation strategy

Synthetic test data files were constructed for:

- **Cs-137** (single daughter Ba-137m, 661.7 keV peak, BR=94.7%, gamma
  emission probability 1/(1+α)=0.901). Reference: 0.853 gammas/decay.
- **K-40** (mixed B-/EC, 1460.8 keV peak from EC branches, total intensity
  10.55-10.66% per literature). Synthetic: 0.1055.
- **Na-22** (B+ to Ne-22 1274.5 keV level, plus annihilation pairs).
  Reference: 1.807 (511 keV from annihilation) + 0.999 (1274.5 keV) = 2.807
  per primary decay.
- **Synthetic U-238 mini-chain** with branched topology (Pa-234m IT to
  Pa-234 in addition to its main β- branch) to exercise the branched-chain
  Bateman.
- **Cs-137 in real Geant4 file format** added after first round of testing
  on the user's actual install — see "Real-data debugging" below.

All synthetic-data tests pass to 6+ decimal places. Bin-edge cases tested:
out-of-range peaks (no spurious counts), log binning, fine binning,
non-uniform bins.

---

## Real-data debugging

After getting `make && python test/run_tests.py` to pass on the test rig, the
user ran the diagnostic against their real `RadioactiveDecay6.1.2` install
and got all-zero output. Several issues surfaced:

1. **String mode names, not integer codes.** The actual Geant4 11.3.0
   `RadioactiveDecay6.1.2` files use string mode names (`"BetaMinus"`,
   `"Alpha"`, `"MshellEC"`, `"IT"`, `"SpFission"`, etc.) in the decay
   channel records, not integer codes as I'd assumed from older versions.
   Geant4's parser uses an `operator>>` overload defined on
   `G4RadioactiveDecayMode` which dispatches on the string. Fixed by
   replacing `modeFromInt` with `modeFromStr` (with integer fallback for
   backward compatibility).

2. **Daughter file with no ground-state P block.** The Ba-137 file in real
   Geant4 has only `P 661.659 ... IT` (the metastable). There is no
   `P 0 ... ` block for the stable ground state. This breaks the file-order
   M assumption (`parents[0]` would be Ba-137m, not Ba-137 ground), so when
   ChainBuilder asked for `IsotopeKey(56, 137, M=1)` it found no match.

3. **Half-life column in the parent line is a placeholder.** Confirmed in
   real data: most parent lines have a number that looks like a half-life
   in seconds but Geant4 actually reads the canonical mean lives from
   `ENSDFSTATE.dat`. (As it happens, the values in the rad files are often
   correct, but this is not guaranteed and the parser does not depend on
   them.)

The fix for (2) and (3) was to make `ChainBuilder` ENSDFSTATE-aware:

- M assignment for daughters uses `EnsdfStateLoader::excitationToM()`,
  which walks the ENSDFSTATE level list for the daughter element and
  returns the M index (M=0 for ground, M=1 for first level above the
  user-configured lifetime threshold, etc.).
- When looking up the decay-data parent block for an `IsotopeKey(Z,A,M)`,
  we first ask ENSDFSTATE for the excitation corresponding to that M, then
  match it to the closest `P`-block in the rad file (within 1 keV).
  This correctly handles missing-ground-state files like Ba-137.
- Mean lives are pulled from ENSDFSTATE (canonical) rather than from the
  rad-file half-life column (placeholder).

Final test suite: 13/13 passing, including a regression test using the
real Geant4 file format and ENSDFSTATE.

---

## Files

```
include/g4gamma/
    Units.hh           Geant4/CLHEP units (energy=MeV, time=ns)
    IsotopeKey.hh      (Z, A, M) tuple
    DataPath.hh        env var / geant4.sh / GEANT4_DATA_DIR resolution
    DecayData.hh       RadioactiveDecay zZ.aA file parser
    PhotonEvap.hh      PhotonEvaporation zZ.aA file parser
    FluorData.hh       G4LEDATA/fluor/fl-tr-pr-Z.dat parser
    EnsdfState.hh      ENSDFSTATE.dat parser (canonical mean lives)
    ChainBuilder.hh    daughter walker, isomer M assignment
    Bateman.hh         analytic generalised Bateman with branching
    GammaSpectrum.hh   main API
src/                   matching .cc files
python/bindings.cc     pybind11 module
test/run_tests.py      self-contained synthetic-data test suite (11 cases)
test/diagnose.py       runtime diagnostic for path/file/format issues
test/validate_against_rdecay01.py   compare against existing rdecay01 CSVs
```

---

## Known limitations

1. **K-shell X-rays only**. L/M not produced (intentional — below NaI/LaBr3
   detection threshold for NORM-relevant cases).
2. **Auger electrons not produced** — they're not gammas anyway. The K-shell
   yield equals ω_K × N_K_vacancies, which is correct because the fluor
   probabilities already include the fluorescence yield.
3. **Isomer M assignment** is by closest excitation match (1 keV tolerance).
   Edge cases involving floating-level conventions in ENSDFSTATE are not
   currently read directly. Validate against rdecay01 for new isotopes.
4. **Top-of-cascade orphan levels** drop probability silently. Geant4's
   behaviour is similar but not identical; cross-validation against rdecay01
   recommended for any heavy nuclide where this matters.
5. **EC vs IC X-ray Z-attribution**: when both EC and IC contribute to the
   K-vacancy count for a decaying node, the daughter Z is used if any EC
   channel exists, else parent Z. Small approximation.

---

## Test results (final, all 13 passing)

```
[1] Cs-137 SE             OK  661 keV peak: 0.852385 = 0.947 × (1/1.111)
[2] Cs-137 finite-t       OK  Ba-137m at T½: 0.4735 = 0.947 × 0.5
[3] K-40 SE               OK  1460 keV peak: 0.1055
[4] K-40 with X-rays      OK  K-α at 2.957 keV, K-β at 3.190 keV
[5] Na-22 with B+         OK  511 keV: 1.8072 = 2 × 0.9036
[6] Na-22 no annihilation OK  511 peak suppressed cleanly
[7] Out-of-range bins     OK  no spurious counts
[8] Cs-137 real format    OK  661 keV peak with ENSDFSTATE-driven M assignment
[8b] Real format chain    OK  Ba-137m identified as M=1 isomer

Plus interactive validation:
- Synthetic U-238 mini-chain (branched topology with Pa-234m IT to Pa-234)
- geant4.sh parsing with commented-out per-dataset env vars
- ENSDFSTATE patching of placeholder half-lives in real Geant4 files
- String mode-name parsing with integer fallback for older datasets
```

---

## Multi-source provider architecture

After validating the Geant4 backend against the user's real install, a second
data source was added: **SandiaDecay** (Sandia Labs, LGPL-2.1, derived from
ENSDF + LBNL ToRI). The motivation was thesis-level cross-validation — a
single-source forward model is fragile to evaluation choices in any one
library, but agreement across two independent ENSDF-derived libraries gives
much stronger confidence.

### Provider abstraction

`IDecayProvider` was introduced with a high-level data model:

```
ParentDecayInfo {
    IsotopeKey, bool stable, double meanLife,
    vector<DecayBranch> branches
}
DecayBranch {
    DecayMode mode,
    double branchingRatio,
    IsotopeKey daughter,
    vector<Emission> emissions    // already aggregated photons per branch
}
Emission { type=Gamma|XRay|AnnihilationPair, energy, intensity }
```

The crucial insight was that Geant4 data is **per-level** (you compute the
cascade yourself by walking the level tree), while SandiaDecay is **per-decay**
(the XML directly tells you "this transition produces these gammas with these
intensities"). The provider abstraction hides this difference: each provider
exposes pre-aggregated photon emissions per branch. `Geant4Provider` runs the
cascade computation internally (moved out of `GammaSpectrum.cc`); 
`SandiaProvider` parses the XML and uses the data directly.

Result: `GammaSpectrum` and `ChainBuilder` are now provider-agnostic — they
work on `IDecayProvider` only and don't care where the data came from. Adding
a third backend (e.g. LARA, ENDF) is a ~150-line `*.cc` implementing 4
methods.

### SandiaDecay XML format gotchas

1. **Two formats coexist**: the verbose `sandia.decay.xml` (~31 MB) uses
   full tag/attribute names (`<nuclide>`, `atomicNumber`, etc.) while the
   minified `sandia.decay.nocoinc.min.xml` (~6 MB) uses abbreviated forms
   (`<n>`, `an`, `mn`, `iso`, `s`, `hl` etc.). `SandiaProvider` auto-detects
   which is in use by counting `<nuclide ` vs `<n ` occurrences and switches
   tag/attribute lookup tables accordingly.

2. **Attribute lookup via substring search is fragile** — if attribute
   `key` is a prefix of another attribute (`isomerNumber` contains `iso`),
   searching for `iso="..."` would match the wrong attribute. Fixed by
   requiring the previous character to be whitespace or `<`.

3. **Photon intensities are per-decay-going-down-this-branch**, not per
   primary decay. Total per-primary contribution is `branch.branchingRatio
   × emission.intensity`. Verified: K-40 → Ar-40 EC branch (BR=0.1086) with
   1460 keV gamma intensity 0.9825 gives 0.1067 per primary, matching DDEP
   reference.

4. **Decay mode strings differ from Geant4**: Sandia uses short codes
   (`b-`, `b+`, `it`, `ec`, `a`, `sf`, etc.). A separate `parseMode` table
   maps these. The IDecayProvider DecayMode enum was extended with a
   generic `EC` value to handle Sandia's shell-unspecified EC, alongside
   Geant4's K/L/M/N-shell-specific variants.

### Integration

- `SpectrumOptions.source` selects backend at builder construction
- `g.DataSource.Geant4` / `g.DataSource.Sandia` exposed via Python
- `result.source_name` returned for diagnostic purposes
- CLI accepts `[geant4|sandia]` as 8th argument
- Cross-validation: `test/showcase.py` produces side-by-side plots
  (`11_geant4_vs_sandia.png`) showing both backends on the same axes;
  agreement is < 1% at all major peaks > 1% intensity

### Bundled data

`data/sandia/sandia.decay.nocoinc.min.xml` (6.2 MB, LGPL-2.1) is bundled
with the project. License attribution and full LGPL text are in
`data/sandia/SANDIA_LICENSE.txt` and `data/sandia/README.md`.

`SandiaProvider::locateXml` resolution order:
1. constructor argument
2. `$SANDIA_DECAY_XML` env var
3. `<repo>/data/sandia/` (bundled)
4. `/usr/local/share/g4gamma/`

### Future providers worth adding

- **LARA / DDEP** (LNHB): ~220 evaluated isotopes, gold-standard for
  metrology. CSV format per-nuclide. Excellent for high-confidence
  validation on headline NORM isotopes. Recommended next addition.
- **ENDF-6 decay sublibrary** (ENDF/B, JEFF, JENDL): ~3800 nuclides,
  standardised format used by transport codes. The format is genuinely
  awful to parse and the data ultimately derives from ENSDF anyway, so
  this is lower priority.
- **NUBASE / GENF**: comprehensive but again ENSDF-derived; lower priority.

---

## Test results (final, all 21 passing)

```
[1] Cs-137 SE             OK  661 keV peak: 0.852385 = 0.947 × (1/1.111)
[2] Cs-137 finite-t       OK  Ba-137m at T½: 0.4735 = 0.947 × 0.5
[3] K-40 SE               OK  1460 keV peak: 0.1055
[4] K-40 with X-rays      OK  K-α at 2.957 keV, K-β at 3.190 keV
[5] Na-22 with B+         OK  511 keV: 1.8072 = 2 × 0.9036
[6] Na-22 no annihilation OK  511 peak suppressed cleanly
[7] Out-of-range bins     OK  no spurious counts
[8] Cs-137 real format    OK  661 keV peak with ENSDFSTATE-driven M assignment
[8b] Real format chain    OK  Ba-137m identified as M=1 isomer
[9] SandiaDecay backend (Cs-137, K-40, Co-60, U-238 chain)   8 sub-checks  OK

Cross-validation:
  Geant4 vs Sandia agreement at major peaks: < 1% across all NORM isotopes
  Cs-137 661 keV:       Geant4=0.8524  Sandia=0.8533  diff=0.1%
  K-40 1460 keV:        Geant4=0.1075  Sandia=0.1067  diff=0.7%
  Co-60 total:          Geant4=1.9988  Sandia=1.9985  diff=0.02%
  U-238 chain Bi-214 609 keV:  Sandia=0.461  (matches DDEP reference)
```
