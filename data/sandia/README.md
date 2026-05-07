# SandiaDecay data

This directory contains `sandia.decay.nocoinc.min.xml.gz`, a redistribution of
the nuclide decay database from the Sandia National Laboratories
[SandiaDecay](https://github.com/sandialabs/SandiaDecay) project.

The file is a minified, no-coincidence-info variant (≈6 MB vs ≈30 MB for the
full file). It contains everything `g4gamma`'s SandiaProvider needs:
parent half-lives, decay branches with branching ratios, and per-decay
gamma + X-ray emissions (energy + intensity).

For the full file with gamma-coincidence information (relevant for ML
applications and detector simulations modelling pile-up), download
`sandia.decay.xml` from upstream:
<https://github.com/sandialabs/SandiaDecay/blob/master/sandia.decay.xml>

## License

LGPL-2.1 — see `SANDIA_LICENSE.txt`. Original copyright notice:

> Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
> (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
> Government retains certain rights in this software.

Underlying nuclear data is sourced from ENSDF (Evaluated Nuclear Structure
Data File, NNDC/Brookhaven) and the LBNL Table of Radioactive Isotopes.
