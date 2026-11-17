# Archiving (Zenodo)

Zenodo (https://zenodo.org) is a free open repository run by CERN. It stores a permanent snapshot of this repository and gives it a DOI,
which is what the paper cites. Metadata come from `.zenodo.json` and `CITATION.cff`.

1. Sign in at https://zenodo.org with GitHub; open https://zenodo.org/account/settings/github/ and switch the repository
   `nonplanarity-binary-graph-rewriting` ON.
2. Publish a GitHub release (tag `v1.0.0`). Zenodo then archives that snapshot and mints two DOIs: one for this version and a
   "concept" DOI that always points to the latest version.
3. Put the DOI into `paper/doi.tex` (`\newcommand{\zenododoi}{10.5281/zenodo.NNNNNNN}`) and rebuild the paper.
The deterministic tarball built by `python3 scripts/make_supplement.py` (SHA-256 in `SHA256SUMS`) is the same tree as the repository.
