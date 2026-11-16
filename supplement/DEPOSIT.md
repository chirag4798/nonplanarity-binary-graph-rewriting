# Depositing the supplement (to do by the authors)

1. Build: `python3 scripts/make_supplement.py` (deterministic; writes `supplement/oligon_supplement_v1.tar.gz` and `supplement/SHA256SUMS`).
2. Choose a licence and fill the author fields in `CITATION.cff`.
3. Upload `oligon_supplement_v1.tar.gz` to a persistent repository (for example Zenodo or OSF), copy the minted DOI into the paper's
   "Data and code availability" section, and keep the SHA-256 of the archive (see `SHA256SUMS`) in Appendix A.
4. Quick check of the headline claim from the archive: `tar xzf oligon_supplement_v1.tar.gz && python3 scripts/verify_certificates.py`
   (Python 3 standard library only; reads `results/necessity/certificates.jsonl` and `necessity_summary.csv`).

## Status of this repository
Public GitHub copy of the supplement tree (code, scripts, results; logs and local paths removed). No licence has been chosen yet:
until the authors add one, all rights are reserved by default. The archive tarball itself (SHA-256 in SHA256SUMS) is not stored here;
rebuild it with `python3 scripts/make_supplement.py`.
