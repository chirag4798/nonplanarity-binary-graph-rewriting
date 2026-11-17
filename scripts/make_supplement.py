#!/usr/bin/env python3
"""Build a deterministic supplement archive (code, scripts, results, README, requirements; NOT the paper sources) and write SHA-256 sums.
Deterministic: sorted members, fixed mtime/uid/gid, gzip header mtime 0.  Output: supplement/oligon_supplement_v1.tar.gz and supplement/SHA256SUMS."""
import gzip, hashlib, io, os, tarfile
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
os.chdir(root)
os.makedirs("supplement", exist_ok=True)
include = ["code", "scripts", "results", "README.md", "LICENSE", "LICENSE-DATA", "CITATION.cff", ".zenodo.json"]
files = []
for top in include:
    if os.path.isfile(top):
        files.append(top)
        continue
    for d, dn, fn in os.walk(top):
        dn[:] = sorted(x for x in dn if x != "__pycache__")
        for f in sorted(fn):
            if f in (".DS_Store",) or f.endswith(".pyc") or f.endswith(".log"):   # logs contain local paths and interpreter warnings
                continue
            files.append(os.path.join(d, f))
files = sorted(files)
buf = io.BytesIO()
with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0, compresslevel=9) as gz:
    with tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as tar:
        for f in files:
            ti = tar.gettarinfo(f, arcname=f)
            ti.mtime = 0; ti.uid = ti.gid = 0; ti.uname = ti.gname = ""
            ti.mode = 0o755 if ti.mode & 0o111 else 0o644
            with open(f, "rb") as fh:
                tar.addfile(ti, fh)
data = buf.getvalue()
open("supplement/oligon_supplement_v1.tar.gz", "wb").write(data)
key = ["results/necessity/certificates.jsonl", "results/necessity/necessity_summary.csv", "results/classification/classification.csv",
       "results/census/rule_census.csv", "results/census_2_4/rule_census.csv", "results/rules_joinable.txt", "results/rules_selected.txt",
       "results/wolfram_rules/selfloop_evolution.csv", "results/wolfram_rules/noninjective_check.txt", "results/incidence_check.txt", "results/selftest_log.txt", "results/verify_independent_log.txt", "results/counterexamples_log.txt",
       "code/oligon_experiment.py", "code/general_rules.py", "scripts/necessity_certificates.py", "scripts/verify_certificates.py",
       "scripts/verify_independent.py", "scripts/general_signature_search.py"]
lines = ["%s  %s" % (hashlib.sha256(open(k, "rb").read()).hexdigest(), k) for k in key if os.path.exists(k)]
lines.append("%s  supplement/oligon_supplement_v1.tar.gz (%d files, %d bytes)" % (hashlib.sha256(data).hexdigest(), len(files), len(data)))
open("supplement/SHA256SUMS", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
