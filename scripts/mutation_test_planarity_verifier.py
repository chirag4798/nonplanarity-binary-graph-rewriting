#!/usr/bin/env python3
"""Mutation test for verify_planarity_certificates.py: corrupt one certificate in four different ways and confirm that the
checker rejects each corrupted file (and accepts the unchanged one).  Writes to a temporary file next to the results."""
import contextlib, io, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify_planarity_certificates as v

lines = [json.loads(l) for l in open(v.CERT)]
tmp = v.CERT + ".mutated"


def run(mut, name, expect):
    with open(tmp, "w") as f:
        for i, c in enumerate(lines):
            c = json.loads(json.dumps(c))
            if i == 5:
                mut(c)
            f.write(json.dumps(c) + "\n")
    orig = v.CERT
    v.CERT = tmp
    with contextlib.redirect_stdout(io.StringIO()):
        rc = v.main()
    v.CERT = orig
    print("%-28s exit %d %s" % (name, rc, "ok" if (rc == 1) == expect else "UNEXPECTED"))
    return (rc == 1) == expect


ok = True
ok &= run(lambda c: c["kuratowski_edges"].pop(), "drop Kuratowski edge", True)
ok &= run(lambda c: c["binding"].reverse(), "permute binding", True)
ok &= run(lambda c: c["output_edges"].pop(), "drop output edge", True)
ok &= run(lambda c: c["host_rotation"].update({"0": c["host_rotation"]["0"][::-1][1:] + c["host_rotation"]["0"][:1]}), "corrupt rotation system", True)
ok &= run(lambda c: None, "unchanged file", False)
os.remove(tmp)
print("MUTATION TEST PASSED" if ok else "MUTATION TEST FAILED")
sys.exit(0 if ok else 1)
