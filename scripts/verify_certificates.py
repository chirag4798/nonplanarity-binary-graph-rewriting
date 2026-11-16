#!/usr/bin/env python3
"""Standalone checker (Python standard library only; no networkx, no code shared with the search) for the certificate file
results/necessity/certificates.jsonl and its coverage of the rule list results/necessity/necessity_summary.csv.

For the statement of Theorem 8, (iii)=>(i) in contrapositive, it checks that
  (1) every certificate uses the host H0 = simple K5 with every edge (a,b), a<b, oriented from lower to higher index ("A5");
  (2) re-applying the recorded single event (rule, match) to that host reproduces the recorded output edge list;
  (3) the recorded rotation system is a valid planar embedding of the simple projection of the output (faces traced, Euler's formula);
  (4) the set of certified rules equals the set of rules that this checker, with its own brute-force decision procedure, finds
      NOT simple-minor-preserving, each rule being certified exactly once, and the summary lists every rule once;
  (5) no simple-minor-preserving rule has a certificate.
Exit status 1 on any failure."""
import csv, itertools, json, os, re, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
NEC = os.path.join(HERE, "..", "results", "necessity")


def parse_rule(s):
    lhs, rhs = s.replace(" ", "").split("->")
    pairs = lambda t: [(int(a) - 1, int(b) - 1) for a, b in re.findall(r"\{(\d+),(\d+)\}", t)]
    L, K = pairs(lhs), pairs(rhs)
    n_old = len({v for e in L for v in e}); n_all = len({v for e in K for v in e} | {v for e in L for v in e})
    return L, K, n_old, n_all - n_old


def host_A5():
    return [(a, b) for a in range(5) for b in range(a + 1, 5)]        # all edges directed low -> high


def apply_event(host, L, K, n_old, n_new, bind):
    vals = [bind[v] for v in range(n_old)]
    if len(set(vals)) != n_old or not all(0 <= x < 5 for x in vals):
        return None
    left, used = [], set()
    for (p, q) in L:                                                 # consume one host edge per left-side edge, same orientation
        want = (bind[p], bind[q])
        idx = next((i for i, e in enumerate(host) if e == want and i not in used), None)
        if idx is None:
            return None
        used.add(idx)
    m = dict(bind)
    for k in range(n_old, n_old + n_new):
        m[k] = 5 + (k - n_old)
    return [e for i, e in enumerate(host) if i not in used] + [(m[u], m[v]) for (u, v) in K]


def valid_rotation_system(edges_after, rotation):
    adj = {}
    for x, y in edges_after:
        if x != y:
            adj.setdefault(x, set()).add(y); adj.setdefault(y, set()).add(x)
    rot = {int(k): v for k, v in rotation.items()}
    if set(rot) != set(adj) or any(sorted(rot[v]) != sorted(adj[v]) for v in adj):
        return False
    pos = {u: {w: i for i, w in enumerate(rot[u])} for u in rot}
    seen, faces = set(), 0
    for u in adj:
        for v in adj[u]:
            if (u, v) in seen:
                continue
            faces += 1; x, y = u, v
            while (x, y) not in seen:
                seen.add((x, y)); lst = rot[y]; z = lst[(pos[y][x] + 1) % len(lst)]; x, y = y, z
    m = sum(len(s) for s in adj.values()) // 2
    comp, todo = 0, set(adj)
    while todo:
        comp += 1; st = [todo.pop()]
        while st:
            w = st.pop()
            for z in adj[w]:
                if z in todo:
                    todo.remove(z); st.append(z)
    return len(adj) - m + faces == 1 + comp


def simple_minor_preserving(L, K, n_old, n_new):
    """Brute force over branch assignments (own implementation)."""
    need = {frozenset(e) for e in L}
    for assign in itertools.product(range(-1, n_old), repeat=n_new):
        g = list(range(n_old)) + list(assign)
        ok = True
        for a in range(n_old):
            S = [v for v in range(len(g)) if g[v] == a]
            seen, todo = {S[0]}, [S[0]]
            while todo:
                v = todo.pop()
                for (x, y) in K:
                    for (p, q) in ((x, y), (y, x)):
                        if p == v and q in S and q not in seen:
                            seen.add(q); todo.append(q)
            if len(seen) != len(S):
                ok = False; break
        if not ok:
            continue
        cov = {frozenset((g[x], g[y])) for (x, y) in K if g[x] != -1 and g[y] != -1 and g[x] != g[y]}
        if need <= cov:
            return True
    return False


def main():
    summary = list(csv.DictReader(open(os.path.join(NEC, "necessity_summary.csv"))))
    rules = [r["rule"] for r in summary]
    certs = [json.loads(l) for l in open(os.path.join(NEC, "certificates.jsonl"))]
    problems = []
    if len(set(rules)) != len(rules):
        problems.append("summary lists a rule more than once")
    cert_rules = Counter(c["rule"] for c in certs)
    if any(v != 1 for v in cert_rules.values()):
        problems.append("some rule has more than one certificate")
    hosts = Counter(c["H0"] for c in certs)
    if set(hosts) != {"A5"}:
        problems.append("certificates on hosts other than A5: %s" % dict(hosts))
    ok_event = ok_rot = 0
    for c in certs:
        L, K, n_old, n_new = parse_rule(c["rule"])
        new = apply_event(host_A5(), L, K, n_old, n_new, {int(k): v for k, v in c["match"].items()})
        if new is not None and Counter(new) == Counter(map(tuple, c["edges_after"])):
            ok_event += 1
        else:
            problems.append("event mismatch: " + c["rule"])
        if valid_rotation_system(c["edges_after"], c["rotation"]):
            ok_rot += 1
        else:
            problems.append("invalid rotation system: " + c["rule"])
    non_sm = set()
    for r in rules:
        L, K, n_old, n_new = parse_rule(r)
        if not simple_minor_preserving(L, K, n_old, n_new):
            non_sm.add(r)
    certified = set(cert_rules)
    if certified != non_sm:
        problems.append("certified rules != non-simple-minor-preserving rules (%d vs %d; %d uncertified, %d certified but simple-minor-preserving)"
                        % (len(certified), len(non_sm), len(non_sm - certified), len(certified - non_sm)))
    # consistency with the summary file flags
    flag_bad = sum(1 for r in summary if (r["simple_minor_preserving"] == "0") != (r["rule"] in non_sm))
    if flag_bad:
        problems.append("%d summary flags disagree with the checker's own decision" % flag_bad)
    n = len(certs)
    print("certificates: %d verified: %d failed: %d %s" % (n, min(ok_event, ok_rot) if not problems else n - len(problems), len(problems), dict(hosts)))
    print("rules in summary: %d | not simple-minor-preserving (checker's own decision): %d | certified rules: %d | sets equal: %s" % (len(rules), len(non_sm), len(certified), certified == non_sm))
    print("hosts used by certificates: %s (A5 = simple K5, edges directed from lower to higher vertex index)" % dict(hosts))
    print("events re-applied and matching: %d/%d | rotation systems valid: %d/%d | summary flags consistent: %s" % (ok_event, n, ok_rot, n, flag_bad == 0))
    for p in problems[:10]:
        print("PROBLEM:", p)
    print("CERTIFICATE VERIFICATION", "PASSED" if not problems else "FAILED")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
