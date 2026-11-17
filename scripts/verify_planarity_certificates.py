#!/usr/bin/env python3
"""Standalone checker (Python standard library only; no networkx; no code shared with the search) for the certificates
results/planarity/certificates.jsonl, which support the computer-assisted direction of the planarity-preservation theorem
(Theorem "Planarity-preserving rules on S").

Statement checked (iii)=>(i) in contrapositive: every rule of the rule list results/necessity/necessity_summary.csv (5225
closed rules) that violates condition P has a planar host and an event whose output has a non-planar simple graph.

For every certificate it checks:
  (1) the host H_s is planar: the recorded rotation system is valid (faces traced, Euler's formula n-m+f=1+c);
  (2) the binding sends the old variables injectively to vertices of H_s and every pair of the left side L is an edge of
      H_s (the host has each of its edges in four copies, both orientations twice, so every left side matches there);
  (3) re-applying the rule (right side with fresh vertices 1000,1001,...) to H_s reproduces the recorded output graph;
  (4) the recorded Kuratowski subgraph lies in the output and is a subdivision of K5 or K3,3 (degrees, path
      contraction, complete or complete-bipartite quotient), so the output is non-planar by Kuratowski's theorem.
Coverage: with its own implementation of condition P (below) the checker finds which rules of the list violate P;
the set of certified rules must equal that set, each rule once; every rule satisfying P must be simple-minor-preserving
under its own brute-force decision procedure (the second claim of the proposition on condition P).
Condition P: let G be the simple graph of K, O the old and N the new variables, pairs(L) the distinct pairs of L.
  (P1) every edge of G between two old variables is in pairs(L);
  (P2) every connected component C of G[N] has at most two old neighbours, and if two, they form a pair of pairs(L)
       and G[C + the two old vertices] + (edge between them) is planar.  The last clause is automatic here: n_new <= 2,
       so the graph has at most 4 vertices, and every graph on at most 4 vertices is planar (asserted below).
Exit status 1 on any failure."""
import csv, itertools, json, os, re, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
NEC = os.path.join(HERE, "..", "results", "necessity", "necessity_summary.csv")
CERT = os.path.join(HERE, "..", "results", "planarity", "certificates.jsonl")


def parse_rule(s):
    lhs, rhs = s.replace(" ", "").split("->")
    pairs = lambda t: [(int(a) - 1, int(b) - 1) for a, b in re.findall(r"\{(\d+),(\d+)\}", t)]
    L, K = pairs(lhs), pairs(rhs)
    n_old = len({v for e in L for v in e}); n_all = len({v for e in K for v in e} | {v for e in L for v in e})
    return L, K, n_old, n_all - n_old


def valid_rotation_system(edges, rotation):
    adj = {}
    for x, y in edges:
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


def is_kuratowski_subdivision(edges):
    """edges: set of frozenset pairs. True iff they form a subdivision of K5 or K3,3."""
    adj = {}
    for e in edges:
        a, b = tuple(e)
        adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    deg = {v: len(s) for v, s in adj.items()}
    branch = [v for v in adj if deg[v] >= 3]
    if any(d not in (2, 3, 4) for d in deg.values()):
        return False
    if sorted(deg[v] for v in branch) == [4] * 5:
        kind = "K5"
    elif sorted(deg[v] for v in branch) == [3] * 6:
        kind = "K33"
    else:
        return False
    # contract paths through degree-2 vertices
    quotient, used = set(), set()
    for b in branch:
        for first in adj[b]:
            prev, cur, length = b, first, 1
            path_edges = [frozenset((b, first))]
            while deg[cur] == 2:
                nxt = [w for w in adj[cur] if w != prev][0]
                path_edges.append(frozenset((cur, nxt)))
                prev, cur = cur, nxt
            if cur == b:
                return False
            quotient.add(frozenset((b, cur)))
            used |= set(path_edges)
    if used != set(edges):
        return False                 # stray cycle component or unvisited edges
    if kind == "K5":
        return len(quotient) == 10
    # K3,3: 9 distinct pairs forming a complete bipartite graph on the six branch vertices
    if len(quotient) != 9:
        return False
    g = {v: set() for v in branch}
    for e in quotient:
        a, b = tuple(e); g[a].add(b); g[b].add(a)
    a0 = branch[0]
    side_a = {a0} | {x for x in branch if x != a0 and x not in g[a0]}
    side_b = set(branch) - side_a
    return len(side_a) == 3 and len(side_b) == 3 and all(g[v] == side_b for v in side_a) and all(g[v] == side_a for v in side_b)


def satisfies_P(L, K, n_old, n_new):
    assert n_new <= 2
    PL = {frozenset(e) for e in L if e[0] != e[1]}
    G = {frozenset(e) for e in K if e[0] != e[1]}
    old = set(range(n_old)); new = list(range(n_old, n_old + n_new))
    for e in G:
        a, b = tuple(e)
        if a in old and b in old and e not in PL:
            return False
    comp = {v: v for v in new}

    def find(x):
        while comp[x] != x:
            x = comp[x]
        return x
    for e in G:
        a, b = tuple(e)
        if a in comp and b in comp:
            comp[find(a)] = find(b)
    att = {}
    for e in G:
        a, b = tuple(e)
        for x, y in ((a, b), (b, a)):
            if x in comp and y in old:
                att.setdefault(find(x), set()).add(y)
    for c, A in att.items():
        if len(A) > 2:
            return False
        if len(A) == 2 and frozenset(A) not in PL:
            return False
    return True


def simple_minor_preserving(L, K, n_old, n_new):
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
    fail = []
    rules = [r["rule"] for r in csv.DictReader(open(NEC))]
    if len(rules) != len(set(rules)):
        fail.append("rule list has duplicates")
    certs = [json.loads(l) for l in open(CERT)]
    cnt = Counter(c["rule"] for c in certs)
    ok_host = ok_bind = ok_out = ok_kur = 0
    for c in certs:
        L, K, n_old, n_new = parse_rule(c["rule"])
        H = [tuple(e) for e in c["host_edges"]]
        if valid_rotation_system(H, c["host_rotation"]):
            ok_host += 1
        else:
            fail.append("host not planar: " + c["rule"])
        sig = c["binding"]
        Hset = {frozenset(e) for e in H}
        if len(sig) == n_old and len(set(sig)) == n_old and all(frozenset((sig[a], sig[b])) in Hset for a, b in L):
            ok_bind += 1
        else:
            fail.append("binding invalid: " + c["rule"])
        m = {i: sig[i] for i in range(n_old)}
        for k in range(n_new):
            m[n_old + k] = 1000 + k
        out = set(Hset) | {frozenset((m[u], m[v])) for u, v in K if m[u] != m[v]}
        if out == {frozenset(e) for e in c["output_edges"]}:
            ok_out += 1
        else:
            fail.append("output mismatch: " + c["rule"])
        kur = {frozenset(e) for e in c["kuratowski_edges"]}
        if kur <= out and is_kuratowski_subdivision(kur):
            ok_kur += 1
        else:
            fail.append("Kuratowski subgraph invalid: " + c["rule"])
    notP, P = [], []
    for s in rules:
        L, K, n_old, n_new = parse_rule(s)
        (P if satisfies_P(L, K, n_old, n_new) else notP).append(s)
    if set(cnt) != set(notP) or any(v != 1 for v in cnt.values()):
        fail.append("certified rules differ from the rules violating P (own decision)")
    P_not_sm = [s for s in P if not simple_minor_preserving(*parse_rule(s))]
    if P_not_sm:
        fail.append("%d rules satisfy P but are not simple-minor-preserving" % len(P_not_sm))
    print("rules in list: %d | satisfy P: %d | violate P: %d | certificates: %d" % (len(rules), len(P), len(notP), len(certs)))
    print("hosts planar (rotation system): %d/%d | bindings valid: %d/%d | outputs reproduced: %d/%d | Kuratowski subdivisions valid: %d/%d"
          % (ok_host, len(certs), ok_bind, len(certs), ok_out, len(certs), ok_kur, len(certs)))
    print("rules satisfying P that are not simple-minor-preserving: %d" % len(P_not_sm))
    for f in fail[:10]:
        print("FAIL:", f)
    print("PLANARITY CERTIFICATE VERIFICATION PASSED" if not fail else "PLANARITY CERTIFICATE VERIFICATION FAILED")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
