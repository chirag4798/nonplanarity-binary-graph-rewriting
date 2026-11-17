#!/usr/bin/env python3
"""
check_noninjective_rigidity.py -- numerical check that the rigidity results (Theorem on rigidity, Proposition on
transported subgraphs) hold when matches are NOT required to be injective.

For every minor-preserving rule of S, on random connected hosts with loops and parallel edges, and for sampled matches in
which distinct pattern variables may bind the same vertex, one event is applied and the following is checked.

 (T) Transport (any delta_R).  For random edge sets D1, D2 of the host (a subgraph is identified with its edge set; its
     vertices are the endpoints of its edges), the transports D1', D2' are built exactly as in the paper (unconsumed edges,
     images of consumed edges, spanning trees of the branch sets rooted at the images of old variables that are vertices
     of D).  Subgraphs are pairs (V,E), isolated vertices allowed; intersection and union act on vertices and edges (a shared vertex without a shared edge
     counts).  Then beta_1 and the number of components of D_i', D1' u D2' and D1' n D2' equal those of D_i, D1 u D2, D1 n D2.
 (R) Rigidity (delta_R = 0).  After contracting the spanning trees, the edges outside X (those with a discarded end)
     form exactly the host H, the edges in X with the discarded vertices form a tree after contracting V(H) to a point,
     and beta_1 and the number of components of H' and H agree.
 (P) Path-preserving rules with delta_R = 0: the 2-core of H' is a subdivision of the 2-core of H (equal after suppressing
     vertices of degree 2).
Fixed seed.  The check exercises the steps of the proofs; it does not replace them.
"""
import argparse, os, random, sys
HERE_ = os.path.dirname(os.path.abspath(__file__))
from collections import Counter, defaultdict
import networkx as nx

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import oligon_experiment as oe
import wolfram_bridge as wb


def invariants(edges, verts=None):
    """edges: list of (u, v); verts: optional vertex set (isolated vertices count). Returns (beta1, components)."""
    parent = {}
    for x in (verts or ()):
        parent[x] = x

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for u, v in edges:
        parent.setdefault(u, u); parent.setdefault(v, v)
    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
    nv = len(parent)
    c = len({find(x) for x in parent}) if nv else 0
    return len(edges) - nv + c, c


def multigraph(edges):
    G = nx.MultiGraph()
    G.add_edges_from(edges)
    return G


def two_core(G):
    G = G.copy()
    while True:
        low = [v for v in G if G.degree(v) <= 1]
        if not low:
            return G
        G.remove_nodes_from(low)


def suppress(G):
    """Suppress vertices of degree 2 that are not carrying a loop; result unique up to isomorphism."""
    G = G.copy()
    changed = True
    while changed:
        changed = False
        for v in list(G):
            if G.degree(v) == 2 and not G.has_edge(v, v):
                nb = [w for _, w in G.edges(v)]
                a, b = nb
                G.remove_node(v)
                G.add_edge(a, b)
                changed = True
                break
    return G


def random_host(rng):
    n = rng.randint(3, 6)
    edges = [(rng.randrange(i), i) for i in range(1, n)]
    for _ in range(rng.randint(2, 7)):
        a = rng.randrange(n); b = rng.randrange(n)
        edges.append((a, b))                     # loops and parallel edges allowed
    edges = [(a, b) if rng.random() < .5 else (b, a) for a, b in edges]
    return edges


def model_data(r, assign):
    br = {a: a for a in range(r.n_old)}
    for k, x in enumerate(assign):
        br[r.n_old + k] = x
    cls = []
    for u, v in r.K:
        if -1 in (br[u], br[v]):
            cls.append("X")
        elif br[u] == br[v]:
            cls.append("C")
        else:
            cls.append("M")
    # spanning tree of each branch set from the C-edges
    tree = {}
    for a in range(r.n_old):
        seen = {a}; tr = []; changed = True
        while changed:
            changed = False
            for ki, (u, v) in enumerate(r.K):
                if cls[ki] == "C" and br[u] == a and ((u in seen) != (v in seen)):
                    tr.append(ki); seen |= {u, v}; changed = True
        tree[a] = tr
    # bijection L edges -> M edges
    used = set(); img = {}
    for li, (a, b) in enumerate(r.L):
        for ki, (u, v) in enumerate(r.K):
            if cls[ki] == "M" and ki not in used and {br[u], br[v]} == {a, b}:
                used.add(ki); img[li] = ki; break
        else:
            return None
    return br, cls, tree, img


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20240601)
    ap.add_argument("--hosts", type=int, default=5)
    ap.add_argument("--cap", type=int, default=12)
    ap.add_argument("--control", action="store_true", help="negative control: transport WITHOUT the branch-set trees (must be detected)")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    S = oe.enumerate_rules(3, 1) + oe.enumerate_rules(4, 2)
    st = Counter(); bad = []
    for r in S:
        mm = oe.minor_models(r)
        if not mm:
            continue
        st["rules"] += 1
        md = model_data(r, mm[0]["assign"])
        if md is None:
            bad.append(("no bijection", oe.rule_str(r))); continue
        br, cls, tree, img = md
        if a.control:
            tree = {k: [] for k in tree}
        dz = oe.delta_R(r)
        pathp = oe.is_path_preserving(r)
        for _ in range(a.hosts):
            E = random_host(rng)
            ms = wb.find_matches_general(E, r)
            if len(ms) > a.cap:
                ms = rng.sample(ms, a.cap)
            for (i, j, bind) in ms:
                nid = 1 + max(max(e) for e in E)
                full = dict(bind)
                for k in range(r.n_new):
                    full[r.n_old + k] = nid + k
                cons = {0: i, 1: j}
                # labelled edges of H'
                lab = {}
                for idx, e in enumerate(E):
                    if idx not in (i, j):
                        lab[("e", idx)] = e
                for ki, (u, v) in enumerate(r.K):
                    lab[("k", ki)] = (full[u], full[v])
                nonin = len(set(bind.values())) < len(bind)
                st["events"] += 1; st["noninjective"] += nonin
                # ---- (T) transport
                allv = sorted({x for e in E for x in e})
                for _t in range(3):
                    D1 = {idx for idx in range(len(E)) if rng.random() < .5}
                    D2 = {idx for idx in range(len(E)) if rng.random() < .5}
                    # general subgraphs (V, E): endpoints plus possibly isolated vertices of the host
                    V1 = {x for idx in D1 for x in E[idx]} | {x for x in allv if rng.random() < .2}
                    V2 = {x for idx in D2 for x in E[idx]} | {x for x in allv if rng.random() < .2}

                    def transport(D, V):
                        out = {("e", idx) for idx in D if idx not in (i, j)}
                        for li, idx in cons.items():
                            if idx in D:
                                out.add(("k", img[li]))
                        for ao in range(r.n_old):
                            if bind[ao] in V:
                                out |= {("k", ki) for ki in tree[ao]}
                        return out

                    def verts_of_T(Tx, V):
                        # old vertices of D persist in H' (closed rule); new vertices are the ends of transported edges
                        return set(V) | {x for l in Tx for x in lab[l]}
                    T1, T2 = transport(D1, V1), transport(D2, V2)
                    W1, W2 = verts_of_T(T1, V1), verts_of_T(T2, V2)
                    pairs = [((D1, V1), (T1, W1)), ((D2, V2), (T2, W2)),
                             ((D1 | D2, V1 | V2), (T1 | T2, W1 | W2)), ((D1 & D2, V1 & V2), (T1 & T2, W1 & W2))]
                    for ((Dx, Vx), (Tx, Wx)) in pairs:
                        before = invariants([E[idx] for idx in Dx], Vx)
                        after = invariants([lab[l] for l in Tx], Wx)
                        st["transport_checks"] += 1
                        if before != after:
                            bad.append(("transport", oe.rule_str(r), E, (i, j, bind), sorted(Dx), before, after))
                # ---- (R) rigidity, delta_R = 0
                if dz == 0:
                    mp = {}
                    for k in range(r.n_new):
                        x = br[r.n_old + k]
                        mp[nid + k] = None if x == -1 else bind[x]
                    G2 = []
                    Xed = []
                    for l, (u, v) in lab.items():
                        u2 = mp.get(u, u) if u >= nid else u
                        v2 = mp.get(v, v) if v >= nid else v
                        if l[0] == "k" and cls[l[1]] == "C":
                            continue                      # contracted tree edge
                        if l[0] == "k" and cls[l[1]] == "X":
                            Xed.append((u2 if u2 is not None else u, v2 if v2 is not None else v)); continue
                        G2.append(tuple(sorted((u2, v2))))
                    ok = Counter(G2) == Counter(tuple(sorted(e)) for e in E)
                    # tree condition for X
                    V0 = {x for e in E for x in e}
                    Q = nx.MultiGraph()
                    for u, v in Xed:
                        Q.add_edge("s" if u in V0 else u, "s" if v in V0 else v)
                    if r.n_new and any(x == -1 for x in br.values()):
                        okQ = Q.number_of_edges() == Q.number_of_nodes() - 1 and nx.is_connected(Q) and not any(u == v for u, v in Q.edges())
                    else:
                        okQ = Q.number_of_edges() == 0
                    Hp = list(lab.values())
                    ok = ok and okQ and invariants(Hp) == invariants(E)
                    st["rigidity_checks"] += 1
                    if not ok:
                        bad.append(("rigidity", oe.rule_str(r), E, (i, j, bind)))
                    if pathp:
                        c0, c1 = two_core(multigraph(E)), two_core(multigraph(Hp))
                        same = nx.is_isomorphic(suppress(c0), suppress(c1))
                        st["two_core_checks"] += 1
                        if not same:
                            bad.append(("2-core", oe.rule_str(r), E, (i, j, bind)))
    print(dict(st))
    gen = os.path.join(HERE_, "..", "paper", "generated")
    if os.path.isdir(gen) and not a.control:
        with open(os.path.join(gen, "noninj_rigidity_nums.tex"), "w") as f:
            for k, v in (("nrRules", st["rules"]), ("nrEvents", st["events"]), ("nrNonInj", st["noninjective"]),
                         ("nrTransport", st["transport_checks"]), ("nrRigid", st["rigidity_checks"]), ("nrTwoCore", st["two_core_checks"])):
                f.write("\\newcommand{\\%s}{%d}\n" % (k, v))
    for b in bad[:8]:
        print("VIOLATION", b)
    if a.control:
        print("negative control: %d violations detected (must be > 0)" % len(bad))
        return 0 if bad else 1
    ok = not bad
    print("NON-INJECTIVE RIGIDITY CHECK PASSED" if ok else "NON-INJECTIVE RIGIDITY CHECK FAILED (%d)" % len(bad))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
