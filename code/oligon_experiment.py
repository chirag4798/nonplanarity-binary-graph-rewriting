#!/usr/bin/env python3
"""
oligon_experiment.py -- reproducible experiments for the paper
"Topological Stability and Discrete Charge Conservation of Minimal Oligons in
2-Way Causally Invariant Hypergraph Rewriting Systems".

Sub-commands
------------
  selftest   Machine-check the paper's lemmas/theorems on small graphs (exact, no sampling
             of physical claims).  Exit code != 0 on any failure.
  census     Enumerate minimal binary rules 2_2 -> k_2 (k = 3, 4), compute rule invariants
             (delta_R, closed, path-preserving) and run a *bounded-depth* joinability filter.
  run        Evolve control / single-defect / two-defect-collision experiments on a lattice
             background and write per-step CSV rows.
  all        census, then `run` on every rule that passes the joinability filter.

IMPORTANT HONESTY NOTES (mirrored in the paper's Limitations section)
---------------------------------------------------------------------
* The confluence filter is a BOUNDED necessary-condition test (finite seeds, finite
  generations, finite joining depth).  Confluence of a non-terminating system is
  undecidable in general; a pass is NOT a proof of causal invariance.
* "dim_est" is a log-log slope of mean ball volume over radii 2..rmax.  It is an
  indicative growth exponent on finite graphs, not a Hausdorff dimension.
* kappa is the Ollivier-Ricci curvature with idleness alpha = 1/2 (kappa_alpha, NOT
  normalised by 1-alpha), computed on the simple undirected projection, on a seeded
  random edge sample (global) and on defect-labelled edges (defect).
* Edge direction matters for MATCHING (rules are on ordered pairs); all topological
  observables use the underlying undirected multigraph.

Dependencies (tested versions): Python 3.13.5, networkx 3.6.1, numpy 2.3.3, scipy 1.16.3.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
import platform
import random
import sys
import time
from collections import Counter, defaultdict, namedtuple
from multiprocessing import Pool

import networkx as nx
import numpy as np
import scipy
from scipy.optimize import linprog

DEFAULT_SEED = 20240601

# --------------------------------------------------------------------------- rules
Rule = namedtuple("Rule", "L K n_old n_new")  # L,K: tuples of (var,var); vars 0..n_old+n_new-1


def rule_str(r: Rule) -> str:
    f = lambda es: "{" + ",".join("{%d,%d}" % (a + 1, b + 1) for a, b in es) + "}"
    return f(r.L) + "->" + f(r.K)


def parse_rule(s: str) -> Rule:
    s = s.replace(" ", "")
    lhs, rhs = s.split("->")

    def edges(t):
        t = t.strip("{}")
        out = []
        for chunk in t.replace("},{", "|").split("|"):
            a, b = chunk.strip("{}").split(",")
            out.append((int(a) - 1, int(b) - 1))
        return tuple(out)

    L, K = edges(lhs), edges(rhs)
    nl = len({v for e in L for v in e})
    nk = len({v for e in K for v in e})
    return Rule(L, K, nl, max(0, nk - nl))


def _relabel_first_appearance(edges):
    m = {}
    out = []
    for a, b in edges:
        for v in (a, b):
            if v not in m:
                m[v] = len(m)
        out.append((m[a], m[b]))
    return tuple(out)


def _lhs_classes():
    """Connected 2-edge, loop-free ordered patterns up to renaming/edge order."""
    seen = set()
    for e1 in itertools.permutations(range(4), 2):
        for e2 in itertools.permutations(range(4), 2):
            if not ({e1[0], e1[1]} & {e2[0], e2[1]}):
                continue
            c = min(_relabel_first_appearance((e1, e2)), _relabel_first_appearance((e2, e1)))
            seen.add(c)
    return sorted(seen)


def _weakly_connected(edges):
    g = nx.Graph()
    g.add_edges_from(edges)
    return g.number_of_nodes() > 0 and nx.is_connected(g)


def enumerate_rules(k_out: int, max_new: int, allow_loops=False, closed_only=True):
    """All rules 2_2 -> k_out_2 up to renaming of variables and edge order.

    closed_only=True keeps rules with: L connected (by construction), every L variable
    occurring in K (no vertex loss), K weakly connected, no K loops unless allow_loops.
    """
    rules = []
    for L in _lhs_classes():
        n_old = len({v for e in L for v in e})
        auts = []
        for perm in itertools.permutations(range(n_old)):
            mp = dict(enumerate(perm))
            if sorted((mp[a], mp[b]) for a, b in L) == sorted(L):
                auts.append(mp)
        for n_new in range(0, max_new + 1):
            nv = n_old + n_new
            pairs = [(a, b) for a in range(nv) for b in range(nv) if allow_loops or a != b]
            seen = set()
            for K in itertools.combinations_with_replacement(pairs, k_out):
                used = {v for e in K for v in e}
                if closed_only:
                    if not set(range(n_old)) <= used:
                        continue
                    if not _weakly_connected(K):
                        continue
                if not set(range(n_old, nv)) <= used:
                    continue
                best = None
                for mp in auts:
                    for newperm in itertools.permutations(range(n_old, nv)):
                        full = dict(mp)
                        full.update(dict(zip(range(n_old, nv), newperm)))
                        img = tuple(sorted((full[a], full[b]) for a, b in K))
                        if best is None or img < best:
                            best = img
                if best in seen:
                    continue
                seen.add(best)
                rules.append(Rule(tuple(L), best, n_old, n_new))
    return rules


def is_closed(r: Rule) -> bool:
    used = {v for e in r.K for v in e}
    return set(range(r.n_old)) <= used and _weakly_connected(r.K) and _weakly_connected(r.L)


def delta_R(r: Rule) -> int:
    """delta_R = |K| - |L| - n_new  (exact Delta beta_1 for closed rules, Lemma 1)."""
    return len(r.K) - len(r.L) - r.n_new


def is_path_preserving(r: Rule) -> bool:
    """Exists edge-disjoint paths P_e in K (undirected multigraph), one per L-edge e,
    joining e's endpoints, with all internal vertices new variables and internal
    vertices pairwise distinct across paths (Definition 'path-preserving')."""
    K = list(r.K)
    new_vars = set(range(r.n_old, r.n_old + r.n_new))

    def paths(a, b, used_edges, used_int):
        # DFS over K-edge indices from a to b
        res = []

        def dfs(cur, visited_int, edge_idx):
            for idx, (u, v) in enumerate(K):
                if idx in used_edges or idx in edge_idx:
                    continue
                if u == cur:
                    nxt = v
                elif v == cur:
                    nxt = u
                else:
                    continue
                if nxt == b:
                    res.append((tuple(edge_idx) + (idx,), tuple(visited_int)))
                elif nxt in new_vars and nxt not in used_int and nxt not in visited_int:
                    dfs(nxt, visited_int + [nxt], edge_idx + [idx])

        dfs(a, [], [])
        return res

    def assign(i, used_edges, used_int):
        if i == len(r.L):
            return True
        a, b = r.L[i]
        for edge_idx, ints in paths(a, b, used_edges, used_int):
            if assign(i + 1, used_edges | set(edge_idx), used_int | set(ints)):
                return True
        return False

    return assign(0, frozenset(), frozenset())


def minor_models(r: Rule, simple: bool = False):
    """All valid 'rooted-minor models' of L inside K: an assignment of each new variable to the
    branch set of an old variable (0..n_old-1) or to -1 (discarded) such that every branch set is
    connected by K-edges inside it, and the K-edges between distinct branch sets contain every
    L-edge (as a multiset, orientation ignored).  Contracting branch sets and deleting discarded
    vertices then turns the image of K into a graph containing the image of L, i.e. H is a minor of H'."""
    new_vars = list(range(r.n_old, r.n_old + r.n_new))
    out = []
    for assign in itertools.product(list(range(r.n_old)) + [-1], repeat=len(new_vars)):
        br = {a: a for a in range(r.n_old)}
        br.update(dict(zip(new_vars, assign)))
        ok = True
        for a in range(r.n_old):
            S = {v for v, b in br.items() if b == a}
            g = nx.Graph()
            g.add_nodes_from(S)
            g.add_edges_from((u, v) for u, v in r.K if u in S and v in S and u != v)
            if not nx.is_connected(g):
                ok = False
                break
        if not ok:
            continue
        C = sum(1 for u, v in r.K if br[u] == br[v] != -1)
        X = sum(1 for u, v in r.K if br[u] == -1 or br[v] == -1)
        M = [(br[u], br[v]) for u, v in r.K if -1 not in (br[u], br[v]) and br[u] != br[v]]
        need = Counter(frozenset(e) for e in r.L)
        if simple:      # planarity only sees the underlying SIMPLE graph: each distinct L pair must be covered once
            need = Counter({k: 1 for k in need})
        have = Counter(frozenset(e) for e in M)
        if all(have[k] >= v for k, v in need.items()):
            out.append(dict(assign=assign, C=C, X=X, M=len(M),
                            A=sum(1 for b in assign if b != -1), Dv=sum(1 for b in assign if b == -1)))
    return out


def is_minor_preserving(r: Rule) -> bool:
    return bool(minor_models(r))


def is_simple_minor_preserving(r: Rule) -> bool:
    """Weaker: the SIMPLE projection of H is a minor of that of H' (enough for planarity; not for beta_1)."""
    return bool(minor_models(r, simple=True))


# --------------------------------------------------------------------------- states
# A state is a list of (a, b, mask) with ints a,b and mask an int bitmask of defect labels.


def graph_to_edges(G: nx.Graph, rng=None):
    out = []
    for u, v in sorted(G.edges()):
        out.append((u, v, 0))
    return out


def find_matches(edges, rule: Rule):
    """All injective matches (i, j, binding): edge i matches rule.L[0], edge j rule.L[1]."""
    inc = defaultdict(list)
    for idx, (a, b, _) in enumerate(edges):
        inc[a].append(idx)
        inc[b].append(idx)
    (p, q), (r, s) = rule.L
    out = []
    for i, (a, b, _) in enumerate(edges):
        if a == b:
            continue
        bind = {p: a, q: b}
        cand = set()
        for v in (r, s):
            if v in bind:
                cand.update(inc[bind[v]])
        for j in sorted(cand):
            if j == i:
                continue
            c, d, _ = edges[j]
            if c == d:
                continue
            b2 = dict(bind)
            ok = True
            for var, val in ((r, c), (s, d)):
                if var in b2:
                    if b2[var] != val:
                        ok = False
                        break
                else:
                    b2[var] = val
            if ok and len(set(b2.values())) == len(b2):
                out.append((i, j, b2))
    return out


def tracking_model(rule: Rule):
    """Label-transport model for minor-preserving rules (None otherwise): the first minor model, the K-edge
    index imaging each L edge, and the K edges inside branch sets (with their branch variable)."""
    key = (rule.L, rule.K)
    if key in _TM_CACHE:
        return _TM_CACHE[key]
    models = minor_models(rule)
    tm = None
    if models:
        assign = models[0]["assign"]
        br = {a: a for a in range(rule.n_old)}
        br.update(dict(zip(range(rule.n_old, rule.n_old + rule.n_new), assign)))
        used, img = set(), []
        for (la, lb) in rule.L:
            for idx, (u, v) in enumerate(rule.K):
                if idx in used or -1 in (br[u], br[v]) or br[u] == br[v]:
                    continue
                if frozenset((br[u], br[v])) == frozenset((la, lb)):
                    used.add(idx)
                    img.append(idx)
                    break
        # only a spanning forest of each branch set's internal edges is labelled (extra internal edges would add cycles)
        parent = {}

        def find(x):
            while parent.setdefault(x, x) != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        C = []
        for idx, (u, v) in enumerate(rule.K):
            if br[u] == br[v] != -1:
                ru, rv = find(u), find(v)
                if ru != rv:
                    parent[ru] = rv
                    C.append((idx, br[u]))
        tm = dict(img=img, C=C)
    _TM_CACHE[key] = tm
    return tm


_TM_CACHE = {}


def apply_matches(edges, chosen, rule: Rule, next_vid: int, vmask=None, tm=None):
    """Apply pairwise edge-disjoint matches simultaneously.
    Default ('cone') labelling: every product edge inherits the bitwise OR of the defect masks of the two
    consumed edges.  With tm (from tracking_model) and vmask (vertex -> OR of incident edge masks, computed
    before the step) the label is transported through the minor model instead: the K-edge imaging L-edge i
    inherits the mask of the consumed edge matched to it; K-edges inside the branch set of old variable a inherit
    vmask of the vertex bound to a (vertex expansion of a vertex lying on the defect); all other K edges are
    unlabelled.  Then the labelled subgraph D' is obtained from D by vertex expansions and edge replacement, so
    D' is homotopy equivalent to D and contains D as a minor."""
    consumed = set()
    produced = []
    for (i, j, bind) in chosen:
        assert i not in consumed and j not in consumed
        mask = edges[i][2] | edges[j][2]
        m = dict(bind)
        for k in range(rule.n_old, rule.n_old + rule.n_new):
            m[k] = next_vid
            next_vid += 1
        cmap = {}
        if tm is not None:
            cmap[tm["img"][0]] = edges[i][2]
            cmap[tm["img"][1]] = edges[j][2]
            for idx, a in tm["C"]:
                cmap[idx] = vmask.get(bind[a], 0)
        for idx, (u, v) in enumerate(rule.K):
            produced.append((m[u], m[v], cmap.get(idx, 0) if tm is not None else mask))
        consumed.update((i, j))
    kept = [e for idx, e in enumerate(edges) if idx not in consumed]
    return kept + produced, next_vid


def next_vertex_id(edges):
    return 1 + max((max(a, b) for a, b, _ in edges), default=-1)


def generation_step(edges, rule, rng, order="random", mode="cone"):
    """One Wolfram-style generation: a maximal set of pairwise disjoint matches."""
    ms = find_matches(edges, rule)
    if order == "random":
        rng.shuffle(ms)
    else:
        ms.sort(key=lambda m: (m[0], m[1]))
    used, chosen = set(), []
    for m in ms:
        if m[0] in used or m[1] in used:
            continue
        used.update((m[0], m[1]))
        chosen.append(m)
    if mode == "model":
        tm = tracking_model(rule)
        if tm is None:
            raise ValueError("rule is not minor-preserving: no label-transport model")
        vmask = defaultdict(int)
        for a_, b_, m_ in edges:
            vmask[a_] |= m_
            vmask[b_] |= m_
        new_edges, _ = apply_matches(edges, chosen, rule, next_vertex_id(edges), vmask=vmask, tm=tm)
    else:
        new_edges, _ = apply_matches(edges, chosen, rule, next_vertex_id(edges))
    return new_edges, len(chosen)


# --------------------------------------------------------------------------- topology
def beta1_of(edge_pairs):
    """(beta1, V, E, c) of the multigraph with the given (a,b) pairs.  beta1 = E - V + c."""
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    E = 0
    for a, b in edge_pairs:
        E += 1
        for v in (a, b):
            parent.setdefault(v, v)
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    V = len(parent)
    c = sum(1 for v in parent if find(v) == v)
    return E - V + c, V, E, c


def simple_graph(edge_pairs):
    G = nx.Graph()
    G.add_edges_from((a, b) for a, b in edge_pairs if a != b)
    return G


def is_planar(G):
    return bool(nx.check_planarity(G)[0]) if G.number_of_nodes() else True


def reduced_core(edges):
    """Multigraph obtained by pruning degree<=1 vertices repeatedly, then suppressing
    degree-2 non-loop vertices (homeomorphism type of the 2-core)."""
    M = nx.MultiGraph()
    for a, b, *_ in edges:
        M.add_edge(a, b)
    changed = True
    while changed:
        changed = False
        for v in list(M.nodes):
            if v in M and M.degree(v) <= 1:
                M.remove_node(v)
                changed = True
    changed = True
    while changed:
        changed = False
        for v in list(M.nodes):
            if v in M and M.degree(v) == 2 and not M.has_edge(v, v):
                nb = [u for u in M[v] for _ in range(M.number_of_edges(v, u))]
                a, b = nb
                M.remove_node(v)
                M.add_edge(a, b)
                changed = True
    return M


# --------------------------------------------------------------------------- curvature
def orc_edge(G, x, y, alpha=0.5):
    """Ollivier-Ricci kappa_alpha(x,y) = 1 - W1(mu_x, mu_y); mu_v = alpha*delta_v +
    (1-alpha)*uniform(N(v)); ground metric = graph distance (all support distances <= 3)."""
    def mu(v):
        nb = list(G[v])
        d = {v: alpha}
        for n in nb:
            d[n] = d.get(n, 0.0) + (1 - alpha) / len(nb)
        return d

    mx, my = mu(x), mu(y)
    Sx, Sy = sorted(mx), sorted(my)
    cost = np.zeros((len(Sx), len(Sy)))
    for i, u in enumerate(Sx):
        dist = nx.single_source_shortest_path_length(G, u, cutoff=3)
        for j, v in enumerate(Sy):
            cost[i, j] = dist[v]
    n, m = len(Sx), len(Sy)
    A = np.zeros((n + m, n * m))
    for i in range(n):
        A[i, i * m:(i + 1) * m] = 1
    for j in range(m):
        A[n + j, j::m] = 1
    b = np.array([mx[u] for u in Sx] + [my[v] for v in Sy])
    res = linprog(cost.ravel(), A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    if not res.success:
        raise RuntimeError("ORC LP failed: " + res.message)
    return 1.0 - res.fun


def mean_orc(G, edge_list, cap, rng):
    edge_list = sorted(edge_list)
    if len(edge_list) > cap:
        edge_list = rng.sample(edge_list, cap)
    if not edge_list:
        return float("nan"), 0
    vals = [orc_edge(G, x, y) for x, y in edge_list]
    return float(np.mean(vals)), len(vals)


def dim_estimate(G, rng, ncent=8, rmax=8):
    nodes = sorted(G.nodes)
    if len(nodes) < 10:
        return float("nan")
    centers = rng.sample(nodes, min(ncent, len(nodes)))
    vol = np.zeros(rmax + 1)
    for c in centers:
        d = nx.single_source_shortest_path_length(G, c, cutoff=rmax)
        cnt = np.bincount(list(d.values()), minlength=rmax + 1)[: rmax + 1].cumsum()
        vol += cnt
    vol /= len(centers)
    r = np.arange(2, rmax + 1)
    return float(np.polyfit(np.log(r), np.log(vol[2:]), 1)[0])


# --------------------------------------------------------------------------- confluence
def _canon(edges2):
    G = nx.DiGraph()
    for a, b in edges2:
        if G.has_edge(a, b):
            G[a][b]["m"] += 1
        else:
            G.add_edge(a, b, m=1)
    for n in G:
        G.nodes[n]["l"] = "%d-%d" % (G.out_degree(n), G.in_degree(n))
    h = nx.weisfeiler_lehman_graph_hash(G, node_attr="l", edge_attr="m", iterations=3)
    return G, (G.number_of_nodes(), G.number_of_edges(), h)


def _iso(g, h):
    return nx.is_isomorphic(g, h, edge_match=lambda a, b: a["m"] == b["m"])


class _Side:
    def __init__(self):
        self.b = defaultdict(list)

    def add(self, edges3):
        g, k = _canon([(a, b) for a, b, _ in edges3])
        for h in self.b[k]:
            if _iso(g, h):
                return False
        self.b[k].append(g)
        return True

    def meets(self, other):
        for k, gs in self.b.items():
            for h in other.b.get(k, []):
                for g in gs:
                    if _iso(g, h):
                        return True
        return False


def _expand(frontier, rule, side, branch_cap):
    out, truncated = [], False
    for st in frontier:
        ms = find_matches(st, rule)
        if len(ms) > branch_cap:
            truncated = True
            ms = ms[:branch_cap]
        nv = next_vertex_id(st)
        for m in ms:
            ns, _ = apply_matches(st, [m], rule, nv)
            if side.add(ns):
                out.append(ns)
    return out, truncated


def joinable(e1, e2, rule, J, branch_cap):
    """Return (joined: bool, truncated: bool): do the two states have isomorphic
    descendants within J further single-event steps each?"""
    s1, s2 = _Side(), _Side()
    s1.add(e1)
    s2.add(e2)
    f1, f2 = [e1], [e2]
    trunc = False
    for depth in range(J + 1):
        if s1.meets(s2):
            return True, trunc
        if depth == J:
            break
        f1, t1 = _expand(f1, rule, s1, branch_cap)
        f2, t2 = _expand(f2, rule, s2, branch_cap)
        trunc = trunc or t1 or t2
    return False, trunc


def _conf_seeds():
    tri = [(0, 1), (1, 2), (2, 0)]
    bowtie = [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)]
    g23 = [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]
    return [[(a, b, 0) for a, b in s] for s in (tri, bowtie, g23)]


def confluence_test(rule: Rule, depth=3, max_edges=14, pairs_cap=12, J=2, branch_cap=30, seed=DEFAULT_SEED):
    """Bounded local-joinability test.  Test states = all states along a canonical-order
    generation evolution (<= `depth` generations, <= `max_edges` edges) from three small
    seeds.  For each test state, up to `pairs_cap` unordered pairs of distinct matches are
    applied separately and required to be joinable within J further steps."""
    rng = random.Random(seed)
    states = []
    for s in _conf_seeds():
        cur = s
        for g in range(depth + 1):
            states.append(cur)
            if len(cur) > max_edges:
                break
            cur, n = generation_step(cur, rule, rng, order="canonical")
            if n == 0:
                break
    tested = failures = incon = 0
    seen_keys = []
    for st in states:
        g, k = _canon([(a, b) for a, b, _ in st])
        if any(k == k2 and _iso(g, g2) for k2, g2 in seen_keys):
            continue
        seen_keys.append((k, g))
        ms = find_matches(st, rule)
        pairs = list(itertools.combinations(range(len(ms)), 2))
        if len(pairs) > pairs_cap:
            pairs = rng.sample(pairs, pairs_cap)
        nv = next_vertex_id(st)
        for (a, b) in pairs:
            s1, _ = apply_matches(st, [ms[a]], rule, nv)
            s2, _ = apply_matches(st, [ms[b]], rule, nv)
            ok, trunc = joinable(s1, s2, rule, J, branch_cap)
            tested += 1
            if not ok:
                if trunc:
                    incon += 1
                else:
                    failures += 1
                    return dict(status="unjoined_within_J", pairs_tested=tested, unjoined=failures, inconclusive=incon)
    status = "inconclusive" if incon else "joinable_all_tested"
    return dict(status=status, pairs_tested=tested, unjoined=failures, inconclusive=incon)


def _census_worker(args):
    r, kw = args
    c = confluence_test(r, **kw)
    return r, c


# --------------------------------------------------------------------------- backgrounds / defects
def make_background(kind: str, n: int):
    dims = {"grid2d": 2, "grid3d": 3}[kind]
    G = nx.grid_graph(dim=[n] * dims)
    idx = {v: i for i, v in enumerate(sorted(G.nodes()))}
    H = nx.relabel_nodes(G, idx)
    return H, idx, dims


def pick_sites(G, center, k, exclude=(), spacing=2):
    dist0 = nx.single_source_shortest_path_length(G, center)
    order = sorted(dist0, key=lambda v: (dist0[v], v))
    chosen = []
    for sp in (spacing, 1):
        chosen = []
        for v in order:
            if v in exclude:
                continue
            if all(nx.shortest_path_length(G, v, s) >= sp for s in chosen):
                chosen.append(v)
            if len(chosen) == k:
                return chosen
    raise ValueError("not enough sites near center for defect; enlarge the background")


def defect_edges(kind: str, sites, mask):
    if kind == "K5":
        pairs = list(itertools.combinations(range(5), 2))
    elif kind == "K33":
        pairs = [(i, 3 + j) for i in range(3) for j in range(3)]
    else:
        raise ValueError(kind)
    return [(min(sites[i], sites[j]), max(sites[i], sites[j]), mask) for i, j in pairs]


# --------------------------------------------------------------------------- evolution runs
FIELDS = [
    "experiment", "tracking", "rule", "background", "defect", "order", "seed", "t", "events",
    "V", "E", "beta1_multi_global", "beta1_simple_global", "planar_simple_global",
    "D1_edges", "beta1_D1", "D1_planar", "D2_edges", "beta1_D2", "D2_planar",
    "union_edges", "beta1_union", "beta1_sum", "iota_union_minus_sum", "Q_sum_initial", "Q_sum_conserved",
    "orc_global_mean", "orc_global_n", "orc_defect_mean", "orc_defect_n", "dim_est", "stopped_reason",
]


def observe(edges, t, events, ctx, rng):
    pairs = [(a, b) for a, b, _ in edges]
    b_multi, V, E, _ = beta1_of(pairs)
    Gs = simple_graph(pairs)
    row = {k: "" for k in FIELDS}
    row.update(experiment=ctx["experiment"], tracking=ctx["tracking"], rule=ctx["rule"], background=ctx["bg"], defect=ctx["defect"],
               order=ctx["order"], seed=ctx["seed"], t=t, events=events, V=V, E=E,
               beta1_multi_global=b_multi, beta1_simple_global=nx.number_of_edges(Gs) - nx.number_of_nodes(Gs) + nx.number_connected_components(Gs))
    if ctx["planarity_global"] and E <= ctx["planarity_max_edges"]:
        row["planar_simple_global"] = int(is_planar(Gs))
    sub = {}
    for bit, name in ((1, "D1"), (2, "D2")):
        if bit in ctx["bits"]:
            es = [(a, b) for a, b, m in edges if m & bit]
            sub[bit] = es
            row[name + "_edges"] = len(es)
            row["beta1_" + name] = beta1_of(es)[0] if es else 0
            row[name + "_planar"] = int(is_planar(simple_graph(es)))
    if ctx["bits"]:
        un = [(a, b) for a, b, m in edges if m & sum(ctx["bits"])]
        row["union_edges"] = len(un)
        row["beta1_union"] = beta1_of(un)[0] if un else 0
        q = sum(row["beta1_D%d" % b] for b in ctx["bits"])
        row["beta1_sum"] = q
        row["iota_union_minus_sum"] = row["beta1_union"] - q
        if t == 0:
            ctx["Q0"] = q
        row["Q_sum_initial"] = ctx["Q0"]
        row["Q_sum_conserved"] = int(q == ctx["Q0"])
    if ctx["orc"]:
        mu, n = mean_orc(Gs, list(Gs.edges()), ctx["orc_sample"], rng)
        row["orc_global_mean"], row["orc_global_n"] = mu, n
        if ctx["bits"]:
            dedges = {(min(a, b), max(a, b)) for a, b, m in edges if m and a != b}
            mu, n = mean_orc(Gs, sorted(dedges), ctx["orc_defect_cap"], rng)
            row["orc_defect_mean"], row["orc_defect_n"] = mu, n
    row["dim_est"] = dim_estimate(Gs, rng, rmax=ctx["dim_rmax"]) if ctx["dim"] else ""
    return row


def run_evolution(rule, experiment, args, bg_kind, defect_kind, sep, seed, order):
    if args.tracking == "model" and tracking_model(rule) is None:
        raise ValueError("model tracking requested but rule is not minor-preserving")
    n = args.n2d if bg_kind == "grid2d" else args.n3d
    G, idx, dims = make_background(bg_kind, n)
    edges = graph_to_edges(G)
    bits = []
    center = tuple(n // 2 for _ in range(dims))
    inv = {v: i for i, v in idx.items()}
    if experiment in ("single", "collision"):
        c1 = idx[center]
        s1 = pick_sites(G, c1, 5 if defect_kind == "K5" else 6, spacing=2)
        edges += defect_edges(defect_kind, s1, 1)
        bits.append(1)
        if experiment == "collision":
            c2c = list(center)
            c2c[0] = center[0] + sep if center[0] + sep < n else center[0] - sep
            if not (0 <= c2c[0] < n):
                raise ValueError("separation too large for grid")
            s2 = pick_sites(G, idx[tuple(c2c)], 5 if defect_kind == "K5" else 6, exclude=set(s1), spacing=2)
            edges += defect_edges(defect_kind, s2, 2)
            bits.append(2)
    rng = random.Random(seed)
    ctx = dict(experiment=experiment, tracking=args.tracking, rule=rule_str(rule), bg=bg_kind, defect=defect_kind if bits else "none",
               order=order, seed=seed, bits=bits, planarity_global=(bg_kind == "grid2d" or args.global_planarity),
               planarity_max_edges=args.planarity_max_edges, orc=not args.no_orc, orc_sample=args.orc_sample,
               orc_defect_cap=args.orc_defect_cap, dim=not args.no_dim, dim_rmax=args.dim_rmax)
    rows = [observe(edges, 0, 0, ctx, random.Random(seed * 1000003))]
    for t in range(1, args.tmax + 1):
        if rows[-1]["E"] > args.max_edges:
            rows[-1]["stopped_reason"] = "edge_cap"
            break
        edges, ev = generation_step(edges, rule, rng, order, args.tracking)
        rows.append(observe(edges, t, ev, ctx, random.Random(seed * 1000003 + t)))
        if ev == 0:
            rows[-1]["stopped_reason"] = "no_matches"
            break
    return rows


def write_rows(path, rows, mode="a"):
    new = (mode == "w") or not os.path.exists(path)
    with open(path, mode if not new else "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new:
            w.writeheader()
        w.writerows(rows)


def manifest(args, out):
    m = dict(python=sys.version, platform=platform.platform(), networkx=nx.__version__, numpy=np.__version__,
             scipy=scipy.__version__, argv=sys.argv, args={k: v for k, v in vars(args).items() if k != "func"},
             time_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    with open(os.path.join(out, "run_manifest.json"), "w") as f:
        json.dump(m, f, indent=2, default=str)


# --------------------------------------------------------------------------- commands
def cmd_census(args):
    os.makedirs(args.out, exist_ok=True)
    manifest(args, args.out)
    rules = []
    for sig in args.sigs:
        k = int(sig.split("_")[1][0]) if "_" in sig else int(sig)
        mx = 1 if k == 3 else 2
        for r in enumerate_rules(k, mx if args.max_new is None else args.max_new):
            rules.append((sig, r))
    kw = dict(depth=args.conf_depth, max_edges=args.conf_max_edges, pairs_cap=args.conf_pairs,
              J=args.conf_join, branch_cap=args.conf_branch, seed=args.seed)
    print("enumerated %d rules" % len(rules), flush=True)
    jobs = [(r, kw) for _, r in rules]
    if args.jobs > 1:
        with Pool(args.jobs) as p:
            results = p.map(_census_worker, jobs, chunksize=4)
    else:
        results = [_census_worker(j) for j in jobs]
    path = os.path.join(args.out, "rule_census.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["signature", "rule", "n_old", "n_new", "closed", "delta_R", "path_preserving",
                    "conf_status", "conf_pairs_tested", "conf_unjoined", "conf_inconclusive"])
        for (sig, _), (r, c) in zip(rules, results):
            w.writerow([sig, rule_str(r), r.n_old, r.n_new, int(is_closed(r)), delta_R(r),
                        int(is_path_preserving(r)), c["status"], c["pairs_tested"], c["unjoined"], c["inconclusive"]])
    print("wrote", path)
    return [(sig, r, c) for (sig, _), (r, c) in zip(rules, results)]


def cmd_run(args, rules=None):
    os.makedirs(args.out, exist_ok=True)
    manifest(args, args.out)
    if rules is None:
        rules = [parse_rule(s) for s in args.rule]
    path = os.path.join(args.out, "evolution.csv")
    summ = os.path.join(args.out, "summary.csv")
    if os.path.exists(path):
        os.remove(path)
    summary = []
    for r in rules:
        for bg in args.backgrounds:
            plan = [("control", "none", 0)]
            for d in args.defects:
                plan.append(("single", d, 0))
                for sep in args.seps:
                    plan.append(("collision", d, sep))
            for (exp, d, sep) in plan:
                for order in args.orders:
                    seed = args.seed + 7919 * len(summary)
                    t0 = time.time()
                    try:
                        rows = run_evolution(r, exp, args, bg, d, sep, seed, order)
                    except ValueError as e:
                        print("SKIP", rule_str(r), bg, exp, d, sep, ":", e)
                        continue
                    write_rows(path, rows)
                    last = rows[-1]
                    summary.append(dict(rule=rule_str(r), background=bg, experiment=exp, defect=d, sep=sep, order=order,
                                        seed=seed, steps=last["t"], final_E=last["E"], final_dim_est=last["dim_est"],
                                        Q_conserved_all_steps=(all(x["Q_sum_conserved"] == 1 for x in rows) if exp != "control" else ""),
                                        beta1_global_constant=int(len({x["beta1_multi_global"] for x in rows}) == 1),
                                        D1_nonplanar_all_steps=(all(x["D1_planar"] == 0 for x in rows) if exp != "control" else ""),
                                        seconds=round(time.time() - t0, 2)))
                    print("done", rule_str(r), bg, exp, d, sep, order, "t=%d E=%d" % (last["t"], last["E"]), flush=True)
    if summary:
        with open(summ, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
            w.writeheader()
            w.writerows(summary)
        print("wrote", path, "and", summ)


def cmd_all(args):
    res = cmd_census(args)
    keep = [r for _, r, c in res if c["status"] == "joinable_all_tested" and is_closed(r)]
    if args.max_rules:
        keep = keep[: args.max_rules]
    print("%d rules pass the bounded joinability filter; running experiments on %d" % (
        sum(1 for _, r, c in res if c["status"] == "joinable_all_tested"), len(keep)))
    cmd_run(args, keep)


# --------------------------------------------------------------------------- selftest
def _random_connected(rng, nmin=6, nmax=10):
    while True:
        n = rng.randint(nmin, nmax)
        m = rng.randint(n, min(n * (n - 1) // 2, 2 * n))
        G = nx.gnm_random_graph(n, m, seed=rng.randint(0, 2**31))
        if nx.is_connected(G):
            return G


def _check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ((" :: " + detail) if detail else ""))
    return ok


def cmd_selftest(args):
    rng = random.Random(args.seed)
    allok = True
    rules = enumerate_rules(3, 1) + enumerate_rules(4, 2)
    print("selftest: %d closed rules enumerated (2_2->3_2 with <=1 new vertex, 2_2->4_2 with <=2)" % len(rules))

    # (1) Lemma 1 (exact change of beta_1 for closed rules, Delta c = 0)
    trials = bad = 0
    for _ in range(600):
        r = rng.choice(rules)
        edges = graph_to_edges(_random_connected(rng))
        ms = find_matches(edges, r)
        if not ms:
            continue
        new, _ = apply_matches(edges, [rng.choice(ms)], r, next_vertex_id(edges))
        b0, V0, E0, c0 = beta1_of([(a, b) for a, b, _ in edges])
        b1, V1, E1, c1 = beta1_of([(a, b) for a, b, _ in new])
        trials += 1
        if not (c1 == c0 and b1 - b0 == delta_R(r)):
            bad += 1
    allok &= _check("Lemma 1 (closed rule: c preserved, Delta beta1 = delta_R)", bad == 0 and trials > 100, "trials=%d violations=%d" % (trials, bad))

    # (1b) general identity incl. non-closed rules
    open_rules = enumerate_rules(3, 1, closed_only=False)
    trials = bad = 0
    for _ in range(600):
        r = rng.choice(open_rules)
        edges = graph_to_edges(_random_connected(rng))
        ms = find_matches(edges, r)
        if not ms:
            continue
        i, j, bind = rng.choice(ms)
        new, _ = apply_matches(edges, [(i, j, bind)], r, next_vertex_id(edges))
        b0, V0, E0, c0 = beta1_of([(a, b) for a, b, _ in edges])
        b1, V1, E1, c1 = beta1_of([(a, b) for a, b, _ in new])
        vb = {v for a, b, _ in edges for v in (a, b)}
        va = {v for a, b, _ in new for v in (a, b)}
        n_lost = len(vb - va)
        n_new_h = len(va - vb)
        trials += 1
        if b1 - b0 != (len(r.K) - len(r.L)) - n_new_h + n_lost + (c1 - c0):
            bad += 1
    allok &= _check("Lemma 1 general Euler identity (open rules too)", bad == 0 and trials > 100, "trials=%d violations=%d" % (trials, bad))

    # (2) Theorems 3/4 : minor-preserving rules
    mp = [r for r in rules if is_minor_preserving(r)]
    pp = [r for r in rules if is_path_preserving(r)]
    allok &= _check("path-preserving subset of minor-preserving", all(is_minor_preserving(r) for r in pp))
    print("   closed rules: %d; minor-preserving: %d; path-preserving: %d; minor-preserving with delta_R=0: %d"
          % (len(rules), len(mp), len(pp), sum(1 for r in mp if delta_R(r) == 0)))
    # Theorem 4 counting identity: for delta_R = 0 every minor model has |C|=A, |X|=Dv, |M|=|L|
    bad_counts = 0
    for r in mp:
        if delta_R(r) == 0:
            for m in minor_models(r):
                if not (m["C"] == m["A"] and m["X"] == m["Dv"] and m["M"] == len(r.L)):
                    bad_counts += 1
    allok &= _check("Theorem 4 counting identity (|C|=A,|X|=Dv,|M|=|L|) on all delta_R=0 minor models", bad_counts == 0, "violations=%d" % bad_counts)
    defects = {"K5": nx.complete_graph(5), "K33": nx.complete_bipartite_graph(3, 3)}
    viol_np = viol_b = viol_core = ev_total = 0
    for dn, D in defects.items():
        base = graph_to_edges(D)
        core0 = reduced_core(base)
        for r in mp:
            for _ in range(3):
                cur = list(base)
                b0 = beta1_of([(a, b) for a, b, _ in cur])[0]
                nev = 0
                for _step in range(6):
                    ms = find_matches(cur, r)
                    if not ms:
                        break
                    cur, _ = apply_matches(cur, [rng.choice(ms)], r, next_vertex_id(cur))
                    nev += 1
                    ev_total += 1
                    if is_planar(simple_graph([(a, b) for a, b, _ in cur])):
                        viol_np += 1
                    bt = beta1_of([(a, b) for a, b, _ in cur])[0]
                    if bt != b0 + nev * delta_R(r):
                        viol_b += 1
                    if is_path_preserving(r) and delta_R(r) == 0 and not nx.is_isomorphic(core0, reduced_core(cur)):
                        viol_core += 1
    allok &= _check("Theorem 3: minor-preserving rules never planarise K5/K33 (random sequences)", viol_np == 0, "events=%d violations=%d" % (ev_total, viol_np))
    allok &= _check("Theorem 1/Cor 1: beta1(H_t) = beta1(H_0) + (#events)*delta_R", viol_b == 0, "violations=%d" % viol_b)
    allok &= _check("Theorem 3(b): path-preserving, delta_R=0 => reduced core of H_t isomorphic to H_0", viol_core == 0, "violations=%d" % viol_core)

    # (2b) exhaustive one-step over ALL injective matches
    viol = 0
    planarising = Counter()
    for dn, D in defects.items():
        base = graph_to_edges(D)
        for r in rules:
            for m in find_matches(base, r):
                new, _ = apply_matches(base, [m], r, next_vertex_id(base))
                pl = is_planar(simple_graph([(a, b) for a, b, _ in new]))
                if pl and is_minor_preserving(r):
                    viol += 1
                if pl and not is_minor_preserving(r):
                    planarising[(dn, rule_str(r))] += 1
    allok &= _check("Theorem 3 one-step exhaustive over all matches on K5,K33", viol == 0, "violations=%d" % viol)
    print("   INFO non-minor-preserving closed rules that planarise K5 or K33 for >=1 match: %d (rule,defect) pairs" % len(planarising))

    # (2c) model-tracked label: beta1(D_t) constant and non-planar for every minor-preserving rule
    bad_track = ntr = 0
    for dn, D in defects.items():
        base = [(a, b, 1) for a, b, _ in graph_to_edges(D)]
        for r in mp:
            cur = list(base)
            for _t in range(3):
                cur, ev = generation_step(cur, r, rng, "random", "model")
                if ev == 0:
                    break
                ntr += 1
                dpairs = [(a, b) for a, b, m in cur if m & 1]
                if beta1_of(dpairs)[0] != beta1_of([(a, b) for a, b, _ in base])[0] or is_planar(simple_graph(dpairs)):
                    bad_track += 1
                if len(cur) > 400:
                    break
    allok &= _check("model-tracked defect label: beta1(D_t) constant, D_t non-planar (isolated K5/K33, 3 gens)", bad_track == 0 and ntr > 0, "steps=%d violations=%d" % (ntr, bad_track))

    # (3) explicit counterexample R_c
    Rc = parse_rule("{{1,2},{1,3}}->{{1,2},{2,3},{2,4}}")
    base = graph_to_edges(defects["K5"])
    ok = is_closed(Rc) and delta_R(Rc) == 0 and not is_path_preserving(Rc)
    found = False
    for m in find_matches(base, Rc):
        new, _ = apply_matches(base, [m], Rc, next_vertex_id(base))
        pr = [(a, b) for a, b, _ in new]
        if beta1_of(pr)[0] == 6 and is_planar(simple_graph(pr)):
            found = True
    allok &= _check("Counterexample R_c: closed, delta_R=0, beta1 stays 6, yet K5 -> planar", ok and found)

    # (4) Mayer-Vietoris additivity defect
    trials = bad = 0
    for _ in range(400):
        G = _random_connected(rng, 6, 10)
        es = sorted(G.edges())
        A = {e for e in es if rng.random() < 0.6}
        B = {e for e in es if rng.random() < 0.6}
        if not A or not B:
            continue
        U = A | B
        bA, bB, bU = (beta1_of(x)[0] for x in (A, B, U))
        VA = {v for e in A for v in e}
        VB = {v for e in B for v in e}
        Ig = nx.Graph()
        Ig.add_nodes_from(VA & VB)          # intersection subgraph: shared VERTICES and shared edges
        Ig.add_edges_from(A & B)
        cI = nx.number_connected_components(Ig)
        bI = Ig.number_of_edges() - Ig.number_of_nodes() + cI
        cA, cB, cU = (beta1_of(x)[3] for x in (A, B, U))
        # beta1(U) = b(A)+b(B)-b(I) + [c(I)+c(U)-c(A)-c(B)]   (derived from chi additivity)
        trials += 1
        if bU != bA + bB - bI + (cI + cU - cA - cB):
            bad += 1
    allok &= _check("Lemma 4 (Mayer-Vietoris identity for graphs)", bad == 0 and trials > 100, "trials=%d violations=%d" % (trials, bad))
    A = nx.complete_graph(5)
    B = nx.relabel_nodes(nx.complete_graph(5), {0: 0, 1: 1, 2: 2, 3: "a", 4: "b"})
    eA = {tuple(sorted(e, key=str)) for e in A.edges()}
    eB = {tuple(sorted(e, key=str)) for e in B.edges()}
    bs = [beta1_of(list(x))[0] for x in (eA, eB, eA | eB)]
    allok &= _check("Counterexample: two K5 sharing a triangle: beta1 = 6,6,11 (not 12)", bs == [6, 6, 11], str(bs))

    print("SELFTEST", "PASSED" if allok else "FAILED")
    return 0 if allok else 1


# --------------------------------------------------------------------------- CLI
def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_census(q):
        q.add_argument("--sigs", nargs="+", default=["2_3"], help="e.g. 2_3 2_4")
        q.add_argument("--max-new", type=int, default=None, help="override max new vertices (default 1 for k=3, 2 for k=4)")
        q.add_argument("--jobs", type=int, default=1)
        q.add_argument("--conf-depth", type=int, default=3, help="generations of the test-state evolution (paper target <=10; cost explodes)")
        q.add_argument("--conf-max-edges", type=int, default=14)
        q.add_argument("--conf-pairs", type=int, default=12)
        q.add_argument("--conf-join", type=int, default=2, help="joining depth J")
        q.add_argument("--conf-branch", type=int, default=30)

    def add_run(q):
        q.add_argument("--backgrounds", nargs="+", default=["grid2d", "grid3d"], choices=["grid2d", "grid3d"])
        q.add_argument("--defects", nargs="+", default=["K5", "K33"], choices=["K5", "K33"])
        q.add_argument("--seps", nargs="+", type=int, default=[3])
        q.add_argument("--orders", nargs="+", default=["random", "canonical"], choices=["random", "canonical"])
        q.add_argument("--tracking", choices=["cone", "model"], default="cone", help="defect-label transport: causal-cone OR (default) or minor-model transport (minor-preserving rules only)")
        q.add_argument("--tmax", type=int, default=10)
        q.add_argument("--max-edges", type=int, default=40000)
        q.add_argument("--n2d", type=int, default=10)
        q.add_argument("--n3d", type=int, default=6)
        q.add_argument("--orc-sample", type=int, default=150)
        q.add_argument("--orc-defect-cap", type=int, default=60)
        q.add_argument("--no-orc", action="store_true")
        q.add_argument("--no-dim", action="store_true")
        q.add_argument("--dim-rmax", type=int, default=8)
        q.add_argument("--global-planarity", action="store_true", help="also test global planarity on 3D backgrounds (always nonplanar for n>=3)")
        q.add_argument("--planarity-max-edges", type=int, default=30000)

    s = sub.add_parser("selftest")
    s.set_defaults(func=cmd_selftest)
    c = sub.add_parser("census")
    c.add_argument("--out", default="results")
    add_census(c)
    c.set_defaults(func=lambda a: (cmd_census(a), 0)[1])
    r = sub.add_parser("run")
    r.add_argument("--out", default="results")
    r.add_argument("--rule", nargs="+", required=True, help='e.g. "{{1,2},{1,3}}->{{1,2},{2,3},{2,4}}"')
    add_run(r)
    r.set_defaults(func=lambda a: (cmd_run(a), 0)[1])
    al = sub.add_parser("all")
    al.add_argument("--out", default="results")
    al.add_argument("--max-rules", type=int, default=0)
    add_census(al)
    add_run(al)
    al.set_defaults(func=lambda a: (cmd_all(a), 0)[1])
    return p


if __name__ == "__main__":
    a = build_parser().parse_args()
    np.random.seed(a.seed)
    random.seed(a.seed)
    sys.exit(a.func(a))
