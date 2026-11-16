"""Generalisation of the rule machinery to left sides with m >= 1 edges (used for the general-signature search).
Independent of the 2-edge code paths in oligon_experiment.py except for the Rule type and minor_models()."""
import itertools
from collections import defaultdict
import oligon_experiment as oe

Rule = oe.Rule


def _connected(edges):
    adj = defaultdict(set)
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    if not adj:
        return False
    start = next(iter(adj)); seen = {start}; st = [start]
    while st:
        v = st.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w); st.append(w)
    return len(seen) == len(adj)


def _relabel(edges):
    m = {}; out = []
    for a, b in edges:
        for v in (a, b):
            if v not in m:
                m[v] = len(m)
        out.append((m[a], m[b]))
    return tuple(out)


def lhs_classes(m):
    """Connected loop-free left sides with m edges, up to renaming and edge order."""
    pairs = list(itertools.permutations(range(m + 1), 2))
    seen = set()
    for L in itertools.combinations_with_replacement(pairs, m):
        if not _connected(L):
            continue
        seen.add(min(_relabel(o) for o in itertools.permutations(L)))
    return sorted(seen)


def enumerate_rules(m, k, max_new):
    """Closed (vertex-preserving, connected-K) loop-free rules m_2 -> k_2 up to renaming/edge order."""
    out = []
    for L in lhs_classes(m):
        n_old = len({v for e in L for v in e})
        auts = []
        for perm in itertools.permutations(range(n_old)):
            mp = dict(enumerate(perm))
            if sorted((mp[a], mp[b]) for a, b in L) == sorted(L):
                auts.append(mp)
        for n_new in range(0, max_new + 1):
            nv = n_old + n_new
            if nv > k + 1:               # connected K has at most k+1 vertices
                break
            pairs = [(a, b) for a in range(nv) for b in range(nv) if a != b]
            newperms = list(itertools.permutations(range(n_old, nv)))
            seen = set()
            for K in itertools.combinations_with_replacement(pairs, k):
                used = {v for e in K for v in e}
                if len(used) != nv or not _connected(K):
                    continue
                best = None
                for mp in auts:
                    for npm in newperms:
                        full = dict(mp); full.update(dict(zip(range(n_old, nv), npm)))
                        img = tuple(sorted((full[a], full[b]) for a, b in K))
                        if best is None or img < best:
                            best = img
                if best in seen:
                    continue
                seen.add(best)
                out.append(Rule(tuple(L), best, n_old, n_new))
    return out


def find_matches(edges, rule):
    """All injective matches of an m-edge left side: list of (edge-index tuple aligned with rule.L, binding)."""
    L = rule.L; m = len(L)
    inc = defaultdict(list)
    for i, (a, b, *_ ) in enumerate(edges):
        inc[a].append(i); inc[b].append(i)
    order = [0]; bound = set(L[0])
    while len(order) < m:
        for j in range(m):
            if j not in order and (L[j][0] in bound or L[j][1] in bound):
                order.append(j); bound |= set(L[j]); break
        else:
            raise ValueError("left side not connected")
    res = []

    def rec(pos, used, bind, idx):
        if pos == m:
            res.append((tuple(idx[j] for j in range(m)), dict(bind)))
            return
        li = order[pos]; p, q = L[li]
        if p in bind:
            cand = inc[bind[p]]
        elif q in bind:
            cand = inc[bind[q]]
        else:
            cand = range(len(edges))
        for e in list(cand):
            if e in used:
                continue
            a, b = edges[e][0], edges[e][1]
            if a == b:
                continue
            nb = dict(bind); ok = True
            for var, val in ((p, a), (q, b)):
                if var in nb:
                    if nb[var] != val:
                        ok = False; break
                else:
                    nb[var] = val
            if not ok or len(set(nb.values())) != len(nb):
                continue
            idx2 = dict(idx); idx2[li] = e
            rec(pos + 1, used | {e}, nb, idx2)

    rec(0, frozenset(), {}, {})
    return res


def apply_event(edges, match, rule, next_vid):
    idxs, bind = match
    m = dict(bind)
    for k in range(rule.n_old, rule.n_old + rule.n_new):
        m[k] = next_vid; next_vid += 1
    kept = [e for i, e in enumerate(edges) if i not in set(idxs)]
    return kept + [(m[u], m[v], 0) for (u, v) in rule.K]


def next_vid(edges):
    return 1 + max((max(a, b) for a, b, *_ in edges), default=-1)
