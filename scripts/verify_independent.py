#!/usr/bin/env python3
"""Independent verification of the two computer-dependent ingredients of the dichotomy theorem.

(A) Completeness of the rule enumeration: count the rules up to renaming by BURNSIDE'S LEMMA (no canonical forms, no
    search), per (n_old, n_new), and compare with the sizes of the enumerated sets.
(B) Correctness of the minor-preservation decision: a second, from-scratch decision procedure (plain loops and an explicit
    contraction of K, no networkx and no code shared with oligon_experiment.minor_models) is compared with the main one on
    ALL closed rules, and its branch assignments are validated SEMANTICALLY: for each rule with a model, random host graphs and
    matches are generated, the event is applied, branch sets are contracted in the actual output graph, and the result is
    checked to contain the host (as a simple graph, and as a multigraph for multigraph-minor-preserving rules)."""
import itertools, os, random, sys
from collections import Counter, defaultdict
from math import factorial
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import oligon_experiment as oe


# ----------------------------------------------------------------------------- (A) Burnside
def connected(edges):
    adj = defaultdict(set)
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    seen = set(); todo = [next(iter(adj))]
    while todo:
        v = todo.pop()
        if v in seen: continue
        seen.add(v); todo += list(adj[v] - seen)
    return len(seen) == len(adj)


def raw_multisets(nv, m, need_all):
    pairs = [(a, b) for a in range(nv) for b in range(nv) if a != b]
    out = []
    for es in itertools.combinations_with_replacement(pairs, m):
        if {v for e in es for v in e} == set(range(nv)) and connected(es):
            out.append(es)
    return out


def orbits(n_old, n_new, m, k):
    nv = n_old + n_new
    Ls = raw_multisets(n_old, m, True)
    Ks = raw_multisets(nv, k, True)
    total = 0; G = 0
    for po in itertools.permutations(range(n_old)):
        for pn in itertools.permutations(range(n_old, nv)):
            g = list(po) + list(pn); G += 1
            fl = sum(1 for L in Ls if sorted((g[a], g[b]) for a, b in L) == sorted(L))
            fk = sum(1 for K in Ks if sorted((g[a], g[b]) for a, b in K) == sorted(K))
            total += fl * fk
    assert total % G == 0
    return total // G


def check_enumeration():
    ok = True
    for m, k, mx in ((2, 3, 1), (2, 4, 2)):
        enum = Counter((r.n_old, r.n_new) for r in oe.enumerate_rules(k, mx))
        cnt = Counter()
        for n_old in (2, 3):
            for n_new in range(0, mx + 1):
                if n_old + n_new > k + 1: continue
                cnt[(n_old, n_new)] = orbits(n_old, n_new, m, k)
        cnt = Counter({a: b for a, b in cnt.items() if b})
        same = cnt == enum
        ok &= same
        print("(A) 2_2->%d_2: Burnside total %d | enumerated total %d | per-(n_old,n_new) equal: %s" % (k, sum(cnt.values()), sum(enum.values()), same))
        if not same: print("    burnside", dict(cnt), "enum", dict(enum))
    return ok


# ----------------------------------------------------------------------------- (B) second decision procedure
def decide(L, K, n_old, n_new, simple):
    """Return a branch assignment (tuple) or None.  Independent of oe.minor_models."""
    for assign in itertools.product(range(-1, n_old), repeat=n_new):
        g = list(range(n_old)) + list(assign)            # g[v] = branch of variable v, -1 = discarded
        good = True
        for a in range(n_old):                           # each branch set connected via K-edges inside it
            S = [v for v in range(len(g)) if g[v] == a]
            seen = {S[0]}; todo = [S[0]]
            while todo:
                v = todo.pop()
                for (x, y) in K:
                    for (p, q) in ((x, y), (y, x)):
                        if p == v and q in S and q not in seen:
                            seen.add(q); todo.append(q)
            if len(seen) != len(S): good = False; break
        if not good: continue
        contracted = [frozenset((g[x], g[y])) for (x, y) in K if g[x] != -1 and g[y] != -1 and g[x] != g[y]]
        cL = Counter(frozenset(e) for e in L)
        cC = Counter(contracted)
        need = Counter({e: 1 for e in cL}) if simple else cL
        if all(cC[e] >= c for e, c in need.items()):
            return assign
    return None


def semantic_test(r, assign, rng, trials, multi):
    n_old, n_new = r.n_old, r.n_new
    for _ in range(trials):
        n = rng.randint(max(n_old + 1, 5), 8)
        while True:
            m_ = rng.randint(n, min(n * (n - 1) // 2, 2 * n))
            pairs = rng.sample([(a, b) for a in range(n) for b in range(a + 1, n)], m_)
            if connected(pairs): break
        edges = [((a, b) if rng.random() < .5 else (b, a)) + (0,) for a, b in pairs]
        ms = oe.find_matches(edges, r)
        if not ms: continue
        mt = rng.choice(ms)
        nid = oe.next_vertex_id(edges)
        new, _ = oe.apply_matches(edges, [mt], r, nid)
        bind = mt[2]
        rep = {}                                         # contraction map on the OUTPUT graph
        for k in range(n_new):
            a = assign[k]
            rep[nid + k] = bind[a] if a >= 0 else None
        out = []
        for (x, y, _) in new:
            x, y = rep.get(x, x), rep.get(y, y)
            if x is None or y is None or x == y: continue
            out.append(frozenset((x, y)))
        host = [frozenset((x, y)) for (x, y, _) in edges]
        if multi:
            ch, co = Counter(host), Counter(out)
            if not all(co[e] >= c for e, c in ch.items()): return False
        else:
            if not set(host) <= set(out): return False
    return True


def check_decision():
    rules = oe.enumerate_rules(3, 1) + oe.enumerate_rules(4, 2)
    rng = random.Random(oe.DEFAULT_SEED)
    mism = 0; sm = mm = 0; sem_bad = 0
    for r in rules:
        a_sm = decide(r.L, r.K, r.n_old, r.n_new, True)
        a_mm = decide(r.L, r.K, r.n_old, r.n_new, False)
        if (a_sm is not None) != oe.is_simple_minor_preserving(r) or (a_mm is not None) != oe.is_minor_preserving(r):
            mism += 1
        sm += a_sm is not None; mm += a_mm is not None
        if a_sm is not None and not semantic_test(r, a_sm, rng, 12, False): sem_bad += 1
        if a_mm is not None and not semantic_test(r, a_mm, rng, 12, True): sem_bad += 1
    print("(B) closed rules %d | second procedure: simple-minor-preserving %d, minor-preserving %d | disagreements with main code: %d | "
          "semantic contraction failures: %d" % (len(rules), sm, mm, mism, sem_bad))
    return mism == 0 and sem_bad == 0


if __name__ == "__main__":
    a = check_enumeration(); b = check_decision()
    print("INDEPENDENT VERIFICATION", "PASSED" if a and b else "FAILED")
    sys.exit(0 if a and b else 1)
