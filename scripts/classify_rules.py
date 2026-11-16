#!/usr/bin/env python3
"""
classify_rules.py -- exhaustive structural classification of minimal binary rules 2_2 -> k_2
acting on an ISOLATED Kuratowski defect H_0 in {K5, K33}.

Buckets (closed rules; delta_R from Theorem 1; 'minor-preserving' = Definition in the paper, generalising path-preserving):
  protected_frozen          minor-preserving, delta_R = 0  (Theorem 3: non-planar forever, all update orders;
                            Theorem 4: cycle space frozen, H_t homotopy equivalent to H_0)
  protected_cycle_growing   minor-preserving, delta_R > 0  (non-planar forever, but beta1 grows by delta_R per event)
  dissolving                NOT minor-preserving and a planar state was OBSERVED (exhaustive one-step over all
                            injective matches, or in sampled multi-generation runs <= tmax).  Certified by that state.
  no_dissolution_observed   NOT minor-preserving, applicable, but no planar state seen in the SAMPLED runs.
                            NOT claimed protected; measures the gap between sufficiency (Thm 3) and necessity.
  inapplicable_to_simple_H0 NOT minor-preserving, but L never matches in H_0 (e.g. L has parallel edges): vacuous.
  open_rule                 not closed (vertex loss and/or disconnected K): beta1 change not a rule constant.
  THEOREM_VIOLATION         must be empty; a planar state observed for a minor-preserving rule.

Usage: python classify_rules.py --out ../results/classification --jobs 8 [--sigs 2_3 2_4] [--tmax 10]
"""
import argparse, csv, os, random, sys
from collections import Counter
from multiprocessing import Pool
import networkx as nx

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import oligon_experiment as oe

DEFECTS = {"K5": nx.complete_graph(5), "K33": nx.complete_bipartite_graph(3, 3)}


def planar_state(edges):
    return oe.is_planar(oe.simple_graph([(a, b) for a, b, _ in edges]))


def connected_state(edges):
    g = oe.simple_graph([(a, b) for a, b, _ in edges])
    return g.number_of_nodes() > 0 and nx.is_connected(g)


def analyse(task):
    r, a = task
    out = dict(rule=oe.rule_str(r), n_old=r.n_old, n_new=r.n_new, closed=int(oe.is_closed(r)),
               delta_R=oe.delta_R(r), path_preserving=int(oe.is_path_preserving(r)), minor_preserving=int(oe.is_minor_preserving(r)))
    out["vertex_preserving"] = int(set(range(r.n_old)) <= {v for e in r.K for v in e})
    out["K_connected"] = int(oe._weakly_connected(r.K))
    first_planar = {}
    onestep_planar = 0
    n_matches0 = 0
    onestep_disconnect = 0
    for dn, D in DEFECTS.items():
        base = oe.graph_to_edges(D)
        for m in oe.find_matches(base, r):
            n_matches0 += 1
            new, _ = oe.apply_matches(base, [m], r, oe.next_vertex_id(base))
            if planar_state(new):
                onestep_planar += 1
            if not connected_state(new):
                onestep_disconnect += 1
        if not out["closed"] and not a.open_dynamics:
            continue
        for order in ("random", "canonical"):
            for s in range(a.seeds if order == "random" else 1):
                rng = random.Random(a.seed + 1000 * s)
                cur = list(base)
                for t in range(1, a.tmax + 1):
                    cur, ev = oe.generation_step(cur, r, rng, order)
                    if ev == 0:
                        break
                    if planar_state(cur):
                        first_planar.setdefault(dn, t)
                        break
                    if len(cur) > a.edge_cap:
                        break
    out["matches_on_H0"] = n_matches0
    out["onestep_planarising_matches"] = onestep_planar
    out["onestep_disconnecting_matches"] = onestep_disconnect
    out["first_planar_K5"] = first_planar.get("K5", "")
    out["first_planar_K33"] = first_planar.get("K33", "")
    observed = onestep_planar > 0 or bool(first_planar)
    if not out["closed"]:
        out["bucket"] = "open_rule"
    elif out["minor_preserving"]:
        if observed:
            out["bucket"] = "THEOREM_VIOLATION"
        else:
            out["bucket"] = "protected_frozen" if out["delta_R"] == 0 else "protected_cycle_growing"
    else:
        if observed:
            out["bucket"] = "dissolving"
        elif n_matches0 == 0:
            out["bucket"] = "inapplicable_to_simple_H0"
        else:
            out["bucket"] = "no_dissolution_observed"
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="../results/classification")
    p.add_argument("--sigs", nargs="+", default=["2_3", "2_4"])
    p.add_argument("--max-new3", type=int, default=1)
    p.add_argument("--max-new4", type=int, default=2)
    p.add_argument("--tmax", type=int, default=10)
    p.add_argument("--edge-cap", type=int, default=1500)
    p.add_argument("--seeds", type=int, default=3)
    p.add_argument("--seed", type=int, default=oe.DEFAULT_SEED)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--open-dynamics", action="store_true", help="also evolve non-closed rules")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tasks = []
    for sig in a.sigs:
        k = int(sig.split("_")[1][0])
        mx = a.max_new3 if k == 3 else a.max_new4
        closed = oe.enumerate_rules(k, mx, closed_only=True)
        allr = oe.enumerate_rules(k, mx, closed_only=False)
        cs = set(closed)
        print(sig, "closed:", len(closed), "all (incl. open):", len(allr), flush=True)
        for r in allr:
            tasks.append((sig, r))
    work = [(r, a) for _, r in tasks]
    if a.jobs > 1:
        with Pool(a.jobs) as pool:
            res = pool.map(analyse, work, chunksize=8)
    else:
        res = [analyse(w) for w in work]
    for (sig, _), o in zip(tasks, res):
        o["signature"] = sig
    fields = ["signature"] + [k for k in res[0] if k != "signature"]
    with open(os.path.join(a.out, "classification.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(res)
    cnt = Counter((o["signature"], o["bucket"], o["delta_R"] if o["closed"] else "n/a") for o in res)
    with open(os.path.join(a.out, "classification_counts.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["signature", "bucket", "delta_R", "count"])
        for k in sorted(cnt, key=str):
            w.writerow(list(k) + [cnt[k]])
    print("%-6s %-26s %-6s %s" % ("sig", "bucket", "dR", "count"))
    for k in sorted(cnt, key=str):
        print("%-6s %-26s %-6s %d" % (k[0], k[1], k[2], cnt[k]))
    bad = [o for o in res if o["bucket"] == "THEOREM_VIOLATION"]
    print("THEOREM VIOLATIONS:", len(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
