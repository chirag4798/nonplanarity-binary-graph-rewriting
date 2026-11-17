#!/usr/bin/env python3
"""Queries the public zbMATH Open API (https://api.zbmath.org) and prints year, authors, title, source and DOI link of the first 8 hits
per query.  Usage: python3 zbmath_search.py 'query one' 'query two' ...   (uses curl; results change over time).
Used for the literature check of 17 November 2026; the output is kept in results/literature/zbmath_queries.txt."""
import json, sys, urllib.parse, urllib.request, time
qs = sys.argv[1:]
for q in qs:
    url = "https://api.zbmath.org/v1/document/_search?search_string=%s&page=0&results_per_page=8" % urllib.parse.quote(q)
    import subprocess
    try:
        d = json.loads(subprocess.run(["curl", "-s", "-m", "40", "-A", "Mozilla/5.0 (research script)", url], capture_output=True, text=True).stdout)
    except Exception as e:
        print("ERR", q, e); continue
    st = d.get("status", {})
    print("== %s  (total %s)" % (q, st.get("nr_total_results")))
    for r in (d.get("result") or []):
        au = "; ".join(a["name"] for a in r["contributors"]["authors"][:3])
        src = r.get("source", {}).get("source", "") if isinstance(r.get("source"), dict) else ""
        print(" - %s | %s | %s | %s | %s" % (r["year"], au, r["title"]["title"][:110], src[:70], r.get("links", [{}])[0].get("url", "") if r.get("links") else ""))
    time.sleep(2)
