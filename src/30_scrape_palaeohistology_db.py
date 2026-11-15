#!/usr/bin/env python3
"""Scrape the Palaeohistology literature database (T. Scheyer, Zurich):
https://www.developmental-palaeontology.net/palaeohistology/db.php
It is a BIBLIOGRAPHY of palaeohistology papers (reference, year, link/DOI,
five taxonomic-group flags, taxa sampled), not a measurement database.
Writes data_raw/external/scheyer_palaeohistology_db/palaeohistology_db.csv
plus the raw HTML pages. The site needs a session cookie; requests are spaced.
Needs internet. Usage: python3 src/30_scrape_palaeohistology_db.py
"""
import csv, html, os, re, subprocess, tempfile, time
BASE = "https://www.developmental-palaeontology.net/palaeohistology/"
OUT = "/Users/jmitchell/Library/CloudStorage/Dropbox/Research/00_current_Histology/data_raw/external/scheyer_palaeohistology_db"
os.makedirs(os.path.join(OUT, "html"), exist_ok=True)
JAR = os.path.join(tempfile.mkdtemp(), "cookies.txt")
def get(url):                                   # curl: system certificate store + cookie jar
    return subprocess.run(["curl", "-sS", "-L", "--max-redirs", "5", "-b", JAR, "-c", JAR,
                           "-A", "Mozilla/5.0 (compatible; research-client)", url],
                          check=True, capture_output=True).stdout.decode("utf-8", "replace")
def clean(x):
    return html.unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", x))).strip()
FLAGS = ["classical_overview", "fossil_synapsida", "parareptilia", "eureptilia_excl_dinos_birds", "eureptilia_dinos_birds"]
rows, total, row = [], None, 0
while True:
    page = get(BASE + "db.php?row=%d" % row)
    open(os.path.join(OUT, "html", "db_row%03d.html" % row), "w").write(page)
    if total is None:
        total = int(re.search(r"<b>(\d+)</b> records found", page).group(1))
    for blk in re.findall(r"<tbody class='body'>(.*?)</tbody>", page, flags=re.S):
        trs = re.findall(r"<tr[^>]*>(.*?)</tr>", blk, flags=re.S)
        tds = re.findall(r"<td([^>]*)>(.*?)</td>\s*(?=<td|$)", trs[1].strip(), flags=re.S)
        if len(tds) != 10:
            raise SystemExit("unexpected row layout at row=%d: %d cells" % (row, len(tds)))
        href = lambda c: (re.findall(r"href='([^']+)'", c) or [""])[0]
        notes = tds[9][1]
        full = re.search(r"class='hover'>\s*<div>(.*?)</div>", notes, flags=re.S)   # full text sits in the hover box
        rec = {"reference": clean(trs[0]), "year": clean(tds[0][1]), "link": href(tds[1][1]),
               "morphobank": href(tds[2][1]), "doi": href(tds[3][1]),
               "taxa_sampled_and_notes": clean(full.group(1) if full else notes)}
        for i, f in enumerate(FLAGS):
            rec[f] = int("background-color" in tds[4 + i][0])
        rows.append(rec)
    row += 50
    if row >= total:
        break
    time.sleep(1.5)
open(os.path.join(OUT, "html", "further_lit.html"), "w").write(get(BASE + "further_lit.php"))
cols = ["reference", "year", "link", "morphobank", "doi"] + FLAGS + ["taxa_sampled_and_notes"]
with open(os.path.join(OUT, "palaeohistology_db.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
print("records on site: %d | scraped: %d" % (total, len(rows)))
