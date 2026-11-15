#!/usr/bin/env python3
"""Extract Griebeler, Klein & Sander 2013 (PLOS ONE e0067012) Table S1 sauropod
growth data from the .doc supplement.

The supplement is legacy binary Word (.doc); textutil mangles the table, so we
convert with wvHtml (preserves <table>) and parse that. Table S1 = per growth
annulus, % final bone length and mass (kg) for 7 sauropodomorph specimens.

We model on the linear dimension (% final bone length); mass kept for reference.
Output: data_raw/external/griebeler/griebeler_sauropod_long.csv
"""
import os, re, subprocess, csv, tempfile
from html.parser import HTMLParser

ROOT = os.environ.get("HISTO_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DOC = os.path.join(ROOT, "data_raw/supplements/pone.0067012.s001.doc")
OUTDIR = os.path.join(ROOT, "data_raw/external/griebeler"); os.makedirs(OUTDIR, exist_ok=True)
HTML = os.path.join(OUTDIR, "griebeler_s1.html")

subprocess.run(["wvHtml", DOC, HTML], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

class T(HTMLParser):
    def __init__(s): super().__init__(); s.rows=[]; s.cur=None; s.cell=None
    def handle_starttag(s,t,a):
        if t=='tr': s.cur=[]
        elif t in('td','th'): s.cell=[]
    def handle_endtag(s,t):
        if t=='tr' and s.cur is not None: s.rows.append(s.cur); s.cur=None
        elif t in('td','th') and s.cell is not None and s.cur is not None:
            s.cur.append(re.sub(r'\s+',' ',''.join(s.cell)).strip()); s.cell=None
    def handle_data(s,d):
        if s.cell is not None: s.cell.append(d)

p=T(); p.feed(open(HTML, encoding='utf-8', errors='ignore').read())
rows=[r for r in p.rows if len(r)>=15]            # the 22 annulus data rows (No + 7 pairs)
data=[r for r in rows if re.match(r'^\d+$', r[0])]

# specimen i (1..7) lives at columns (pct=2i-1, kg=2i); metadata per specimen:
specs=[
 dict(species="Plateosaurus_engelhardti", specimen="Plateosaurus_engelhardti", genus="Plateosaurus", family="Plateosauridae", clade="Sauropodomorpha"),
 dict(species="Mamenchisauridae_indet",   specimen="mamenchisaurid_indet",     genus="Mamenchisaurus", family="Mamenchisauridae", clade="Sauropoda"),
 dict(species="Apatosaurus_sp",           specimen="Apatosaurus_SMA_0014",     genus="Apatosaurus",  family="Diplodocidae", clade="Sauropoda"),
 dict(species="Apatosaurus_sp",           specimen="Apatosaurus_BYU_601-17328",genus="Apatosaurus",  family="Diplodocidae", clade="Sauropoda"),
 dict(species="Diplodocidae_indet",       specimen="diplodocid_MfN_R2625",     genus="Diplodocus",   family="Diplodocidae", clade="Sauropoda"),
 dict(species="Diplodocidae_indet",       specimen="diplodocid_MfN_RNW4",      genus="Diplodocus",   family="Diplodocidae", clade="Sauropoda"),
 dict(species="Camarasaurus_sp",          specimen="Camarasaurus_sp",          genus="Camarasaurus", family="Camarasauridae", clade="Macronaria"),
]
def num(x):
    x=(x or "").strip().replace("%","");  m=re.search(r"-?\d+\.?\d*", x); return float(m.group()) if m else None

out=[]
for r in data:
    ann=int(r[0])
    for i,meta in enumerate(specs):
        pct=num(r[2*i+1]) if 2*i+1 < len(r) else None
        kg =num(r[2*i+2]) if 2*i+2 < len(r) else None
        if pct is None: continue
        out.append({**meta, "element":"long_bone", "annulus":ann, "age":ann,
                    "pct_bone_length":pct, "mass_kg":kg,
                    "source":"Griebeler, Klein & Sander 2013 (PLOS ONE e0067012)"})

cols=["clade","family","species","genus","specimen","element","annulus","age","pct_bone_length","mass_kg","source"]
with open(os.path.join(OUTDIR,"griebeler_sauropod_long.csv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(out)

import collections
bs=collections.Counter((o["species"],o["specimen"]) for o in out)
print(f"{len(out)} annulus points | {len(set(o['species'] for o in out))} species | {len(bs)} specimens")
for (sp,sc),n in bs.items(): print(f"  {sp:28s} {sc:28s} {n} annuli")
