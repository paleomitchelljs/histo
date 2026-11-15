#!/usr/bin/env python3
"""Extract the Woodward et al. 2015 Maiasaura growth data from the supplement docx.

The Dryad deposit (doi:10.5061/dryad.7vf45) ships these as Word documents, not
CSVs, so we parse the Word tables directly (<w:tbl>/<w:tr>/<w:tc>, namespace-aware).

Outputs:
  data_raw/external/maiasaura/suppl_table1_crosssectional.csv  (faithful)
  data_raw/external/maiasaura/suppl_table2_longitudinal.csv    (faithful)
  data_clean/maiasaura_lag_long.csv     modeling-ready longitudinal (circ in mm)
  data_clean/lag_circumference_long.csv UNION of theropod femora + Maiasaura tibiae

Run with a single integer arg to instead dump that table index as TSV (debug).
"""
import sys, os, zipfile, re, csv
import xml.etree.ElementTree as ET

ROOT = "/Users/jmitchell/Library/CloudStorage/Dropbox/Research/00_current_Histology"
DOCX = os.path.join(ROOT, "data_raw/supplements/Supplementary+Figures+and+Tables.docx")
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def cell_text(tc):
    return re.sub(r"\s+", " ", "".join(t.text or "" for t in tc.iter(W + "t"))).strip()


def tables(docx):
    with zipfile.ZipFile(docx) as z:
        root = ET.fromstring(z.read("word/document.xml"))
    out = []
    for tbl in root.iter(W + "tbl"):
        out.append([[cell_text(tc) for tc in tr.findall(W + "tc")]
                    for tr in tbl.iter(W + "tr")])
    return out


def num(s):
    """Float from a messy cell: strip flag glyphs/units; non-numeric -> None."""
    if s is None:
        return None
    s = s.replace("*", "").replace("^", "").replace("□", "").replace("◊", "")
    m = re.search(r"-?\d+\.?\d*", s)
    return float(m.group()) if m else None


def base_code(s):
    """Strip trailing flag glyphs from a specimen code."""
    return re.sub(r"[\*\^□◊\s]+$", "", s).strip()


tbls = tables(DOCX)

if len(sys.argv) > 1:                      # debug dump
    for r in tbls[int(sys.argv[1])]:
        print("\t".join(r))
    sys.exit()

os.makedirs(os.path.join(ROOT, "data_raw/external/maiasaura"), exist_ok=True)
os.makedirs(os.path.join(ROOT, "data_clean"), exist_ok=True)

# ---- Suppl Table 1: cross-sectional (one row per individual) -------
t1 = tbls[0][2:]                            # drop caption + header rows
cross = []
for r in t1:
    r = (r + [""] * 6)[:6]
    code = base_code(r[0])
    if not code:
        continue
    cross.append({
        "specimen": code,
        "flags": re.sub(r"[A-Za-z0-9\s]", "", r[0]),          # the * ^ symbols
        "tibia_length_cm": num(r[1]),
        "min_diaph_circ_cm": num(r[2]),
        "calc_femur_length_cm": num(r[3]),
        "lag_count": num(r[4]),
        "efs_lags": num(r[5]),
    })
with open(os.path.join(ROOT, "data_raw/external/maiasaura/suppl_table1_crosssectional.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(cross[0].keys())); w.writeheader(); w.writerows(cross)

# ---- Suppl Table 2: longitudinal (perinate -> LAGs -> periosteal) --
t2 = tbls[1][2:]
long_raw, cur = [], None
for r in t2:
    r = (r + [""] * 7)[:7]
    if base_code(r[0]):
        cur = base_code(r[0])
    loc = r[1].strip()
    if not loc:
        continue
    long_raw.append({"specimen": cur, "location": loc,
                     "total_area_cm2": num(r[2]), "circ_cm": num(r[3]),
                     "calc_tibia_len_cm": num(r[4]), "calc_femur_len_cm": num(r[5]),
                     "mass_dme_kg": num(r[6])})
with open(os.path.join(ROOT, "data_raw/external/maiasaura/suppl_table2_longitudinal.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(long_raw[0].keys())); w.writeheader(); w.writerows(long_raw)

# ---- modeling-ready longitudinal: assign age, convert cm -> mm -----
# age convention: perinate = 0; LAG n = n; periosteal surface = (max LAG in
# specimen) + 0.5 (size at death, just past the last ring). Medullary cavity is
# an internal remodeling boundary, not an external size-at-age -> excluded.
by_spec = {}
for row in long_raw:
    by_spec.setdefault(row["specimen"], []).append(row)

maia = []
for spec, rows in by_spec.items():
    lag_idxs = [int(re.search(r"\d+", x["location"]).group())
                for x in rows if x["location"].lower().startswith("lag")]
    maxlag = max(lag_idxs) if lag_idxs else 0
    for x in rows:
        loc = x["location"].lower()
        if x["circ_cm"] is None:
            continue
        if loc.startswith("perinate"):
            pos, age = "perinate", 0.0
        elif loc.startswith("lag"):
            pos, age = "lag", float(re.search(r"\d+", loc).group())
        elif loc.startswith("periosteal"):
            pos, age = "periosteal", maxlag + 0.5
        else:                                # medullary cavity -> skip
            continue
        maia.append({
            "clade": "Ornithischia", "family": "Hadrosauridae",
            "species": "Maiasaura_peeblesorum", "genus": "Maiasaura",
            "specimen": spec, "element": "tibia", "zone": pos,
            "lag_index": "" if pos != "lag" else int(age),
            "age": age, "circ_mm": round(x["circ_cm"] * 10, 1),
            "is_efs": False, "source": "Woodward et al. 2015 (Paleobiology)"})
maia.sort(key=lambda d: (d["specimen"], d["age"]))
COLS = ["clade", "family", "species", "genus", "specimen", "element", "zone",
        "lag_index", "age", "circ_mm", "is_efs", "source"]
with open(os.path.join(ROOT, "data_clean/maiasaura_lag_long.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS); w.writeheader(); w.writerows(maia)

# ---- unified long file: theropod femora + Maiasaura tibiae ---------
femur = []
with open(os.path.join(ROOT, "data_clean/lag_femur_long.csv")) as f:
    for d in csv.DictReader(f):
        d["source"] = "D'Emic et al. 2023 (Science)"
        femur.append({c: d.get(c, "") for c in COLS})
with open(os.path.join(ROOT, "data_clean/lag_circumference_long.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS); w.writeheader()
    w.writerows(femur); w.writerows(maia)

# ---- summary -------------------------------------------------------
ages = [d["age"] for d in maia]
circs = [d["circ_mm"] for d in maia]
print(f"Suppl Table 1 (cross-sectional): {len(cross)} individuals")
print(f"Suppl Table 2 (longitudinal):    {len(by_spec)} specimens, {len(maia)} usable size-at-age points")
print(f"  age range {min(ages):.1f}-{max(ages):.1f} yr | circ {min(circs):.0f}-{max(circs):.0f} mm")
print(f"Unified file: {len(femur)} theropod-femur + {len(maia)} Maiasaura-tibia = {len(femur)+len(maia)} rows")
print("Wrote data_clean/maiasaura_lag_long.csv and data_clean/lag_circumference_long.csv")
