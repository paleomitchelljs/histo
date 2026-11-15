"""
Reshape measures.csv from wide hierarchical format to long form.

Input:  measures.csv  (species/specimen blocks with LAG circumferences)
Output: measures_long.csv (one row per observation)

Columns: species, specimen, bone, lag_index, circumference_mm, est_mass_kg, source, notes
"""

import csv
import re

INPUT = 'measures.csv'
OUTPUT = 'measures_long.csv'


def is_source(text):
    """Heuristic: does this column-B entry look like a literature citation?"""
    if not text:
        return False
    if re.search(r'\(\d{4}', text):
        return True
    if 'this study' in text.lower():
        return True
    if 'pers. comm.' in text.lower():
        return True
    return False


def parse_float(s):
    """Parse a numeric string, returning None for blanks / N/A / 'maximum'."""
    if not s:
        return None
    s = s.strip().replace('%', '').replace(',', '')
    if not s or s.lower() in ('n/a', 'maximum'):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def main():
    with open(INPUT, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    # --- Split rows into blocks separated by empty rows ---
    blocks = []
    current_block = []
    for row in rows:
        while len(row) < 13:
            row.append('')
        if not any(c.strip() for c in row[:13]):
            if current_block:
                blocks.append(current_block)
                current_block = []
        else:
            current_block.append(row)
    if current_block:
        blocks.append(current_block)

    # --- Process each block ---
    records = []
    current_species = None

    for block in blocks:

        # ---- Identify species ----
        species = None
        for row in block:
            col_a = row[0].strip()
            # Skip notes that ended up in the species column (e.g. line 46)
            if col_a and not col_a.lower().startswith('note'):
                species = col_a
                break
        if species:
            current_species = species
        else:
            species = current_species

        # ---- Identify bone element ----
        bone = None
        for row in block:
            col_c = row[2].strip().lower()
            if col_c in ('femur', 'tibia'):
                bone = col_c
                break

        # ---- Histological vs composite ----
        # Histological blocks have circumference values in column 5;
        # composite blocks (Albertosaurus-type) have only mass in column 6.
        has_circ = any(
            parse_float(r[5]) is not None
            for r in block
            if r[5].strip().lower() != 'maximum'
        )

        if has_circ:
            # ============ HISTOLOGICAL BLOCK ============
            # Scan column B for specimen ID(s) and source(s)
            specimen = None
            sources = []
            for row in block:
                col_b = row[1].strip()
                if not col_b:
                    continue
                if is_source(col_b):
                    sources.append(col_b)
                elif specimen is None:
                    specimen = col_b
            source = '; '.join(sources) if sources else ''

            # Collect orphan notes (notes on rows without data)
            orphan_notes = []
            for row in block:
                col_f = row[5].strip()
                col_g = row[6].strip()
                col_m = row[12].strip()
                col_a_note = row[0].strip() if row[0].strip().lower().startswith('note') else ''
                note_parts = [p for p in [col_m, col_a_note] if p]
                circ = parse_float(col_f)
                mass = parse_float(col_g)
                if note_parts and circ is None and mass is None:
                    orphan_notes.extend(note_parts)

            # Extract data rows
            lag_idx = 0
            first_data = True
            for row in block:
                col_f = row[5].strip()
                col_g = row[6].strip()
                col_m = row[12].strip()
                col_a_note = row[0].strip() if row[0].strip().lower().startswith('note') else ''

                # Skip maximum summary rows
                if col_f.lower() == 'maximum' or row[9].strip().lower() == 'maximum':
                    continue

                circ = parse_float(col_f)
                mass = parse_float(col_g)

                if circ is not None or mass is not None:
                    lag_idx += 1
                    # Build notes for this record
                    row_notes = [p for p in [col_m, col_a_note] if p]
                    if first_data and orphan_notes:
                        row_notes = orphan_notes + row_notes
                        first_data = False
                    else:
                        first_data = False
                    notes_str = '; '.join(row_notes)

                    records.append({
                        'species': species,
                        'specimen': specimen or '',
                        'bone': bone or '',
                        'lag_index': lag_idx,
                        'circumference_mm': circ if circ is not None else '',
                        'est_mass_kg': mass if mass is not None else '',
                        'source': source,
                        'notes': notes_str,
                    })
        else:
            # ============ COMPOSITE BLOCK ============
            # Each row may represent a different specimen; source on last citation row.
            block_source = None
            for row in block:
                col_b = row[1].strip()
                if is_source(col_b):
                    block_source = col_b

            lag_idx = 0
            for row in block:
                col_b = row[1].strip()
                col_g = row[6].strip()
                col_m = row[12].strip()

                if row[9].strip().lower() == 'maximum':
                    continue
                if col_g and col_g.lower() == 'maximum':
                    continue

                mass = parse_float(col_g)
                if mass is not None:
                    lag_idx += 1
                    row_specimen = col_b if (col_b and not is_source(col_b)) else ''
                    records.append({
                        'species': species,
                        'specimen': row_specimen,
                        'bone': '',
                        'lag_index': lag_idx,
                        'circumference_mm': '',
                        'est_mass_kg': mass,
                        'source': block_source or '',
                        'notes': col_m,
                    })

    # --- Write output ---
    fieldnames = ['species', 'specimen', 'bone', 'lag_index',
                  'circumference_mm', 'est_mass_kg', 'source', 'notes']
    with open(OUTPUT, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print(f'Wrote {len(records)} records to {OUTPUT}')

    # Summary stats
    species_set = set(r['species'] for r in records)
    specimens = set((r['species'], r['specimen']) for r in records if r['specimen'])
    print(f'  {len(species_set)} species, {len(specimens)} unique species-specimen pairs')
    histo = sum(1 for r in records if r['circumference_mm'] != '')
    comp = sum(1 for r in records if r['circumference_mm'] == '')
    print(f'  {histo} histological observations, {comp} composite observations')


if __name__ == '__main__':
    main()
