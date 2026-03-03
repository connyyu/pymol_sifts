from pymol import cmd
import requests
from collections import defaultdict

def pdb_sort_key(pdb_id):
    number_part = int(pdb_id[0])
    letter_part = pdb_id[1:]
    return (number_part, letter_part)

def get_uniprot_mappings_from_pdb(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"

    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
    except requests.RequestException as e:
        print(f"API request for {pdb_id} failed: {e}")
        return [], []

    pdb_entry = data.get(pdb_id.lower())
    if not pdb_entry:
        return [], []

    dbref_entries = []

    for ac, details in pdb_entry.get("UniProt", {}).items():
        identifier = details.get("identifier", "")
        for mapping in details.get("mappings", []):
            chain_id = mapping.get("chain_id")
            start = mapping.get("unp_start")
            end = mapping.get("unp_end")

            if start is None or end is None:
                continue

            dbref_entries.append({
                "model_identifier": pdb_id.lower(),
                "chain_id": chain_id,
                "uniprot_start": int(start),
                "uniprot_end": int(end),
                "ac": ac,
                "id": identifier,
            })

    return [pdb_id.lower()], dbref_entries

def merge_dbref_lines(dbref_lines):
    if not dbref_lines:
        return []

    dbref_lines = sorted(
        dbref_lines,
        key=lambda x: (
            pdb_sort_key(x["model_identifier"]),
            x["chain_id"],
            x["ac"],
            x["id"],
            x["uniprot_start"]
        ),
    )

    merged = []
    current = dbref_lines[0]

    for nxt in dbref_lines[1:]:
        same_entry = (
            current["model_identifier"] == nxt["model_identifier"] and
            current["chain_id"] == nxt["chain_id"] and
            current["ac"] == nxt["ac"] and
            current["id"] == nxt["id"] and
            current["uniprot_end"] + 1 == nxt["uniprot_start"]
        )

        if same_entry:
            current["uniprot_end"] = nxt["uniprot_end"]
        else:
            merged.append(current)
            current = nxt

    merged.append(current)
    return merged

def sifts(pdb_list="", _self=None):
    obj_names = pdb_list.split()
    
    if not obj_names:  # Also change this check
        print("Usage: sifts <object_name> [<object_name2> ...]")
        print("<object_name> should match with PDB identifier.")
        return

    # Delete old selections for UniProt IDs to avoid conflicts
    existing_selections = cmd.get_names("selections")
    for sel in existing_selections:
        if sel.isupper():  # crude filter for previous UniProt selections
            cmd.delete(sel)

    # Reset colors of all involved objects
    for obj in obj_names:
        if obj in cmd.get_object_list():
            cmd.color("white", f"model {obj}")

    all_dbrefs = []

    for obj_name in obj_names:
        if obj_name not in cmd.get_object_list():
            print(f"Object '{obj_name}' not found in PyMOL, skipping.")
            continue

        pdb_code = obj_name.lower()
        identifiers, dbref_data = get_uniprot_mappings_from_pdb(pdb_code)

        if not identifiers or not dbref_data:
            print(f"No UniProt mappings found for {pdb_code}.")
            continue

        merged_dbref = merge_dbref_lines(dbref_data)
        all_dbrefs.extend(merged_dbref)

    if not all_dbrefs:
        print("No DBREF entries found across the provided objects.")
        return

    # Group by UniProt ID
    uniprot_groups = defaultdict(list)
    for entry in all_dbrefs:
        uniprot_groups[entry["id"]].append(entry)

    # Predefined colors
    colors = [
        "carbon", "cyan", "lightmagenta", "yellow", "salmon", "slate", "orange",
        "deepteal", "violetpurple", "marine", "olive", "smudge",
        "teal", "wheat", "lightpink", "skyblue"
    ]
    color_index = 0

    # Apply base coloring individually
    for obj in obj_names:
        cmd.util.cbaw(f"model {obj}")

    for uniprot_id, entries in uniprot_groups.items():
        selection_parts = []
        model_chains = defaultdict(list)

        for e in entries:
            model_chains[e["model_identifier"]].append(e["chain_id"])

        for model, chains in model_chains.items():
            chains_str = "+".join(sorted(set(chains)))
            selection_parts.append(f'(chain {chains_str} and model "{model}")')

        selection_expr = " or ".join(selection_parts)
        selection_name = uniprot_id.upper()  # ensure uppercase for consistency
        cmd.select(selection_name, selection_expr)

        color = colors[color_index % len(colors)]
        cmd.color(color, selection_name)
        print(f"sele {selection_name}, {selection_expr}; color {color}, {selection_name};")

        color_index += 1

    cmd.deselect()
    cmd.util.cnc

    print(f"\n{', '.join(obj_names)} annotated with SIFTS mappings from UniProt and PDBe.")

cmd.extend("sifts", sifts)
