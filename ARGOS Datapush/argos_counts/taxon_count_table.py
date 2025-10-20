#!/usr/bin/env python3
# Christie Woodside
"""
Count unique taxon IDs by kingdom from a TSV/CSV table.

Prefers in-file columns:
- taxonomy_id (string or int) -> unique key
- lineage (NCBI-style semicolon-separated) -> kingdom classification
- organism_name -> only used if taxid/lineage missing; we fallback to NCBI

Usage:
  python3 taxon_count_table.py argos.tsv
  python3 taxon_count_table.py argos.csv --sep , --organism-col organism_name --taxid-col taxonomy_id --lineage-col lineage
"""

import argparse
import csv
import json
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

try:
    import requests
except ImportError:
    requests = None  # We'll error only if we actually need the web fallback.

COMMON_ORG_COL_NAMES = {
    "organism", "organism_name", "scientific_name", "species",
    "taxonomy_name", "tax_name"
}
COMMON_TAXID_COL_NAMES = {"taxid", "tax_id", "taxonomy_id", "ncbi_tax_id"}
COMMON_LINEAGE_COL_NAMES = {"lineage", "ncbi_lineage"}

CACHE_DEFAULT_PATH = ".ncbi_tax_cache.json"


def load_cache(path: Path):
    if path.exists():
        try:
            return json.loads(path.read_text())
        except Exception:
            return {}
    return {}


def save_cache(path: Path, cache: dict):
    try:
        path.write_text(json.dumps(cache))
    except Exception:
        pass


def sniff_delimiter(file_path, user_sep=None):
    if user_sep:
        return "\t" if user_sep == "\\t" else user_sep
    with open(file_path, "r", newline="") as fh:
        sample = fh.read(4096)
    tabs = sample.count("\t")
    commas = sample.count(",")
    if tabs >= commas:
        return "\t"
    return ","


def find_col(header, user_col, common_names):
    # exact case-insensitive match first
    if user_col:
        for i, name in enumerate(header):
            if name.strip().lower() == user_col.strip().lower():
                return i, header[i]
    # common names
    lowered = [h.strip().lower() for h in header]
    for i, name in enumerate(lowered):
        if name in common_names:
            return i, header[i]
    # heuristic match
    for i, name in enumerate(lowered):
        if any(k in name for k in common_names):
            return i, header[i]
    return None, None


def classify_from_lineage(lineage_text: str):
    if not lineage_text:
        return "Unknown"
    # NCBI lineage contains capitalized clades/kingdoms
    if "Viruses" in lineage_text:
        return "virus"
    if "Fungi" in lineage_text:
        return "fungi"
    if "Bacteria" in lineage_text:
        return "bacteria"
    return "other"


def ncbi_esearch_taxid(organism, api_key=None, session=None):
    if requests is None:
        raise RuntimeError("requests is required for NCBI fallback but is not installed.")
    s = session or requests.Session()
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": "taxonomy", "term": organism, "retmode": "json"}
    if api_key:
        params["api_key"] = api_key
    resp = s.get(url, params=params, timeout=30)
    if resp.status_code == 429:
        time.sleep(1.0)
        resp = s.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    ids = data.get("esearchresult", {}).get("idlist") or []
    return ids[0] if ids else None


def ncbi_fetch_lineage(tax_id, api_key=None, session=None):
    if requests is None:
        raise RuntimeError("requests is required for NCBI fallback but is not installed.")
    s = session or requests.Session()
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {"db": "taxonomy", "id": tax_id, "retmode": "xml"}
    if api_key:
        params["api_key"] = api_key
    resp = s.get(url, params=params, timeout=30)
    if resp.status_code == 429:
        time.sleep(1.0)
        resp = s.get(url, params=params, timeout=30)
    resp.raise_for_status()
    try:
        root = ET.fromstring(resp.content)
        return root.findtext("Taxon/Lineage") or ""
    except Exception:
        return ""


def get_taxid_and_kingdom_from_row(row, org_idx, taxid_idx, lineage_idx,
                                   api_key, cache, session):
    """
    Strategy:
      1) If taxonomy_id present -> use it
         - If lineage present -> classify from lineage
         - Else -> fetch lineage from NCBI (fallback)
      2) Else if organism present -> resolve via NCBI
      3) Else -> unknown
    """
    org = (row[org_idx].strip() if (org_idx is not None and org_idx < len(row)) else "")
    raw_taxid = (row[taxid_idx].strip() if (taxid_idx is not None and taxid_idx < len(row)) else "")
    lineage = (row[lineage_idx].strip() if (lineage_idx is not None and lineage_idx < len(row)) else "")

    # Normalize taxid to string
    taxid = None
    if raw_taxid:
        # keep as string but ensure numeric-only content if present
        m = re.search(r"\d+", raw_taxid)
        if m:
            taxid = m.group(0)

    # If we have taxid:
    if taxid:
        if lineage:
            return taxid, classify_from_lineage(lineage)
        # need lineage → NCBI fallback (cached)
        cache_key = f"taxid:{taxid}"
        if cache_key in cache:
            return taxid, cache[cache_key].get("kingdom", "Unknown")
        try:
            fetched_lineage = ncbi_fetch_lineage(taxid, api_key=api_key, session=session)
        except Exception as e:
            print(f"--- EFetch error for taxid={taxid}: {e}")
            fetched_lineage = ""
        kingdom = classify_from_lineage(fetched_lineage)
        cache[cache_key] = {"kingdom": kingdom}
        return taxid, kingdom

    # Else no taxid: try organism via NCBI
    if org:
        cache_key = f"org:{org.lower()}"
        if cache_key in cache:
            entry = cache[cache_key]
            return entry.get("tax_id"), entry.get("kingdom")
        tax_from_name = None
        try:
            tax_from_name = ncbi_esearch_taxid(org, api_key=api_key, session=session)
        except Exception as e:
            print(f"--- ESearch error for {org!r}: {e}")
        if not tax_from_name:
            cache[cache_key] = {"tax_id": None, "kingdom": "Unknown"}
            return None, "Unknown"
        # fetch lineage to classify
        lineage2 = ""
        try:
            lineage2 = ncbi_fetch_lineage(tax_from_name, api_key=api_key, session=session)
        except Exception as e:
            print(f"--- EFetch error for taxid={tax_from_name} ({org!r}): {e}")
        kingdom = classify_from_lineage(lineage2)
        cache[cache_key] = {"tax_id": tax_from_name, "kingdom": kingdom}
        return tax_from_name, kingdom

    return None, "Unknown"


def main():
    p = argparse.ArgumentParser(description="Count unique taxa by kingdom from a TSV/CSV.")
    p.add_argument("input_file", help="Path to input TSV/CSV with a header row")
    p.add_argument("--sep", choices=[",", "\\t"], help="Field delimiter (default: auto-detect)")
    p.add_argument("--organism-col", help="Header name of organism (default: auto-detect)")
    p.add_argument("--taxid-col", help="Header name of taxonomy id (default: auto-detect)")
    p.add_argument("--lineage-col", help="Header name of lineage (default: auto-detect)")
    p.add_argument("--cache", default=CACHE_DEFAULT_PATH, help=f"Path to JSON cache (default: {CACHE_DEFAULT_PATH})")
    args = p.parse_args()

    sep = sniff_delimiter(args.input_file, user_sep=args.sep)

    with open(args.input_file, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter=sep)
        try:
            header = next(reader)
        except StopIteration:
            print("❌ Empty file.", file=sys.stderr)
            sys.exit(1)

        org_idx, org_name = find_col(header, args.organism_col, COMMON_ORG_COL_NAMES)
        taxid_idx, taxid_name = find_col(header, args.taxid_col, COMMON_TAXID_COL_NAMES)
        lineage_idx, lineage_name = find_col(header, args.lineage_col, COMMON_LINEAGE_COL_NAMES)

        if org_idx is None and taxid_idx is None:
            print("❌ Need at least an organism column or a taxonomy_id column.", file=sys.stderr)
            print("   Header columns detected:", header, file=sys.stderr)
            sys.exit(2)

        api_key = os.getenv("NCBI_API_KEY")
        cache_path = Path(args.cache)
        cache = load_cache(cache_path)
        session = requests.Session() if requests else None

        seen_taxa = {}
        kingdom_counts = defaultdict(set)
        unknown_organisms = []

        # Read all rows into memory to compute total unique progress denominator nicely
        rows = list(reader)
        total = len(rows)
        for i, row in enumerate(rows, 1):
            tax_id, kingdom = get_taxid_and_kingdom_from_row(
                row, org_idx, taxid_idx, lineage_idx, api_key, cache, session
            )

            # Progress print (organism if available, else taxid)
            label = ""
            if org_idx is not None and org_idx < len(row) and row[org_idx].strip():
                label = row[org_idx].strip()
            elif tax_id:
                label = f"TaxID:{tax_id}"
            else:
                label = "<no-organism>"

            print(f"[{i}/{total}] Processed: {label} → {kingdom} (TaxID: {tax_id})")

            if tax_id and kingdom in {"bacteria", "fungi", "virus", "other", "Unknown"}:
                if tax_id not in seen_taxa:
                    seen_taxa[tax_id] = kingdom
                    kingdom_counts[kingdom].add(tax_id)
            elif kingdom == "Unknown":
                if org_idx is not None and org_idx < len(row):
                    org = row[org_idx].strip()
                    if org:
                        unknown_organisms.append(org)

            # be polite if we are likely hitting NCBI (no API key and needed fallback)
            # when we had to use the web (detected via None taxid + org present and no lineage)
            if not api_key and ((taxid_idx is None or not (row[taxid_idx].strip() if taxid_idx < len(row) else "")) and not (row[lineage_idx].strip() if lineage_idx is not None and lineage_idx < len(row) else "")):
                time.sleep(0.34)  # ~3 req/s

        save_cache(cache_path, cache)

    # Summary / Output
    if unknown_organisms:
        print("\nOrganisms that could not be classified:")
        for org in sorted(set(unknown_organisms)):
            print(f"- {org}")

    print("\nUnique taxa counts by kingdom:")
    print(f"fungi\t{len(kingdom_counts['fungi'])}")
    print(f"bacteria\t{len(kingdom_counts['bacteria'])}")
    print(f"virus\t{len(kingdom_counts['virus'])}")
    print(f"other\t{len(kingdom_counts['other'])}")
    print(f"unknown\t{len(kingdom_counts['Unknown'])}")


if __name__ == "__main__":
    main()
