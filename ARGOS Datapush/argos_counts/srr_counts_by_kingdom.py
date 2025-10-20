#!/usr/bin/env python3
# Christie Woodside
"""
Counts the number of SRR IDs in each kingdom from an ARGOS table (TSV/CSV).

Prefers in-file lineage/taxonomy_id to avoid web calls:
- If `lineage` present: classify via lineage (Viruses/Fungi/Bacteria -> virus/fungi/bacteria; else other)
- Else if `taxonomy_id` present: fetch lineage once (cached) and classify
- Else: fallback to resolving from organism name via NCBI (cached)

Only counts a row if it has a non-empty SRR value.

Output format matches the original script:
  - fungi, bacteria, virus counts
  - Summary block with totals and failure/skip tallies
"""

import csv
import argparse
import os
import time
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

try:
    import requests
except ImportError:
    requests = None  # Only needed if we must fallback to NCBI

# ---------- Config / heuristics ----------

COMMON_ORG_COL_NAMES = {
    "organism", "organism_name", "scientific_name", "species",
    "taxonomy_name", "tax_name"
}
COMMON_TAXID_COL_NAMES = {"taxid", "tax_id", "taxonomy_id", "ncbi_tax_id"}
COMMON_LINEAGE_COL_NAMES = {"lineage", "ncbi_lineage"}
COMMON_SRR_COL_NAMES = {"sra_run_id", "sra", "srr", "run", "run_accession"}

CACHE_DEFAULT_PATH = ".ncbi_tax_cache.json"

# ---------- Utilities ----------

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
                return name
    # common names
    lowered = [h.strip().lower() for h in header]
    for i, name in enumerate(lowered):
        if name in common_names:
            return header[i]
    # heuristic substring
    for i, name in enumerate(lowered):
        if any(k in name for k in common_names):
            return header[i]
    return None

def classify_from_lineage(lineage_text: str):
    if not lineage_text:
        return None
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

def normalize_taxid(raw_taxid: str):
    if not raw_taxid:
        return None
    m = re.search(r"\d+", raw_taxid)
    return m.group(0) if m else None

# ---------- Main logic ----------

def get_kingdom_from_row(row, cols, api_key, cache, session, be_polite=False):
    """
    Determine kingdom for a row using lineage/taxid/organism fallbacks.
    Returns "fungi"/"bacteria"/"virus"/"other" or None on failure.
    """
    org = row.get(cols["organism"], "").strip() if cols["organism"] else ""
    raw_taxid = row.get(cols["taxid"], "").strip() if cols["taxid"] else ""
    lineage = row.get(cols["lineage"], "").strip() if cols["lineage"] else ""
    taxid = normalize_taxid(raw_taxid)

    # 1) lineage present → classify directly
    if lineage:
        return classify_from_lineage(lineage)

    # 2) taxid present → fetch lineage (cached)
    if taxid:
        cache_key = f"taxid:{taxid}"
        if cache_key in cache:
            return cache[cache_key].get("kingdom")
        fetched_lineage = ncbi_fetch_lineage(taxid, api_key=api_key, session=session)
        k = classify_from_lineage(fetched_lineage) or "other"
        cache[cache_key] = {"kingdom": k}
        return k

    # 3) organism present → resolve via NCBI (cached)
    if org:
        cache_key = f"org:{org.lower()}"
        if cache_key in cache:
            return cache[cache_key].get("kingdom")
        tax_from_name = ncbi_esearch_taxid(org, api_key=api_key, session=session)
        if not tax_from_name:
            cache[cache_key] = {"kingdom": None}
            return None
        lineage2 = ncbi_fetch_lineage(tax_from_name, api_key=api_key, session=session)
        k = classify_from_lineage(lineage2)
        cache[cache_key] = {"kingdom": k}
        if be_polite and not api_key:
            time.sleep(0.34)  # ~3 req/s if unauthenticated
        return k

    return None

def main():
    parser = argparse.ArgumentParser(description="Count SRR entries by kingdom from a TSV/CSV file")
    parser.add_argument("input_file", help="Path to the TSV/CSV file with header")
    parser.add_argument("--sep", choices=[",", "\\t"], help="Field delimiter (default: auto-detect)")
    parser.add_argument("--organism-col", help="Header name for organism (e.g., organism_name)")
    parser.add_argument("--taxid-col", help="Header name for taxonomy id (e.g., taxonomy_id)")
    parser.add_argument("--lineage-col", help="Header name for lineage (e.g., lineage)")
    parser.add_argument("--srr-col", help="Header name for SRR/run id (e.g., sra_run_id)")
    parser.add_argument("--cache", default=CACHE_DEFAULT_PATH, help=f"Path to JSON cache (default: {CACHE_DEFAULT_PATH})")
    args = parser.parse_args()

    sep = sniff_delimiter(args.input_file, user_sep=args.sep)

    with open(args.input_file, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=sep)
        if reader.fieldnames is None:
            print("❌ Empty file or missing header.", file=sys.stderr)
            sys.exit(1)

        header = reader.fieldnames

        org_col = find_col(header, args.organism_col, COMMON_ORG_COL_NAMES)
        taxid_col = find_col(header, args.taxid_col, COMMON_TAXID_COL_NAMES)
        lineage_col = find_col(header, args.lineage_col, COMMON_LINEAGE_COL_NAMES)
        srr_col = find_col(header, args.srr_col, COMMON_SRR_COL_NAMES)

        if srr_col is None:
            print("❌ Could not find SRR column. Supply one with --srr-col", file=sys.stderr)
            print("   Header columns detected:", header, file=sys.stderr)
            sys.exit(2)

        api_key = os.getenv("NCBI_API_KEY")
        cache_path = Path(args.cache)
        cache = load_cache(cache_path)
        session = requests.Session() if requests else None

        # tallies (match original output keys)
        kingdom_counts = defaultdict(int)
        organism_to_kingdom = {}  # cache at organism-level to reduce fallbacks
        skipped_blank = 0
        failed_taxonomy = 0
        counted_total = 0
        other_kingdom = 0

        entries = list(reader)
        total_rows = len(entries)

        cols = {"organism": org_col, "taxid": taxid_col, "lineage": lineage_col, "srr": srr_col}

        for i, row in enumerate(entries, 1):
            srr = (row.get(srr_col, "") or "").strip()
            org = (row.get(org_col, "") or "").strip() if org_col else ""

            # Require an SRR to count this row
            if not srr:
                # silently skip rows without SRR since the goal is SRR counts
                continue

            if not org and not (taxid_col or lineage_col):
                # original script tracked "blank organism" — treat rows with no organism
                # (and no tax/lineage to classify) as skipped_blank
                skipped_blank += 1
                continue

            # Try to reuse a cached kingdom by organism label when available
            kingdom = None
            if org:
                kingdom = organism_to_kingdom.get(org)

            if kingdom is None:
                try:
                    kingdom = get_kingdom_from_row(
                        row, cols, api_key, cache, session, be_polite=True
                    )
                except Exception as e:
                    kingdom = None

                if kingdom is None:
                    failed_taxonomy += 1
                    continue

                if org:
                    organism_to_kingdom[org] = kingdom

            if kingdom in {"fungi", "bacteria", "virus"}:
                kingdom_counts[kingdom] += 1
                counted_total += 1
            else:
                # mirror original behavior for untracked
                other_kingdom += 1

            # Optional progress (commented to keep output identical)
            # print(f"[{i}/{total_rows}] {org or '—'} (SRR:{srr}) → {kingdom}")

        save_cache(cache_path, cache)

    # Final outputs — keep the same shape/labels as your original script
    print("\nTotal SRR counts by kingdom:")
    print(f"fungi\t{kingdom_counts['fungi']}")
    print(f"bacteria\t{kingdom_counts['bacteria']}")
    print(f"virus\t{kingdom_counts['virus']}")

    print(f"\nSummary:")
    print(f"Total rows in CSV: {total_rows}")
    print(f"Counted entries: {counted_total}")
    print(f"Skipped (blank organism): {skipped_blank}")
    print(f"Failed taxonomy lookup: {failed_taxonomy}")
    print(f"Other/untracked kingdom: {other_kingdom}")

if __name__ == "__main__":
    main()
