#!/usr/bin/env python3
# Christie Woodside
"""
Counts the number of BioSamples associated with each kingdom for the FDA-ARGOS BioProject.
Input: biosampleMeta_ARGOS_extended.tsv (or CSV) with columns like:
  - organism_name
  - lineage
  - taxonomy_id
  - biosample   (e.g., SAMNxxxxxx)

We dedupe by BioSample accession, classify each unique BioSample into
virus/fungi/bacteria/Other/Unknown, and print totals matching the original script.
"""

import argparse
import csv
import json
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

try:
    import requests
except ImportError:
    requests = None  # Only needed if we must fallback to NCBI

# ---------- Column name heuristics ----------

COMMON_ORG_COL_NAMES = {
    "organism", "organism_name", "scientific_name", "species",
    "taxonomy_name", "tax_name"
}
COMMON_TAXID_COL_NAMES = {"taxid", "tax_id", "taxonomy_id", "ncbi_tax_id"}
COMMON_LINEAGE_COL_NAMES = {"lineage", "ncbi_lineage"}
COMMON_BIOSAMPLE_COL_NAMES = {
    "biosample", "biosample_acc", "biosample_accession", "sample_accession",
    "sample", "biosample_id", "biosampleid", "biosample accession", "biosample id",
    "BioSample", "BioSample Accession"
}

CACHE_DEFAULT_PATH = ".ncbi_tax_cache.json"

# ---------- Utilities ----------

def sniff_delimiter(file_path, user_sep=None):
    if user_sep:
        return "\t" if user_sep == "\\t" else user_sep
    with open(file_path, "r", newline="") as fh:
        sample = fh.read(4096)
    tabs = sample.count("\t")
    commas = sample.count(",")
    return "\t" if tabs >= commas else ","

def find_col(header, user_col, common_names):
    # exact case-insensitive match first
    if user_col:
        for name in header:
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
        return "Unknown"
    if "Viruses" in lineage_text:
        return "virus"
    if "Fungi" in lineage_text:
        return "fungi"
    if "Bacteria" in lineage_text:
        return "bacteria"
    return "Other"

def normalize_taxid(raw_taxid: str):
    if not raw_taxid:
        return None
    m = re.search(r"\d+", str(raw_taxid))
    return m.group(0) if m else None

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
    ids = resp.json().get("esearchresult", {}).get("idlist") or []
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

# ---------- Core classification per row ----------

def determine_kingdom(row, cols, api_key, cache, session, polite=False):
    """
    Determine kingdom for a BioSample row using lineage/taxid/organism fallbacks.
    Returns one of: 'fungi', 'bacteria', 'virus', 'Other', 'Unknown'
    """
    org = row.get(cols["organism"], "").strip() if cols["organism"] else ""
    lineage = row.get(cols["lineage"], "").strip() if cols["lineage"] else ""
    raw_taxid = row.get(cols["taxid"], "").strip() if cols["taxid"] else ""
    taxid = normalize_taxid(raw_taxid)

    # 1) lineage present → classify directly
    if lineage:
        return classify_from_lineage(lineage)

    # 2) taxid present → fetch lineage (cached)
    if taxid:
        ck = f"taxid:{taxid}"
        if ck in cache:
            return cache[ck].get("kingdom", "Unknown")
        fetched = ncbi_fetch_lineage(taxid, api_key=api_key, session=session)
        k = classify_from_lineage(fetched)
        cache[ck] = {"kingdom": k}
        return k

    # 3) organism present → resolve via NCBI (cached)
    if org:
        ck = f"org:{org.lower()}"
        if ck in cache:
            return cache[ck].get("kingdom", "Unknown")
        t = ncbi_esearch_taxid(org, api_key=api_key, session=session)
        if not t:
            cache[ck] = {"kingdom": "Unknown"}
            return "Unknown"
        fetched = ncbi_fetch_lineage(t, api_key=api_key, session=session)
        k = classify_from_lineage(fetched)
        cache[ck] = {"kingdom": k}
        if polite and not api_key:
            time.sleep(0.34)  # ~3 req/s when unauthenticated
        return k

    return "Unknown"

# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(description="Count BioSamples by Kingdom from a TSV/CSV file")
    parser.add_argument("input_file", help="Path to biosampleMeta_ARGOS_extended.tsv (or CSV) with header")
    parser.add_argument("--sep", choices=[",", "\\t"], help="Field delimiter (default: auto-detect)")
    parser.add_argument("--organism-col", help="Header name for organism (e.g., organism_name)")
    parser.add_argument("--taxid-col", help="Header name for taxonomy id (e.g., taxonomy_id)")
    parser.add_argument("--lineage-col", help="Header name for lineage (e.g., lineage)")
    parser.add_argument("--biosample-col", help="Header name for biosample accession (e.g., biosample)")
    parser.add_argument("--cache", default=CACHE_DEFAULT_PATH, help=f"Path to JSON cache (default: {CACHE_DEFAULT_PATH})")
    args = parser.parse_args()

    sep = sniff_delimiter(args.input_file, user_sep=args.sep)

    with open(args.input_file, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=sep)
        if not reader.fieldnames:
            print("❌ Empty file or missing header.", file=sys.stderr)
            sys.exit(1)

        header = reader.fieldnames
        org_col = find_col(header, args.organism_col, COMMON_ORG_COL_NAMES)
        taxid_col = find_col(header, args.taxid_col, COMMON_TAXID_COL_NAMES)
        lineage_col = find_col(header, args.lineage_col, COMMON_LINEAGE_COL_NAMES)
        biosample_col = find_col(header, args.biosample_col, COMMON_BIOSAMPLE_COL_NAMES)

        if biosample_col is None:
            print("❌ Could not find BioSample column. Supply one with --biosample-col", file=sys.stderr)
            print("   Header columns detected:", header, file=sys.stderr)
            sys.exit(2)

        api_key = os.getenv("NCBI_API_KEY")
        cache_path = Path(args.cache)
        cache = load_cache(cache_path)
        session = requests.Session() if requests else None

        # Read all rows and dedupe by BioSample accession
        rows = list(reader)
        seen_biosamples = {}
        for row in rows:
            bs = (row.get(biosample_col, "") or "").strip()
            if not bs:
                continue
            if bs not in seen_biosamples:
                seen_biosamples[bs] = row

        # Tally
        kingdom_counts = Counter()
        unknown_organisms = []

        cols = {"organism": org_col, "taxid": taxid_col, "lineage": lineage_col, "biosample": biosample_col}

        total_unique = len(seen_biosamples)
        for idx, (bs, row) in enumerate(seen_biosamples.items(), 1):
            k = determine_kingdom(row, cols, api_key, cache, session, polite=True)
            kingdom_counts[k] += 1
            if k == "Unknown":
                # For continuity with your original script which listed unknown organisms,
                # record the organism name if available; otherwise the BioSample ID.
                org = (row.get(org_col, "") or "").strip() if org_col else ""
                unknown_organisms.append(org if org else bs)

            # Optional progress:
            # if idx % 100 == 0 or idx == total_unique:
            #     print(f"[{idx}/{total_unique}] {bs} → {k}")

        save_cache(cache_path, cache)

    # Output (same shape/labels as your original)
    if unknown_organisms:
        print("\nOrganisms that could not be classified:")
        for org in sorted(set(unknown_organisms)):
            print(f"- {org}")

    print("\nBiosample counts by kingdom:")
    print(f"fungi\t{kingdom_counts.get('fungi', 0)}")
    print(f"bacteria\t{kingdom_counts.get('bacteria', 0)}")
    print(f"virus\t{kingdom_counts.get('virus', 0)}")
    print(f"other\t{kingdom_counts.get('Other', 0)}")
    print(f"unknown\t{kingdom_counts.get('Unknown', 0)}")
    print(f"total: \t{sum(kingdom_counts.values())}")
    print(f"expected total organisms: {total_unique}")

if __name__ == "__main__":
    main()
