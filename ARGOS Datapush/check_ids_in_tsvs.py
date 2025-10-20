#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path
from typing import Union, Optional
import re
import os


'''Takes in the file genome_assembly_id-breakdown.tsv which stephen provided. It then takes in the biosampleMeta tsv and ngsQC tsv to see if ids in assemblyQC are not found in the other two tsvs.'''
# to activate environment (because python won't really update): source ~/venvs/py313/bin/activate
# --- keep your existing imports, read_tsv, col_exists, to_id_set ---

'''command below'''
# python3 /Users/christiewoodside/Desktop/ARGOS/code/HIVE3/current/check_ids_in_tsvs.py \
  # --input <(awk 'NR==20{print; exit}' /Users/christiewoodside/Desktop/genome_assembly_id-breakdown.tsv | tr '\t' '\n' | sed '/^$/d') \
  # --no-header \
  # --id-col 0 \
  # --tsv1 /Users/christiewoodside/Desktop/ngsQC_ARGOS.tsv --tsv1-col genome_assembly_id \
  # --tsv2 /Users/christiewoodside/Desktop/biosampleMeta_ARGOS.tsv --tsv2-col genome_assembly_id \
  # --out /Users/christiewoodside/Desktop/missing_ids_unique_to_assembly.tsv


def read_tsv(path: str, header: bool = True) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"File not found: {path}")
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, header=(0 if header else None))

def col_exists(df: pd.DataFrame, col: Union[str, int]) -> bool:
    if isinstance(col, str):
        return col in df.columns
    elif isinstance(col, int):
        return 0 <= col < df.shape[1]
    else:
        return False

def to_id_set_from_df(df: pd.DataFrame, col: str) -> set:
    if not col_exists(df, col):
        raise KeyError(f"Column '{col}' not found. Columns: {list(df.columns)}")
    return set(x.strip() for x in df[col].astype(str) if x.strip() != "")

# GCF/GCA with optional .version (missing version => 1)
_ACC_RE = re.compile(r'^(?P<prefix>GC[AF])_(?P<digits>\d+)(?:\.(?P<ver>\d+))?$')

def parse_acc(acc: str) -> Optional[tuple[str, str, int]]:
    m = _ACC_RE.match(acc)
    if not m:
        return None
    prefix = m.group("prefix")
    digits = m.group("digits")
    ver_s = m.group("ver")
    ver = int(ver_s) if ver_s and ver_s.isdigit() else 1
    return prefix, digits, ver

def find_gca_for_one_table(gcf_id: str, target_set: set, max_increments: int) -> tuple[str, bool]:
    """
    For a given GCF id and ONE table's id set, try GCA variants in that table only.
    Returns (gca_candidate, fixed_bool).
    """
    parsed = parse_acc(gcf_id)
    if not parsed:
        return ("", False)
    prefix, digits, ver = parsed
    if prefix != "GCF":
        # per requirement, we do not try to 'fix' GCA inputs
        return ("", False)

    base = f"GCA_{digits}"
    # try same version, then increments up to max_increments
    for i in range(0, max_increments + 1):
        cand = f"{base}.{ver + i}"
        if cand in target_set:
            return (cand, True)
    return (f"{base}.{ver}", False)  # initial candidate for reference only

def default_updated_path(orig_path: str) -> str:
    p = Path(orig_path)
    suffixes = ''.join(p.suffixes)
    if suffixes.endswith(".tsv.gz"):
        return str(p.with_name(p.name.replace(".tsv.gz", "_updated.tsv.gz")))
    elif suffixes.endswith(".tsv"):
        return str(p.with_name(p.name.replace(".tsv", "_updated.tsv")))
    else:
        return str(p.with_name(p.name + "_updated.tsv"))

def main():
    ap = argparse.ArgumentParser(
        description=("Check IDs from an input TSV in two tables. If a table lacks the GCF, "
                     "search that table for a GCA variant (same digits+version, then increment). "
                     "Write per-table subset files containing only rows that were fixed (GCA→GCF). "
                     "Report file lists only unresolved NOs (any table still missing with no GCA match).")
    )
    ap.add_argument("--input", required=True)
    ap.add_argument("--no-header", action="store_true")
    ap.add_argument("--id-col", default="genome_assembly_id")
    ap.add_argument("--org-col", default=None)
    ap.add_argument("--tsv1", required=True, help="First TSV (e.g., ngsQC_ARGOS.tsv)")
    ap.add_argument("--tsv1-col", required=True)
    ap.add_argument("--tsv2", required=True, help="Second TSV (e.g., biosampleMeta_ARGOS.tsv)")
    ap.add_argument("--tsv2-col", required=True)
    ap.add_argument("--max-increments", type=int, default=10)
    ap.add_argument("--out", default=None, help="Report TSV path (NOs only)")
    ap.add_argument("--tsv1-updated-out", default=None)
    ap.add_argument("--tsv2-updated-out", default=None)
    args = ap.parse_args()

    # Parse input column types
    id_col = int(args.id_col) if args.no_header and str(args.id_col).isdigit() else args.id_col
    org_col = None
    if args.org_col is not None:
        org_col = int(args.org_col) if args.no_header and str(args.org_col).isdigit() else args.org_col

    # Load input & tables
    df_in  = read_tsv(args.input, header=not args.no_header)
    if not col_exists(df_in, id_col):
        raise KeyError(f"ID column '{id_col}' not in INPUT: {list(df_in.columns)}")
    if org_col is not None and not col_exists(df_in, org_col):
        raise KeyError(f"Organism column '{org_col}' not in INPUT")

    df1 = read_tsv(args.tsv1, header=True)
    df2 = read_tsv(args.tsv2, header=True)
    if not col_exists(df1, args.tsv1_col):
        raise KeyError(f"Column '{args.tsv1_col}' not in {args.tsv1}")
    if not col_exists(df2, args.tsv2_col):
        raise KeyError(f"Column '{args.tsv2_col}' not in {args.tsv2}")

    set1 = to_id_set_from_df(df1, args.tsv1_col)
    set2 = to_id_set_from_df(df2, args.tsv2_col)

    ids  = df_in[id_col].astype(str).str.strip()
    orgs = df_in[org_col].astype(str).str.strip() if org_col is not None else pd.Series([""] * len(df_in))

    # Outputs
    report_rows = []  # unresolved NOs
    gca_to_gcf_for_tsv1: dict[str, str] = {}
    gca_to_gcf_for_tsv2: dict[str, str] = {}

    for idv, org in zip(ids, orgs):
        if idv == "":
            continue

        parsed = parse_acc(idv)
        # If input is already GCA, we don't fix—just consider unresolved for both tables that lack it.
        if parsed and parsed[0] == "GCA":
            # unresolved NO (user wants GCA left alone)
            report_rows.append({
                "missing_id": idv,
                "organism": org,
                "gca_candidate": "",
                "gca_present": "no"
            })
            continue

        # Original is (likely) GCF — evaluate per table
        present1 = idv in set1
        present2 = idv in set2

        fixed1 = fixed2 = False
        cand1 = cand2 = ""

        # Table 1: if missing GCF, try to find a GCA variant ONLY in table1
        if not present1:
            cand1, fixed1 = find_gca_for_one_table(idv, set1, args.max_increments)
            if fixed1 and cand1:
                gca_to_gcf_for_tsv1[cand1] = idv

        # Table 2: if missing GCF, try to find a GCA variant ONLY in table2
        if not present2:
            cand2, fixed2 = find_gca_for_one_table(idv, set2, args.max_increments)
            if fixed2 and cand2:
                gca_to_gcf_for_tsv2[cand2] = idv

        # If either table remains missing and unfixed, log a NO for manual review
        unresolved = ((not present1 and not fixed1) or (not present2 and not fixed2))
        if unresolved:
            # 'gca_present' is 'yes' if we fixed at least one table, 'no' if none
            gca_present = "yes" if (fixed1 or fixed2) else "no"
            # For reference, show one candidate if any (prefer cand1 then cand2)
            gca_candidate = cand1 or cand2
            report_rows.append({
                "missing_id": idv,
                "organism": org,
                "gca_candidate": gca_candidate,
                "gca_present": gca_present
            })

    # ----- Report: ONLY NOs -----
    report_df = pd.DataFrame(report_rows, columns=["missing_id", "organism", "gca_candidate", "gca_present"])
    # Keep only rows where at least one table is still unresolved (that's what we collected)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        report_df.to_csv(args.out, sep="\t", index=False)
        print(f"Wrote report (unresolved NOs and partials): {args.out} ({len(report_df)} rows)")
    else:
        print(report_df.to_csv(sep="\t", index=False).rstrip("\n"))

    # ----- Updated subset for TSV1: ONLY rows that needed fixing in TSV1 -----
    if gca_to_gcf_for_tsv1:
        out1_path = args.tsv1_updated_out or default_updated_path(args.tsv1)
        subset1 = df1[df1[args.tsv1_col].isin(gca_to_gcf_for_tsv1.keys())].copy()
        if not subset1.empty:
            subset1[args.tsv1_col] = subset1[args.tsv1_col].map(lambda x: gca_to_gcf_for_tsv1.get(x, x))
            Path(out1_path).parent.mkdir(parents=True, exist_ok=True)
            subset1.to_csv(out1_path, sep="\t", index=False)
            print(f"Wrote TSV1 updated subset: {out1_path} ({len(subset1)} rows)")
        else:
            print("No matching rows found in TSV1 to update.")
    else:
        print("No fixes required for TSV1.")

    # ----- Updated subset for TSV2: ONLY rows that needed fixing in TSV2 -----
    if gca_to_gcf_for_tsv2:
        out2_path = args.tsv2_updated_out or default_updated_path(args.tsv2)
        subset2 = df2[df2[args.tsv2_col].isin(gca_to_gcf_for_tsv2.keys())].copy()
        if not subset2.empty:
            subset2[args.tsv2_col] = subset2[args.tsv2_col].map(lambda x: gca_to_gcf_for_tsv2.get(x, x))
            Path(out2_path).parent.mkdir(parents=True, exist_ok=True)
            subset2.to_csv(out2_path, sep="\t", index=False)
            print(f"Wrote TSV2 updated subset: {out2_path} ({len(subset2)} rows)")
        else:
            print("No matching rows found in TSV2 to update.")
    else:
        print("No fixes required for TSV2.")

if __name__ == "__main__":
    main()