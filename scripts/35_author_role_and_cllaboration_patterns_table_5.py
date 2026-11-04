#!/usr/bin/env python3
"""Independent reconstruction of author_ABC_profile.csv from author_program_publication_ABC.csv.

This is NOT the original NIH-MAP upstream script. The publication-level input is
already disambiguated and annotated with program, PI status, and A/B/C/Other.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

KEY = ["canonical_author_id", "atlas_program"]
COLS = [
    "canonical_author_id", "canonical_author_name", "atlas_program", "PI_status",
    "multi_atlas_PI_status", "n_all_program_publications", "n_IC_publications",
    "n_A", "n_B", "n_C", "participates_A", "participates_B",
    "participates_C", "C_share", "repeated_C", "pattern_breadth",
    "multi_pattern", "dominant_pattern", "pattern_entropy",
    "effective_number_of_patterns", "activity_category",
]


def reconstruct(path):
    df = pd.read_csv(path, dtype={"pmid": str})
    required = KEY + ["pmid", "canonical_author_name", "PI_status", "multi_atlas_PI_status", "publication_pattern"]
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    if not set(df.publication_pattern.unique()).issubset({"A", "B", "C", "Other"}):
        raise ValueError("Unexpected collaboration patterns")
    if df[KEY + ["pmid"]].duplicated().any():
        raise ValueError("Duplicate author/program/PMID records detected")
    base = df[KEY + ["canonical_author_name", "PI_status", "multi_atlas_PI_status"]].drop_duplicates()
    if base.duplicated(KEY).any():
        raise ValueError("Person metadata not constant within author/program")
    counts = (df.groupby(KEY).pmid.nunique()
              .rename("n_all_program_publications").reset_index())
    tab = df.groupby(KEY + ["publication_pattern"]).pmid.nunique().unstack("publication_pattern", fill_value=0)
    tab = tab.reindex(columns=["A", "B", "C"], fill_value=0).rename(columns={"A":"n_A", "B":"n_B", "C":"n_C"}).reset_index()
    prof = base.merge(counts, on=KEY, validate="one_to_one").merge(tab,on=KEY,validate="one_to_one")
    prof["n_IC_publications"] = prof[["n_A", "n_B", "n_C"]].sum(axis=1)
    for p in "ABC": prof[f"participates_{p}"] = (prof[f"n_{p}"] > 0).astype(int)
    prof["C_share"] = np.where(prof.n_IC_publications > 0, prof.n_C / prof.n_IC_publications, np.nan)
    prof["repeated_C"] = (prof.n_C >= 2).astype(int)
    prof["pattern_breadth"] = prof[[f"participates_{p}" for p in "ABC"]].sum(axis=1)
    prof["multi_pattern"] = (prof.pattern_breadth >= 2).astype(int)
    arr = prof[["n_A", "n_B", "n_C"]].to_numpy()
    top = arr.max(axis=1)
    num_top = (arr == top[:,None]).sum(axis=1)
    prof["dominant_pattern"] = np.where(top==0, "none", np.where(num_top>1, "tie", np.asarray(["A","B","C"])[arr.argmax(axis=1)]))
    # Entropy is defined only for those with at least two internal-collaboration publications.
    total = prof.n_IC_publications.to_numpy(dtype=float)
    proportions = np.divide(arr, total[:, None], out=np.zeros_like(arr,dtype=float), where=total[:,None]>0)
    with np.errstate(divide="ignore",invalid="ignore"):
        logprob = np.where(proportions>0,np.log(proportions),0.)
    entropy = -(proportions*logprob).sum(axis=1)
    prof["pattern_entropy"] = np.where(total>=2,entropy,np.nan)
    prof["effective_number_of_patterns"] = np.exp(prof.pattern_entropy)
    n = prof.n_IC_publications
    prof["activity_category"] = np.select(
        [n>=10,n>=5,n>=2,n>=1],
        ["Highly Active","Very Active","Moderately Active","Slightly Active"],
        default="Inactive",
    )
    return prof[COLS].sort_values(KEY).reset_index(drop=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input",help="Path to author_program_publication_ABC.csv")
    ap.add_argument(
        "--output",
        default="outputs/tables/table_5_author_role_and_collaboration_pattern.csv",
    )
    ap.add_argument("--reference",help="Optional original author_ABC_profile.csv for validation")
    args=ap.parse_args()
    result=reconstruct(args.input)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output,index=False)
    print(f"Saved {len(result):,} author-program profiles to {args.output}")
    if args.reference:
        reference=pd.read_csv(args.reference).sort_values(KEY).reset_index(drop=True)
        pd.testing.assert_frame_equal(result,reference,check_dtype=False,check_exact=False,atol=1e-12,rtol=1e-12)
        print("PASS: all 21 fields match the reference author_ABC_profile.csv")

if __name__=="__main__": main()
