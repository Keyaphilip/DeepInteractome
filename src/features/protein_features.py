"""
Fetch protein structural confidence features from the EBI AlphaFold API.

Per-residue pLDDT scores are summarized into per-protein statistics:
  - mean_pLDDT   : mean confidence score (0–100)
  - min_pLDDT    : minimum confidence score
  - pct_confident: % of residues with pLDDT ≥ 70 (structured)
  - pct_disordered: % of residues with pLDDT < 50 (likely disordered)

Usage:
    python src/features/protein_features.py --accession P04637
    python src/features/protein_features.py --accession-list data/interim/uniprot_map.csv

Results cached to data/interim/protein_features_cache.csv.
"""
import os
import sys
import csv
import time
import json
import argparse
import requests
from typing import Optional

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
CACHE_PATH = os.path.join("data", "interim", "protein_features_cache.csv")
CACHE_FIELDS = [
    "accession", "mean_pLDDT", "min_pLDDT", "pct_confident", "pct_disordered",
    "residue_count",
]


# ── Cache helpers ─────────────────────────────────────────────────────────────

def load_cache() -> dict:
    cache = {}
    if os.path.exists(CACHE_PATH):
        with open(CACHE_PATH, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                acc = row.pop("accession")
                cache[acc] = {k: float(v) for k, v in row.items()}
    return cache


def save_cache(cache: dict) -> None:
    os.makedirs(os.path.dirname(CACHE_PATH), exist_ok=True)
    with open(CACHE_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CACHE_FIELDS)
        writer.writeheader()
        for acc, stats in cache.items():
            writer.writerow({"accession": acc, **stats})


# ── API fetch ─────────────────────────────────────────────────────────────────

def fetch_plddt(accession: str, retries: int = 3) -> Optional[list]:
    """
    Fetch per-residue pLDDT scores for a UniProt accession from AlphaFold EBI.
    Returns a list of float scores, or None if the protein is not in AlphaFold.
    """
    url = ALPHAFOLD_API.format(accession=accession)
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=20)
            if r.status_code == 404:
                return None  # not in AlphaFold DB
            r.raise_for_status()
            entries = r.json()
            if not entries:
                return None
            # The pLDDT is stored in the cif/mmCIF file; the JSON has a pdbUrl.
            # We use the 'plddt' key available in the summary endpoint.
            entry = entries[0]
            plddt_scores = entry.get("plddt", [])
            if plddt_scores:
                return [float(s) for s in plddt_scores]
            return None
        except requests.RequestException as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
            else:
                print(f"[warn] AlphaFold fetch failed for {accession}: {e}")
                return None
    return None


def summarize_plddt(scores: list) -> dict:
    if not scores:
        return {
            "mean_pLDDT": 0.0, "min_pLDDT": 0.0,
            "pct_confident": 0.0, "pct_disordered": 0.0,
            "residue_count": 0,
        }
    n = len(scores)
    mean_p = sum(scores) / n
    min_p = min(scores)
    pct_conf = sum(1 for s in scores if s >= 70) / n * 100
    pct_dis = sum(1 for s in scores if s < 50) / n * 100
    return {
        "mean_pLDDT": round(mean_p, 2),
        "min_pLDDT": round(min_p, 2),
        "pct_confident": round(pct_conf, 2),
        "pct_disordered": round(pct_dis, 2),
        "residue_count": n,
    }


# ── Public API ────────────────────────────────────────────────────────────────

def get_protein_features(accessions: list, verbose: bool = True) -> dict:
    """
    Fetch structural features for a list of UniProt accessions.
    Returns dict {accession: {mean_pLDDT, min_pLDDT, pct_confident, pct_disordered}}.
    Caches results locally.
    """
    cache = load_cache()
    to_fetch = [a for a in accessions if a and a not in cache]

    if verbose and to_fetch:
        print(f"Fetching AlphaFold pLDDT for {len(to_fetch)} accession(s)…")

    for acc in to_fetch:
        scores = fetch_plddt(acc)
        stats = summarize_plddt(scores or [])
        cache[acc] = stats
        if verbose:
            if scores:
                print(f"  {acc}: mean_pLDDT={stats['mean_pLDDT']}, residues={stats['residue_count']}")
            else:
                print(f"  {acc}: not found in AlphaFold")
        time.sleep(0.3)

    save_cache(cache)
    return {a: cache[a] for a in accessions if a in cache}


def build_protein_df(gene_to_acc: dict) -> "pd.DataFrame":
    """
    Given a {gene: accession} dict, return a DataFrame with structural features
    keyed by gene. Suitable for merging with the variant feature DataFrame.
    """
    import pandas as pd

    accessions = [v for v in gene_to_acc.values() if v]
    acc_features = get_protein_features(accessions)

    rows = []
    for gene, acc in gene_to_acc.items():
        row = {"gene": gene}
        if acc and acc in acc_features:
            row.update(acc_features[acc])
        else:
            row.update({
                "mean_pLDDT": None, "min_pLDDT": None,
                "pct_confident": None, "pct_disordered": None,
                "residue_count": None,
            })
        rows.append(row)
    return pd.DataFrame(rows)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Fetch AlphaFold pLDDT structural features.",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--accession", help="Single UniProt accession (e.g. P04637)")
    group.add_argument(
        "--accession-list",
        help="CSV file with 'gene' and 'accession' columns (output of fetch_uniprot.py)",
    )
    args = parser.parse_args()

    if args.accession:
        accessions = [args.accession]
    else:
        import csv as _csv
        accessions = []
        with open(args.accession_list, encoding="utf-8") as f:
            for row in _csv.DictReader(f):
                if row.get("accession"):
                    accessions.append(row["accession"])

    results = get_protein_features(accessions)
    print(f"\nFetched structural features for {len(results)} protein(s).")
    for acc, stats in results.items():
        print(f"  {acc}: {stats}")


if __name__ == "__main__":
    main()
