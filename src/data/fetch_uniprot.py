"""
Fetch UniProt accession IDs from gene symbols using the UniProt REST API.

Usage:
    python src/data/fetch_uniprot.py --genes BRCA1 TP53 CFTR
    python src/data/fetch_uniprot.py --genes-file data/interim/genes.txt

Results are cached to data/interim/uniprot_map.csv.
"""
import os
import sys
import csv
import time
import argparse
import requests

CACHE_PATH = os.path.join("data", "interim", "uniprot_map.csv")
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


def load_cache() -> dict:
    """Load existing gene → accession mapping from cache."""
    mapping = {}
    if os.path.exists(CACHE_PATH):
        with open(CACHE_PATH, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                mapping[row["gene"]] = row["accession"]
    return mapping


def save_cache(mapping: dict) -> None:
    os.makedirs(os.path.dirname(CACHE_PATH), exist_ok=True)
    with open(CACHE_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "accession"])
        writer.writeheader()
        for gene, acc in mapping.items():
            writer.writerow({"gene": gene, "accession": acc})


def query_uniprot(gene: str, organism: str = "9606", retries: int = 3) -> str:
    """
    Query UniProt for the canonical human accession for a gene symbol.
    Returns the accession string or '' if not found.
    """
    params = {
        "query": f"gene_exact:{gene} AND organism_id:{organism} AND reviewed:true",
        "fields": "accession,gene_names",
        "format": "json",
        "size": 1,
    }
    for attempt in range(retries):
        try:
            r = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=15)
            r.raise_for_status()
            results = r.json().get("results", [])
            if results:
                return results[0]["primaryAccession"]
            return ""
        except requests.RequestException as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
            else:
                print(f"[warn] UniProt query failed for {gene}: {e}")
                return ""
    return ""


def resolve_genes(genes: list, verbose: bool = True) -> dict:
    """
    Resolve a list of gene symbols → UniProt accessions.
    Uses on-disk cache; only queries the API for new genes.
    Returns dict {gene: accession}.
    """
    mapping = load_cache()
    new_genes = [g for g in genes if g not in mapping]

    if verbose and new_genes:
        print(f"Resolving {len(new_genes)} gene(s) via UniProt API…")

    for gene in new_genes:
        acc = query_uniprot(gene)
        mapping[gene] = acc
        if verbose:
            status = acc if acc else "not found"
            print(f"  {gene} → {status}")
        time.sleep(0.2)  # polite rate-limiting

    save_cache(mapping)
    return mapping


def main():
    parser = argparse.ArgumentParser(description="Map gene symbols to UniProt accessions.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--genes", nargs="+", help="Gene symbols (space-separated)")
    group.add_argument("--genes-file", help="Plain-text file with one gene symbol per line")
    args = parser.parse_args()

    if args.genes_file:
        with open(args.genes_file, encoding="utf-8") as f:
            genes = [line.strip() for line in f if line.strip()]
    else:
        genes = args.genes

    mapping = resolve_genes(genes)
    found = sum(1 for v in mapping.values() if v)
    print(f"\nResolved {found}/{len(genes)} gene(s). Cache saved to {CACHE_PATH}")


if __name__ == "__main__":
    main()
