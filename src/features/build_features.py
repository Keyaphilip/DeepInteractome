"""
Feature engineering for DeepInteractome.

Reads a ClinVar VCF (plain-text or gzip) and produces a CSV of numeric
features suitable for training or inference.

Usage:
    python src/features/build_features.py                         # sample VCF
    python src/features/build_features.py --vcf data/raw/clinvar.vcf.gz  --limit 50000
    python src/features/build_features.py --vcf data/raw/clinvar.vcf.gz  # full run
"""
import os
import sys
import gzip
import argparse
import random

import pandas as pd
import numpy as np

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
from src.utils.vcf_parser import SimpleVCFParser


# ── Nucleotide one-hot encoding ───────────────────────────────────────────────

_BASE_MAP = {
    "A": [1, 0, 0, 0],
    "C": [0, 1, 0, 0],
    "G": [0, 0, 1, 0],
    "T": [0, 0, 0, 1],
    "N": [0, 0, 0, 0],
}

_CLNVC_MAP = {
    "single_nucleotide_variant": 0,
    "Deletion": 1,
    "Insertion": 2,
    "Duplication": 3,
    "Indel": 4,
    "Inversion": 5,
    "Microsatellite": 6,
}

_ORIGIN_MAP = {
    "0": 0, "1": 1, "2": 2, "4": 4, "8": 8, "16": 16, "32": 32,
}


def one_hot_base(base: str) -> list:
    return _BASE_MAP.get(base.upper(), [0, 0, 0, 0])


def one_hot_sequence(sequence: str, max_len: int = 1) -> list:
    """One-hot encode up to max_len bases, padding with zeros."""
    enc = []
    for b in (sequence.upper() + "N" * max_len)[:max_len]:
        enc.extend(_BASE_MAP.get(b, [0, 0, 0, 0]))
    return enc


def random_flanking(size: int = 5) -> tuple:
    """Placeholder flanking context (random). Overridden when a genome FASTA is available."""
    bases = "ACGT"
    up = "".join(random.choices(bases, k=size))
    dn = "".join(random.choices(bases, k=size))
    return up, dn


# ── Consequence (MC field) helpers ────────────────────────────────────────────

_HIGH_IMPACT = {
    "frameshift_variant", "stop_gained", "stop_lost", "start_lost",
    "splice_acceptor_variant", "splice_donor_variant",
}
_MOD_IMPACT = {"missense_variant", "inframe_deletion", "inframe_insertion"}


def parse_mc(mc_field: str) -> dict:
    """
    Parse the ClinVar MC (molecular consequence) INFO field.
    Returns dict with flags: is_frameshift, is_nonsense, is_missense, is_synonymous.
    Example MC value: 'SO:0001583|missense_variant'
    """
    terms = set()
    for part in mc_field.split(","):
        if "|" in part:
            terms.add(part.split("|")[1].strip())
        else:
            terms.add(part.strip())
    return {
        "is_frameshift": int(bool(terms & {"frameshift_variant"})),
        "is_nonsense": int(bool(terms & {"stop_gained", "stop_lost", "start_lost"})),
        "is_missense": int(bool(terms & {"missense_variant"})),
        "is_synonymous": int(bool(terms & {"synonymous_variant"})),
        "is_high_impact": int(bool(terms & _HIGH_IMPACT)),
        "is_moderate_impact": int(bool(terms & _MOD_IMPACT)),
    }


# ── Core feature extraction ───────────────────────────────────────────────────

class GenomicFeatureEngineer:
    """Converts a ClinVar VCF into a labelled feature DataFrame."""

    def __init__(self, vcf_path: str):
        self.vcf_path = vcf_path
        self.parser = SimpleVCFParser(vcf_path)

    # kept for backwards-compat with existing tests
    def _one_hot_encode_dna(self, sequence: str, max_len: int = 1) -> list:
        return one_hot_sequence(sequence, max_len)

    def _get_flanking_sequence(self, chrom, pos, ref, size: int = 5):
        return random_flanking(size)

    def _extract_row(self, v: dict, flank_size: int = 5) -> dict:
        info = v["INFO"]
        ref = v["REF"]
        alt = v["ALT"][0] if v["ALT"] else "N"

        # ── Label ─────────────────────────────────────────────────────────────
        clnsig = info.get("CLNSIG", "Unknown")
        is_pathogenic = int("Pathogenic" in clnsig or "Likely_pathogenic" in clnsig)

        # ── Allele frequency ──────────────────────────────────────────────────
        try:
            af = float(info.get("AF_ESP", info.get("AF", 0.0)))
        except (ValueError, TypeError):
            af = 0.0

        # ── Allele length features ────────────────────────────────────────────
        ref_len = len(ref)
        alt_len = len(alt)
        variant_size = abs(alt_len - ref_len)
        is_snv = int(ref_len == 1 and alt_len == 1)

        # ── CLNVC (variant class) ─────────────────────────────────────────────
        clnvc = info.get("CLNVC", "unknown")
        clnvc_encoded = _CLNVC_MAP.get(clnvc, -1)

        # ── Origin (germline / somatic …) ─────────────────────────────────────
        origin_raw = info.get("ORIGIN", "0")
        try:
            origin = int(origin_raw.split(",")[0])
        except Exception:
            origin = 0

        # ── dbSNP RS id ───────────────────────────────────────────────────────
        has_rs = int(info.get("RS", "") != "")

        # ── Molecular consequence ─────────────────────────────────────────────
        mc_features = parse_mc(info.get("MC", ""))

        # ── HGVS length (proxy for annotation complexity) ─────────────────────
        hgvs = info.get("CLNHGVS", "")
        hgvs_len = len(hgvs)

        # ── One-hot base encoding ─────────────────────────────────────────────
        ref_oh = one_hot_base(ref[0] if ref else "N")
        alt_oh = one_hot_base(alt[0] if alt else "N")

        # ── Flanking context ──────────────────────────────────────────────────
        upstream, downstream = self._get_flanking_sequence(
            v["CHROM"], v["POS"], ref, size=flank_size
        )
        up_oh = one_hot_sequence(upstream, flank_size)
        dn_oh = one_hot_sequence(downstream, flank_size)

        row = {
            "CHROM": v["CHROM"],
            "POS": v["POS"],
            # Length features
            "REF_len": ref_len,
            "ALT_len": alt_len,
            "variant_size": variant_size,
            "is_snv": is_snv,
            # Population frequency
            "AF": af,
            # Variant class
            "CLNVC_code": clnvc_encoded,
            # Germline/somatic origin
            "ORIGIN": origin,
            # dbSNP
            "has_rs": has_rs,
            # HGVS annotation length
            "hgvs_len": hgvs_len,
            # Base one-hot
            "Ref_A": ref_oh[0], "Ref_C": ref_oh[1], "Ref_G": ref_oh[2], "Ref_T": ref_oh[3],
            "Alt_A": alt_oh[0], "Alt_C": alt_oh[1], "Alt_G": alt_oh[2], "Alt_T": alt_oh[3],
            # Label
            "Target_Label": is_pathogenic,
            "ClinSig_Raw": clnsig,
        }

        # Molecular consequence flags
        row.update(mc_features)

        # Flanking sequence features
        for i in range(flank_size):
            bi = i * 4
            row[f"Up_{i}_A"] = up_oh[bi];   row[f"Up_{i}_C"] = up_oh[bi + 1]
            row[f"Up_{i}_G"] = up_oh[bi + 2]; row[f"Up_{i}_T"] = up_oh[bi + 3]
            row[f"Down_{i}_A"] = dn_oh[bi]; row[f"Down_{i}_C"] = dn_oh[bi + 1]
            row[f"Down_{i}_G"] = dn_oh[bi + 2]; row[f"Down_{i}_T"] = dn_oh[bi + 3]

        return row

    def process_data(
        self,
        output_path: str = None,
        limit: int = None,
        protein_df: pd.DataFrame = None,
    ) -> pd.DataFrame:
        """
        Parse VCF, extract features, optionally merge protein structural data.

        Args:
            output_path: If given, save CSV here.
            limit: Cap the number of variants processed (useful for fast iteration).
            protein_df: Optional DataFrame with protein structural features keyed
                        by 'gene' column (from protein_features.py).
        """
        print(f"Processing {self.vcf_path} {'(limit=' + str(limit) + ')' if limit else '(full)'}…")
        self.parser.parse_header()

        data = []
        for i, v in enumerate(self.parser.parse_variants()):
            if limit and i >= limit:
                break
            try:
                row = self._extract_row(v)
                data.append(row)
            except Exception as exc:
                print(f"[warn] skipping variant at line ~{i}: {exc}")

        df = pd.DataFrame(data)

        # ── Merge protein structural features (Phase 2) ───────────────────────
        if protein_df is not None and "gene" in df.columns:
            df = df.merge(protein_df, on="gene", how="left")
            # Fill NA for variants without a protein mapping
            struct_cols = [c for c in protein_df.columns if c != "gene"]
            df[struct_cols] = df[struct_cols].fillna(0.0)

        if output_path:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            df.to_csv(output_path, index=False)
            print(f"Saved → {output_path}  ({len(df):,} rows, {df.shape[1]} features)")

        return df


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Build genomic feature CSV from a ClinVar VCF.")
    parser.add_argument("--vcf", default=os.path.join("data", "raw", "sample.vcf"),
                        help="Input VCF path (plain or .gz)")
    parser.add_argument("--out", default=os.path.join("data", "processed", "genomic_features.csv"),
                        help="Output CSV path")
    parser.add_argument("--limit", type=int, default=None,
                        help="Maximum number of variants to process")
    args = parser.parse_args()

    if not os.path.exists(args.vcf):
        print(f"VCF not found: {args.vcf}")
        sys.exit(1)

    engineer = GenomicFeatureEngineer(args.vcf)
    df = engineer.process_data(output_path=args.out, limit=args.limit)
    print(f"\nDataset summary:")
    print(f"  Rows      : {len(df):,}")
    print(f"  Features  : {df.shape[1]}")
    print(f"  Pathogenic: {df['Target_Label'].sum():,} ({df['Target_Label'].mean()*100:.1f}%)")
    print(f"  Benign    : {(df['Target_Label'] == 0).sum():,}")


if __name__ == "__main__":
    main()
