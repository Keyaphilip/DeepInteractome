"""
Download ClinVar VCF data from NCBI FTP.

Usage:
    python src/data/download_data.py               # full download
    python src/data/download_data.py --check-only  # verify URL only
"""
import os
import sys
import gzip
import shutil
import argparse
import requests
from tqdm import tqdm


CLINVAR_VCF_GZ_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
)
CLINVAR_TBI_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"
)


def download_file(url: str, target_path: str, chunk_size: int = 1024 * 64) -> None:
    """Stream-download a file with a progress bar (memory-efficient)."""
    response = requests.get(url, stream=True, timeout=60)
    response.raise_for_status()
    total = int(response.headers.get("content-length", 0))

    os.makedirs(os.path.dirname(target_path), exist_ok=True)
    print(f"Downloading {url}\n  → {target_path}")
    with open(target_path, "wb") as f, tqdm(
        total=total, unit="iB", unit_scale=True, unit_divisor=1024
    ) as bar:
        for chunk in response.iter_content(chunk_size):
            f.write(chunk)
            bar.update(len(chunk))


def decompress_gz(gz_path: str) -> str:
    """Decompress a .gz file without loading into RAM. Returns the output path."""
    out_path = gz_path[:-3]  # strip .gz
    if os.path.exists(out_path):
        print(f"Already decompressed: {out_path}")
        return out_path
    print(f"Decompressing {gz_path} …")
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out, length=1024 * 1024)  # 1 MB chunks
    print(f"Done → {out_path}")
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Download ClinVar VCF (GRCh38).")
    parser.add_argument(
        "--raw-dir",
        default=os.path.join("data", "raw"),
        help="Target directory for raw downloads (default: data/raw)",
    )
    parser.add_argument(
        "--skip-decompress",
        action="store_true",
        help="Keep only the .gz file; skip decompression.",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Only check that the URL is reachable (HEAD request), then exit.",
    )
    args = parser.parse_args()

    gz_path = os.path.join(args.raw_dir, "clinvar.vcf.gz")
    tbi_path = gz_path + ".tbi"

    if args.check_only:
        r = requests.head(CLINVAR_VCF_GZ_URL, timeout=10)
        size_mb = int(r.headers.get("content-length", 0)) / 1e6
        print(f"URL reachable ✓  |  size ≈ {size_mb:.1f} MB  |  status {r.status_code}")
        sys.exit(0)

    # ── Download compressed VCF ───────────────────────────────────────────────
    if not os.path.exists(gz_path):
        download_file(CLINVAR_VCF_GZ_URL, gz_path)
    else:
        print(f"Already downloaded: {gz_path}")

    # ── Download index (.tbi) — useful for tabix-enabled tools ───────────────
    if not os.path.exists(tbi_path):
        try:
            download_file(CLINVAR_TBI_URL, tbi_path)
        except Exception as e:
            print(f"[warn] Could not download .tbi index: {e}")

    # ── Decompress ────────────────────────────────────────────────────────────
    if not args.skip_decompress:
        decompress_gz(gz_path)


if __name__ == "__main__":
    main()
