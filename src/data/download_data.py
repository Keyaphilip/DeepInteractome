import os
import requests
import gzip
import shutil
from tqdm import tqdm

def download_file(url, target_path):
    """Downloads a file from a URL with a progress bar."""
    response = requests.get(url, stream=True)
    total_size_in_bytes = int(response.headers.get('content-length', 0))
    block_size = 1024 # 1 Kibibyte
    
    print(f"Downloading {url} to {target_path}...")
    
    with open(target_path, 'wb') as file, tqdm(
        desc=target_path,
        total=total_size_in_bytes,
        unit='iB',
        unit_scale=True,
        unit_divisor=1024,
    ) as bar:
        for data in response.iter_content(block_size):
            size = file.write(data)
            bar.update(size)

def main():
    # URL for ClinVar VCF (GRCh38) - Using a specific small summary file or the main one
    # The main ClinVar VCF is ~50MB compressed, which is manageable
    vcf_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
    
    raw_dir = os.path.join("data", "raw")
    os.makedirs(raw_dir, exist_ok=True)
    
    vcf_filename = "clinvar.vcf.gz"
    vcf_path = os.path.join(raw_dir, vcf_filename)
    
    # Download
    if not os.path.exists(vcf_path):
        try:
            download_file(vcf_url, vcf_path)
            print("Download complete!")
        except Exception as e:
            print(f"Error downloading file: {e}")
            return
    else:
        print(f"File {vcf_path} already exists. Skipping download.")

    # Decompress for easier inspection (optional, but good for learning)
    extracted_path = vcf_path.replace(".gz", "")
    if not os.path.exists(extracted_path):
        print(f"Unzipping to {extracted_path}...")
        try:
            with gzip.open(vcf_path, 'rb') as f_in:
                with open(extracted_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print("Unzip complete!")
        except Exception as e:
            print(f"Error unzipping: {e}")

if __name__ == "__main__":
    main()
