import gzip
import os

class SimpleVCFParser:
    """
    A simple VCF parser to understand the format from scratch.
    Does not depend on complex libraries like pysam.
    """
    def __init__(self, filepath):
        self.filepath = filepath
        self.header = []
        self.samples = []
        
    def _open_file(self):
        if self.filepath.endswith('.gz'):
            return gzip.open(self.filepath, 'rt', encoding='utf-8')
        else:
            return open(self.filepath, 'r', encoding='utf-8')

    def parse_header(self):
        """Metadata lines and header line."""
        with self._open_file() as f:
            for line in f:
                if line.startswith('##'):
                    self.header.append(line.strip())
                elif line.startswith('#CHROM'):
                    # This is the header line with sample names
                    parts = line.strip().split('\t')
                    self.samples = parts[9:]
                    return
                else:
                    return # Data started

    def parse_variants(self):
        """Yields variants as dictionaries."""
        with self._open_file() as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                
                # Standard VCF fields
                variant = {
                    'CHROM': parts[0],
                    'POS': int(parts[1]),
                    'ID': parts[2],
                    'REF': parts[3],
                    'ALT': parts[4].split(','),
                    'QUAL': parts[5],
                    'FILTER': parts[6],
                    'INFO': self._parse_info(parts[7]),
                    'FORMAT': parts[8].split(':') if len(parts) > 8 else [],
                    'SAMPLES': parts[9:]  # List of sample data strings
                }
                yield variant

    def _parse_info(self, info_str):
        """Parses the INFO column ID=Value;ID2=Value..."""
        info_dict = {}
        for item in info_str.split(';'):
            if '=' in item:
                key, val = item.split('=', 1)
                info_dict[key] = val
            else:
                info_dict[item] = True
        return info_dict

if __name__ == "__main__":
    # Test on a file if it exists
    test_file = "data/raw/clinvar.vcf.gz"
    if os.path.exists(test_file):
        print(f"Parsing {test_file}...")
        parser = SimpleVCFParser(test_file)
        parser.parse_header()
        print(f"Found {len(parser.header)} header lines.")
        
        count = 0
        for v in parser.parse_variants():
            print(v)
            count += 1
            if count >= 3: break
    else:
        print("Test file not found.")
