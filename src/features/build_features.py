import pandas as pd
import numpy as np
import os
import sys

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from src.utils.vcf_parser import SimpleVCFParser

class GenomicFeatureEngineer:
    def __init__(self, vcf_path):
        self.vcf_path = vcf_path
        self.parser = SimpleVCFParser(vcf_path)
        
    def _one_hot_encode_dna(self, sequence, max_len=1):
        """
        One-hot encodes a DNA sequence.
        A: [1,0,0,0], C: [0,1,0,0], G: [0,0,1,0], T: [0,0,0,1]
        """
        mapping = {
            'A': [1, 0, 0, 0],
            'C': [0, 1, 0, 0],
            'G': [0, 0, 1, 0],
            'T': [0, 0, 0, 1],
            'N': [0, 0, 0, 0] # Unknown
        }
        
        encoding = []
        for base in sequence[:max_len]:
            encoding.extend(mapping.get(base, [0, 0, 0, 0]))
            
        # Pad if shorter than max_len
        current_len = len(sequence)
        if current_len < max_len:
            padding = [0, 0, 0, 0] * (max_len - current_len)
            encoding.extend(padding)
            
        return encoding

    def process_data(self, output_path=None):
        """
        Reads VCF, extracts features, and returns a processed DataFrame.
        """
        print(f"Processing {self.vcf_path}...")
        self.parser.parse_header()
        
        data = []
        for v in self.parser.parse_variants():
            # Target Label: Clinical Significance
            clnsig = v['INFO'].get('CLNSIG', 'Unknown')
            # Binary target: 1 if Pathogenic/Likely_pathogenic, 0 otherwise
            is_pathogenic = 1 if 'Pathogenic' in clnsig or 'Likely_pathogenic' in clnsig else 0
            
            # Feature: Allele Frequency
            try:
                af = float(v['INFO'].get('AF', 0.0))
            except:
                af = 0.0
                
            # Feature: Variant Type length
            ref_len = len(v['REF'])
            alt_len = len(v['ALT'][0])
            
            # Feature: One-Hot Encoded REF and ALT (just first base for simplicity in basic model)
            ref_onehot = self._one_hot_encode_dna(v['REF'], max_len=1)
            alt_onehot = self._one_hot_encode_dna(v['ALT'][0], max_len=1)
            
            row = {
                'CHROM': v['CHROM'],
                'POS': v['POS'],
                'REF_len': ref_len,
                'ALT_len': alt_len,
                'AF': af,
                'Ref_A': ref_onehot[0], 'Ref_C': ref_onehot[1], 'Ref_G': ref_onehot[2], 'Ref_T': ref_onehot[3],
                'Alt_A': alt_onehot[0], 'Alt_C': alt_onehot[1], 'Alt_G': alt_onehot[2], 'Alt_T': alt_onehot[3],
                'Target_Label': is_pathogenic,
                'ClinSig_Raw': clnsig
            }
            data.append(row)
            
        df = pd.DataFrame(data)
        
        if output_path:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            df.to_csv(output_path, index=False)
            print(f"Saved processed data to {output_path}")
            
        return df

if __name__ == "__main__":
    # Test with sample data
    input_vcf = "data/raw/sample.vcf"
    output_csv = "data/processed/genomic_features.csv"
    
    if os.path.exists(input_vcf):
        engineer = GenomicFeatureEngineer(input_vcf)
        df_processed = engineer.process_data(output_path=output_csv)
        print("\nFirst 5 processed rows:")
        print(df_processed[['POS', 'Ref_A', 'Alt_A', 'Target_Label']].head())
    else:
        print(f"Input file {input_vcf} not found.")
