import pytest
import os
import sys
import pandas as pd

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.features.build_features import GenomicFeatureEngineer

@pytest.fixture
def mock_vcf(tmp_path):
    """Creates a temporary VCF file for testing features."""
    vcf_content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tC\t.\t.\tAF=0.5;CLNSIG=Pathogenic
"""
    p = tmp_path / "test_features.vcf"
    p.write_text(vcf_content, encoding='utf-8')
    return str(p)

def test_one_hot_encoding(mock_vcf):
    engineer = GenomicFeatureEngineer(mock_vcf)
    
    # Test A
    assert engineer._one_hot_encode_dna('A') == [1, 0, 0, 0]
    # Test C
    assert engineer._one_hot_encode_dna('C') == [0, 1, 0, 0]
    # Test G
    assert engineer._one_hot_encode_dna('G') == [0, 0, 1, 0]
    # Test T
    assert engineer._one_hot_encode_dna('T') == [0, 0, 0, 1]
    # Test N (unknown)
    assert engineer._one_hot_encode_dna('N') == [0, 0, 0, 0]

def test_process_data(mock_vcf):
    engineer = GenomicFeatureEngineer(mock_vcf)
    df = engineer.process_data()
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert df.iloc[0]['POS'] == 100
    assert df.iloc[0]['Ref_A'] == 1 # REF was A
    assert df.iloc[0]['Alt_C'] == 1 # ALT was C (index 1 in 1-hot depends on mapping: A, C, G, T)
    # Mapping is A:0, C:1, G:2, T:3.
    # So if ALT is C, Alt_onehot should be [0, 1, 0, 0].
    # So Alt_C corresponds to index 1.
    # Let's check the column name logic in build_features.py:
    # 'Alt_A': alt_onehot[0], 'Alt_C': alt_onehot[1]...
    assert df.iloc[0]['Alt_C'] == 1
    assert df.iloc[0]['Target_Label'] == 1 # Pathogenic
