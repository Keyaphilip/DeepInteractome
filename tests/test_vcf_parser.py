import pytest
import os
import sys

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.utils.vcf_parser import SimpleVCFParser

@pytest.fixture
def sample_vcf(tmp_path):
    """Creates a temporary VCF file for testing."""
    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t100\trs1\tA\tC\t.\tPASS\tAF=0.5;CLNSIG=Pathogenic\tGT\t0/1
1\t200\trs2\tG\tT\t.\tPASS\tAF=0.01;CLNSIG=Benign\tGT\t0/0
"""
    p = tmp_path / "test.vcf"
    p.write_text(vcf_content, encoding='utf-8')
    return str(p)

def test_parse_header(sample_vcf):
    parser = SimpleVCFParser(sample_vcf)
    parser.parse_header()
    assert len(parser.header) == 3
    assert parser.samples == ['SAMPLE1']

def test_parse_variants(sample_vcf):
    parser = SimpleVCFParser(sample_vcf)
    variants = list(parser.parse_variants())
    
    assert len(variants) == 2
    
    v1 = variants[0]
    assert v1['CHROM'] == '1'
    assert v1['POS'] == 100
    assert v1['REF'] == 'A'
    assert v1['ALT'] == ['C']
    assert v1['INFO']['AF'] == '0.5'
    assert v1['INFO']['CLNSIG'] == 'Pathogenic'

    v2 = variants[1]
    assert v2['POS'] == 200
    assert v2['INFO']['CLNSIG'] == 'Benign'
