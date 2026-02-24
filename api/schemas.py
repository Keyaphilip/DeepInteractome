"""
Pydantic schemas for the DeepInteractome prediction API.
"""
from pydantic import BaseModel, ConfigDict, Field
from typing import List


class VariantInput(BaseModel):
    """A single genomic variant to score."""
    model_config = ConfigDict(
        json_schema_extra={
            "example": {"chrom": "chr1", "pos": 10177, "ref": "A", "alt": "AC", "af": 0.425}
        }
    )
    chrom: str = Field(..., description="Chromosome (e.g. chr1, 1, chrX)")
    pos: int = Field(..., description="1-based genomic position")
    ref: str = Field(..., description="Reference allele")
    alt: str = Field(..., description="Alternate allele")
    af: float = Field(0.0, ge=0.0, le=1.0, description="Allele frequency [0,1]")


class PredictionRequest(BaseModel):
    """Request body containing one or more variants to predict."""
    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "variants": [
                    {"chrom": "chr1", "pos": 10177, "ref": "A", "alt": "AC", "af": 0.425},
                    {"chrom": "chr1", "pos": 10235, "ref": "T", "alt": "TA", "af": 0.001},
                ]
            }
        }
    )
    variants: List[VariantInput] = Field(..., min_length=1)


class VariantPrediction(BaseModel):
    """Prediction result for a single variant."""
    chrom: str
    pos: int
    ref: str
    alt: str
    result: str = Field(..., description="PATHOGENIC or BENIGN")
    pathogenic_probability: float = Field(..., description="Probability of being pathogenic [0,1]")


class PredictionResponse(BaseModel):
    """Response body with predictions for all submitted variants."""
    model_used: str = Field(..., description="Name of the model that made the predictions")
    predictions: List[VariantPrediction]
