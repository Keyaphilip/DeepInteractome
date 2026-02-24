"""
Integration tests for the DeepInteractome FastAPI application.
Uses FastAPI's built-in TestClient (based on httpx / requests).
"""
import sys
import os
import pytest

# Ensure project root is importable
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from fastapi.testclient import TestClient
from api.main import app

client = TestClient(app)


# ── Health check ─────────────────────────────────────────────────────────────

def test_health_check():
    """GET / should return 200 with status ok."""
    response = client.get("/")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert "message" in data


# ── Predict endpoint ─────────────────────────────────────────────────────────

SAMPLE_PAYLOAD = {
    "variants": [
        {"chrom": "chr1", "pos": 10177, "ref": "A", "alt": "AC", "af": 0.425},
        {"chrom": "chr1", "pos": 10235, "ref": "T", "alt": "TA", "af": 0.001},
        {"chrom": "chr1", "pos": 10505, "ref": "G", "alt": "C",  "af": 0.15},
    ]
}


def test_predict_returns_200():
    """POST /predict with valid data should return 200."""
    response = client.post("/predict", json=SAMPLE_PAYLOAD)
    assert response.status_code == 200


def test_predict_response_schema():
    """Response must include model_used and predictions list."""
    response = client.post("/predict", json=SAMPLE_PAYLOAD)
    data = response.json()
    assert "model_used" in data
    assert "predictions" in data
    assert isinstance(data["predictions"], list)
    assert len(data["predictions"]) == 3


def test_predict_each_result_has_required_fields():
    """Each prediction entry must have all expected fields."""
    response = client.post("/predict", json=SAMPLE_PAYLOAD)
    for pred in response.json()["predictions"]:
        assert "chrom" in pred
        assert "pos" in pred
        assert "ref" in pred
        assert "alt" in pred
        assert pred["result"] in ("PATHOGENIC", "BENIGN")
        assert 0.0 <= pred["pathogenic_probability"] <= 1.0


def test_predict_empty_variants_rejected():
    """POST /predict with empty variants list should fail validation (422)."""
    response = client.post("/predict", json={"variants": []})
    assert response.status_code == 422


def test_predict_missing_body_rejected():
    """POST /predict with no body should return 422."""
    response = client.post("/predict")
    assert response.status_code == 422


def test_predict_af_out_of_range_rejected():
    """AF > 1.0 should fail Pydantic validation."""
    response = client.post("/predict", json={
        "variants": [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "C", "af": 1.5}]
    })
    assert response.status_code == 422
