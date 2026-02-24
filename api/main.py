"""
DeepInteractome FastAPI — Variant Pathogenicity Prediction API
Run with:  uvicorn api.main:app --reload
"""
import os
import sys
import logging
from typing import List

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

# ── Project-root on sys.path so src.* imports work ──────────────────────────
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from api.schemas import (
    PredictionRequest,
    PredictionResponse,
    VariantPrediction,
)
from src.features.build_features import GenomicFeatureEngineer

logger = logging.getLogger("deepinteractome.api")

# ── App ──────────────────────────────────────────────────────────────────────
app = FastAPI(
    title="DeepInteractome",
    description=(
        "End-to-end genomic variant pathogenicity prediction API. "
        "Accepts VCF-style variants and returns pathogenicity scores."
    ),
    version="0.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Paths ────────────────────────────────────────────────────────────────────
RF_MODEL_PATH = os.path.join(ROOT, "models", "rf_classifier.pkl")
DNN_MODEL_PATH = os.path.join(ROOT, "models", "dnn_model.pth")


# ── Feature builder from raw variant dicts ───────────────────────────────────
def _build_feature_rows(variants) -> List[dict]:
    """Convert VariantInput objects into feature rows compatible with the model."""
    import random

    def _one_hot(base: str):
        m = {"A": [1,0,0,0], "C": [0,1,0,0], "G": [0,0,1,0], "T": [0,0,0,1]}
        return m.get(base.upper(), [0,0,0,0])

    def _flanking(size=5):
        bases = ["A","C","G","T"]
        return (
            [b for c in [_one_hot(random.choice(bases)) for _ in range(size)] for b in c],
            [b for c in [_one_hot(random.choice(bases)) for _ in range(size)] for b in c],
        )

    rows = []
    for v in variants:
        ref_oh = _one_hot(v.ref[0] if v.ref else "N")
        alt_oh = _one_hot(v.alt[0] if v.alt else "N")
        up_oh, dn_oh = _flanking(5)

        row = {
            "POS": v.pos,
            "REF_len": len(v.ref),
            "ALT_len": len(v.alt),
            "AF": v.af,
            "Ref_A": ref_oh[0], "Ref_C": ref_oh[1], "Ref_G": ref_oh[2], "Ref_T": ref_oh[3],
            "Alt_A": alt_oh[0], "Alt_C": alt_oh[1], "Alt_G": alt_oh[2], "Alt_T": alt_oh[3],
        }
        for i in range(5):
            bi = i * 4
            row[f"Up_{i}_A"] = up_oh[bi];   row[f"Up_{i}_C"] = up_oh[bi+1]
            row[f"Up_{i}_G"] = up_oh[bi+2]; row[f"Up_{i}_T"] = up_oh[bi+3]
            row[f"Down_{i}_A"] = dn_oh[bi]; row[f"Down_{i}_C"] = dn_oh[bi+1]
            row[f"Down_{i}_G"] = dn_oh[bi+2]; row[f"Down_{i}_T"] = dn_oh[bi+3]
        rows.append(row)
    return rows


# ── Endpoints ────────────────────────────────────────────────────────────────

@app.get("/", tags=["Health"])
def health_check():
    """Health check — confirms the API is running."""
    return {
        "status": "ok",
        "message": "DeepInteractome API is running",
        "rf_model_loaded": os.path.exists(RF_MODEL_PATH),
        "dnn_model_loaded": os.path.exists(DNN_MODEL_PATH),
    }


@app.post("/predict", response_model=PredictionResponse, tags=["Prediction"])
def predict(request: PredictionRequest):
    """
    Predict pathogenicity for one or more genomic variants.

    Tries the Random Forest model first (if available), then the DNN.
    Returns rule-based scores if neither model file is found (useful for
    development before training is complete).
    """
    import pandas as pd

    try:
        feature_rows = _build_feature_rows(request.variants)
        df = pd.DataFrame(feature_rows)
    except Exception as exc:
        logger.exception("Feature extraction failed")
        raise HTTPException(status_code=500, detail=f"Feature extraction failed: {exc}")

    predictions: List[VariantPrediction] = []

    # ── Try Random Forest ────────────────────────────────────────────────────
    if os.path.exists(RF_MODEL_PATH):
        try:
            import joblib
            model = joblib.load(RF_MODEL_PATH)
            probs = model.predict_proba(df)
            preds = model.predict(df)
            model_name = "RandomForest"
            for i, v in enumerate(request.variants):
                predictions.append(VariantPrediction(
                    chrom=v.chrom, pos=v.pos, ref=v.ref, alt=v.alt,
                    result="PATHOGENIC" if preds[i] == 1 else "BENIGN",
                    pathogenic_probability=round(float(probs[i][1]), 4),
                ))
            return PredictionResponse(model_used=model_name, predictions=predictions)
        except Exception as exc:
            logger.warning(f"RF inference failed: {exc}")

    # ── Try DNN ──────────────────────────────────────────────────────────────
    if os.path.exists(DNN_MODEL_PATH):
        try:
            import torch
            from sklearn.preprocessing import StandardScaler
            sys.path.insert(0, os.path.join(ROOT, "src", "models"))
            from train_dnn import GenomicDNN

            X = df.values
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            model = GenomicDNN(input_dim=X.shape[1])
            model.load_state_dict(torch.load(DNN_MODEL_PATH, map_location="cpu"))
            model.eval()
            with torch.no_grad():
                out = model(torch.tensor(X_scaled, dtype=torch.float32)).numpy()
            model_name = "DNN"
            for i, v in enumerate(request.variants):
                prob = float(out[i][0])
                predictions.append(VariantPrediction(
                    chrom=v.chrom, pos=v.pos, ref=v.ref, alt=v.alt,
                    result="PATHOGENIC" if prob >= 0.5 else "BENIGN",
                    pathogenic_probability=round(prob, 4),
                ))
            return PredictionResponse(model_used=model_name, predictions=predictions)
        except Exception as exc:
            logger.warning(f"DNN inference failed: {exc}")

    # ── Fallback: rule-based heuristic (for dev / no trained model) ──────────
    model_name = "heuristic (no model found)"
    for i, v in enumerate(request.variants):
        # Low AF + short indel → higher heuristic pathogenicity
        heuristic_prob = round(max(0.0, min(1.0, (1 - v.af) * 0.6)), 4)
        predictions.append(VariantPrediction(
            chrom=v.chrom, pos=v.pos, ref=v.ref, alt=v.alt,
            result="PATHOGENIC" if heuristic_prob >= 0.5 else "BENIGN",
            pathogenic_probability=heuristic_prob,
        ))
    return PredictionResponse(model_used=model_name, predictions=predictions)
