import torch
import joblib
import pandas as pd
import numpy as np
import os
import sys
from sklearn.preprocessing import StandardScaler

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from src.features.build_features import GenomicFeatureEngineer

# Import the DNN model architecture
from train_dnn import GenomicDNN

def predict_with_random_forest(vcf_path, model_path):
    """
    Make predictions using the Random Forest model.
    """
    print("=== Random Forest Prediction ===")
    print(f"Loading model from {model_path}...")
    
    if not os.path.exists(model_path):
        print(f"Error: Model not found at {model_path}")
        return
    
    # Load model
    rf_model = joblib.load(model_path)
    
    # Process VCF to features
    print(f"Processing VCF file: {vcf_path}")
    engineer = GenomicFeatureEngineer(vcf_path)
    df = engineer.process_data()
    
    # Prepare features (drop non-feature columns)
    # Keep POS as it was included during training
    drop_cols = ['CHROM', 'ClinSig_Raw', 'Target_Label']
    drop_cols = [c for c in drop_cols if c in df.columns]
    X = df.drop(columns=drop_cols)
    
    # Make predictions
    predictions = rf_model.predict(X)
    probabilities = rf_model.predict_proba(X)
    
    # Display results
    print("\n--- Predictions ---")
    for i, (pred, prob) in enumerate(zip(predictions, probabilities)):
        pos = df.iloc[i]['POS'] if 'POS' in df.columns else i
        pathogenic_prob = prob[1] * 100  # Probability of class 1 (pathogenic)
        result = "PATHOGENIC" if pred == 1 else "BENIGN"
        print(f"Variant {i+1} (POS: {pos}): {result} (Pathogenic probability: {pathogenic_prob:.2f}%)")
    
    return predictions, probabilities

def predict_with_dnn(vcf_path, model_path):
    """
    Make predictions using the Deep Neural Network model.
    """
    print("\n=== Deep Neural Network Prediction ===")
    print(f"Loading model from {model_path}...")
    
    if not os.path.exists(model_path):
        print(f"Error: Model not found at {model_path}")
        return
    
    # Process VCF to features
    print(f"Processing VCF file: {vcf_path}")
    engineer = GenomicFeatureEngineer(vcf_path)
    df = engineer.process_data()
    
    # Prepare features (same as Random Forest)
    drop_cols = ['CHROM', 'ClinSig_Raw', 'Target_Label']
    drop_cols = [c for c in drop_cols if c in df.columns]
    X = df.drop(columns=drop_cols).values
    
    # Scale features (important for NN)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Load model
    model = GenomicDNN(input_dim=X.shape[1])
    model.load_state_dict(torch.load(model_path))
    model.eval()
    
    # Make predictions
    with torch.no_grad():
        X_tensor = torch.tensor(X_scaled, dtype=torch.float32)
        predictions_prob = model(X_tensor).numpy()
    
    predictions = (predictions_prob > 0.5).astype(int).flatten()
    
    # Display results
    print("\n--- Predictions ---")
    for i, (pred, prob) in enumerate(zip(predictions, predictions_prob)):
        pos = df.iloc[i]['POS'] if 'POS' in df.columns else i
        pathogenic_prob = prob[0] * 100
        result = "PATHOGENIC" if pred == 1 else "BENIGN"
        print(f"Variant {i+1} (POS: {pos}): {result} (Pathogenic probability: {pathogenic_prob:.2f}%)")
    
    return predictions, predictions_prob

if __name__ == "__main__":
    # Paths
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))
    vcf_file = os.path.join(base_dir, "data/raw/sample.vcf")
    rf_model_path = os.path.join(base_dir, "models/rf_classifier.pkl")
    dnn_model_path = os.path.join(base_dir, "models/dnn_model.pth")
    
    print("DeepInteractome - Variant Pathogenicity Prediction\n")
    
    # Check if VCF exists
    if not os.path.exists(vcf_file):
        print(f"Error: VCF file not found at {vcf_file}")
        print("Please provide a valid VCF file path.")
        sys.exit(1)
    
    # Run predictions with both models
    print(f"Input VCF: {vcf_file}\n")
    
    # Random Forest
    if os.path.exists(rf_model_path):
        rf_preds, rf_probs = predict_with_random_forest(vcf_file, rf_model_path)
    else:
        print(f"Random Forest model not found. Skipping...")
    
    # Deep Neural Network
    if os.path.exists(dnn_model_path):
        dnn_preds, dnn_probs = predict_with_dnn(vcf_file, dnn_model_path)
    else:
        print(f"DNN model not found. Skipping...")
    
    print("\n=== Prediction Complete ===")
