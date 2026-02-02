import pandas as pd
import numpy as np
import os
import sys
import joblib
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score, confusion_matrix

def train_model(data_path, model_output_path):
    """
    Trains a Random Forest classifier on the processed genomic data.
    """
    print(f"Loading data from {data_path}...")
    try:
        df = pd.read_csv(data_path)
    except FileNotFoundError:
        print(f"Error: File not found at {data_path}")
        return

    # Prepare features and target
    # Dropping non-feature columns. Keep numeric features.
    # Based on build_features.py, features are Ref_A...T, Alt_A...T, variants length, AF
    
    target_col = 'Target_Label'
    drop_cols = ['CHROM', 'ClinSig_Raw', target_col]
    
    # Ensure we only drop columns that exist
    drop_cols = [c for c in drop_cols if c in df.columns]
    
    X = df.drop(columns=drop_cols)
    y = df[target_col]
    
    print(f"Features: {list(X.columns)}")
    print(f"Target: {target_col}")
    
    # Split data
    # Handle small datasets gracefully
    if len(X) < 10:
        print("Warning: Dataset is very small. Skipping stratification.")
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    else:
        try:
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
        except ValueError as e:
            print(f"Warning: Stratification failed ({e}). Falling back to random split.")
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Initialize and train model
    print("Training RandomForestClassifier...")
    clf = RandomForestClassifier(n_estimators=100, random_state=42)
    clf.fit(X_train, y_train)
    
    # Evaluate
    print("Evaluating model...")
    y_pred = clf.predict(X_test)
    
    acc = accuracy_score(y_test, y_pred)
    print(f"\nAccuracy: {acc:.4f}")
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred))
    print("\nConfusion Matrix:")
    print(confusion_matrix(y_test, y_pred))
    
    # Save model
    os.makedirs(os.path.dirname(model_output_path), exist_ok=True)
    joblib.dump(clf, model_output_path)
    print(f"\nModel saved to {model_output_path}")

if __name__ == "__main__":
    # Define paths relative to project root assuming script is run from project root
    # or handle absolute paths.
    # Let's assume we run from project root.
    
    # Check if we are in src/models or root
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))
    
    input_csv = os.path.join(base_dir, "data/processed/genomic_features.csv")
    output_model = os.path.join(base_dir, "models/rf_classifier.pkl")
    
    if os.path.exists(input_csv):
        train_model(input_csv, output_model)
    else:
        print(f"Data file not found at {input_csv}. Please run src/features/build_features.py first.")
