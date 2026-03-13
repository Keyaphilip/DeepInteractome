import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as plt_sns
import joblib
import torch
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score, confusion_matrix, ConfusionMatrixDisplay
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from src.models.train_dnn import GenomicDNN

def evaluate_models(data_path, rf_model_path, dnn_model_path, output_dir):
    """Evaluates both the RF and DNN models and saves plots."""
    print("Loading data for evaluation...")
    if not os.path.exists(data_path):
        print(f"Error: {data_path} not found.")
        return

    os.makedirs(output_dir, exist_ok=True)
    
    df = pd.read_csv(data_path)
    target_col = 'Target_Label'
    drop_cols = ['CHROM', 'ClinSig_Raw', target_col]
    drop_cols = [c for c in drop_cols if c in df.columns]
    
    X = df.drop(columns=drop_cols)
    y = df[target_col].values
    feature_names = X.columns.tolist()
    
    # Same split as training for consistent validation set
    X_train, X_test, y_train, y_test = train_test_split(X.values, y, test_size=0.2, random_state=42)
    
    # Scale Data
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    results = {}

    plt.style.use('seaborn-v0_8-whitegrid')
    
    # --- Random Forest Evaluation ---
    if os.path.exists(rf_model_path):
        print("Evaluating Random Forest...")
        rf_model = joblib.load(rf_model_path)
        y_prob_rf = rf_model.predict_proba(X_test_scaled)[:, 1]
        results['Random Forest'] = y_prob_rf
        
        # Feature Importance Plot
        importance = rf_model.feature_importances_
        indices = np.argsort(importance)[::-1][:20] # Top 20
        plt.figure(figsize=(10, 8))
        plt.title("Random Forest Top 20 Feature Importances")
        plt.bar(range(len(indices)), importance[indices], align="center", color='skyblue')
        plt.xticks(range(len(indices)), [feature_names[i] for i in indices], rotation=90)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'rf_feature_importance.png'), dpi=300)
        plt.close()
    
    # --- DNN Evaluation ---
    if os.path.exists(dnn_model_path):
        print("Evaluating Deep Neural Network...")
        model = GenomicDNN(input_dim=X_test_scaled.shape[1])
        model.load_state_dict(torch.load(dnn_model_path))
        model.eval()
        with torch.no_grad():
            y_prob_dnn = model(torch.tensor(X_test_scaled, dtype=torch.float32)).numpy().ravel()
        results['DNN'] = y_prob_dnn

    # --- Plot ROC Curves ---
    if results:
        plt.figure(figsize=(8, 8))
        for name, y_prob in results.items():
            fpr, tpr, _ = roc_curve(y_test, y_prob)
            roc_auc = auc(fpr, tpr)
            plt.plot(fpr, tpr, lw=2, label=f'{name} (area = {roc_auc:.3f})')
        
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic (ROC)')
        plt.legend(loc="lower right")
        plt.savefig(os.path.join(output_dir, 'roc_curve.png'), dpi=300)
        plt.close()

    # --- Plot Precision-Recall Curves ---
    if results:
        plt.figure(figsize=(8, 8))
        for name, y_prob in results.items():
            precision, recall, _ = precision_recall_curve(y_test, y_prob)
            pr_auc = average_precision_score(y_test, y_prob)
            plt.plot(recall, precision, lw=2, label=f'{name} (AP = {pr_auc:.3f})')
            
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('Precision-Recall Curve')
        plt.legend(loc="lower left")
        plt.savefig(os.path.join(output_dir, 'pr_curve.png'), dpi=300)
        plt.close()
        
    print(f"Evaluation complete. Plots saved to {output_dir}")

if __name__ == "__main__":
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))
    data_path = os.path.join(base_dir, "data/processed/genomic_features.csv")
    rf_model_path = os.path.join(base_dir, "models/rf_classifier.pkl")
    dnn_model_path = os.path.join(base_dir, "models/dnn_model.pth")
    output_dir = os.path.join(base_dir, "reports/figures")
    
    evaluate_models(data_path, rf_model_path, dnn_model_path, output_dir)
