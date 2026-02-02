import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import numpy as np
import os
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# Define the Model
class GenomicDNN(nn.Module):
    def __init__(self, input_dim):
        super(GenomicDNN, self).__init__()
        self.fc1 = nn.Linear(input_dim, 64)
        self.relu1 = nn.ReLU()
        self.dropout1 = nn.Dropout(0.3)
        self.fc2 = nn.Linear(64, 32)
        self.relu2 = nn.ReLU()
        self.dropout2 = nn.Dropout(0.2)
        self.fc3 = nn.Linear(32, 1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.fc1(x)
        x = self.relu1(x)
        x = self.dropout1(x)
        x = self.fc2(x)
        x = self.relu2(x)
        x = self.dropout2(x)
        x = self.fc3(x)
        x = self.sigmoid(x)
        return x

# Custom Dataset
class GenomicDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.float32).unsqueeze(1) # [batch, 1]

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

def train_dnn(data_path, model_output_path, epochs=20, batch_size=32):
    print(f"Loading data from {data_path}...")
    if not os.path.exists(data_path):
        print(f"Error: {data_path} not found.")
        return

    df = pd.read_csv(data_path)
    
    target_col = 'Target_Label'
    drop_cols = ['CHROM', 'ClinSig_Raw', target_col]
    drop_cols = [c for c in drop_cols if c in df.columns]
    
    X = df.drop(columns=drop_cols).values
    y = df[target_col].values
    
    print(f"Input Feature Dimension: {X.shape[1]}")
    
    # Split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Scale Data (Important for NN)
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)
    
    # create Datasets
    train_dataset = GenomicDataset(X_train, y_train)
    test_dataset = GenomicDataset(X_test, y_test)
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    
    # Init Model
    model = GenomicDNN(input_dim=X.shape[1])
    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    
    best_loss = float('inf')
    
    print("Starting training...")
    for epoch in range(epochs):
        model.train()
        running_loss = 0.0
        for inputs, labels in train_loader:
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
            
        avg_loss = running_loss / len(train_loader)
        
        # Validation (simple, on test set for now)
        model.eval()
        with torch.no_grad():
            outputs = model(torch.tensor(X_test, dtype=torch.float32))
            val_loss = criterion(outputs, torch.tensor(y_test, dtype=torch.float32).unsqueeze(1)).item()
        
        print(f"Epoch {epoch+1}/{epochs}, Loss: {avg_loss:.4f}, Val Loss: {val_loss:.4f}")
        
        if val_loss < best_loss:
            best_loss = val_loss
            # Save best model
            os.makedirs(os.path.dirname(model_output_path), exist_ok=True)
            torch.save(model.state_dict(), model_output_path)
            
    print(f"Best model saved to {model_output_path} with Val Loss: {best_loss:.4f}")

if __name__ == "__main__":
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))
    input_csv = os.path.join(base_dir, "data/processed/genomic_features.csv")
    output_model = os.path.join(base_dir, "models/dnn_model.pth")
    
    train_dnn(input_csv, output_model)
