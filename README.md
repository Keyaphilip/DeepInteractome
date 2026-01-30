# DeepInteractome: Genomic Disease Prediction

Welcome to your ML journey! This project aims to build an end-to-end machine learning pipeline for predicting diseases from genomic data.

## 🧬 Project Structure

- `data/`: Contains raw and processed genomic data.
  - `raw/`: Original VCF/FASTA files.
  - `processed/`: Cleaned CSV/Parquet files for ML.
- `src/`: Source code.
  - `data/`: Scripts for data downloading and loading.
  - `features/`: Feature engineering.
  - `models/`: Model definitions and training scripts.
- `notebooks/`: Jupyter notebooks for experimentation.
- `docs/`: References and documentation.

## 🚀 Quick Start

### 1. Setup Environment
```bash
# Create virtual environment
python -m venv venv

# Activate (Windows)
.\venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### 2. Download Data
 We use the ClinVar dataset for disease-related variants.
```bash
python src/data/download_data.py
```

### 3. Next Steps
- Explore the data in `notebooks/01_data_exploration.ipynb` (Coming soon!)
- Build your first Logistic Regression model.

## 📚 Resources
- Check `docs/references/` for papers on cardiovascular genetics and pathways.
