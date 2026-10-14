# 🧬 PCA on Gene Expression: Tumor vs Normal

This project uses PCA to compress high-dimensional gene expression data (500 genes) into 2D, so that **Tumor** and **Normal** samples can be told apart. It then uses the **PC1 loadings** to find candidate biomarker genes.

## 🎯 What it does
1. Generates a simulated log2 expression matrix (100 samples × 500 genes, with 20 true biomarker genes). You can also load your own CSV with `--input`.
2. Scales every gene with `StandardScaler`.
3. Fits `PCA` (10 components).
4. Creates plots: a scree plot, a PC1 vs PC2 scatter plot, and the top loading genes.
5. Saves the top genes to a CSV file.

## 📂 Structure
```
pca-gene-expression-biomarkers/
├── pca_gene_expression.py
├── requirements.txt
├── .gitignore
└── results/   (plots + CSV, created automatically)
```

## ▶️ How to run
```bash
pip install -r requirements.txt
python pca_gene_expression.py
# with your own data:
python pca_gene_expression.py --input my_expression.csv   # rows = samples, columns = genes + 'label'
```

## 📊 Output

### 1. Scree plot + sample scatter plot
**Left:** how much variance each PC explains. **Right:** every dot is one sample. Tumor (red) and Normal (green) separate along PC1.

![PCA overview](results/pca_overview.png)

### 2. Top 15 genes driving PC1
Longer bar = bigger effect on PC1. Red = positive loading, green = negative loading. These are the candidate biomarker genes.

![Top PC1 genes](results/top_pc1_genes.png)

### 3. CSV file
`results/top_pc1_genes.csv` has the same top genes with their loading values.

## 🔬 Key insight
PC1 explains only about 3.5% of the variance, yet Tumor and Normal samples still separate cleanly. This is because the signal sits in just a few genes. The top PC1 loading genes (GENE_001–020) are exactly the ones that were truly different in the simulation.

## 🚀 Next steps
- Try a real dataset (GEO / TCGA / Golub leukemia)
- Add a Logistic Regression or Random Forest classifier after PCA
- Compare with t-SNE / UMAP

## 🛠️ Tech stack
Python · NumPy · pandas · scikit-learn · Matplotlib
## Complete Workflow

![PCA Gene Expression Workflow](images/pca_workflow.png)
