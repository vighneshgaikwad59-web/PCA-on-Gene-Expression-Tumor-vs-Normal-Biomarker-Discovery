import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

OUT = Path("pca_output")


def simulate_expression(n_per_group=50, n_genes=500, n_markers=20, seed=42):
    """Simulated log2 expression: first n_markers genes are up/down in tumor."""
    rng = np.random.default_rng(seed)
    genes = [f"GENE_{i+1:03d}" for i in range(n_genes)]
    base = rng.normal(8, 1.5, n_genes)                      # per-gene baseline
    normal = base + rng.normal(0, 1, (n_per_group, n_genes))
    tumor = base + rng.normal(0, 1, (n_per_group, n_genes))
    shift = rng.choice([-1, 1], n_markers) * rng.uniform(1.5, 3, n_markers)
    tumor[:, :n_markers] += shift                           # true biomarkers
    X = np.vstack([normal, tumor])
    labels = ["Normal"] * n_per_group + ["Tumor"] * n_per_group
    df = pd.DataFrame(X, columns=genes)
    df["label"] = labels
    return df


def load_data(path):
    return simulate_expression() if path is None else pd.read_csv(path)


def run_pca(df, n_components=10):
    X = df.drop(columns="label")
    Xs = StandardScaler().fit_transform(X)                  # genes -> mean 0, std 1
    pca = PCA(n_components=n_components, random_state=42)
    scores = pca.fit_transform(Xs)
    return X.columns, pca, scores


def plot_all(df, genes, pca, scores):
    OUT.mkdir(exist_ok=True)
    labels = df["label"].to_numpy()
    colors = {"Normal": "#2a9d8f", "Tumor": "#e63946"}
    evr = pca.explained_variance_ratio_ * 100

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ax[0].bar(range(1, len(evr) + 1), evr, color="#457b9d", label="Per PC")
    ax[0].plot(range(1, len(evr) + 1), np.cumsum(evr), "o-", color="#e76f51", label="Cumulative")
    ax[0].set(xlabel="Principal Component", ylabel="Variance explained (%)", title="Scree Plot")
    ax[0].legend()
    for g, c in colors.items():
        m = labels == g
        ax[1].scatter(scores[m, 0], scores[m, 1], c=c, label=g, alpha=0.8, edgecolor="k")
    ax[1].set(xlabel=f"PC1 ({evr[0]:.1f}%)", ylabel=f"PC2 ({evr[1]:.1f}%)", title="Samples in PC space")
    ax[1].legend()
    plt.tight_layout()
    plt.savefig(OUT / "pca_overview.png", dpi=150)
    plt.close()

    top = pd.Series(pca.components_[0], index=genes).sort_values(key=abs, ascending=False).head(15)
    plt.figure(figsize=(8, 5))
    top[::-1].plot.barh(color=["#e63946" if v > 0 else "#2a9d8f" for v in top[::-1]])
    plt.xlabel("PC1 loading")
    plt.title("Top 15 genes driving PC1 (candidate biomarkers)")
    plt.tight_layout()
    plt.savefig(OUT / "top_pc1_genes.png", dpi=150)
    plt.close()
    top.rename("PC1_loading").to_csv(OUT / "top_pc1_genes.csv")
    return top


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", help="CSV: samples x genes + 'label' column")
    args = ap.parse_args()

    df = load_data(args.input)
    genes, pca, scores = run_pca(df)
    top = plot_all(df, genes, pca, scores)

    evr = pca.explained_variance_ratio_ * 100
    print(f"Samples: {len(df)} | Genes: {len(genes)}")
    print(f"PC1 = {evr[0]:.1f}% | PC2 = {evr[1]:.1f}% | PC1+PC2 = {evr[:2].sum():.1f}%")
    print("\nTop 10 genes driving PC1:")
    print(top.head(10).round(3).to_string())
    print(f"\nSaved plots + CSV in {OUT}/")


if __name__ == "__main__":
    main()