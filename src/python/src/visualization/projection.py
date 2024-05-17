from typing import Optional

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import umap
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio


def plot_explained_variance(df: pd.DataFrame):
    """
    Plot the cumulative explained variance.

    Parameters:
        df (pd.DataFrame): The input data frame.
    """
    # Fit PCA and calculate the cumulative explained variance (CEV)
    pca = PCA()
    pca.fit(df)
    cev = pca.explained_variance_ratio_.cumsum()

    # Plot the cumulative explained variance
    plt.figure(figsize=(8, 6))
    plt.plot(range(1, len(cev) + 1), cev, marker="o")
    plt.xlabel("# Principal Components", fontsize=16)
    plt.ylabel("CEV", fontsize=16)
    plt.title("Explained Variance by # of Principal Components", fontsize=18)
    plt.grid(True)
    plt.axhline(y=0.90, color="r", linestyle="--", label="90% Variance")
    plt.axhline(y=0.95, color="g", linestyle="--", label="95% Variance")
    plt.legend(fontsize=14)
    plt.show()


def interactive_scatterplot(data: pd.DataFrame, attributes: pd.DataFrame, title: str,
                            num_dims: int, technique: str = "PCA") -> go.Figure:
    """
    Generate an interactive scatterplot with plotly.

    Parameters:
        data (pd.DataFrame): The input data, 1 row per point.
        attributes (pd.DataFrame): The attributes/labels for the points.
        title (str): The title of the plot.
        technique (str): Dimensionality reduction technique, PCA or UMAP. Default="PCA".

    Returns:
        go.Figure: The interactive scatterplot ready to be displayed or saved.
    """


