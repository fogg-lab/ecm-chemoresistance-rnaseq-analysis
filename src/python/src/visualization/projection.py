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

def make_fixed_len_color_sequence(original_colors: list, n: int):
    """
    Create a new color sequence of length `n` from `original_colors` via interpolation.
    """
    colors_rgb = np.array(
        [
            color.replace("rgb(", "").replace(")", "").split(",")
            for color in original_colors
        ]
    ).astype(int)
    new_colors = []
    color_index_scale_factor = (len(colors_rgb) - 1) / (n - 1)

    for i in range(n):
        position = i * color_index_scale_factor
        index_1 = int(position)
        index_2 = int(np.ceil(position))
        t = position - index_1
        new_color = (1 - t) * colors_rgb[index_1] + t * colors_rgb[index_2]
        new_color = np.round(new_color).astype(int)
        new_colors.append(f"rgb{tuple(new_color)}")

    return new_colors


def interactive_scatterplot(
    components_df: pd.DataFrame,
    color_col: str,
    shape_col: str,
    n_dims: int,
    explained_variance: Optional[np.ndarray] = None,
    title: str = "",
) -> go.Figure:
    """
    Plot an interactive scatter plot with samples colored by a categorical attribute and shaped by another attribute.
    This function works for both PCA and UMAP plots.

    Adapted from: https://stackoverflow.com/a/69733822/21293703

    Parameters:
        components_df (pd.DataFrame): DataFrame containing components and metadata.
        color_col (str): The attribute to color the points by.
        shape_col (str): The attribute to shape the points by.
        n_dims (int): Number of dimensions/components.
        explained_variance (Optional[np.ndarray]): Array of explained variance ratios for each principal component. If None, assume UMAP.

    Returns:
        go.Figure: The plotly figure object.
    """
    components_df = components_df.copy()
    components_df["ID"] = "ID=" + components_df.index.astype(str)

    components_df[color_col] = f"{color_col}=" + components_df[color_col].astype(str)
    color_col_values = np.unique(components_df[color_col])
    n_colors = len(color_col_values)
    palette = px.colors.sequential.Jet[1:-1]
    palette = make_fixed_len_color_sequence(palette, n_colors)
    assert len(color_col_values) == len(palette)
    color_map = dict(zip(color_col_values, palette))

    # Map the shape attribute to shapes
    shapes = ["circle", "diamond", "square", "triangle-up"]
    components_df[shape_col] = f"{shape_col}=" + components_df[shape_col].astype(str)
    shape_col_values = np.unique(components_df[shape_col])
    shape_map = {
        shape_col_values[i]: shapes[i % len(shapes)]
        for i in range(len(shape_col_values))
    }

    if explained_variance is not None:
        dims = [f"PC-{i+1}" for i in range(n_dims)]
        dim_labels = {
            f"PC-{i+1}": f"PC-{i+1} ({explained_variance[i]*100:.2f}% variance)"
            for i in range(n_dims)
        }
    else:
        dims = [f"UMAP-{i+1}" for i in range(n_dims)]
        dim_labels = dict(zip(dims, dims))

    fig = go.Figure(
        go.Scatter3d(
            x=components_df[dims[0]],
            y=components_df[dims[1]],
            z=components_df[dims[2]],
            customdata=components_df.loc[:, ["ID", color_col, shape_col]],
            marker=dict(
                color=components_df[color_col].map(color_map),
                symbol=components_df[shape_col].map(shape_map),
                size=5,
            ),
            mode="markers",
            hovertemplate=f"%{{customdata[0]}}<br>%{{customdata[1]}}<br>%{{customdata[2]}}<extra></extra>",
            showlegend=False,
        )
    ).update_layout(
        template="presentation",
        scene=dict(
            xaxis_title_text=dim_labels[dims[0]],
            yaxis_title_text=dim_labels[dims[1]],
            zaxis_title_text=dim_labels[dims[2]],
        ),
        height=700,
    )

    if n_dims > 3:
        axes = ["x", "y", "z"]
        fig.update_layout(
            updatemenus=[
                {
                    "active": axes.index(ax),
                    "buttons": [
                        {
                            "label": f"{dim_labels[dim]}",
                            "method": "update",
                            "args": [
                                {f"{ax}": [components_df[dim]]},
                                {f"{ax}axis": {"title": {"text": dim_labels[dim]}}},
                                [0],
                            ],
                        }
                        for dim in dims
                    ],
                    "y": 1.1 - (.1 * (axes.index(ax) + 1)),
                }
                for ax in axes
            ]
        ).update_traces(showlegend=False)

    # Add dummy traces for color legend and shape legend
    for value, shape in shape_map.items():
        fig.add_trace(
            go.Scatter3d(
                x=[None],
                y=[None],
                z=[None],
                mode="markers",
                marker=dict(symbol=shape, color="black", size=12),
                showlegend=True,
                name=value,
            )
        )
    for value, color in color_map.items():
        fig.add_trace(
            go.Scatter3d(
                x=[None],
                y=[None],
                z=[None],
                mode="markers",
                marker=dict(color=color, size=5),
                showlegend=True,
                name=value,
            )
        )

    fig.update_layout(
        legend_itemclick=False,
        legend_itemdoubleclick=False,
        title=title,
    )

    fig.show()
    return fig
