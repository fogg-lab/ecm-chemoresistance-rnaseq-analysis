import argparse
import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go


def read_data(counts_file: str, coldata_file: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read the counts and coldata files.

    Args:
        counts_file (str): Path to the VST-transformed counts CSV file.
        coldata_file (str): Path to the coldata CSV file.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: A tuple containing the counts and coldata DataFrames.
    """
    counts = pd.read_csv(counts_file, index_col=0)
    counts.index.name = "Ensembl_gene_id"

    coldata = pd.read_csv(coldata_file, index_col=0)
    coldata.index.name = "sample_id"

    return counts, coldata


def calculate_correlation(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the sample-to-sample correlation matrix using np.corrcoef.

    Args:
        counts (pd.DataFrame): The VST-transformed counts DataFrame.

    Returns:
        pd.DataFrame: The correlation matrix.
    """
    corr_matrix = np.corrcoef(counts.T)
    return pd.DataFrame(corr_matrix, index=counts.columns, columns=counts.columns)


def create_heatmap(
    corr_matrix: pd.DataFrame, coldata: pd.DataFrame, output_dir: str
) -> None:
    """
    Create and save an interactive heatmap of sample-to-sample correlation using Plotly.

    Args:
        corr_matrix (pd.DataFrame): The correlation matrix.
        coldata (pd.DataFrame): The coldata DataFrame.
        output_dir (str): Directory to save the output HTML file.

    """
    hover_text = []
    for sample in corr_matrix.index:
        sample_info = [f"{sample}"]
        for trait, value in coldata.loc[sample].items():
            sample_info.append(f"{trait}: {value}")
        hover_text.append("<br>".join(sample_info))

    fig = go.Figure(
        data=go.Heatmap(
            z=corr_matrix.values,
            x=corr_matrix.columns,
            y=corr_matrix.index,
            hoverongaps=False,
            text=[
                [hover_text[i]] * len(corr_matrix.columns)
                for i in range(len(corr_matrix.index))
            ],
            hoverinfo="text+z",
            colorscale="Viridis",
        )
    )

    fig.update_layout(
        title="Sample-to-Sample Correlation Heatmap",
        xaxis_title="Samples",
        yaxis_title="Samples",
        width=800,
        height=800,
    )

    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "sample_correlation_heatmap.html")
    fig.write_html(output_file)
    print(f"Saved sample correlation heatmap to {output_file}")


def main(counts_file: str, coldata_file: str, output_dir: str) -> None:
    """
    Main function to create a sample-to-sample correlation heatmap.

    Args:
        counts_file (str): Path to the VST-transformed counts CSV file.
        coldata_file (str): Path to the coldata CSV file.
        output_dir (str): Directory to save the output HTML file.

    """
    counts, coldata = read_data(counts_file, coldata_file)
    corr_matrix = calculate_correlation(counts)
    create_heatmap(corr_matrix, coldata, output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create a sample-to-sample correlation heatmap using Plotly"
    )
    parser.add_argument(
        "counts_file", help="Path to the VST-transformed counts CSV file"
    )
    parser.add_argument("coldata_file", help="Path to the coldata CSV file")
    parser.add_argument("output_dir", help="Directory to save the output HTML file")

    args = parser.parse_args()

    main(args.counts_file, args.coldata_file, args.output_dir)
