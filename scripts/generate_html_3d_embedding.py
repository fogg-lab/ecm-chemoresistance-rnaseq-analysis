"""Generate an interactive 3D scatter plot using scatter-gl.

it'll look like this: https://pair-code.github.io/scatter-gl/,
but with fewer toggles, and more informative point info on hover or click.

Example usage:
python generate_html_3d_embedding.py "Thing" "3D PCA" html_3d_data/pca_ex.csv html_3d_data/metadata_ex.csv html_3d_data/template.html PCA_3D.html
"""

import argparse

def generate_html(title1, title2, data_csv_path, metadata_csv_path, template_path, output_path):
    # Read the template
    with open(template_path, 'r') as template_file:
        template = template_file.read()

    # Read in data and metadata CSV files as strings
    with open(data_csv_path, 'r') as data_csv_file:
        data_csv_content = data_csv_file.read().strip()
    with open(metadata_csv_path, 'r') as metadata_csv_file:
        metadata_csv_content = metadata_csv_file.read().strip()

    # Replace placeholders in the template
    html_content = template.replace('{{TITLE1}}', title1)
    html_content = html_content.replace('{{TITLE2}}', title2)
    html_content = html_content.replace('ID,x,y,z\na,1,2,3\nb,2,3,4', data_csv_content)
    html_content = html_content.replace('ID,ch1,ch2,ch3\na,A,B,C\nb,B,C,D', metadata_csv_content)

    # Write the generated HTML to the output file
    with open(output_path, 'w') as output_file:
        output_file.write(html_content)

def main():
    parser = argparse.ArgumentParser(description='Generate a standalone HTML file for an interactive ScatterGL plot.')
    parser.add_argument('title1', help='Window title')
    parser.add_argument('title2', help='Plot title')
    parser.add_argument('data_csv', help='Path to the data CSV file')
    parser.add_argument('metadata_csv', help='Path to the metadata CSV file')
    parser.add_argument('template', help='Path to the HTML template file')
    parser.add_argument('output', help='Path for the output HTML file')

    args = parser.parse_args()

    generate_html(args.title1, args.title2, args.data_csv, args.metadata_csv, args.template, args.output)
    print(f"HTML file generated successfully: {args.output}")

if __name__ == "__main__":
    main()
