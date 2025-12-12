#!/usr/bin/env python3
"""
Create an HTML report summarizing all modification detection results.
This script generates a comprehensive report with statistics and visualizations.
"""

import pandas as pd
import os
import sys
from datetime import datetime


def create_html_template():
    """
    Create the HTML template for the report.
    
    Returns:
    --------
    str
        HTML template string
    """
    
    template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>RNAModBench Analysis Report</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <style>
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }
        .header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 2rem 0; }
        .section { margin: 2rem 0; padding: 1.5rem; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
        .stats-card { background: #f8f9fa; border-left: 4px solid #007bff; padding: 1rem; margin: 0.5rem 0; }
        .tool-card { border: 1px solid #dee2e6; border-radius: 8px; padding: 1rem; margin: 0.5rem 0; }
        .tool-card:hover { box-shadow: 0 4px 8px rgba(0,0,0,0.1); }
        .table-responsive { max-height: 400px; overflow-y: auto; }
        .plot-container { text-align: center; margin: 1rem 0; }
        .plot-container img { max-width: 100%; height: auto; border-radius: 8px; }
        .footer { background: #343a40; color: white; padding: 1rem 0; margin-top: 2rem; }
        .download-btn { margin: 0.25rem; }
    </style>
</head>
<body>
    <div class="header">
        <div class="container">
            <div class="row align-items-center">
                <div class="col-md-8">
                    <h1><i class="fas fa-dna"></i> RNAModBench Analysis Report</h1>
                    <p class="lead">Comprehensive RNA Modification Detection Analysis</p>
                    <p>Generated on: {generation_date}</p>
                </div>
                <div class="col-md-4 text-end">
                    <i class="fas fa-microscope" style="font-size: 4rem; opacity: 0.7;"></i>
                </div>
            </div>
        </div>
    </div>

    <div class="container">
        <!-- Executive Summary -->
        <div class="section">
            <h2><i class="fas fa-chart-line"></i> Executive Summary</h2>
            <div class="row">
                <div class="col-md-3">
                    <div class="stats-card">
                        <h4>Total Tools</h4>
                        <h2>{total_tools}</h2>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stats-card">
                        <h4>Total Modifications</h4>
                        <h2>{total_modifications}</h2>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stats-card">
                        <h4>Unique Positions</h4>
                        <h2>{unique_positions}</h2>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="stats-card">
                        <h4>Consensus Sites</h4>
                        <h2>{consensus_sites}</h2>
                    </div>
                </div>
            </div>
        </div>

        <!-- Tool Performance -->
        <div class="section">
            <h2><i class="fas fa-tools"></i> Tool Performance</h2>
            <div class="table-responsive">
                {tool_table}
            </div>
        </div>

        <!-- Visualizations -->
        <div class="section">
            <h2><i class="fas fa-chart-bar"></i> Visualizations</h2>
            <div class="row">
                <div class="col-md-6">
                    <div class="plot-container">
                        <h4>Modification Counts by Tool</h4>
                        <img src="modification_counts.png" alt="Modification Counts" class="img-fluid">
                    </div>
                </div>
                <div class="col-md-6">
                    <div class="plot-container">
                        <h4>Chromosome Distribution</h4>
                        <img src="chromosome_distribution.png" alt="Chromosome Distribution" class="img-fluid">
                    </div>
                </div>
            </div>
        </div>

        <!-- Tool Comparison -->
        <div class="section">
            <h2><i class="fas fa-balance-scale"></i> Tool Comparison</h2>
            <p>This section shows the overlap between different modification detection tools.</p>
            <div class="row">
                {tool_comparison_plots}
            </div>
        </div>

        <!-- Method Details -->
        <div class="section">
            <h2><i class="fas fa-info-circle"></i> Method Details</h2>
            <div class="row">
                {method_details}
            </div>
        </div>

        <!-- Download Section -->
        <div class="section">
            <h2><i class="fas fa-download"></i> Download Results</h2>
            <div class="row">
                <div class="col-md-6">
                    <h5>Summary Files</h5>
                    <a href="modification_summary.tsv" class="btn btn-primary download-btn">
                        <i class="fas fa-file-csv"></i> Summary Statistics
                    </a>
                    <a href="tool_comparison.tsv" class="btn btn-primary download-btn">
                        <i class="fas fa-file-csv"></i> Tool Comparison
                    </a>
                </div>
                <div class="col-md-6">
                    <h5>Raw Results</h5>
                    <a href="../" class="btn btn-secondary download-btn">
                        <i class="fas fa-folder-open"></i> All Results
                    </a>
                </div>
            </div>
        </div>

        <!-- Quality Metrics -->
        <div class="section">
            <h2><i class="fas fa-award"></i> Quality Metrics</h2>
            <div class="row">
                <div class="col-md-4">
                    <div class="stats-card">
                        <h5>Average Coverage</h5>
                        <p>{avg_coverage}</p>
                    </div>
                </div>
                <div class="col-md-4">
                    <div class="stats-card">
                        <h5>Modification Rate</h5>
                        <p>{mod_rate}%</p>
                    </div>
                </div>
                <div class="col-md-4">
                    <div class="stats-card">
                        <h5>Strand Bias</h5>
                        <p>{strand_bias}</p>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <div class="footer">
        <div class="container">
            <div class="row">
                <div class="col-md-6">
                    <p>&copy; 2024 RNAModBench Pipeline. Generated by academic research software.</p>
                </div>
                <div class="col-md-6 text-end">
                    <p>For research use only. Please cite the original publication.</p>
                </div>
            </div>
        </div>
    </div>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>
    """
    
    return template

def read_summary_data(summary_file):
    """
    Read summary statistics from file.
    
    Parameters:
    -----------
    summary_file : str
        Path to summary file
    
    Returns:
    --------
    pd.DataFrame
        Summary statistics DataFrame
    """
    
    try:
        return pd.read_csv(summary_file, sep='\t')
    except Exception as e:
        print(f"Error reading summary file: {e}", file=sys.stderr)
        return pd.DataFrame()


def read_comparison_data(comparison_file):
    """
    Read tool comparison data from file.
    
    Parameters:
    -----------
    comparison_file : str
        Path to comparison file
    
    Returns:
    --------
    pd.DataFrame
        Tool comparison DataFrame
    """
    
    try:
        return pd.read_csv(comparison_file, sep='\t')
    except Exception as e:
        print(f"Error reading comparison file: {e}", file=sys.stderr)
        return pd.DataFrame()

def generate_tool_table(summary_df):
    """
    Generate HTML table for tool performance.
    
    Parameters:
    -----------
    summary_df : pd.DataFrame
        Summary statistics DataFrame
    
    Returns:
    --------
    str
        HTML table string
    """
    
    if summary_df.empty:
        return '<p>No summary data available.</p>'
    
    # Convert DataFrame to HTML table
    html_table = summary_df.to_html(
        classes='table table-striped table-hover',
        index=False,
        table_id='toolTable'
    )
    
    return html_table

def generate_method_details():
    """
    Generate method details section.
    
    Returns:
    --------
    str
        HTML string with method details
    """
    
    methods = [
        {
            'name': 'CHEUI',
            'description': 'Deep learning-based m6A detection using signal-level features',
            'features': ['Signal intensity', 'K-mer models', 'Two-stage prediction'],
            'thresholds': 'Probability > 0.999, Ratio > 0.1'
        },
        {
            'name': 'm6Anet',
            'description': 'Neural network model for m6A detection from nanopore signals',
            'features': ['Multiple features', 'Probabilistic output', 'Transcript-level'],
            'thresholds': 'Probability > 0.5, Ratio > 0.1'
        },
        {
            'name': 'Nanocompore',
            'description': 'Comparative analysis of signal intensities between samples',
            'features': ['Statistical testing', 'Gaussian mixture models', 'Signal comparison'],
            'thresholds': 'P-value < 0.05, |Log OR| > 0.5'
        },
        {
            'name': 'ELIGOS2',
            'description': 'Statistical framework for RNA modification detection',
            'features': ['Odds ratio testing', 'Multiple testing correction', 'Strand-specific'],
            'thresholds': 'Adj P-value < 0.0001, OR > 1.2'
        },
        {
            'name': 'Epinano',
            'description': 'Error-based approach for modification detection',
            'features': ['Basecalling errors', 'Sliding window', 'SVM classification'],
            'thresholds': 'Delta error > 0.1'
        },
        {
            'name': 'MINES',
            'description': 'Tombo-based modification detection with machine learning',
            'features': ['Fraction modified', 'Coverage filtering', 'K-mer analysis'],
            'thresholds': 'Coverage > 20, Ratio > 0.1'
        }
    ]
    
    html_content = ""
    for method in methods:
        html_content += f"""
        <div class="col-md-4">
            <div class="tool-card">
                <h5>{method['name']}</h5>
                <p>{method['description']}</p>
                <h6>Key Features:</h6>
                <ul>
        """
        for feature in method['features']:
            html_content += f"<li>{feature}</li>"
        
        html_content += f"""
                </ul>
                <h6>Thresholds:</h6>
                <p><small>{method['thresholds']}</small></p>
            </div>
        </div>
        """
    
    return html_content

def main(summary_file, comparison_file, output_file):
    """
    Main function to create the HTML report.
    
    Parameters:
    -----------
    summary_file : str
        Path to summary statistics file
    comparison_file : str
        Path to tool comparison file
    output_file : str
        Path to output HTML file
    """
    
    # Read data
    summary_df = read_summary_data(summary_file)
    comparison_df = read_comparison_data(comparison_file)
    
    # Calculate statistics
    total_tools = len(summary_df) if not summary_df.empty else 0
    total_modifications = summary_df['Total_Modifications'].sum() if not summary_df.empty else 0
    unique_positions = comparison_df['Position_Key'].nunique() if not comparison_df.empty else 0
    consensus_sites = len(comparison_df[comparison_df['Tool_Count'] > 1]) if not comparison_df.empty else 0
    
    # Generate table
    tool_table = generate_tool_table(summary_df)
    
    # Generate method details
    method_details = generate_method_details()
    
    # Generate tool comparison plots (placeholder for now)
    tool_comparison_plots = ""
    if not comparison_df.empty:
        # Find pairs of tools for Venn diagrams
        tools = comparison_df['Tool'].unique()
        for i, tool1 in enumerate(tools):
            for tool2 in tools[i+1:]:
                plot_file = f"venn_{tool1}_vs_{tool2}.png"
                if os.path.exists(os.path.join(os.path.dirname(output_file), plot_file)):
                    tool_comparison_plots += f"""
                    <div class="col-md-6">
                        <div class="plot-container">
                            <h4>Overlap: {tool1} vs {tool2}</h4>
                            <img src="{plot_file}" alt="Venn Diagram" class="img-fluid">
                        </div>
                    </div>
                    """
    
    if not tool_comparison_plots:
        tool_comparison_plots = "<p>No comparison plots available.</p>"
    
    # Quality metrics (placeholder values)
    avg_coverage = "N/A"
    mod_rate = "N/A"
    strand_bias = "N/A"
    
    # Fill template
    template = create_html_template()
    html_content = template.format(
        generation_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        total_tools=total_tools,
        total_modifications=total_modifications,
        unique_positions=unique_positions,
        consensus_sites=consensus_sites,
        tool_table=tool_table,
        tool_comparison_plots=tool_comparison_plots,
        method_details=method_details,
        avg_coverage=avg_coverage,
        mod_rate=mod_rate,
        strand_bias=strand_bias
    )
    
    # Write HTML file
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"HTML report generated: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python create_report.py <summary_file> <comparison_file> <output_html>")
        sys.exit(1)
    
    summary_file = sys.argv[1]
    comparison_file = sys.argv[2]
    output_file = sys.argv[3]
    
    main(summary_file, comparison_file, output_file)
