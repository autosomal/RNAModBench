#!/usr/bin/env python3
"""
Generate summary statistics and comparison across all modification detection tools.
This script aggregates results from all tools and creates comparative analyses.
"""

import pandas as pd
import numpy as np
import os
import sys
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn3


def read_tool_results(input_files):
    """
    Read and combine results from all tools.
    
    Parameters:
    -----------
    input_files : list
        List of processed result files from all tools
    
    Returns:
    --------
    dict
        Dictionary with tool names as keys and DataFrames as values
    """
    
    tool_results = {}
    
    for file_path in input_files:
        # Extract tool name from file path
        tool_name = None
        if 'CHEUI' in file_path:
            tool_name = 'CHEUI'
        elif 'ELIGOS2' in file_path:
            tool_name = 'ELIGOS2'
        elif 'm6Anet' in file_path:
            tool_name = 'm6Anet'
        elif 'Nanocompore' in file_path:
            tool_name = 'Nanocompore'
        elif 'DENA' in file_path:
            tool_name = 'DENA'
        elif 'Epinano' in file_path:
            tool_name = 'Epinano'
        elif 'MINES' in file_path:
            tool_name = 'MINES'
        elif 'xPore' in file_path:
            tool_name = 'xPore'
        elif 'yanocomp' in file_path:
            tool_name = 'yanocomp'
        elif 'NanoSPA' in file_path:
            tool_name = 'NanoSPA'
        
        if tool_name:
            try:
                df = pd.read_csv(file_path, sep='\t')
                df['Tool'] = tool_name
                df['File'] = file_path
                tool_results[tool_name] = df
                print(f"Loaded {len(df)} modifications from {tool_name}")
            except Exception as e:
                print(f"Error reading {file_path}: {e}", file=sys.stderr)
    
    return tool_results

def read_tool_results(input_files):
    """
    Read and combine results from all tools.
    
    Parameters:
    -----------
    input_files : list
        List of processed result files from all tools
    
    Returns:
    --------
    dict
        Dictionary with tool names as keys and DataFrames as values
    """
    
    tool_results = {}
    
    for file_path in input_files:
        # Extract tool name from file path
        tool_name = None
        if 'CHEUI' in file_path:
            tool_name = 'CHEUI'
        elif 'ELIGOS2' in file_path:
            tool_name = 'ELIGOS2'
        elif 'm6Anet' in file_path:
            tool_name = 'm6Anet'
        elif 'Nanocompore' in file_path:
            tool_name = 'Nanocompore'
        elif 'DENA' in file_path:
            tool_name = 'DENA'
        elif 'Epinano' in file_path:
            tool_name = 'Epinano'
        elif 'MINES' in file_path:
            tool_name = 'MINES'
        elif 'xPore' in file_path:
            tool_name = 'xPore'
        elif 'yanocomp' in file_path:
            tool_name = 'yanocomp'
        elif 'NanoSPA' in file_path:
            tool_name = 'NanoSPA'
        
        if tool_name:
            try:
                df = pd.read_csv(file_path, sep='\t')
                df['Tool'] = tool_name
                df['File'] = file_path
                tool_results[tool_name] = df
                print(f"Loaded {len(df)} modifications from {tool_name}")
            except Exception as e:
                print(f"Error reading {file_path}: {e}", file=sys.stderr)
    
    return tool_results

def create_modification_summary(tool_results, output_file):
    """
    Create summary statistics for all modifications.
    
    Parameters:
    -----------
    tool_results : dict
        Dictionary with tool results
    output_file : str
        Path to output summary file
    """
    
    summary_stats = []
    
    for tool_name, df in tool_results.items():
        stats = {
            'Tool': tool_name,
            'Total_Modifications': len(df),
            'Unique_Chromosomes': df['Chr'].nunique() if 'Chr' in df.columns else 'N/A',
            'Average_Probability': df['Prob'].mean() if 'Prob' in df.columns else 'N/A',
            'Average_Modification_Ratio': df['mod_ratio'].mean() if 'mod_ratio' in df.columns else 'N/A',
            'Strand_Distribution': df['Strand'].value_counts().to_dict() if 'Strand' in df.columns else 'N/A'
        }
        summary_stats.append(stats)
    
    # Create summary DataFrame
    summary_df = pd.DataFrame(summary_stats)
    
    # Save summary
    summary_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Modification summary saved to {output_file}")
    return summary_df

def create_tool_comparison(tool_results, output_file):
    """
    Create comparison of modifications across tools.
    
    Parameters:
    -----------
    tool_results : dict
        Dictionary with tool results
    output_file : str
        Path to output comparison file
    """
    
    # Find common chromosomes and positions across tools
    all_modifications = []
    
    for tool_name, df in tool_results.items():
        if 'Chr' in df.columns and 'Start' in df.columns:
            for _, row in df.iterrows():
                all_modifications.append({
                    'Chr': row['Chr'],
                    'Start': row['Start'],
                    'Tool': tool_name,
                    'Status': row.get('Status', 'Mod'),
                    'Prob': row.get('Prob', np.nan),
                    'Strand': row.get('Strand', '*')
                })
    
    # Create comparison DataFrame
    comparison_df = pd.DataFrame(all_modifications)
    
    if len(comparison_df) > 0:
        # Find overlapping modifications
        comparison_df['Position_Key'] = comparison_df['Chr'].astype(str) + ':' + comparison_df['Start'].astype(str)
        
        # Count tools per position
        position_counts = comparison_df.groupby('Position_Key')['Tool'].nunique().reset_index()
        position_counts.columns = ['Position_Key', 'Tool_Count']
        
        # Merge with original data
        comparison_df = comparison_df.merge(position_counts, on='Position_Key')
        
        # Sort by tool agreement (positions detected by more tools first)
        comparison_df = comparison_df.sort_values(['Tool_Count', 'Position_Key'], ascending=[False, True])
        
        # Save comparison
        comparison_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"Tool comparison saved to {output_file}")
        print(f"Total unique positions: {comparison_df['Position_Key'].nunique()}")
        print(f"Positions detected by multiple tools: {len(position_counts[position_counts['Tool_Count'] > 1])}")
    else:
        print("No modifications found for comparison")
        # Create empty file with headers
        pd.DataFrame(columns=['Chr', 'Start', 'Tool', 'Status', 'Prob', 'Strand', 'Position_Key', 'Tool_Count']).to_csv(output_file, sep='\t', index=False)
    
    return comparison_df

def create_visualizations(tool_results, output_dir):
    """
    Create visualization plots for the results.
    
    Parameters:
    -----------
    tool_results : dict
        Dictionary with tool results
    output_dir : str
        Directory to save plots
    """
    
    # Set up plotting style
    plt.style.use('seaborn-v0_8')
    sns.set_palette("husl")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Bar plot of modification counts per tool
    fig, ax = plt.subplots(figsize=(12, 8))
    tool_counts = {tool: len(df) for tool, df in tool_results.items()}
    
    bars = ax.bar(tool_counts.keys(), tool_counts.values())
    ax.set_xlabel('Tool')
    ax.set_ylabel('Number of Modifications')
    ax.set_title('RNA Modifications Detected by Each Tool')
    ax.tick_params(axis='x', rotation=45)
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                f'{int(height)}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'modification_counts.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Venn diagram of tool overlaps (if applicable)
    if len(tool_results) >= 2:
        from matplotlib_venn import venn2, venn3
        
        # Get positions for each tool
        tool_positions = {}
        for tool_name, df in tool_results.items():
            if 'Chr' in df.columns and 'Start' in df.columns:
                positions = set(df['Chr'].astype(str) + ':' + df['Start'].astype(str))
                tool_positions[tool_name] = positions
        
        # Create Venn diagrams for pairs of tools
        tools = list(tool_positions.keys())
        for i in range(len(tools)):
            for j in range(i+1, len(tools)):
                if len(tool_positions[tools[i]]) > 0 and len(tool_positions[tools[j]]) > 0:
                    fig, ax = plt.subplots(figsize=(10, 8))
                    venn2([tool_positions[tools[i]], tool_positions[tools[j]]], 
                          set_labels=(tools[i], tools[j]), ax=ax)
                    ax.set_title(f'Overlap between {tools[i]} and {tools[j]}')
                    plt.tight_layout()
                    plt.savefig(os.path.join(output_dir, f'venn_{tools[i]}_vs_{tools[j]}.png'), 
                               dpi=300, bbox_inches='tight')
                    plt.close()
    
    # 3. Chromosome distribution
    fig, ax = plt.subplots(figsize=(15, 8))
    chromosome_data = []
    
    for tool_name, df in tool_results.items():
        if 'Chr' in df.columns:
            chr_counts = df['Chr'].value_counts()
            for chr_name, count in chr_counts.items():
                chromosome_data.append({
                    'Tool': tool_name,
                    'Chromosome': str(chr_name),
                    'Count': count
                })
    
    if chromosome_data:
        chr_df = pd.DataFrame(chromosome_data)
        chr_pivot = chr_df.pivot(index='Chromosome', columns='Tool', values='Count').fillna(0)
        
        chr_pivot.plot(kind='bar', ax=ax, width=0.8)
        ax.set_xlabel('Chromosome')
        ax.set_ylabel('Number of Modifications')
        ax.set_title('Modification Distribution Across Chromosomes')
        ax.tick_params(axis='x', rotation=45)
        ax.legend(title='Tool', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'chromosome_distribution.png'), dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"Visualizations saved to {output_dir}")

def main(input_files, summary_file, comparison_file, plot_dir):
    """
    Main function to generate all summaries and comparisons.
    
    Parameters:
    -----------
    input_files : list
        List of processed result files
    summary_file : str
        Path to output summary file
    comparison_file : str
        Path to output comparison file
    plot_dir : str
        Directory to save plots
    """
    
    print(f"Processing {len(input_files)} result files...")
    
    # Read all tool results
    tool_results = read_tool_results(input_files)
    
    if not tool_results:
        print("No valid results found. Exiting.")
        return
    
    # Create summary statistics
    summary_df = create_modification_summary(tool_results, summary_file)
    print("\nSummary Statistics:")
    print(summary_df.to_string(index=False))
    
    # Create tool comparison
    comparison_df = create_tool_comparison(tool_results, comparison_file)
    
    # Create visualizations
    create_visualizations(tool_results, plot_dir)
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python generate_summary.py <input_files> <summary_output> <comparison_output> [plot_dir]")
        print("Example: python generate_summary.py \"results/*/*_processed.txt\" summary.tsv comparison.tsv plots")
        sys.exit(1)
    
    # Handle file patterns
    input_pattern = sys.argv[1]
    summary_file = sys.argv[2]
    comparison_file = sys.argv[3]
    plot_dir = sys.argv[4] if len(sys.argv) > 4 else "plots"
    
    # Expand file pattern
    input_files = glob.glob(input_pattern)
    
    if not input_files:
        print(f"No files found matching pattern: {input_pattern}")
        sys.exit(1)
    
    main(input_files, summary_file, comparison_file, plot_dir)