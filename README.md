# RNAModBench - Comprehensive Nanopore RNA Modification Detection Pipeline

## Overview

RNAModBench is a comprehensive, cross-species pipeline for detecting RNA modifications from nanopore direct RNA sequencing data. This pipeline integrates multiple state-of-the-art tools and provides a standardized workflow for comparative analysis.

## Features

- **Multi-tool Integration**: Supports 11 different RNA modification detection tools
- **Cross-species Compatibility**: Works with any species with available reference genome
- **Standardized Output**: Converts all tool outputs to unified BED-like format
- **Comparative Analysis**: Generates tool overlap statistics and consensus calls
- **Quality Control**: Includes comprehensive QC metrics and filtering
- **Scalable**: Supports high-throughput analysis with parallel processing
- **Modular Design**: Run individual tools or complete pipeline

## Installation

### Prerequisites

- Python 3.8+
- Conda or Miniconda
- Snakemake
- Git

### Quick Setup

```bash
# Clone the repository
git clone https://github.com/autosomal/RNAModBench.git
cd RNAModBench

# Create conda environments
conda env create -f envs/cheui.yaml
conda env create -f envs/m6anet.yaml
conda env create -f envs/nanocompore.yaml
# ... create other environments as needed

# Install pipeline dependencies
conda install -c conda-forge -c bioconda snakemake minimap2 samtools nanopolish
```

### Reference Data Setup

```bash
# Prepare reference files
mkdir -p reference
cd reference

# Download your reference genome and annotation
# Example for human:
wget ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa

# Download gene annotation
wget ftp://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
gunzip Homo_sapiens.GRCh38.112.gtf.gz
mv Homo_sapiens.GRCh38.112.gtf genes.gtf

# Create transcriptome (this will be done automatically by the pipeline)
```

## Usage

### Basic Usage

```bash
# Edit configuration file
nano config/config.yaml

# Run complete pipeline
snakemake --use-conda --cores 40

# Run specific tools only
snakemake --use-conda --cores 40 results/CHEUI/sample1/sample1_CHEUI_processed.txt

# Dry run to check workflow
snakemake --dry-run
```

### Configuration

Edit `config/config.yaml` to specify:

- Sample names
- Tools to run
- Reference paths
- Filtering thresholds
- Computational resources

### Sample Configuration

```yaml
# Sample configuration
samples:
  - sample1
  - sample2
  - control1

# Tools to run
tools:
  - CHEUI
  - m6Anet
  - Nanocompore
  - ELIGOS2

# Reference paths
reference_dir: "reference"
data_dir: "data"
results_dir: "results"

# Filtering thresholds
cheui:
  prob_threshold: 0.999
  ratio_threshold: 0.1

m6anet:
  prob_threshold: 0.5
  ratio_threshold: 0.1
```

## Input Data Structure

```
data/
├── sample1/
│   └── fast5/
│       ├── read1.fast5
│       ├── read2.fast5
│       └── ...
├── sample2/
│   └── fast5/
└── control1/
        └── fast5/
```

## Output Structure

```
results/
├── basecalling/
├── alignment/
├── nanopolish/
├── CHEUI/
│   └── sample1/
│       └── sample1_CHEUI_processed.txt
├── m6Anet/
├── Nanocompore/
├── summary/
│   ├── modification_summary.tsv
│   └── tool_comparison.tsv
└── report/
    └── RNAModBench_report.html
```

## Output Format

All tools produce standardized BED-like output with the following columns:

1. **Chr**: Chromosome or transcript ID
2. **Start**: Start position (0-based)
3. **End**: End position
4. **Status**: Modification status (Mod/Unmod)
5. **Prob**: Probability or score
6. **Strand**: Strand information (+, -, or *)
7. **mod_ratio**: Modification ratio/stoichiometry

## Quality Control

The pipeline includes several QC steps:

- Basecalling quality filtering
- Read alignment quality control
- Coverage-based filtering
- Statistical significance testing
- Tool-specific quality metrics

## Performance

- **Parallel Processing**: Utilizes multiple cores for speed
- **Memory Efficient**: Optimized for large datasets
- **GPU Support**: Accelerated computation for deep learning tools
- **Modular Execution**: Run only required tools



## Troubleshooting

### Common Issues

1. **Memory Issues**: Reduce thread counts in config.yaml
2. **Missing Dependencies**: Ensure all conda environments are created
3. **Reference Issues**: Check genome and annotation file formats
4. **Tool Errors**: Check individual tool documentation

### Debug Mode

```bash
# Run with verbose output
snakemake --use-conda --cores 40 --verbose

# Run specific rule with debug
snakemake --use-conda --cores 1 results/CHEUI/sample1/sample1_CHEUI_processed.txt --debug
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## License

This pipeline is released under the MIT License. See LICENSE file for details.

## Support

For issues and questions:

1. Check existing GitHub issues
2. Create a new issue with detailed description
3. Include error logs and configuration files

## Acknowledgments

This pipeline integrates tools developed by multiple research groups. We acknowledge the developers of all integrated tools for their contributions to the field.

## References

1. CHEUI: [Publication details]
2. m6Anet: [Publication details]
3. Nanocompore: [Publication details]
4. ELIGOS2: [Publication details]
5. Epinano: [Publication details]
6. MINES: [Publication details]
7. xPore: [Publication details]
8. yanocomp: [Publication details]
9. NanoSPA: [Publication details]
10. DENA: [Publication details]
11. DRUMMER: [Publication details]

---

For more information, please visit the [GitHub repository](https://github.com/autosomal/RNAModBench).