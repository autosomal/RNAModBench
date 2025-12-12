#!/usr/bin/env python3
"""
Setup script for RNAModBench pipeline.
This script helps with initial setup and configuration.
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path

def check_dependencies():
    """Check if required dependencies are installed."""
    
    print("Checking dependencies...")
    
    dependencies = {
        'conda': 'conda --version',
        'snakemake': 'snakemake --version',
        'python': 'python --version'
    }
    
    missing_deps = []
    
    for dep, cmd in dependencies.items():
        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            print(f"✓ {dep} is installed")
        except subprocess.CalledProcessError:
            print(f"✗ {dep} is not installed")
            missing_deps.append(dep)
    
    return missing_deps

def create_directories():
    """Create necessary directory structure."""
    
    print("Creating directory structure...")
    
    directories = [
        'data',
        'reference',
        'results',
        'results/summary',
        'results/report',
        'results/logs',
        'resources'
    ]
    
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        print(f"✓ Created {directory}")

def setup_conda_envs():
    """Set up conda environments for all tools."""
    
    print("Setting up conda environments...")
    
    env_files = [
        'envs/cheui.yaml',
        'envs/m6anet.yaml',
        'envs/nanocompore.yaml',
        'envs/eligos2.yaml',
        'envs/epinano.yaml',
        'envs/dena.yaml',
        'envs/tombo.yaml',
        'envs/mines.yaml',
        'envs/xpore.yaml',
        'envs/yanocomp.yaml',
        'envs/nanospa.yaml'
    ]
    
    for env_file in env_files:
        if os.path.exists(env_file):
            env_name = os.path.splitext(os.path.basename(env_file))[0]
            print(f"Creating environment: {env_name}")
            try:
                subprocess.run(f"conda env create -f {env_file}", shell=True, check=True)
                print(f"✓ Created {env_name} environment")
            except subprocess.CalledProcessError:
                print(f"✗ Failed to create {env_name} environment")
        else:
            print(f"✗ Environment file not found: {env_file}")

def download_test_data():
    """Download test data for pipeline validation."""
    
    print("Downloading test data...")
    
    # This is a placeholder - in a real implementation, you would
    # download actual test data from a public repository
    
    test_data_dir = "data/test"
    os.makedirs(test_data_dir, exist_ok=True)
    
    # Create placeholder files
    placeholder_files = [
        "README.md",
        "sample_list.txt"
    ]
    
    for file in placeholder_files:
        with open(os.path.join(test_data_dir, file), 'w') as f:
            f.write("# Test data placeholder\n")
            f.write("Replace with actual fast5 files for testing\n")
    
    print("✓ Created test data directory structure")
    print("  Please replace placeholder files with actual fast5 data")

def validate_config():
    """Validate the configuration file."""
    
    print("Validating configuration...")
    
    config_file = "config/config.yaml"
    
    if not os.path.exists(config_file):
        print(f"✗ Configuration file not found: {config_file}")
        return False
    
    # Basic validation - in a real implementation, you would
    # use a YAML parser to validate the configuration
    
    print(f"✓ Configuration file exists: {config_file}")
    print("  Please review and edit the configuration file as needed")
    return True

def main():
    """Main setup function."""
    
    parser = argparse.ArgumentParser(description='Setup RNAModBench pipeline')
    parser.add_argument('--skip-conda', action='store_true', 
                       help='Skip conda environment setup')
    parser.add_argument('--skip-test-data', action='store_true',
                       help='Skip test data download')
    parser.add_argument('--full-setup', action='store_true',
                       help='Run full setup including conda environments')
    
    args = parser.parse_args()
    
    print("=== RNAModBench Pipeline Setup ===\n")
    
    # Check dependencies
    missing_deps = check_dependencies()
    if missing_deps:
        print(f"\nMissing dependencies: {', '.join(missing_deps)}")
        print("Please install missing dependencies before continuing.")
        sys.exit(1)
    
    # Create directories
    create_directories()
    
    # Setup conda environments
    if args.full_setup or not args.skip_conda:
        setup_conda_envs()
    
    # Download test data
    if not args.skip_test_data:
        download_test_data()
    
    # Validate configuration
    validate_config()
    
    print("\n=== Setup Complete ===")
    print("\nNext steps:")
    print("1. Edit config/config.yaml with your sample names and settings")
    print("2. Place your fast5 files in the data/ directory")
    print("3. Download reference genome and annotation files")
    print("4. Run the pipeline: snakemake --use-conda --cores 40")
    print("\nFor more information, see README.md")

if __name__ == "__main__":
    main()