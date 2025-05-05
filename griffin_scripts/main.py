"""
Entry point
"""

import os
import argparse
import json
import sys

from model_generation import full_model
from spectrum_processing import compare_and_remake

def main():
    # Setting up command-line argument parsing
    parser = argparse.ArgumentParser(description="Supernova Model Generation")
    parser.add_argument("--config", type=str, default="config.json", help="Path to the configuration file")
    args = parser.parse_args()

    # Load configuration file
    try:
        base_path = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(base_path, args.config)
        with open(config_path, 'r') as f:
            config = json.load(f)
    except FileNotFoundError:
        print(f"Config file not found: {args.config}")
        sys.exit(1)

    # Validate configuration format (TODO: Add validation logic here)

    # Process the configuration for spectra and filter files
    project_root = os.path.abspath(os.path.join(base_path, ".."))
    spectra_files = config.get("input", {}).get("spectra_files")
    if isinstance(spectra_files, dict) and "path" in spectra_files and "filenames" in spectra_files:
        spectra_dir = os.path.join(project_root, spectra_files["path"])
        config["input"]["spectra_files"] = [os.path.join(spectra_dir, fname) for fname in spectra_files["filenames"]]
    
    filter_files = config.get("input", {}).get("filter_files")
    if isinstance(filter_files, dict) and "path" in filter_files and "filenames" in filter_files:
        filter_dir = os.path.join(project_root, filter_files["path"])
        config["input"]["filter_files"]["filenames"] = [os.path.join(filter_dir, fname) for fname in filter_files["filenames"]]

    # Set up model dictionary. The key should match the experiment's name.
    model_dict = {
        "full_model": full_model,
        "compare_and_remake": compare_and_remake
    }

    # Run the model experiments using the experiment name defined in each experiment block.
    experiments = config.get("experiments", [])
    skip = []
    for idx, exp in enumerate(experiments):
        # Get list of experiments to skip if any
        if idx == 0:
            skip = exp.get("skip", [])
            continue

        # Run or skip experiement
        experiment_name = exp.get("name")
        if experiment_name in skip:
            continue
        try:
            model_dict[experiment_name](config, experiment_name=experiment_name)
        except KeyError:
            print(f"Model '{experiment_name}' not found in model dictionary.")
            continue
        except Exception as e:
            print(f"Error running model '{experiment_name}': {e}")
            continue

if __name__ == '__main__':
    main()
