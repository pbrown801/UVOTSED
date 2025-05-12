import os
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.interpolate import interp1d
from utils import (
    kfun,
    get_experiment_param_grid,
    lazy_parameter_grid,
    blanket_spectrum,
    apply_reddening,
    apply_redshift,
    get_filters,
    apply_foleymodel
)


#def generate_spectrum(dmb, f11_scale=1.0):
#    return f11_scale * f_11_fun(slambda_wave) + slambda_fun(slambda_wave) * (dmb - 1.1)
#def generate_f11spectrum(dmb, spec_wave, spec_flux):
#    return f11_scale * f_11_fun(spec_wave) + slambda_fun(spec_wave) * (dmb - 1.1)

def full_model(config: dict, experiment_name: str = "full_model"):
    # Get parameter grid configuration for the experiment.
    param_grid_config = get_experiment_param_grid(config, experiment_name)
    
    # Process spectra files.
    input_config = config.get("input", {})
    spectra_files = input_config.get("fewspectra_files", [])
    if isinstance(spectra_files, dict):
        base_path = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.abspath(os.path.join(base_path, ".."))
        spectra_dir = os.path.join(project_root, spectra_files.get("path", ""))
        spectra_files = [os.path.join(spectra_dir, fname) for fname in spectra_files.get("filenames", [])]
    
    # Get output directory.
    output_config = config.get("output", {})
    if "output_dir" not in output_config:
        raise ValueError("Output directory not specified in the configuration.")
    output_dir = output_config["output_dir"]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    filters = get_filters(
        path=input_config.get("filter_files", {}).get("path", ""),
        filenames=input_config.get("filter_files", {}).get("filenames", [])
    )
    
    results = []
    total_to_process = len(spectra_files) * len(param_grid_config)
    print(f"Total combinations to process: {total_to_process}")
    for spectrum_file in tqdm(spectra_files, desc="Processing spectra files"):
        for params in tqdm(lazy_parameter_grid(param_grid_config), desc="Processing parameter grid", leave=False):
            lb = params.get("lb")
            av26 = params.get("av26")
            av31 = params.get("av31")
            z = params.get("z")
            dmb = params.get("dmb")
            f11_param = params.get("f11")
        
        
            wavelength, _, flux_modified = blanket_spectrum(spectrum_file, lb, uv_cutoff=4000, alpha=0.2)
            if wavelength is None:
                continue

            sort_idx = np.argsort(wavelength)
            wavelength = wavelength[sort_idx]
            flux_modified = flux_modified[sort_idx]


            flux_dmb=apply_foleymodel(wavelength, flux_modified, dmb)

            mwreddened_flux = apply_reddening(wavelength, flux_dmb, av31, 3.1)
            hostreddened_flux = apply_reddening(wavelength, mwreddened_flux, av31, 2.6)
            wave_redshifted, flux_redshifted = apply_redshift(wavelength, hostreddened_flux, z)
            final_flux = flux_redshifted * f11_param * (dmb / 1.1)

            
            mags = kfun(wave_redshifted, final_flux, filters)
            
            spectrum_name = os.path.splitext(os.path.basename(spectrum_file))[0]
            result = {
                "Spectrum": spectrum_name,
                "lb": round(lb, 3),
                "av26": round(av26,3 ),
                "av31": round(av31, 3),
                "redshift": round(z, 3),
                "dmb": round(dmb, 3),
            }
            for key, mag in mags.items():
                result[key] = round(mag, 6)
            results.append(result)
    
    results_df = pd.DataFrame(results)
    summary_path = os.path.join(output_dir, "summary.csv")
    results_df.to_csv(summary_path, index=False)
    print(f"Model results saved to {summary_path}")
    return results_df
