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
    get_filters
)

uv_model = pd.read_csv('../spectra/uvmodel.data', sep='\\s+', comment='#')
uv_model.columns = ['wavelength', 'f_11', 'Slambda']
uv_model_wave = np.asarray(uv_model.wavelength).astype(float)
uv_model_flux = np.asarray(uv_model.f_11)
uv_model_dif  = np.asarray(uv_model.Slambda)


uv_model_wave=np.insert(uv_model_wave,0,1400, axis=None)
uv_model_wave=np.append(uv_model_wave,8000, axis=None)

uv_model_flux=np.insert(uv_model_flux,0,uv_model_flux[0], axis=None)
uv_model_flux=np.append(uv_model_flux,uv_model_flux[len(uv_model_flux)-1], axis=None)

uv_model_dif=np.insert(uv_model_dif,0,uv_model_dif[0], axis=None)
uv_model_dif=np.append(uv_model_dif,uv_model_dif[len(uv_model_dif)-1], axis=None)


f_11_fun = interp1d(uv_model_wave,uv_model_flux, kind='cubic')
slambda_fun = interp1d(uv_model_wave, uv_model_dif, kind='cubic')

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
            rv = params.get("rv")
            av = params.get("av")
            z = params.get("z")
            dmb = params.get("dmb")
            f11_param = params.get("f11")
        
        
            wavelength, _, flux_modified = blanket_spectrum(spectrum_file, lb, uv_cutoff=4000, alpha=0.2)
            if wavelength is None:
                continue

            sort_idx = np.argsort(wavelength)
            wavelength = wavelength[sort_idx]
            flux_modified = flux_modified[sort_idx]

            if rv == 1.6:
                av16 = av
                av31 = None
            else:
                av16 = None
                av31 = av

            reddened_flux = apply_reddening(wavelength, flux_modified, av, rv)
            wave_redshifted, flux_redshifted = apply_redshift(wavelength, reddened_flux, z)
            final_flux = flux_redshifted * f11_param * (dmb / 1.1)

            
            mags = kfun(wave_redshifted, final_flux, filters)
            
            spectrum_name = os.path.splitext(os.path.basename(spectrum_file))[0]
            result = {
                "Spectrum": spectrum_name,
                "lb": round(lb, 3),
                "rv": round(rv, 3),
                "av16": round(av16,3 ) if av16 is not None else None,
                "av31": round(av31, 3) if av31 is not None else None,
                "redshift": round(z, 3),
                "dmb": round(dmb, 3),
                "f11": round(f11_param, 3)
            }
            for key, mag in mags.items():
                result[key] = round(mag, 6)
            results.append(result)
    
    results_df = pd.DataFrame(results)
    summary_path = os.path.join(output_dir, "summary.csv")
    results_df.to_csv(summary_path, index=False)
    print(f"Model results saved to {summary_path}")
    return results_df
