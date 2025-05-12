"""

"""

from pyphot import Filter
from pyphot import unit
import numpy as np
import pandas as pd
import time
import os
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import seaborn as sns
import extinction
from scipy.interpolate import interp1d


def resolve_parameter(param_config: dict, global_params: dict = None) -> dict:
    """
    Resolves a parameter definition. It supports three formats:
    
    1. Reference: { "ref": "some_key" }
    2. Inline explicit list: { "values": [...] }
    3. Inline range: { "start": X, "end": Y, "steps": Z }
    """
    # Assume inline parameter definition
    if not (isinstance(param_config, dict) and "ref" in param_config):
        return param_config
    
    # Resolve reference parameter
    if global_params is None:
        raise ValueError("Global parameters must be provided for reference resolution.")
    ref_key = param_config["ref"]
    if ref_key not in global_params:
        raise ValueError(f"Reference parameter '{ref_key}' not found in global parameters.")
    return global_params[ref_key]



def get_experiment_param_grid(config: dict, experiment_name: str) -> dict:
    """
    Retrieves the parameter grid for a specified experiment from the configuration.

    Example input:

    """
    experiments = config.get("experiments", {})
    global_params = config.get("model_parameters", {})

    for exp in experiments:
        if exp.get("name") != experiment_name:
            continue
        param_grid = exp.get("parameter_grid", {})
        resolved_grid = {}
        for param_name, param_config in param_grid.items():
            resolved_grid[param_name] = resolve_parameter(param_config, global_params)
        return resolved_grid
    
    raise ValueError(f"Experiment '{experiment_name}' not found in configuration.")





def generate_values_lazy(param_config):
    """
    Yields values for a parameter defined by either an explicit "values" list,
    or by a range with "start", "end", and "steps".
    """
    if "values" in param_config:
        for v in param_config["values"]:
            yield v
    elif all(key in param_config for key in ["start", "end", "steps"]):
        start = param_config["start"]
        end = param_config["end"]
        steps = param_config["steps"]
        if steps < 1:
            raise ValueError("Number of steps must be at least 1.")
        if steps == 1:
            yield start
        else:
            step = (end - start) / (steps - 1)
            for i in range(steps):
                yield start + i * step
    else:
        raise ValueError("Parameter must have 'values' or 'start', 'end', and 'steps' defined.")

def lazy_parameter_grid(param_dict):
    """
    Yields all combinations of parameters defined in param_dict.

    Example input:
    ```
    {
        "param1": {"values": [1, 2, 3]},
        "param2": {"start": 0, "end": 10, "steps": 5},
        "param3": {"values": [0.1, 0.2]}
    }
    ```
    """
    keys = list(param_dict.keys())
    def recursive_grid(index, current):
        if index == len(keys):
            yield current.copy()
        else:
            current_key = keys[index]
            for value in generate_values_lazy(param_dict[current_key]):
                current[current_key] = value
                yield from recursive_grid(index + 1, current)
    yield from recursive_grid(0, {})


'''
slambda_data = pd.read_csv('../spectra/uvmodel.data', sep='\\s+', comment='#')
slambda_data.columns = ['wavelength', 'f_11', 'Slambda']
slambda_wave = np.asarray(slambda_data.wavelength).astype(float)
print("when does the slambda in utils run?")
f_11 = np.asarray(slambda_data.f_11)
f_11_fun = interp1d(slambda_wave, f_11, kind='cubic')
slambda_fun = interp1d(slambda_wave, slambda_data.Slambda, kind='cubic')
'''


#def generate_spectrum(dmb: float, use_base: bool = False):
#    """Generate a spectrum with a given DMB."""
#    if use_base:
#        return f_11_fun(slambda_wave)
#
#    spectrum_flux = f_11_fun(slambda_wave) + slambda_fun(slambda_wave) * (dmb - 1.1)
#    return spectrum_flux



def get_filters(path: str, filenames: list):
    path = os.path.abspath(path)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Path {path} does not exist.")
    
    filters = {}
    for name in filenames:
        if not os.path.exists(os.path.join(path, name)):
            print(f"Filter file {name} not found in {path}.")
            continue
        try:
            data = pd.read_csv(os.path.join(path, name), sep="\\s+", comment="#")
            data.columns = ["Wavelength", "Area"]
            filter = Filter(
                data["Wavelength"].values * unit["Angstrom"],
                data["Area"].values,
                name=name.split("_")[0].split('/')[-1],
                dtype="photon",
                unit="Angstrom",
            )
            filter.transmission = data["Area"].values
            filters[name.split("_")[0].split('/')[-1]] = filter
        except Exception as e:
            print(f"Error reading filter file {name}: {e}")
            continue
    return filters


def setup_filters(
    filenames: list = None,
    path: str = "",
    filter_names: list = None,
    columns: list = None,
) -> dict:
    """Read filter files and return a dictionary of Filter objects."""
    filter_dict = {}
    if not os.path.exists(path):
        return filter_dict
    if path != "" and not path.endswith("/"):
        path += "/"
    for fname, name in zip(filenames, filter_names):
        data = pd.read_csv(f"{path}{fname}", sep="\\s+", comment="#")
        data.columns = columns
        filt = Filter(
            data["Wavelength"].values * unit["Angstrom"],
            data["Area"].values,
            name=name,
            dtype="photon",
            unit="Angstrom",
        )
        filt.transmission = data["Area"].values
        filter_dict[name] = filt
    return filter_dict



def kfun(wave, flux, filters, subtractor_key=None) -> dict:
    """Calculate the k-correction for a given spectrum and filters."""
    my_dict = {}
    for key, filter_ in filters.items():
        f = filter_.get_flux(wave * unit["Angstrom"], flux)
        my_dict[key] = -2.5 * np.log10(f) - filter_.Vega_zero_mag
    if subtractor_key is not None:
        for key in my_dict:
            if key != subtractor_key:
                my_dict[key] -= my_dict[subtractor_key]
    return my_dict


# Rename metallicity to line blanketing (lb)
def blanket_spectrum(spectrum_file: str, lb: float, uv_cutoff: float = 4000, alpha: float = 0.2):
    try:
        data = pd.read_csv(spectrum_file, sep='\\s+', comment='#', header=None)
    except Exception as ex:
        print(f"Error reading {spectrum_file}: {ex}")
        return None, None, None

    if data.shape[1] < 2:
        print(f"File {spectrum_file} does not have enough columns.")
        return None, None, None

    if data.shape[1] == 2:
        data.columns = ['wavelength', 'flux']
    else:
        data.columns = ['wavelength', 'flux', 'ferr']

    data['wavelength'] = pd.to_numeric(data['wavelength'], errors='coerce')
    data['flux'] = pd.to_numeric(data['flux'], errors='coerce')

    data.dropna(subset=['wavelength', 'flux'], inplace=True)

    wavelength = np.array(data['wavelength'], dtype=float)
    flux = np.array(data['flux'], dtype=float)

    flux_original = flux / np.max(flux)
    flux_modified = flux_original.copy()

    factor = max(1 - alpha * (lb - 1), 0.1)
    uv_mask = wavelength < uv_cutoff
    flux_modified[uv_mask] *= factor

    return wavelength, flux_original, flux_modified




def apply_reddening(wavelength: np.ndarray, flux: np.ndarray, AV: float, RV: float):
    """
        AV: Total V-band extinction (use negative for dereddening)
        RV: Ratio of total to selective extinction
    """
    extinction_curve = extinction.ccm89(wavelength, AV, RV)
    reddened_flux = extinction.apply(extinction_curve, flux)
    return reddened_flux

def apply_redshift(wavelength: np.ndarray, flux: np.ndarray, z: float):
    if z < 0:
        raise ValueError("Redshift must be non-negative.")
    wavelength_observed = wavelength * (1 + z)
    flux_observed = flux / (1 + z)
    return wavelength_observed, flux_observed
