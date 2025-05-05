import numpy as np
import pandas as pd
import os
from scipy.interpolate import interp1d
from pyphot import Filter
import matplotlib.pyplot as plt

from utils import apply_reddening, apply_redshift, get_filters, kfun

def process_spectrum(spectrum_file: str, lb: float, uv_cutoff: float = 4000, alpha: float = 0.2):
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



def compute_color(m1, m2):
    """Compute the color between two magnitudes."""
    return m1 - m2

def compute_color_error(err1, err2):
    """Compute the error in color."""
    return np.sqrt(err1**2 + err2**2)

def filter_models_by_photometry(models_df, observed, errors):
    # Color
    obs_color1 = compute_color(observed["UVW2"], observed["V"])
    obs_color2 = compute_color(observed["U"], observed["B"])
    obs_color3 = compute_color(observed["B"], observed["V"])
    obs_color4 = compute_color(observed["UVW2"], observed["UVM2"])
    obs_color5 = compute_color(observed["UVM2"], observed["UVW1"])

    # Tolerance
    tol_color1 = compute_color_error(errors["UVW2"], errors["V"])
    tol_color2 = compute_color_error(errors["U"], errors["B"])
    tol_color3 = compute_color_error(errors["B"], errors["V"])
    tol_color4 = compute_color_error(errors["UVW2"], errors["UVM2"])
    tol_color5 = compute_color_error(errors["UVM2"], errors["UVW1"])

    rows = []
    floor = 0.1



    # w1mv = []
    # bmv = []




    for idx, row in models_df.iterrows():

        # bv = compute_color(row["B"], row["V"])
        # w1v = compute_color(row["UVW1"], row["V"])


        # bmv.append(bv)
        # w1mv.append(w1v)
        # continue


        m_color1 = compute_color(row["UVW2"], row["V"])
        if abs(m_color1 - obs_color1) > max(tol_color1, floor) and not np.isnan(obs_color1):
            # print(abs(m_color1 - obs_color1), tol_color1)
            # print("1")
            continue
        m_color2 = compute_color(row["U"], row["B"])
        if abs(m_color2 - obs_color2) > max(tol_color2, floor)  and not np.isnan(obs_color2):
            # print("2")
            continue
        m_color3 = compute_color(row["B"], row["V"])
        if abs(m_color3 - obs_color3) > max(tol_color3, floor)  and not np.isnan(obs_color3):
            # print("3")
            continue
        m_color4 = compute_color(row["UVW2"], row["UVM2"])
        if abs(m_color4 - obs_color4) > max(tol_color4, floor)  and not np.isnan(obs_color4):
            # print("4")
            continue
        m_color5 = compute_color(row["UVM2"], row["UVW1"])
        if abs(m_color5 - obs_color5) > max(tol_color5, floor)  and not np.isnan(obs_color5):
            # print("5")
            continue

        # Within tolerance
        rows.append(row)
    
    
    # plt.scatter(np.array(bmv), np.array(w1mv))
    # plt.xlabel('B-V')
    # plt.ylabel('UVW1-V')
    # plt.title("TEst")
    # plt.show
    # plt.savefig("src/")

    # print(f"Filtered models: {len(rows)}")

    if not rows:
        # print("No models within tolerance.")
        return pd.DataFrame()
    
    filtered_df = pd.DataFrame(rows)
    return filtered_df

"""{
    "input": {
        "spectra_files": {
            "path": "spectra",
            "filenames": [
                "SN1992A_UV.dat",
                "SN2011by_peak_11fe_appended.dat"
            ]
        },
        "filter_files": {
            "path": "filters",
            "filenames": [
                "U_P08.txt",
                "B_P08.txt",
                "V_P08.txt",
                "UVW2_B11.txt",
                "UVM2_B11.txt",
                "UVW1_B11.txt"
            ]
        }
    },"""
def remake_spectrum(row, config):
    spectrum_basename = row["Spectrum"]
    # spectrum_basename: SN1992A_UV (example)

    # config:
    """{'input': {'spectra_files': ['/Users/griffinbeaudreau/Desktop/School/UVOTSED/spectra/SN1992A_UV.dat', '/Users/griffinbeaudreau/Desktop/School/UVOTSED/spectra/SN2011by_peak_11fe_appended.dat'], """

    # match basename to spectrum_file
    files = config["input"]["spectra_files"]
    for file in files:
        if spectrum_basename in file:
            spectrum_file = file
            break
    else:
        print(f"File not found for spectrum: {spectrum_basename}")
        return None
    
    # print(f"Found spectrum file: {spectrum_file}")


    lb = row["lb"]
    rv = row["rv"]
    av = row["av16"] if rv == 1.6 else row["av31"]
    redshift = row["redshift"]
    dmb = row["dmb"]
    f11 = row["f11"]

    # U,B,V,UVW2,UVM2,UVW1
    magnitudes = {
        "U": row["U"],
        "B": row["B"],
        "V": row["V"],
        "UVW2": row["UVW2"],
        "UVM2": row["UVM2"],
        "UVW1": row["UVW1"]
    }

    # print(f"Processing spectrum: {spectrum_file}")

    """SN2011by_peak_11fe_appended,2.0,3.1,,0.5,0.0,1.8,0.8,-19.728388,-19.60335,-19.805873,-16.472109,-15.396212,-18.16598

    
        This done by next week, and then next monday work on trasferring it                                     !!
        Then have a week where I'm playing with it (test using command line, make sure python versions match)   !!
        One page write up summarizing what I did and email peter                                                !!


        1. From row, take magnitudes
        2. De-redden (apply the reddening function using the MW_ebv with a negative sign)
        3. Divide wavelength by (1+redshift (the redshift specified in the experiemnt config))
        4. Multiply the flux by (1 + redshift)
        5. Apply reddening function with negative host_ebv (host galaxy redenning)
        6. Calculate magnitudes then subtract from row to get correction factor
        7. Instead of saving, add row of these corrections to a list (don't need to save as csv files yay) (use current files and plot wavelength, flux/flux@wave=6000)
        8. After, calcualte mean correction and standard deviation (final product of this project)
    """

    wavelength, _, flux_modified = process_spectrum(spectrum_file, lb, uv_cutoff=4000, alpha=0.2)
    if wavelength is None:
        print(f"Failed to process spectrum: {spectrum_file}")
        return None
    
    sort_idx = np.argsort(wavelength)
    wavelength = wavelength[sort_idx]
    flux_modified = flux_modified[sort_idx]

    reddened_flux = apply_reddening(wavelength, flux_modified, av, rv)
    wave_redshifted, flux_redshifted = apply_redshift(wavelength, reddened_flux, redshift)
    final_flux = flux_redshifted * f11 * (dmb / 1.1)

    return wavelength, final_flux






def compare_and_remake(config, experiment_name=None):
    import os
    import pandas as pd
    import numpy as np  # for NaN checks and imputation
    from tqdm import tqdm

    print(f"DEBUG: Starting compare_and_remake for experiment '{experiment_name}'")
    exp_cfg = next((e for e in config.get("experiments", []) if e.get("name") == experiment_name), None)
    if exp_cfg is None:
        raise ValueError(f"No experiment named '{experiment_name}' in config")

    # Load data
    data_csv_path = exp_cfg.get("data_csv")
    if not os.path.exists(data_csv_path):
        print(f"Data CSV file not found: {data_csv_path}")
        return
    data_df = pd.read_csv(data_csv_path)
    if data_df.empty:
        print(f"Data CSV file is empty: {data_csv_path}")
        return

    # Load and clean model summary
    target_csv = exp_cfg.get("target_csv")
    models_df = pd.read_csv(target_csv)
    numeric_cols = ['U','B','V','UVW2','UVM2','UVW1','lb']
    for col in numeric_cols:
        if col in models_df:
            models_df[col] = pd.to_numeric(models_df[col], errors='coerce')
    models_df.dropna(subset=numeric_cols, inplace=True)

    # Load filters
    filters = get_filters(
        path=config["input"]["filter_files"]["path"],
        filenames=config["input"]["filter_files"]["filenames"]
    )

    # Column fallback
    cols = [c.strip() for c in data_df.columns]
    col_idx = {col: f"_{i+1}" for i, col in enumerate(cols)}

    # Photometric mapping
    phot_bands = [
        ("B_Mmax","B"),("V_MBmax","V"),("U_MBmax","U"),
        ("W2_MBmax","UVW2"),("M2_MBmax","UVM2"),("W1_MBmax","UVW1")
    ]
    err_bands  = [
        ("B_eMmax","B"),("V_eBMmax","V"),("U_eMBmax","U"),
        ("W2_eMBmax","UVW2"),("M2_eMBmax","UVM2"),("W1_eMBmax","UVW1")
    ]

    # Pre-compute model medians for imputation
    model_medians = {band: models_df[band].median() for _, band in phot_bands}
    # print(model_medians)
    # Use a large default error tolerance
    default_err = max(models_df[col].std() for col in numeric_cols if col in models_df)

    # Threshold for missing mags
    missing_threshold = 1

    all_stats = []
    models_passed = 0
    for idx, row in enumerate(data_df.itertuples()):
        # Check raw magnitude completeness
        raw_mags = []
        for col_name, band in phot_bands:
            val = getattr(row, col_name, getattr(row, col_idx[col_name])).strip()
            raw_mags.append(float(val))
        # print(raw_mags)
        missing_count = sum(pd.isna(v) for v in raw_mags)
        if missing_count > missing_threshold:
            # print(f"DEBUG: Skipping row {idx}: {missing_count} missing mags > {missing_threshold}")
            continue

        # Extract mags and errs
        abs_mags, abs_errs = {}, {}
        for col_name, band in phot_bands:
            val = getattr(row, col_name, getattr(row, col_idx[col_name]))
            abs_mags[band] = float(val) if pd.notna(val) else model_medians[band]
        for col_name, band in err_bands:
            val = getattr(row, col_name, getattr(row, col_idx[col_name]))
            abs_errs[band] = float(val) if pd.notna(val) else default_err
        
        # print(abs_mags)
        if any(np.isnan(mag) for mag in abs_mags.values()):
            # print(f"DEBUG: NaN value found in magnitudes for row {idx}")
            continue
        if any(np.isnan(err) for err in abs_errs.values()):
            # print(f"DEBUG: NaN value found in errors for row {idx}")
            continue

        # Extract redshift and MW reddening
        z_global = float(getattr(row, "zhel", getattr(row, col_idx.get("zhel"))))
        ebv_mw   = float(getattr(row, "EBVmw", getattr(row, col_idx.get("EBVmw"))))
        # Host galaxy E(B-V)
        dmb = float(getattr(row, "dm15b", getattr(row, col_idx["dm15b"])))
        bv  = float(getattr(row, "BV", getattr(row, col_idx["BV"])))
        ebv_host = bv - (0.114 * (dmb - 1.1) - 0.07)

        # print(ebv_host)

        if pd.isna(ebv_host):
            # print(f"DEBUG: NaN value found in host E(B-V) for row {idx}")
            continue

        # Filter models by photometry
        filtered = filter_models_by_photometry(models_df, abs_mags, abs_errs)
        if filtered.empty:
            continue

        # print(filtered)

        models_passed += 1

        # print(abs_mags)

        # Compute corrections
        corr_rows = []
        bands = [band for _, band in phot_bands]
        for _, mod in tqdm(filtered.iterrows(), total=len(filtered), desc=f"Row {idx}"):
            spec_file = next(
                f for f in config["input"]["spectra_files"]
                if os.path.basename(f).startswith(mod["Spectrum"]))
            wave, _, flux_mod = process_spectrum(spec_file, lb=mod["lb"])
            if wave is None:
                continue
            synth1 = kfun(wave, flux_mod, filters)
            flux1  = apply_reddening(wave, flux_mod, -ebv_mw*3.1, RV=3.1)
            wave_r = wave / (1 + z_global)
            flux2  = flux1 * (1 + z_global)
            flux_f = apply_reddening(wave_r, flux2, -ebv_host*2.6, RV=2.6)
            synth2 = kfun(wave_r, flux_f, filters)
            corr_rows.append({b: synth1[b] - synth2[b] for b in bands})

        if not corr_rows:
            print(f"DEBUG: No corrections computed for row {idx}")
            continue

        df_corr = pd.DataFrame(corr_rows)
        mean_c, std_c = df_corr.mean(), df_corr.std()
        stats = {"Spectrum": getattr(row, "snname2", getattr(row, "Spectrum", idx))}
        for b in bands:
            stats[f"mean_{b}"] = mean_c[b]
            stats[f"std_{b}"]  = std_c[b]
        all_stats.append(stats)

    print(f"{models_passed} models passed")

    # Write summary CSV
    summary_df = pd.DataFrame(all_stats)
    if summary_df.empty:
        return
    out_dir = config.get("output", {}).get("output_dir", ".")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "per_sn_corrections.csv")
    summary_df.to_csv(out_path, index=False)
    print(f"Wrote per-supernova corrections to {out_path}")
