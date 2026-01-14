import numpy as np
import matplotlib.pyplot as plt
import warnings
import uproot
from scipy.optimize import curve_fit
import os
import re

warnings.filterwarnings('ignore')

# ==========================================
# CONFIGURATION
# ==========================================

CONFIG = {
    'file_path': '/Users/yucao/Desktop/WORKINGDIR/Trigger/DY_madgraph_passingHLTSeeded.root',
    'mass_z': 91.1876,
    'width_z': 2.4952,
    'mass_range': (60, 120),
    'convolution_bins': 2000,
}

# Updated parameters for Double-Sided Crystal Ball with Asymmetric Sigma
# sigma   : Width of the Gaussian core on the Left side
# sigma_2 : Width of the Gaussian core on the Right side
# alpha, n: Tail parameters for the Left side
# alpha_2, n_2: Tail parameters for the Right side
FIT_PARAMETERS = {
    'norm':     {'val': 5000, 'min': 0,      'max': np.inf},
    'mean':     {'val': 0.0,  'min': -10,    'max': 10},
    'sigma':    {'val': 1.5,  'min': 0.01,   'max': 10},   # Left Width
    'sos':      {'val': 0.1,  'min': 0.0,    'max': 1.0},  # Smearing
    'alpha':    {'val': 1.5,  'min': 0.01,   'max': 10},   # Left transition
    'n':        {'val': 2.0,  'min': 0.1,    'max': 20},   # Left power
    'alpha_2':  {'val': 1.5,  'min': 0.01,   'max': 10},   # Right transition
    'n_2':      {'val': 2.0,  'min': 0.1,    'max': 20},   # Right power
    'sigma_2':  {'val': 1.5,  'min': 0.01,   'max': 10}    # Right Width
}

Norm_Settings = {
    'norm':     {'val': 5000, 'min': 0,      'max': np.inf},
}

# Define the order of parameters passed to the fit function
PARAM_ORDER = ['norm', 'mean', 'sigma', 'sos', 'alpha', 'n', 'alpha_2', 'n_2', 'sigma_2']

def fix_param(config, param_name, fix_value=None):
    if param_name not in config:
        # If the parameter is missing (e.g. alpha_2 in old configs), add it temporarily
        config[param_name] = {'val': 0, 'min': 0, 'max': 0}
    
    val = config[param_name]['val'] if fix_value is None else fix_value
    config[param_name]['val'] = val
    eps = abs(val) * 1e-6 if val != 0 else 1e-6
    config[param_name]['min'] = val - eps
    config[param_name]['max'] = val + eps

    return config

def update_parameters_from_config(config, param_order=PARAM_ORDER):
    p0_base = []
    bounds_min = []
    bounds_max = []
    
    for k in param_order:
        if k in config:
            p0_base.append(config[k]['val'])
            bounds_min.append(config[k]['min'])
            bounds_max.append(config[k]['max'])
        else:
            # Fallback for missing new parameters, should not happen if config is correct
            p0_base.append(1.0)
            bounds_min.append(0.01)
            bounds_max.append(10.0)
            
    return p0_base, (bounds_min, bounds_max)

global p0_base, bounds
p0_base, bounds = update_parameters_from_config(FIT_PARAMETERS)

# ==========================================
# PHYSICS FUNCTIONS
# ==========================================

def breit_wigner(m):
    """Standard relativistic Breit-Wigner distribution."""
    mz = CONFIG['mass_z']
    wz = CONFIG['width_z']
    gamma = np.sqrt(mz**2 * (mz**2 + wz**2))
    k = (2 * np.sqrt(2) * mz * wz * gamma) / (np.pi * np.sqrt(mz**2 + gamma))
    denom = (m**2 - mz**2)**2 + (mz * wz)**2
    return k / denom

def asym_dscb(delta_x, sigma_L, alpha_L, n_L, alpha_R, n_R, sigma_R):
    """
    Double-Sided Crystal Ball with Asymmetric Sigma.
    
    Parameters:
    -----------
    delta_x : x - mean
    sigma_L : Width for x < 0
    sigma_R : Width for x > 0
    alpha_L : Transition point for left tail (positive value)
    n_L     : Power for left tail
    alpha_R : Transition point for right tail (positive value)
    n_R     : Power for right tail
    """
    # Ensure positive widths and tail parameters
    sigma_L = np.maximum(sigma_L, 1e-5)
    sigma_R = np.maximum(sigma_R, 1e-5)
    abs_alpha_L = np.abs(alpha_L)
    abs_n_L = np.abs(n_L)
    abs_alpha_R = np.abs(alpha_R)
    abs_n_R = np.abs(n_R)

    # Define t based on asymmetric sigma
    # t = (x - mean) / sigma_L  if x < mean
    # t = (x - mean) / sigma_R  if x > mean
    t = np.where(delta_x < 0, delta_x / sigma_L, delta_x / sigma_R)
    
    result = np.zeros_like(delta_x, dtype=np.float64)

    # 1. Gaussian Core: -alpha_L < t < alpha_R
    mask_core = (t > -abs_alpha_L) & (t < abs_alpha_R)
    result[mask_core] = np.exp(-0.5 * t[mask_core]**2)

    # 2. Left Tail: t <= -alpha_L
    if np.any(t <= -abs_alpha_L):
        mask_L = t <= -abs_alpha_L
        A_L = np.power(abs_n_L / abs_alpha_L, abs_n_L) * np.exp(-0.5 * abs_alpha_L**2)
        B_L = abs_n_L / abs_alpha_L - abs_alpha_L
        result[mask_L] = A_L * np.power(B_L - t[mask_L], -abs_n_L)

    # 3. Right Tail: t >= alpha_R
    if np.any(t >= abs_alpha_R):
        mask_R = t >= abs_alpha_R
        A_R = np.power(abs_n_R / abs_alpha_R, abs_n_R) * np.exp(-0.5 * abs_alpha_R**2)
        B_R = abs_n_R / abs_alpha_R - abs_alpha_R
        result[mask_R] = A_R * np.power(B_R + t[mask_R], -abs_n_R)

    return result

def convolution_model(x, norm, mean, sigma, sos, alpha, n, alpha_2, n_2, sigma_2,
                      t_grid, bw_values, dt):
    """
    Convolves Breit-Wigner with Asymmetric Double-Sided Crystal Ball.
    sigma   -> Left Width
    sigma_2 -> Right Width
    """
    # Apply smearing (sos) to both left and right sigmas
    s_L = np.sqrt(sigma**2 + sos**2)
    s_R = np.sqrt(sigma_2**2 + sos**2)
    
    # Grid for convolution
    delta_matrix = x[:, None] - t_grid[None, :] - mean
    
    # Calculate Resolution Model
    res_matrix = asym_dscb(delta_matrix.ravel(), s_L, alpha, n, alpha_2, n_2, s_R)
    res_matrix = res_matrix.reshape(delta_matrix.shape)
    
    # Numerical Convolution
    conv_result = np.sum(bw_values[None, :] * res_matrix, axis=1) * dt
    
    return norm * conv_result

# ==========================================
# DATA LOADING
# ==========================================

def load_histogram(fpath, bin_name):
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"ROOT file not found: {fpath}")
    
    with uproot.open(fpath) as f:
        keys = []
        for k in f.keys():
            key_str = k.decode('utf-8') if isinstance(k, bytes) else k
            keys.append(key_str.split(';')[0])
        
        hist_key = None
        for k in keys:
            if k == bin_name or k.startswith(bin_name + ';'):
                hist_key = k
                break
        if hist_key is None:
            candidates = [k for k in keys if bin_name in k and ('Pass' in k or 'Fail' in k)]
            if candidates:
                hist_key = candidates[0]
            else:
                raise KeyError(f"Histogram '{bin_name}' not found in ROOT file.")
        
        h = f[hist_key]
        try: 
            values, edges = h.to_numpy()
            centers = (edges[:-1] + edges[1:]) / 2
            errors = h.errors()
        except Exception:
            values = h.values
            edges = h.edges
            centers = (edges[:-1] + edges[1:]) / 2
            if hasattr(h, 'variances') and h.variances is not None:
                errors = np.sqrt(h.variances)
            else:
                errors = np.sqrt(np.maximum(values, 0))
        return centers, values, errors

# ==========================================
# FITTING
# ==========================================

def fit_data(x, y, yerr, t_grid, bw_values, dt):
    current_p0 = p0_base.copy()
    # Estimate normalization
    current_p0[0] = np.max(y) * 5 
    
    best_res = None
    
    try:
        def fit_func(x, *p):
            return convolution_model(x, *p, t_grid, bw_values, dt)
        
        # Single fit attempt using Double-Sided model (no tail loop needed)
        popt, pcov = curve_fit(
            fit_func, x, y,
            p0=current_p0, sigma=yerr, absolute_sigma=True,
            bounds=bounds, maxfev=2000
        )
        
        chi2 = np.sum(((y - fit_func(x, *popt)) / yerr)**2) / (len(x) - len(popt))
        best_res = {'popt': popt, 'chi2': chi2}
        
    except Exception as e:
        print(f"Fit failed: {e}")
        return None

    return best_res

# ==========================================
# PLOTTING
# ==========================================

def plot_result(x_data, y_data, y_err, fit_res, t_grid, bw_values, dt, bin_name, output_path=None):
    full_min, full_max = 50, 130
    x_smooth = np.linspace(full_min, full_max, 1000)
    
    if fit_res:
        popt = fit_res['popt']
        y_fit = convolution_model(x_smooth, *popt, t_grid, bw_values, dt)
        status = "SUCCESS"
        chi2_str = f"χ² = {fit_res['chi2']:.4f}"
        param_dict = dict(zip(PARAM_ORDER, popt))
    else:
        popt = p0_base
        y_fit = convolution_model(x_smooth, *p0_base, t_grid, bw_values, dt)
        status = "FAILED"
        chi2_str = "N/A"
        param_dict = dict(zip(PARAM_ORDER, p0_base))

    fig = plt.figure(figsize=(13, 6))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.5, 1], hspace=0.35, wspace=0.2)
    ax_main = fig.add_subplot(gs[:, 0])
    ax_left = fig.add_subplot(gs[0, 1])
    ax_right = fig.add_subplot(gs[1, 1])

    # Main plot
    ax_main.errorbar(x_data, y_data, yerr=y_err, fmt='ko', ms=3, label='Data', alpha=0.6)
    ax_main.plot(x_smooth, y_fit, 'r--', lw=2, label='Fit' if fit_res else 'Initial Guess')
    ax_main.set_title(f"{bin_name}\n(DSCB Asym)", fontsize=12, fontweight='bold')
    ax_main.set_xlabel('Mass (GeV)')
    ax_main.set_ylabel('Events')
    ax_main.legend()
    ax_main.grid(True, alpha=0.3)

    # Parameter box
    param_text = f"✓ Fit Status: {status}\n{chi2_str}\n\nFit Parameters:\n" + "─" * 30 + "\n"
    for k, v in param_dict.items():
        param_text += f"{k:10s} = {v:10.4f}\n"
        
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax_main.text(0.01, 0.97, param_text, transform=ax_main.transAxes, fontsize=8,
                 verticalalignment='top', fontfamily='monospace', bbox=props)

    # Left zoom (Left Tail region)
    ax_left.errorbar(x_data, y_data, yerr=y_err, fmt='ko', ms=2, alpha=0.6)
    ax_left.plot(x_smooth, y_fit, 'r--', lw=1.5)
    ax_left.set_xlim(50, 80)
    ax_left.set_title("Left Tail Region [50, 80]", fontsize=10)
    mask = (x_data >= 50) & (x_data <= 80)
    if np.any(mask):
        ax_left.set_ylim(0, np.max(y_data[mask] + y_err[mask]) * 1.2)
    ax_left.grid(True, alpha=0.3)

    # Right zoom (Right Tail region)
    ax_right.errorbar(x_data, y_data, yerr=y_err, fmt='ko', ms=2, alpha=0.6)
    ax_right.plot(x_smooth, y_fit, 'r--', lw=1.5)
    ax_right.set_xlim(100, 130)
    ax_right.set_title("Right Tail Region [100, 130]", fontsize=10)
    ax_right.set_xlabel("Mass (GeV)")
    mask = (x_data >= 100) & (x_data <= 130)
    if np.any(mask):
        ax_right.set_ylim(0, np.max(y_data[mask] + y_err[mask]) * 1.2)
    ax_right.grid(True, alpha=0.3)

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path)
    plt.close()

def fit_func(fpath, histogram_name, fit_low, fit_high, ngrid):
    t_grid = np.linspace(fit_low, fit_high, ngrid)
    dt = (fit_high - fit_low) / (ngrid - 1)
    bw_values = breit_wigner(t_grid)

    x, y, yerr = load_histogram(fpath, histogram_name)
    fit_result = fit_data(x, y, yerr, t_grid, bw_values, dt)

    plot_parameters = {
        'x_data': x,
        'y_data': y,
        'y_err': yerr,
        'fit_res': fit_result,
        't_grid': t_grid,
        'bw_values': bw_values,
        'dt': dt,
        'bin_name': histogram_name
    }

    return fit_result, plot_parameters

# ==========================================
# PARSING HELPERS
# ==========================================

def convert_root_params_to_dict(root_param_list, is_pass=True):
    """Convert ROOT-style param strings to a clean dict for fitting."""
    suffix = "P" if is_pass else "F"
    
    # Mapping ROOT parameter names to script internal names
    full_to_simple = {
        f"mean{suffix}": "mean",
        f"sigma{suffix}": "sigma",      # Left width
        f"alpha{suffix}": "alpha",      # Left transition
        f"n{suffix}": "n",              # Left power
        f"alpha{suffix}_2": "alpha_2",  # Right transition
        f"n{suffix}_2": "n_2",          # Right power
        f"sigma{suffix}_2": "sigma_2",  # Right width
        f"sos{suffix}": "sos"
    }

    result = {}
    pattern = re.compile(r'^([a-zA-Z_]\w*)\[\s*([^,\]]+)\s*,\s*([^,\]]+)\s*,\s*([^,\]]+)\s*\]$')

    for item in root_param_list:
        match = pattern.match(item.strip())
        if not match:
            continue
        name = match.group(1)
        try:
            val = float(match.group(2))
            min_val = float(match.group(3))
            max_val = float(match.group(4))
        except ValueError:
            continue
            
        if name in full_to_simple:
            simple_name = full_to_simple[name]
            result[simple_name] = {'val': val, 'min': min_val, 'max': max_val}

    # Ensure all expected keys exist
    return {k: result[k] for k in PARAM_ORDER if k in result}

def update_root_params_with_fit(root_param_list, fit_result, is_pass=True, bound_margin=0.5, param_order=PARAM_ORDER, keep_original_range=True):
    """
    Replace initial values in root_param_list with fitted values.
    """
    if fit_result is None:
        return root_param_list[:] 

    suffix = "P" if is_pass else "F"
    
    # Mapping script internal names back to ROOT parameter names
    simple_to_full = {
        "mean": f"mean{suffix}",
        "sigma": f"sigma{suffix}",
        "alpha": f"alpha{suffix}",
        "n": f"n{suffix}",
        "alpha_2": f"alpha{suffix}_2",
        "n_2": f"n{suffix}_2",
        "sigma_2": f"sigma{suffix}_2",
        "sos": f"sos{suffix}"
    }

    updated = {}
    for simple_name, value in zip(param_order, fit_result['popt']):
        if simple_name not in simple_to_full:
            continue
        full_name = simple_to_full[simple_name]
        low = value * (1 - bound_margin)
        high = value * (1 + bound_margin)
        
        # Handle cases where value might be negative or zero (unlikely for sigma/alpha but safe to check)
        if simple_name in ['mean']:
            low = value - abs(value)*bound_margin - 0.1
            high = value + abs(value)*bound_margin + 0.1
        elif low >= high or low < 0: # Ensure positive bounds for widths/tails
             low = value * 0.5
             high = value * 1.5
        if keep_original_range:
            for item in root_param_list:
                if item.startswith(full_name + "["):
                    match = re.match(r'^' + re.escape(full_name) + r'\[\s*([^,\]]+)\s*,\s*([^,\]]+)\s*,\s*([^,\]]+)\s*\]$', item.strip())
                    if match:
                        orig_low = float(match.group(2))
                        orig_high = float(match.group(3))
                        low = orig_low
                        high = orig_high
                    break
        updated[full_name] = f"{full_name}[{value:.6g},{low:.6g},{high:.6g}]"

    new_list = []
    found_params = set()
    
    for item in root_param_list:
        if '[' in item:
            name = item.split('[')[0]
            if name in updated:
                new_list.append(updated[name])
                found_params.add(name)
                continue
        new_list.append(item)
        
    # Append new parameters (like alpha_2) if they were fitted but not in the original list
    for simple_name, full_name in simple_to_full.items():
        if full_name in updated and full_name not in found_params:
            new_list.append(updated[full_name])

    return new_list

# ==========================================
# EXTERNAL CALL INTERFACE
# ==========================================

def extract_fit_range_from_bin(bin_name, tag_cut=35, fit_low_default=60, fit_low_max=80):
    if 'et_' in bin_name:
        pattern = r'et_(\d+p\d+)To(\d+p\d+)'
    elif 'pt_' in bin_name:
        pattern = r'pt_(\d+p\d+)To(\d+p\d+)'
    else:
        print(f"unable to extract et range from bin name: {bin_name}")
        return fit_low_default
    match = re.search(pattern, bin_name)
    
    if not match:
        print(f"unable to extract et range from bin name: {bin_name}")
        return fit_low_default
    
    probe_min_str = match.group(1).replace('p', '.')
    
    try:
        probe_min = float(probe_min_str)
    except ValueError:
        print(f"unable to convert probe min '{probe_min_str}' to float")
        return fit_low_default
    
    fit_low = max(fit_low_default, tag_cut + probe_min)
    fit_low = min(fit_low, fit_low_max)
    return fit_low

def perform_pre_AltSigFit(sample, tnpBin, tnpWorkspaceParam, do_plot=False):
    tnpBin_name = tnpBin['name']
    low, high = CONFIG['mass_range']
    low = extract_fit_range_from_bin(tnpBin_name) or low
    ngrid = CONFIG['convolution_bins']

    fpath = sample.histFile
    work_dir = os.path.dirname(fpath)
    work_dir = os.path.join(work_dir, 'python_prefit_dcb', sample.name)
    os.makedirs(work_dir, exist_ok=True)

    for isPass in [False, True]:
        histogram_name = f"{tnpBin_name}_{'Pass' if isPass else 'Fail'}"
        print(f"--- PreFitting histogram for MC: {histogram_name} ---")
        
        config = convert_root_params_to_dict(tnpWorkspaceParam, is_pass=isPass)
        config.update(Norm_Settings)
        
        # Fix sos (resolution smearing) usually to a small value for MC
        config = fix_param(config, 'sos', fix_value=0.01)

        global p0_base, bounds
        p0_base, bounds = update_parameters_from_config(config)

        fit_result, plot_parameters = fit_func(fpath, histogram_name, low, high, ngrid)
        
        if fit_result:
            tnpWorkspaceParam = update_root_params_with_fit(tnpWorkspaceParam, fit_result, is_pass=isPass, bound_margin=0.2)

        if do_plot:
            os.makedirs(os.path.join(work_dir, 'plots'), exist_ok=True)
            plot_name = os.path.join(work_dir, 'plots', f"PreAltSigFit_{tnpBin_name}_{'Pass' if isPass else 'Fail'}.png")
            plot_result(**plot_parameters, output_path=plot_name)

    # Return params and dummy info (tail logic removed as it's now double-sided)
    more_to_return = {'tailLeft': 0} 
    return tnpWorkspaceParam, more_to_return

# ==========================================
# MAIN EXECUTION
# ==========================================

if __name__ == "__main__":

    # Updated parameter list including alpha_2 and n_2 for the double-sided fit
    tnpParAltSigFit = [
        "meanP[-0.0,-5.0,5.0]", "sigmaP[1,1,6.0]", "sigmaP_2[1.5,1,6.0]",
        "alphaP[2.0,0.5,5.0]", "nP[2.0,0.1,10]",
        "alphaP_2[2.0,0.5,5.0]", "nP_2[2.0,0.1,10]", # New Right Tail Params
        "sosP[1,0.0,5.0]",

        "meanF[-0.0,-5.0,5.0]", "sigmaF[2,1,15.0]", "sigmaF_2[2.0,1,15.0]",
        "alphaF[2.0,0.5,5.0]", "nF[2.0,0.1,10]",
        "alphaF_2[2.0,0.5,5.0]", "nF_2[2.0,0.1,10]", # New Right Tail Params
        "sosF[1,0.0,5.0]",
        
        # Background parameters (unused in signal fit but kept for format)
        "acmsP[60.,50.,75.]", "betaP[0.04,0.01,0.06]", "gammaP[0.1, 0.005, 0.5]", "peakP[89.0,82.0,90.0]",
        "acmsF[60.,50.,75.]", "betaF[0.04,0.01,0.06]", "gammaF[0.1, 0.005, 0.5]", "peakF[89.0,82.0,90.0]",
    ]

    class Sample:
        def __init__(self, name, histFile):
            self.name = name
            self.histFile = histFile
            
    sample = Sample('DY_madgraph', CONFIG['file_path'])
    
    tnpbins = {}
    if os.path.exists(CONFIG['file_path']):
        with uproot.open(CONFIG['file_path']) as f:
            keys = [k.split(';')[0] for k in f.keys()]
            for k in keys:
                if 'Pass' in k or 'Fail' in k:
                    bin_name = k.replace('_Pass', '').replace('_Fail', '')
                    tnpbins[bin_name] = {'name': bin_name}
    
    for tnpBin_name, tnpBin in tnpbins.items():
        print(f"Processing TnP bin: {tnpBin_name}")
        updated_params, more_info = perform_pre_AltSigFit(sample, tnpBin, tnpParAltSigFit, do_plot=True)
        # print("Updated Parameters:")
        # for p in updated_params:
        #     print(p)