import numpy as np
import matplotlib.pyplot as plt
import warnings
import uproot
from scipy.optimize import curve_fit
import os

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

FIT_PARAMETERS = {
    'norm':     {'val': 5000, 'min': 0,      'max': np.inf},
    'mean':    {'val': 0.0,  'min': -10,    'max': 10},
    'sigma':   {'val': 1.5,  'min': 0.01,   'max': 10},
    'sos':     {'val': 0.1,  'min': 0.0999, 'max': 0.101},
    'alpha':   {'val': 1.5,  'min': 0.01,   'max': 10},
    'n':       {'val': 2.0,  'min': 0.1,    'max': 10},
    'sigma_2': {'val': 2.0,  'min': 0.01,   'max': 10}
}

Norm_Settings = {
    'norm':     {'val': 5000, 'min': 0,      'max': np.inf},
}

PARAM_ORDER = ['norm', 'mean', 'sigma', 'sos', 'alpha', 'n', 'sigma_2']

def fix_param(config, param_name, fix_value=None):
    if param_name not in config:
        raise KeyError(f"Parameter '{param_name}' not found in config.")
    
    val = config[param_name]['val'] if fix_value is None else fix_value
    config[param_name]['val'] = val
    # Use a tiny epsilon relative to the value; fallback to absolute if val is zero
    eps = abs(val) * 1e-6 if val != 0 else 1e-6
    config[param_name]['min'] = val - eps
    config[param_name]['max'] = val + eps

    return config

def update_parameters_from_config(config, param_order=PARAM_ORDER):
    p0_base = [config[k]['val'] for k in param_order]
    bounds = ([config[k]['min'] for k in param_order],
              [config[k]['max'] for k in param_order])
    return p0_base, bounds

global p0_base, bounds
p0_base, bounds = update_parameters_from_config(FIT_PARAMETERS)

# ==========================================
# PHYSICS FUNCTIONS
# ==========================================

def breit_wigner(m):
    mz = CONFIG['mass_z']
    wz = CONFIG['width_z']
    gamma = np.sqrt(mz**2 * (mz**2 + wz**2))
    k = (2 * np.sqrt(2) * mz * wz * gamma) / (np.pi * np.sqrt(mz**2 + gamma))
    denom = (m**2 - mz**2)**2 + (mz * wz)**2
    return k / denom

def roo_cb_ex_gauss(delta_x, sigma, alpha, n, sigma_2, tail_left):
    delta_x = np.atleast_1d(delta_x)
    sigma = np.maximum(sigma, 1e-5)
    sigma_2 = np.maximum(sigma_2, 1e-5)
    
    t = delta_x / sigma
    t0 = delta_x / sigma_2
    abs_alpha = np.abs(alpha)
    abs_n = np.abs(n)
    
    result = np.zeros_like(delta_x, dtype=np.float64)
    
    if tail_left >= 0:  # Left tail
        mask1 = t > 0
        result[mask1] = np.exp(-0.5 * t0[mask1]**2)
        mask2 = (t <= 0) & (t > -abs_alpha)
        result[mask2] = np.exp(-0.5 * t[mask2]**2)
        mask3 = t <= -abs_alpha
        a = np.exp(-0.5 * abs_alpha**2)
        arg = np.clip(n * (t[mask3] + abs_alpha), -100, 100)
        result[mask3] = a * np.exp(arg)
    else:  # Right tail
        mask1 = t0 < 0
        result[mask1] = np.exp(-0.5 * t[mask1]**2)
        mask2 = (t0 >= 0) & (t0 < abs_alpha)
        result[mask2] = np.exp(-0.5 * t0[mask2]**2)
        mask3 = t0 >= abs_alpha
        const_a = np.power(abs_n / abs_alpha, abs_n) * np.exp(-0.5 * abs_alpha**2)
        const_b = abs_n / abs_alpha - abs_alpha
        denom = np.maximum(const_b + t0[mask3], 1e-7)
        result[mask3] = const_a / np.power(denom, abs_n)
        
    return result

def convolution_model(x, norm, mean, sigma, sos, alpha, n, sigma_2,
                      t_grid, bw_values, dt, tail_mode):
    s1 = np.sqrt(sigma**2 + sos**2)
    s2 = np.sqrt(sigma_2**2 + sos**2)
    delta_matrix = x[:, None] - t_grid[None, :] - mean
    res_matrix = roo_cb_ex_gauss(delta_matrix.ravel(), s1, alpha, n, s2, tail_mode)
    res_matrix = res_matrix.reshape(delta_matrix.shape)
    conv_result = np.sum(bw_values[None, :] * res_matrix, axis=1) * dt
    return norm * conv_result

# ==========================================
# DATA LOADING
# ==========================================

def load_histogram(fpath, bin_name):
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"ROOT file not found: {fpath}")
    
    with uproot.open(fpath) as f:
        # Handle both bytes and str keys from uproot
        keys = []
        for k in f.keys():
            key_str = k.decode('utf-8') if isinstance(k, bytes) else k
            keys.append(key_str.split(';')[0])
        
        # Find exact match or best guess
        hist_key = None
        for k in keys:
            if k == bin_name or k.startswith(bin_name + ';'):
                hist_key = k
                break
        if hist_key is None:
            # Try to find by pattern
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
    current_p0[0] = np.max(y) * 5  # scale norm
    
    best_res = None
    best_chi2 = float('inf')
    
    for tail in [1, -1]:
        try:
            def fit_func(x, *p):
                return convolution_model(x, *p, t_grid, bw_values, dt, tail)
            
            popt, pcov = curve_fit(
                fit_func, x, y,
                p0=current_p0, sigma=yerr, absolute_sigma=True,
                bounds=bounds, maxfev=1500
            )
            chi2 = np.sum(((y - fit_func(x, *popt)) / yerr)**2) / (len(x) - len(popt))
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_res = {'popt': popt, 'chi2': chi2, 'tail_mode': tail}
        except Exception:
            continue
    return best_res

# ==========================================
# PLOTTING
# ==========================================

def plot_result(x_data, y_data, y_err, fit_res, t_grid, bw_values, dt, bin_name, output_path=None):
    full_min, full_max = 50, 130
    x_smooth = np.linspace(full_min, full_max, 1000)
    
    if fit_res:
        popt = fit_res['popt']
        tail = fit_res['tail_mode']
        y_fit = convolution_model(x_smooth, *popt, t_grid, bw_values, dt, tail)
        status = "SUCCESS"
        chi2_str = f"χ² = {fit_res['chi2']:.4f}"
        param_dict = dict(zip(PARAM_ORDER, popt))
        param_dict['tail_mode'] = tail
    else:
        popt = p0_base
        tail = 1
        y_fit = convolution_model(x_smooth, *p0_base, t_grid, bw_values, dt, tail)
        status = "FAILED"
        chi2_str = "N/A"
        param_dict = dict(zip(PARAM_ORDER, p0_base))
        param_dict['tail_mode'] = tail

    fig = plt.figure(figsize=(13, 6))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.5, 1], hspace=0.35, wspace=0.2)
    ax_main = fig.add_subplot(gs[:, 0])
    ax_left = fig.add_subplot(gs[0, 1])
    ax_right = fig.add_subplot(gs[1, 1])

    # Main plot
    ax_main.errorbar(x_data, y_data, yerr=y_err, fmt='ko', ms=3, label='Data', alpha=0.6)
    ax_main.plot(x_smooth, y_fit, 'r--', lw=2, label='Fit' if fit_res else 'Initial Guess')
    ax_main.set_title(f"{bin_name}", fontsize=14, fontweight='bold')
    ax_main.set_xlabel('Mass (GeV)')
    ax_main.set_ylabel('Events')
    ax_main.legend()
    ax_main.grid(True, alpha=0.3)

    # Parameter box
    param_text = f"✓ Fit Status: {status}\n{chi2_str}\n\nFit Parameters:\n" + "─" * 30 + "\n"
    for k, v in param_dict.items():
        if k != 'tail_mode':
            param_text += f"{k:10s} = {v:10.4f}\n"
        else:
            param_text += f"{k:10s} = {v:10d}\n"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax_main.text(0.01, 0.97, param_text, transform=ax_main.transAxes, fontsize=8,
                 verticalalignment='top', fontfamily='monospace', bbox=props)

    # Left zoom
    ax_left.errorbar(x_data, y_data, yerr=y_err, fmt='ko', ms=2, alpha=0.6)
    ax_left.plot(x_smooth, y_fit, 'r--', lw=1.5)
    ax_left.set_xlim(50, 80)
    ax_left.set_title("Left Tail [50, 80]", fontsize=10)
    mask = (x_data >= 50) & (x_data <= 80)
    if np.any(mask):
        ax_left.set_ylim(0, np.max(y_data[mask] + y_err[mask]) * 1.2)
    ax_left.grid(True, alpha=0.3)

    # Right zoom
    ax_right.errorbar(x_data, y_data, yerr=y_err, fmt='ko', ms=2, alpha=0.6)
    ax_right.plot(x_smooth, y_fit, 'r--', lw=1.5)
    ax_right.set_xlim(100, 130)
    ax_right.set_title("Right Tail [100, 130]", fontsize=10)
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
    # Precompute convolution grid
    t_grid = np.linspace(fit_low, fit_high, ngrid)
    dt = (fit_high - fit_low) / (ngrid - 1)
    bw_values = breit_wigner(t_grid)

    # Load data
    x, y, yerr = load_histogram(fpath, histogram_name)

    # Fit
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

def convert_root_params_to_dict(root_param_list, is_pass=True):
    """Convert ROOT-style param strings to a clean dict for fitting."""
    import re

    suffix = "P" if is_pass else "F"
    full_to_simple = {
        f"mean{suffix}": "mean",
        f"sigma{suffix}": "sigma",
        f"alpha{suffix}": "alpha",
        f"n{suffix}": "n",
        f"sigma{suffix}_2": "sigma_2",
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

    # Return only expected keys in fixed order
    expected = ['mean', 'sigma', 'alpha', 'n', 'sigma_2', 'sos']
    return {k: result[k] for k in expected if k in result}

def update_root_params_with_fit(root_param_list, fit_result, is_pass=True, bound_margin=0.5, param_order=PARAM_ORDER):
    """
    Replace initial values in root_param_list with fitted values (if fit succeeded).
    New bounds = [val*(1-bound_margin), val*(1+bound_margin)].
    If fit failed (fit_result is None), return original list unchanged.
    """
    if fit_result is None:
        return root_param_list[:]  # no change on fit failure

    suffix = "P" if is_pass else "F"
    simple_to_full = {
        "mean": f"mean{suffix}",
        "sigma": f"sigma{suffix}",
        "alpha": f"alpha{suffix}",
        "n": f"n{suffix}",
        "sigma_2": f"sigma{suffix}_2",
        "sos": f"sos{suffix}"
    }

    # Build mapping from full param name to formatted string
    updated = {}
    for simple_name, value in zip(param_order, fit_result['popt']):
        if simple_name not in simple_to_full:
            continue
        full_name = simple_to_full[simple_name]
        low = value * (1 - bound_margin)
        high = value * (1 + bound_margin)
        # Ensure low < high and avoid zero/negative for scale-like params
        if low >= high:
            low = value * 0.9
            high = value * 1.1
        updated[full_name] = f"{full_name}[{value:.6g},{low:.6g},{high:.6g}]"

    # Replace matching entries in original list
    new_list = []
    for item in root_param_list:
        # Extract param name (before '[')
        if '[' in item:
            name = item.split('[')[0]
            if name in updated:
                new_list.append(updated[name])
                continue
        new_list.append(item)
    return new_list

# for outter usage
def perform_pre_AltSigFit(sample, tnpBin, tnpWorkspaceParam, do_plot=False):
    tnpBin_name = tnpBin['name']

    low, high = CONFIG['mass_range']
    ngrid = CONFIG['convolution_bins']

    fpath = sample.histFile
    work_dir = os.path.dirname(fpath)
    work_dir = os.path.join(work_dir, 'python_prefit', sample.name)
    os.makedirs(work_dir, exist_ok=True)

    more_to_return = {'tailLeft': 0}

    for isPass in [False, True]:
        histogram_name = f"{tnpBin_name}_{'Pass' if isPass else 'Fail'}"
        print(f"--- PreFitting histogram for MC: {histogram_name} ---")
        config = convert_root_params_to_dict(tnpWorkspaceParam, is_pass=isPass)
        config.update(Norm_Settings)

        config = fix_param(config, 'sos', fix_value=0.01)

        global p0_base, bounds
        p0_base, bounds = update_parameters_from_config(config)

        fit_result, plot_parameters = fit_func(fpath, histogram_name, low, high, ngrid)
        more_to_return = {'tailLeft': fit_result['tail_mode'] if fit_result else more_to_return['tailLeft']}

        tnpWorkspaceParam = update_root_params_with_fit(tnpWorkspaceParam, fit_result, is_pass=isPass, bound_margin=0.2)

        if do_plot:
            os.makedirs(os.path.join(work_dir, 'plots'), exist_ok=True)
            plot_name = os.path.join(work_dir, 'plots', f"PreAltSigFit_{tnpBin_name}_{'Pass' if isPass else 'Fail'}.png")
            plot_result(**plot_parameters, output_path=plot_name)
    tail_mode = more_to_return['tailLeft']  # should be 1, -1, or fallback value
    for suffix in ['_tailLeft', '_tailRight']:
        candidate = os.path.join(work_dir, tnpBin_name + suffix)
        if os.path.exists(candidate):
            os.remove(candidate)

    # Create new signal file based on tail mode
    write_dir = os.path.join(work_dir, 'result')
    os.makedirs(write_dir, exist_ok=True)
    if tail_mode == 1:
        with open(os.path.join(write_dir, f"{tnpBin_name}_tailLeft"), 'w') as f:
            pass  # create empty file
    elif tail_mode == -1:
        with open(os.path.join(write_dir, f"{tnpBin_name}_tailRight"), 'w') as f:
            pass  # create empty file
    return tnpWorkspaceParam, more_to_return



# ==========================================
# MAIN EXECUTION
# ==========================================

if __name__ == "__main__":

    tnpParAltSigFit = [
        "meanP[-0.0,-5.0,5.0]","sigmaP[1,1,6.0]","alphaP[2.0,0.5,3.5]" ,'nP[1,0.1,5]',"sigmaP_2[1.5,1,6.0]","sosP[1,0.1,5.0]",
        "meanF[-0.0,-5.0,5.0]","sigmaF[2,1,15.0]","alphaF[2.0,0.5,3.5]",'nF[1,0.1,5]',"sigmaF_2[2.0,1,6.0]","sosF[1,0.1,5.0]",
        "acmsP[60.,50.,75.]","betaP[0.04,0.01,0.06]","gammaP[0.1, 0.005, 0.5]","peakP[89.0,82.0,90.0]",
        "acmsF[60.,50.,75.]","betaF[0.04,0.01,0.06]","gammaF[0.1, 0.005, 0.5]","peakF[89.0,82.0,90.0]",
        ]

    # test:
    class Sample:
        def __init__(self, name, histFile):
            self.name = name
            self.histFile = histFile
    sample = Sample('DY_madgraph', CONFIG['file_path'])
    tnpbins = {}
    with uproot.open(CONFIG['file_path']) as f:
        keys = [k.split(';')[0] for k in f.keys()]
        for k in keys:
            if 'Pass' in k or 'Fail' in k:
                bin_name = k.replace('_Pass', '').replace('_Fail', '')
                tnpbins[bin_name] = {'name': bin_name}
    for tnpBin_name, tnpBin in tnpbins.items():
        print(f"Processing TnP bin: {tnpBin_name}")
        updated_params, more_info = perform_pre_AltSigFit(sample, tnpBin, tnpParAltSigFit, do_plot=True)