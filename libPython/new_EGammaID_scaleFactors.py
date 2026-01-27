#!/usr/bin/env python

import sys,os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from math import sqrt
import ROOT as rt
import CMS_lumi, tdrstyle
from new_efficiencyUtils import *

# --- Style and Configuration ---
tdrstyle.setTDRStyle()


effiMin = 0.68
effiMax = 1.2

sfMin = 0.78
sfMax = 1.12

# --- Utility Functions ---
def isFloat( myFloat ):
    try:
        float(myFloat)
        return True
    except:
        return False



graphColors = [rt.kBlack, rt.kGray+1, rt.kRed +1, rt.kRed-2, rt.kAzure+2, rt.kAzure-1, 
               rt.kSpring-1, rt.kYellow -2 , rt.kYellow+1,
               rt.kBlue-7, rt.kGreen+3, rt.kOrange-5, 
               rt.kMagenta+2, rt.kCyan-4, rt.kPink+1, rt.kTeal-2, rt.kViolet+3, rt.kOrange-3, rt.kOrange+1 ]
def get_axis_title(var_name):
    """Generates a readable axis title from a variable name."""
    name_lower = var_name.lower()
    if 'eta' in name_lower: return "SuperCluster #eta"
    if 'et' in name_lower or 'pt' in name_lower: return "p_{T} [GeV]"
    if 'r9' in name_lower: return "R9"
    if 'vtx' in name_lower or 'pv' in name_lower: return "N_{vtx}"
    return var_name

def find_min_max_y(graphs, force_min=None, force_max=None):
    """Finds the min and max Y values across multiple TGraphs for setting plot ranges."""
    min_y, max_y = 999, -999
    if not graphs: return 0.8, 1.2
    
    for g in graphs:
        for i in range(g.GetN()):
            val = g.GetY()[i]
            err = g.GetEY()[i]
            min_y = min(min_y, val - err)
            max_y = max(max_y, val + err)
    
    if force_min is not None: min_y = min(min_y, force_min)
    if force_max is not None: max_y = max(max_y, force_max)
    hard_coded_min, hard_coded_max = 0.0, 2.0
    min_y = max(min_y, hard_coded_min)
    max_y = min(max_y, hard_coded_max)

    range_y = max(max_y - min_y, 0.1) # Avoid zero range
    return min_y - 0.1 * range_y, max_y + 0.25 * range_y

# --- Plotting Functions ---
def better_1D_legend(legend, mapping_func=None, sep=';'):
    def default_mapping(label):
        import re
        # if has abs -> ||
        if 'abs(' in label:
            label = label.replace('abs(', '|').replace(')', '|')
        if 'abs_el_sc_eta' in label:
            label = label.replace('abs_el_sc_eta', '|SC eta|')

        # Remove common particle prefixes if they exist (generic cleanup)
        label = label.replace('el_', '')
        label = label.replace('pho_', '')

        label = label.replace('_', ' ')

        if 'sc' in label.lower():
            label = re.sub(r'sc', 'SC', label, flags=re.IGNORECASE)
        if 'eta' in label.lower():
            label = re.sub(r'eta', '#eta', label, flags=re.IGNORECASE)
        if 'pt' in label.lower():
            label = re.sub(r'pt', 'p_{T}', label, flags=re.IGNORECASE)
        if 'et' in label.lower() and 'eta' not in label.lower():
            label = re.sub(r'et', 'E_{T}', label, flags=re.IGNORECASE)
        if 'r9' in label.lower():
            label = re.sub(r'r9', 'R9', label, flags=re.IGNORECASE)

        # remove 0 if not needed (e.g. 20.0 -> 20, but keep 2.5)
        label = re.sub(r'(\d)\.0+(?!\d)', r'\1', label)
        return label.strip()

    if mapping_func is None:
        mapping_func = default_mapping

    for entry in legend.GetListOfPrimitives():
        label = entry.GetLabel()
        if sep in label:
            parts = label.split(sep)
            new_parts = [mapping_func(part.strip()) for part in parts]
            entry.SetLabel(sep.join(new_parts))
        else:
            entry.SetLabel(mapping_func(entry.GetLabel()))


def EffiGraph1D(data_graphs, mc_graphs, sf_graphs, nameout, plot_var, slice_vars_info, eff_range=None, sf_range=None):

    has_sf = (len(sf_graphs) > 0)

    W = 1200
    H = 1200
    yUp = 0.45
    c = rt.TCanvas(f"c_{plot_var}", f"c_{plot_var}", W, H)
    c.SetTopMargin(0.055)
    c.SetBottomMargin(0.10)
    c.SetLeftMargin(0.12)
    
    
    if has_sf:
        p1 = rt.TPad( 'up', 'up', 0, yUp, 1,   1, 0,0,0)
        p2 = rt.TPad( 'do', 'do', 0,   0, 1, yUp, 0,0,0)
        p1.SetBottomMargin(0.0075)
        p1.SetTopMargin(   c.GetTopMargin()*1/(1-yUp))
        p2.SetTopMargin(   0.0075)
        p2.SetBottomMargin( c.GetBottomMargin()*1/yUp)
        p1.SetLeftMargin( c.GetLeftMargin() )
        p2.SetLeftMargin( c.GetLeftMargin() )
    else:
        p1 = rt.TPad( 'up', 'up', 0, 0, 1,   1, 0,0,0)
        p1.SetBottomMargin(c.GetBottomMargin()+0.05)
        p1.SetTopMargin(c.GetTopMargin())
        p1.SetLeftMargin(c.GetLeftMargin()+0.02)
        p2 = None

    # if 'et' in plot_var.lower() or 'pt' in plot_var.lower():
    #     p1.SetLogx(); p2.SetLogx()

    n_entries = len(data_graphs)
    n_cols = 1
    if n_entries > 15:
        n_cols = 3
    elif n_entries > 10:
        n_cols = 2

    leg = rt.TLegend(0.3,0.80,0.95 ,0.92)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetNColumns(n_cols)
    leg.SetTextSize(0.015 if n_cols >= 2 else 0.02)

    # Determine plot ranges
    if eff_range:
        effi_min, effi_max = eff_range
    else:
        effi_min, effi_max = find_min_max_y(list(data_graphs.values()) + list(mc_graphs.values()), force_min=effiMin, force_max=effiMax)
    
    if has_sf:
        if sf_range:
            sf_min, sf_max = sf_range
        else:
            sf_min, sf_max = find_min_max_y(list(sf_graphs.values()), force_min=sfMin, force_max=sfMax)

    # Style and add to legend
    for i, key in enumerate(sorted(data_graphs.keys())):
        color = graphColors[i % len(graphColors)]
        
        data_graphs[key].SetMarkerColor(color); data_graphs[key].SetLineColor(color); data_graphs[key].SetLineWidth(2)
        if key in sf_graphs:
            sf_graphs[key].SetMarkerColor(color); sf_graphs[key].SetLineColor(color); sf_graphs[key].SetLineWidth(2)
        if key in mc_graphs:
            mc_graphs[key].SetLineColor(color); mc_graphs[key].SetLineStyle(rt.kDashed); mc_graphs[key].SetLineWidth(2); mc_graphs[key].SetMarkerSize(0)
        
        leg.AddEntry(data_graphs[key], key, "PL")
    
    better_1D_legend(leg)

    # Draw Efficiency Panel
    p1.cd()
    # collect global x from data/mc/sf
    all_x = []
    all_x_err = []
    for gdict in (data_graphs, mc_graphs, sf_graphs):
        for g in gdict.values():
            all_x.extend([g.GetX()[i] for i in range(g.GetN())])
            all_x_err.extend([g.GetEX()[i] for i in range(g.GetN())])
    all_x = [x - ex for x, ex in zip(all_x, all_x_err)] + [x + ex for x, ex in zip(all_x, all_x_err)]
    if not all_x:
        xmin, xmax = 0.0, 1.0
    else:
        xmin = min(all_x); xmax = max(all_x)
        if xmin == xmax:
            xmin -= 0.5; xmax += 0.5

    # Create a dummy histogram to set axes so the full x-range is visible
    nbins_dummy = max(1, int(min(100, max(1, round((xmax - xmin) / (0.1 if xmax>xmin else 1.0))))))
    h_dummy = rt.TH1F(f"h_dummy_eff_{plot_var}", "", nbins_dummy, xmin, xmax)
    h_dummy.SetMinimum(effi_min)
    h_dummy.SetMaximum(effi_max)
    h_dummy.GetYaxis().SetTitle("Data efficiency")
    h_dummy.GetYaxis().SetTitleOffset(1.0)
    h_dummy.GetXaxis().SetTitle(get_axis_title(plot_var))
    h_dummy.Draw()  # draw axes first

    lineAtOne = rt.TLine(xmin,1,xmax,1)
    lineAtOne.SetLineStyle(rt.kDashed)
    lineAtOne.SetLineWidth(2)
    lineAtOne.Draw()

    # now draw all graphs on top
    first_drawn = False
    for key in sorted(data_graphs.keys()):
        gr = data_graphs[key]
        if not first_drawn:
            gr.Draw("P same")
            first_drawn = True
        else:
            gr.Draw("P same")
        if key in mc_graphs:
            mc_graphs[key].Draw("L same")

    # Draw Scale Factor Panel
    if has_sf:
        p2.cd()

        nbins_dummy_sf = max(1, int(min(100, max(1, round((xmax - xmin) / (0.1 if xmax>xmin else 1.0))))))
        h_dummy_sf = rt.TH1F(f"h_dummy_sf_{plot_var}", "", nbins_dummy_sf, xmin, xmax)
        h_dummy_sf.SetMinimum(sf_min)
        h_dummy_sf.SetMaximum(sf_max)
        h_dummy_sf.GetYaxis().SetTitle("Data / MC")
        h_dummy_sf.GetXaxis().SetTitle(get_axis_title(plot_var))
        h_dummy_sf.GetYaxis().SetTitleOffset(1.0)
        h_dummy_sf.GetXaxis().SetTitleOffset(1.0)
        h_dummy_sf.Draw()

        first_drawn = False
        for key in sorted(sf_graphs.keys()):
            gr = sf_graphs[key]
            gr.Draw("P same")

        lineAtOne = rt.TLine(xmin,1,xmax,1)
        lineAtOne.SetLineStyle(rt.kDashed)
        lineAtOne.SetLineWidth(2)
        lineAtOne.Draw()

    c.cd(); p1.Draw(); 
    if has_sf: p2.Draw(); 
    leg.Draw()
    CMS_lumi.CMS_lumi(c, 5, 10)
    
    # Save plots
    slice_str = '_'.join([f"{v}{m}-{M}" for v, (m, M) in slice_vars_info.items()]).replace('.','p')
    base_savename = f"{nameout}_SFvs{plot_var}"
    if slice_str: base_savename += f"_slice_{slice_str}"

    for ext in ['png']:#["pdf", "png", "C"]:
        c.SaveAs(f"{base_savename}.{ext}")

def make_uniform_2d_from_variable_bins(h2):
    import numpy as _np
    nx = h2.GetXaxis().GetNbins()
    ny = h2.GetYaxis().GetNbins()
    if nx == 0 or ny == 0:
        return None, [], []


    x_edges = _np.arange(0, nx+1, dtype=_np.float64)
    y_edges = _np.arange(0, ny+1, dtype=_np.float64)

    hname = h2.GetName() + "_uniform"
    htitle = h2.GetTitle() + " (uniform index axis)"
    h_uni = rt.TH2F(hname, htitle, nx, x_edges, ny, y_edges)

    x_labels = []
    y_labels = []

    for ix in range(1, nx+1):
        xlow = h2.GetXaxis().GetBinLowEdge(ix)
        xup  = h2.GetXaxis().GetBinUpEdge(ix)

        if ix < nx:
            lbl = f"{xlow:.2f}"
        else:
            lbl = f"{xup:.2f}"
        x_labels.append(lbl)
        h_uni.GetXaxis().SetBinLabel(ix, lbl)

        for iy in range(1, ny+1):
            if ix == 1:
                ylow = h2.GetYaxis().GetBinLowEdge(iy)
                yup  = h2.GetYaxis().GetBinUpEdge(iy)

                if iy < ny:
                    ylbl = f"{ylow:.2f}"
                else:
                    ylbl = f"{yup:.2f}"
                y_labels.append(ylbl)
                h_uni.GetYaxis().SetBinLabel(iy, ylbl)

            val = h2.GetBinContent(ix, iy)
            err = h2.GetBinError(ix, iy)
            h_uni.SetBinContent(ix, iy, val)
            h_uni.SetBinError(ix, iy, err)


    try:
        h_uni.SetMinimum(h2.GetMinimum())
        h_uni.SetMaximum(h2.GetMaximum())
    except Exception:
        pass

    return h_uni, x_labels, y_labels

def plot_2d(h2, nameout, title, z_min=None, z_max=None, is_sf=False, x_var=None, y_var=None, use_uniform=False):
    """
    Generic 2D plotting function.

    If use_uniform is True, convert a variable-binned TH2 to an index-uniform TH2
    (each bin shown with equal width) using make_uniform_2d_from_variable_bins.
    """
    rt.gStyle.SetPaintTextFormat('1.3f')
    rt.gStyle.SetOptTitle(1)
    rt.gStyle.SetPalette(1)

    # optionally convert variable-binned histogram to uniform index grid for even visual spacing
    h2_plot = h2
    if use_uniform:
        try:
            h_uni, x_labels, y_labels = make_uniform_2d_from_variable_bins(h2)
            if h_uni is not None:
                h2_plot = h_uni
        except Exception:
            h2_plot = h2

    c = rt.TCanvas(f"c_{h2_plot.GetName()}", title, 1350, 1200)
    c.SetRightMargin(0.18); c.SetLeftMargin(0.12); c.SetTopMargin(0.08); c.SetBottomMargin(0.12)

    # preserve legacy behaviour: set logx if variable looks like pt/et (only when relevant)
    try:
        if h2_plot.GetXaxis().GetTitle().lower().startswith('p'):
            c.SetLogx()
    except Exception:
        pass

    # scale z-range
    if is_sf:
        try:
            dmin = 1.0 - h2_plot.GetMinimum()
            dmax = h2_plot.GetMaximum() - 1.0
            dall = max(dmin, dmax, 0.01)
            h2_plot.SetMinimum(1 - dall); h2_plot.SetMaximum(1 + dall)
        except Exception:
            pass
    else:
        if z_min is not None: h2_plot.SetMinimum(z_min)
        if z_max is not None: h2_plot.SetMaximum(z_max)

    h2_plot.SetTitle(title)
    # Set readable axis titles using provided variable names (if given) or fall back to existing axis titles
    if x_var:
        h2_plot.GetXaxis().SetTitle(get_axis_title(x_var))
    else:
        h2_plot.GetXaxis().SetTitle(get_axis_title(h2_plot.GetXaxis().GetTitle() or "X"))
    if y_var:
        h2_plot.GetYaxis().SetTitle(get_axis_title(y_var))
    else:
        h2_plot.GetYaxis().SetTitle(get_axis_title(h2_plot.GetYaxis().GetTitle() or "Y"))

    # if using uniform mode and labels are long, rotate / shrink labels preserved in make_uniform_2d_from_variable_bins
    if use_uniform:
        try:
            h2_plot.GetXaxis().LabelsOption("v")
            h2_plot.GetXaxis().SetLabelSize(0.04)
            h2_plot.GetYaxis().SetLabelSize(0.04)
        except Exception:
            pass

    h2_plot.Draw("colz TEXT45")

    CMS_lumi.CMS_lumi(c, 5, 10)
    for ext in ['png']:#["pdf", "png", "C"]:
        try:
            c.SaveAs(f"{nameout}.{ext}")
        except Exception:
            pass
    return c

def plot_systematic_diagnostic(h2_abs, h2_rel, nameout_base, syst_name, x_var, y_var, pdfout, use_uniform=False):
    """
    Draw a 2-pad canvas: left = absolute uncertainty map, right = relative uncertainty map.
    If use_uniform is True, convert each input histogram to uniform-index axes before drawing.
    Appends to pdfout and saves individual files named with nameout_base and syst_name.
    """
    c2D_Err = rt.TCanvas(f'canScaleFactor_{syst_name}', f'canScaleFactor: {syst_name}', 1000, 600)
    c2D_Err.Divide(2, 1)

    for ip in (1, 2):
        pad = c2D_Err.cd(ip)
        pad.SetRightMargin(0.15)
        pad.SetLeftMargin(0.15)
        pad.SetTopMargin(0.10)
        # do not force logy globally; keep previous behaviour conservative
        # pad.SetLogy()

    # helper to optionally convert to uniform
    def _maybe_uniform(h):
        if h is None: return None
        if not use_uniform: return h
        try:
            h_uni, _, _ = make_uniform_2d_from_variable_bins(h)
            return h_uni if h_uni is not None else h
        except Exception:
            return h

    h2a = _maybe_uniform(h2_abs)
    h2r = _maybe_uniform(h2_rel)

    if h2a is not None:
        h2a.SetMinimum(0)
        try:
            mx = h2a.GetMaximum()
            # cap display max to reasonable value for absolute uncertainties
            h2a.SetMaximum(min(mx, 0.2) if mx > 0 else 0.2)
        except Exception:
            pass
        h2a.SetTitle(f"e/#gamma absolute SF syst: {syst_name}")
        if x_var: h2a.GetXaxis().SetTitle(get_axis_title(x_var))
        if y_var: h2a.GetYaxis().SetTitle(get_axis_title(y_var))
        c2D_Err.cd(1)
        h2a.DrawCopy("colz TEXT45")

    if h2r is not None:
        h2r.SetMinimum(0)
        try:
            h2r.SetMaximum(1)
        except Exception:
            pass
        h2r.SetTitle(f"e/#gamma relative SF syst: {syst_name}")
        if x_var: h2r.GetXaxis().SetTitle(get_axis_title(x_var))
        if y_var: h2r.GetYaxis().SetTitle(get_axis_title(y_var))
        c2D_Err.cd(2)
        h2r.DrawCopy("colz TEXT45")

    CMS_lumi.CMS_lumi(c2D_Err, 5, 10)
    # append to multipage pdf
    try:
        c2D_Err.Print(pdfout)
    except Exception:
        pass

    # save individual files
    for iext in ['png']:#["pdf", "png", "C"]:
        try:
            c2D_Err.SaveAs(f"{nameout_base}_Syst_{syst_name}_2D_{y_var}_vs_{x_var}.{iext}")
        except Exception:
            pass

    return c2D_Err

# --- Main Logic ---

def doEGM_SFs(filein, lumi, plot_vars_1d=None, plot_vars_2d=None, slices=None):
    print(f"Opening file: {filein} (lumi: {lumi:.1f} fb^-1)")
    CMS_lumi.lumi_13TeV = f"{lumi:.1f} fb^{{-1}}"

    if not os.path.exists(filein):
        print(f"Error: File '{filein}' not found."); sys.exit(1)

    # --- File Parsing & Default Configuration ---
    var_names = []
    is_json = filein.endswith('.json')
    json_data = None

    if is_json:
        import json
        try:
            with open(filein, 'r') as f:
                # first replace replace('ph_', 'el_').replace('Pho_', 'Ele_')
                tmp_content = f.read().replace('ph_', 'el_').replace('Pho_', 'Ele_')
                json_data = json.loads(tmp_content)
            if not json_data:
                print("Error: JSON file is empty."); sys.exit(1)
            # Infer variables from the first bin
            first_bin = next(iter(json_data.values()))
            if 'def' in first_bin:
                var_names = list(first_bin['def'].keys())
            else:
                print("Error: JSON format incorrect (missing 'def' in entries)."); sys.exit(1)
        except Exception as e:
            print(f"Error reading JSON file: {e}"); sys.exit(1)
    else:
        with open(filein, 'r') as f:
            for line in f:
                if line.strip().startswith("### var"):
                    try: var_names.append(line.split(":")[1].strip())
                    except IndexError: print(f"Warning: Could not parse: {line.strip()}")
    
    if not var_names:
        print("Error: No variable definitions found in file."); sys.exit(1)
    print(f"Found variables: {var_names}")

    # Set defaults if arguments are not provided
    if slices is None:
        slices = {}
    if plot_vars_1d is None:
        plot_vars_1d = var_names
        print(f"INFO: No 1D plot variables provided. Defaulting to all found variables: {plot_vars_1d}")
    if plot_vars_2d is None:
        if len(var_names) >= 2:
            plot_vars_2d = [var_names[-1], var_names[-2]] # Common convention: X=et, Y=eta
            print(f"INFO: No 2D plot variables provided. Defaulting to X='{plot_vars_2d[0]}', Y='{plot_vars_2d[1]}'")
        else:
            plot_vars_2d = [] # Not enough vars for 2D plots

    num_vars = len(var_names)
    eff_list = efficiencyList(var_names)

    def bin_range_translation(low, high, var_name=None):
        import fnmatch
        # convert 9999 and 999999 to reasonable values
        _translation_dict = {
            '*_et': {'global_min': 20.0, 'global_max': 130.0},
            '*_pt': {'global_min': 20.0, 'global_max': 130.0},
            # '*_r9': {'global_min': 0.75, 'global_max': 1.1},
        }
        
        if var_name is not None and var_name != '':
            for pattern, limits in _translation_dict.items():
                if fnmatch.fnmatch(var_name, pattern):
                    high = min(high, limits['global_max'])
                    low = max(low, limits['global_min'])
                    return (low, high)
        
        if high > 998:
            high = 2 * low
        
        return (low, high)

    if is_json:
        for bin_key, bin_data in json_data.items():
            # Apply bin translation logic to JSON data for consistency
            if 'def' in bin_data:
                for v in var_names:
                    if v in bin_data['def']:
                        low = float(bin_data['def'][v]['min'])
                        high = float(bin_data['def'][v]['max'])
                        new_low, new_high = bin_range_translation(low, high, v)
                        bin_data['def'][v]['min'] = new_low
                        bin_data['def'][v]['max'] = new_high
            bins = []
            for v in var_names:
                if 'def' in bin_data and v in bin_data['def']:
                    bins.append((float(bin_data['def'][v]['min']), float(bin_data['def'][v]['max'])))
                else:
                    bins.append((-999.0, -999.0))
            effi_info = bin_data.get('info', {})
            eff_list.addEfficiency(efficiencyJson(tuple(bins), effi_info))
        print(f"Successfully loaded {len(eff_list.effs)} efficiency points from JSON.")

    else:
        with open(filein, 'r') as f:
            for line in f:
                if line.strip().startswith('#') or not line.strip(): continue
                parts = line.strip().split()
                if len(parts) < 2 * num_vars + 8 or not all(isFloat(p) for p in parts[:2*num_vars+8]):
                    print(f"Warning: Skipping malformed line: {line.strip()}"); continue
                
                numbers = [float(p) for p in parts]
                bins = tuple(bin_range_translation(numbers[2*i], numbers[2*i+1], var_names[i]) for i in range(num_vars))
                vals = numbers[2*num_vars:]
                eff_list.addEfficiency(efficiency(bins, *vals[:8]))

    print(f"Successfully loaded {len(eff_list.effs)} efficiency points.")
    nameout_base = os.path.splitext(filein)[0]
    pdfout = f"{nameout_base}_egammaPlots.pdf"
    cDummy = rt.TCanvas(); cDummy.Print(pdfout + "[")

    EEEB_bin_name = {
        # 'abs(el_sc_eta)' : {
        #     (0.0, 1.444) : "EB",
        #     (1.566, 2.5)   : "EE"
        # }
        'abs(el_sc_eta)' : {
            (0.0, 1.5) : "EB",
            (1.5, 3.0)   : "EE"
        },
        'abs_el_sc_eta' : {
            (0.0, 1.5) : "EB",
            (1.5, 3.0)   : "EE"
        }
    }

    # --- 1D Plotting ---
    if plot_vars_1d:
        for plot_var in plot_vars_1d:
            print(f"\n--- Generating 1D plots for '{plot_var}' ---")
            data_points = eff_list.get_1d_projection(plot_var, slices, do_sf=False, do_mc=False, bin_name_dict=EEEB_bin_name)
            mc_points = eff_list.get_1d_projection(plot_var, slices, do_sf=False, do_mc=True, bin_name_dict=EEEB_bin_name)
            sf_points = eff_list.get_1d_projection(plot_var, slices, do_sf=True, do_mc=False, bin_name_dict=EEEB_bin_name)
            
            if not data_points:
                print(f"Warning: No data found for '{plot_var}' with current slices. Skipping."); continue

            data_graphs = {k: makeTGraphFromList(v, 'min', 'max') for k, v in data_points.items()}
            mc_graphs = {k: makeTGraphFromList(v, 'min', 'max') for k, v in mc_points.items()}
            sf_graphs = {k: makeTGraphFromList(v, 'min', 'max') for k, v in sf_points.items()}
            
            mc_graphs = {}
            sf_graphs = {}
            
            # Standard plot
            EffiGraph1D(data_graphs, mc_graphs, sf_graphs, nameout_base, plot_var, slices)

            # Zoomed plot
            # EffiGraph1D(data_graphs, mc_graphs, sf_graphs, nameout_base + "_zoom", plot_var, slices, eff_range=(0.5, 1.2), sf_range=(0.9, 1.1))
            # break
    
    do_1D_sum = True
    if do_1D_sum:
        for plot_var in plot_vars_1d:
            print(f"\n--- Generating summed 1D plots for '{plot_var}' ---")
            data_points = eff_list.get_1d_summation(plot_var, {}, do_sf=False, do_mc=False, bin_name_dict=EEEB_bin_name)
            mc_points = eff_list.get_1d_summation(plot_var, {}, do_sf=False, do_mc=True, bin_name_dict=EEEB_bin_name)
            sf_points = eff_list.get_1d_summation(plot_var, {}, do_sf=True, do_mc=False, bin_name_dict=EEEB_bin_name)
            
            if not data_points:
                print(f"Warning: No data found for '{plot_var}' with current slices. Skipping."); continue

            data_graphs = {k: makeTGraphFromList(v, 'min', 'max') for k, v in data_points.items()}
            mc_graphs = {k: makeTGraphFromList(v, 'min', 'max') for k, v in mc_points.items()}
            sf_graphs = {k: makeTGraphFromList(v, 'min', 'max') for k, v in sf_points.items()}

            mc_graphs = {}
            sf_graphs = {}
            
            # Standard plot
            EffiGraph1D(data_graphs, mc_graphs, sf_graphs, nameout_base + "_sum", plot_var, {})
            # Zoomed plot
            EffiGraph1D(data_graphs, mc_graphs, sf_graphs, nameout_base + "_sum_zoom", plot_var, {}, eff_range=(0.5, 1.2), sf_range=(0.9, 1.1))
            # break

    # --- 2D Plotting ---
    exit()
    if len(plot_vars_2d) == 2:
        x_2d, y_2d = plot_vars_2d
        print(f"\n--- Generating 2D plots for '{y_2d}' vs '{x_2d}' ---")
        h2_sf = eff_list.get_2d_histogram(x_2d, y_2d, slices, 'sf')
        h2_err_abs = eff_list.get_2d_histogram(x_2d, y_2d, slices, 'unc_abs')
        
        if h2_sf:
            plot_2d(h2_sf, f"{nameout_base}_SF_2D_{y_2d}_vs_{x_2d}", "Data/MC Scale Factor", is_sf=True, x_var=x_2d, y_var=y_2d, use_uniform=True)
            plot_2d(h2_err_abs, f"{nameout_base}_SFError_2D_{y_2d}_vs_{x_2d}", "SF Absolute Uncertainty", z_min=0, z_max=min(h2_err_abs.GetMaximum(), 0.2), x_var=x_2d, y_var=y_2d, use_uniform=True)
        else:
            print("Warning: No data for 2D plots with current slices. Skipping.")

        # --- Diagnostic Plots ---
        print("\n--- Generating diagnostic plots for systematics ---")
        for syst_name in efficiency.getSystematicNames():
            h2_abs = eff_list.get_2d_histogram(x_2d, y_2d, slices, f'unc_abs_{syst_name}')
            h2_rel = eff_list.get_2d_histogram(x_2d, y_2d, slices, f'unc_rel_{syst_name}')
            if h2_abs is None and h2_rel is None:
                continue
            plot_systematic_diagnostic(h2_abs, h2_rel, nameout_base, syst_name, x_2d, y_2d, pdfout, use_uniform=True)
    else:
        print("\nINFO: Skipping 2D plots (not requested or not enough variables).")

    cDummy.Print(pdfout + "]")
    print(f"\nAll plots saved. Main PDF: {pdfout}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Generic EGamma scale factor plotter.')
    parser.add_argument('txtFile', help='EGM formatted txt file with variable headers.')
    parser.add_argument('--lumi', type=float, default=35.9, help='Luminosity for plot label (fb^-1).')
    parser.add_argument('--plot-var-1d', action='append', help='Variable to plot on X-axis for 1D plots. Can be used multiple times.')
    parser.add_argument('--plot-var-2d', nargs=2, metavar=('X_VAR', 'Y_VAR'), help='Variables for 2D plots (X and Y).')
    parser.add_argument('--slice', nargs=3, action='append', metavar=('VAR', 'MIN', 'MAX'),
                        help="Slice a variable to a range. E.g., --slice el_r9 0.9 999")

    args = parser.parse_args()

    CMS_lumi.writeExtraText = 1
    CMS_lumi.lumi_sqrtS = "13.6 TeV"

    slices = {}
    if args.slice:
        for var, min_val, max_val in args.slice:
            if not isFloat(min_val) or not isFloat(max_val):
                print(f"Error: Slice values for '{var}' must be numbers."); sys.exit(1)
            slices[var] = (float(min_val), float(max_val))

    # Default plot variables if not specified
    plot_vars_1d = None#args.plot_var_1d if args.plot_var_1d else ['el_et', 'abs(el_sc_eta)']
    plot_vars_2d = None#args.plot_var_2d if args.plot_var_2d else ['el_et', 'abs(el_sc_eta)']

    doEGM_SFs(args.txtFile, args.lumi, plot_vars_1d, plot_vars_2d, slices)