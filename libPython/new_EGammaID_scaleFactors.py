#!/usr/bin/env python

import sys,os
from math import sqrt
import ROOT as rt
from . import CMS_lumi, tdrstyle
from .new_efficiencyUtils import *

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
               rt.kBlack, rt.kBlack, rt.kBlack, 
               rt.kBlack, rt.kBlack, rt.kBlack, rt.kBlack, rt.kBlack, rt.kBlack, rt.kBlack ]
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


def EffiGraph1D(data_graphs, mc_graphs, sf_graphs, nameout, plot_var, slice_vars_info):

    W = 800
    H = 800
    yUp = 0.45
    c = rt.TCanvas(f"c_{plot_var}", f"c_{plot_var}", H, W)
    c.SetTopMargin(0.055)
    c.SetBottomMargin(0.10)
    c.SetLeftMargin(0.12)
    
    
    p1 = rt.TPad( 'up', 'up', 0, yUp, 1,   1, 0,0,0)
    p2 = rt.TPad( 'do', 'do', 0,   0, 1, yUp, 0,0,0)
    p1.SetBottomMargin(0.0075)
    p1.SetTopMargin(   c.GetTopMargin()*1/(1-yUp))
    p2.SetTopMargin(   0.0075)
    p2.SetBottomMargin( c.GetBottomMargin()*1/yUp)
    p1.SetLeftMargin( c.GetLeftMargin() )
    p2.SetLeftMargin( c.GetLeftMargin() )

    # if 'et' in plot_var.lower() or 'pt' in plot_var.lower():
    #     p1.SetLogx(); p2.SetLogx()

    leg = rt.TLegend(0.5,0.80,0.95 ,0.92)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)

    # Determine plot ranges
    effi_min, effi_max = find_min_max_y(list(data_graphs.values()) + list(mc_graphs.values()), force_min=effiMin, force_max=effiMax)
    sf_min, sf_max = find_min_max_y(list(sf_graphs.values()), force_min=sfMin, force_max=sfMax)

    # Style and add to legend
    for i, key in enumerate(sorted(data_graphs.keys())):
        color = graphColors[i % len(graphColors)]
        
        data_graphs[key].SetMarkerColor(color); data_graphs[key].SetLineColor(color); data_graphs[key].SetLineWidth(2)
        sf_graphs[key].SetMarkerColor(color); sf_graphs[key].SetLineColor(color); sf_graphs[key].SetLineWidth(2)
        if key in mc_graphs:
            mc_graphs[key].SetLineColor(color); mc_graphs[key].SetLineStyle(rt.kDashed); mc_graphs[key].SetLineWidth(2); mc_graphs[key].SetMarkerSize(0)
        
        leg.AddEntry(data_graphs[key], key, "PL")

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
    # ...existing code...

    c.cd(); p1.Draw(); p2.Draw(); leg.Draw()
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
    with open(filein, 'r') as f:
        for line in f:
            if line.strip().startswith("### var"):
                try: var_names.append(line.split(":")[1].strip())
                except IndexError: print(f"Warning: Could not parse: {line.strip()}")
    
    if not var_names:
        print("Error: No variable definitions ('### varX : name') found in file header."); sys.exit(1)
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

    def bin_valid_translation(low, high):
        # convert 9999 and 999999 to reasonable values
        if high > 999:
            high = 2*low 
        return (low, high)

    with open(filein, 'r') as f:
        for line in f:
            if line.strip().startswith('#') or not line.strip(): continue
            parts = line.strip().split()
            if len(parts) < 2 * num_vars + 8 or not all(isFloat(p) for p in parts[:2*num_vars+8]):
                print(f"Warning: Skipping malformed line: {line.strip()}"); continue
            
            numbers = [float(p) for p in parts]
            bins = tuple(bin_valid_translation(numbers[2*i], numbers[2*i+1]) for i in range(num_vars))
            vals = numbers[2*num_vars:]
            eff_list.addEfficiency(efficiency(bins, *vals[:8]))

    print(f"Successfully loaded {len(eff_list.effs)} efficiency points.")
    nameout_base = os.path.splitext(filein)[0]
    pdfout = f"{nameout_base}_egammaPlots.pdf"
    cDummy = rt.TCanvas(); cDummy.Print(pdfout + "[")

    EEEB_bin_name = {
        'abs(el_sc_eta)' : {
            (0.0, 1.444) : "EB",
            (1.566, 2.5)   : "EE"
        }
    }

    # --- 1D Plotting ---
    if plot_vars_1d:
        for plot_var in plot_vars_1d:
            print(f"\n--- Generating 1D plots for '{plot_var}' ---")
            data_points = eff_list.get_1d_projection(plot_var, slices, do_sf=False, do_mc=False, bin_name_dict=EEEB_bin_name)
            print(data_points)
            mc_points = eff_list.get_1d_projection(plot_var, slices, do_sf=False, do_mc=True, bin_name_dict=EEEB_bin_name)
            sf_points = eff_list.get_1d_projection(plot_var, slices, do_sf=True, do_mc=False, bin_name_dict=EEEB_bin_name)
            
            if not data_points:
                print(f"Warning: No data found for '{plot_var}' with current slices. Skipping."); continue

            data_graphs = {k: makeTGraphFromList(v, 'min', 'max') for k, v in data_points.items()}
            mc_graphs = {k: makeTGraphFromList(v, 'min', 'max') for k, v in mc_points.items()}
            sf_graphs = {k: makeTGraphFromList(v, 'min', 'max') for k, v in sf_points.items()}
            
            EffiGraph1D(data_graphs, mc_graphs, sf_graphs, nameout_base, plot_var, slices)
            # break

    # --- 2D Plotting ---
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
    plot_vars_1d = args.plot_var_1d if args.plot_var_1d else ['el_et', 'abs(el_sc_eta)']
    plot_vars_2d = args.plot_var_2d if args.plot_var_2d else ['el_et', 'abs(el_sc_eta)']

    doEGM_SFs(args.txtFile, args.lumi, plot_vars_1d, plot_vars_2d, slices)