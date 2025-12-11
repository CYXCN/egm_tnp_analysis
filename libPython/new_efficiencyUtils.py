import math
from collections import defaultdict

class efficiency:
    """
    A class to hold efficiency data for a single N-dimensional bin.
    """
    iAltBkgModel = 0
    iAltSigModel = 1
    iAltMCSignal = 2
    iAltTagSelec = 3
    iPUup        = 4
    iPUdown      = 5
    iAltFitRange = 6

    def __init__(self, bins, effData, errEffData, effMC, errEffMC, effAltBkgModel, effAltSigModel, effAltMCSignal, effAltTagSel):
        """
        bins: A tuple of tuples, e.g., ((var1_min, var1_max), (var2_min, var2_max), ...)
        """
        self.bins       = bins
        self.effData    = effData
        self.effMC      = effMC
        self.errEffData = errEffData        
        self.errEffMC   = errEffMC
        self.altEff = [-1]*7
        self.syst   = [-1]*9
        self.altEff[self.iAltBkgModel] = effAltBkgModel
        self.altEff[self.iAltSigModel] = effAltSigModel
        self.altEff[self.iAltMCSignal] = effAltMCSignal
        self.altEff[self.iAltTagSelec] = effAltTagSel
        self.systCombined = 0

    def __str__(self):
        bin_str = "\t".join([f"{b[0]:.3f}\t{b[1]:.3f}" for b in self.bins])
        val_str = f"{self.effData:.4f}\t{self.errEffData:.4f}\t{self.effMC:.4f}\t{self.errEffMC:.4f}"
        alt_str = "\t".join([f"{x:.4f}" for x in self.altEff[:4]])
        return f"{bin_str}\t{val_str}\t{alt_str}"

    @staticmethod
    def getSystematicNames():
        return [ 'statData', 'statMC', 'altBkgModel', 'altSignalModel', 'altMCEff', 'altTagSelection' ]

    def combineSyst(self):
        """
        Calculates systematic uncertainties and their combination.
        This version uses the object's own efficiency values as the central value.
        """
        systAltBkg      = self.altEff[self.iAltBkgModel] - self.effData if self.altEff[self.iAltBkgModel] >= 0 else 0
        systAltSig      = self.altEff[self.iAltSigModel] - self.effData if self.altEff[self.iAltSigModel] >= 0 else 0
        systAltMC       = self.altEff[self.iAltMCSignal] - self.effMC   if self.altEff[self.iAltMCSignal] >= 0 else 0
        systAltTagSelec = self.altEff[self.iAltTagSelec] - self.effMC   if self.altEff[self.iAltTagSelec] >= 0 else 0

        self.syst[0] = self.errEffData
        self.syst[1] = self.errEffMC
        self.syst[self.iAltBkgModel+2] = systAltBkg
        self.syst[self.iAltSigModel+2] = systAltSig
        self.syst[self.iAltMCSignal+2] = systAltMC
        self.syst[self.iAltTagSelec+2] = systAltTagSelec
        
        self.systCombined = math.sqrt(sum(s**2 for s in self.syst if s is not None and s >= 0))
        

    def __add__(self,eff):
        if self.effData < 0 :
            return eff.deepcopy()
        if eff.effData < 0 :
            return self.deepcopy()
        
        bins = self.bins
        errData2 = 1.0 / (1.0/(self.errEffData*self.errEffData)+1.0/(eff.errEffData*eff.errEffData))
        wData1   = 1.0 / (self.errEffData * self.errEffData) * errData2
        wData2   = 1.0 / (eff .errEffData * eff .errEffData) * errData2
        newEffData      = wData1 * self.effData + wData2 * eff.effData;
        newErrEffData   = math.sqrt(errData2)
        
        #        errMC2 = 1.0 / (1.0/(self.errEffMC*self.errEffMC)+1.0/(eff.errEffMC*eff.errEffMC))
        #wMC1   = 1.0 / (self.errEffMC * self.errEffMC) * errMC2
        #wMC2   = 1.0 / (eff .errEffMC * eff .errEffMC) * errMC2
        newEffMC      = wData1 * self.effMC + wData2 * eff.effMC;
        newErrEffMC   = 0.00001#math.sqrt(errMC2)

        newEffAltBkgModel = wData1 * self.altEff[self.iAltBkgModel] + wData2 * eff.altEff[self.iAltBkgModel]
        newEffAltSigModel = wData1 * self.altEff[self.iAltSigModel] + wData2 * eff.altEff[self.iAltSigModel]
        newEffAltMCSignal = wData1 * self.altEff[self.iAltMCSignal] + wData2 * eff.altEff[self.iAltMCSignal]
        newEffAltTagSelec = wData1 * self.altEff[self.iAltTagSelec] + wData2 * eff.altEff[self.iAltTagSelec]

        effout = efficiency(self.bins,newEffData,newErrEffData,newEffMC,newErrEffMC,newEffAltBkgModel,newEffAltSigModel,newEffAltMCSignal,newEffAltTagSelec)
        return effout
    


import ROOT as rt
import numpy as np

def makeTGraphFromList( listOfEfficiencies , keyMin, keyMax ):
    import math as _math
    # filter and validate entries
    valid = []
    for p in listOfEfficiencies:
        try:
            xmin = float(p[keyMin])
            xmax = float(p[keyMax])
            y = float(p['val'])
            yerr = float(p.get('err', 0.0))
        except Exception:
            continue
        if not (_math.isfinite(xmin) and _math.isfinite(xmax) and _math.isfinite(y) and _math.isfinite(yerr)):
            continue
        # normalize negative errors to non-negative
        if yerr < 0:
            yerr = abs(yerr)
        valid.append({'min': xmin, 'max': xmax, 'val': y, 'err': yerr})

    if not valid:
        return rt.TGraphErrors(0)

    # sort by bin min to ensure monotonic x for nicer lines
    valid.sort(key=lambda q: q['min'])
    
    grOut = rt.TGraphErrors(len(valid))
    for ip, point in enumerate(valid):
        x = (point['min'] + point['max']) / 2.0
        xerr = (point['max'] - point['min']) / 2.0
        grOut.SetPoint(ip, x, point['val'])
        grOut.SetPointError(ip, xerr, point['err'])
        
    return grOut

def bin_name_conversion(bindef, bin_name_dict, threshold=0.01):
    """
    Converts a bin definition dictionary to use custom bin names if provided.
    bindef: {var_name: (min, max), ...}
    bin_name_dict: {var_name: { (min, max): "custom_name", ... }, ... }
    Returns a new dictionary with converted names.
    """
    def default_bin_name_maker(var_name, bmin, bmax):
        return f"{bmin:.2f}<{var_name}<{bmax:.2f}"
    if not bin_name_dict:
        return [default_bin_name_maker(var_name, bmin, bmax) for var_name, (bmin, bmax) in bindef.items()]

    converted = []
    for var_name, (bmin, bmax) in bindef.items():
        if var_name in bin_name_dict:
            name_map_dict = bin_name_dict[var_name]
            # Try exact match first
            if (bmin, bmax) in name_map_dict:
                converted.append(name_map_dict[(bmin, bmax)])
            else:
                # Try approximate match within threshold
                found = False
                for (min_key, max_key), custom_name in name_map_dict.items():
                    if abs(min_key - bmin) < threshold and abs(max_key - bmax) < threshold:
                        converted.append(custom_name)
                        found = True
                        break
                if not found:
                    converted.append(default_bin_name_maker(var_name, bmin, bmax))
        else:
            converted.append(default_bin_name_maker(var_name, bmin, bmax))
    return converted

class efficiencyList: 
    """
    A class to hold a collection of efficiency objects, indexed by their N-dimensional bins.
    """
    def __init__(self, var_names):
        self.effs = {} # Key: N-dim bin tuple, Value: efficiency object
        self.var_names = var_names
        self.var_map = {name: i for i, name in enumerate(var_names)}
        self._syst_combined = False
        self.efficiency_dict = {}

    def _ensure_eff_dict(self):
        """
        Make sure efficiency_dict is populated from self.effs if empty.
        """
        if not self.efficiency_dict and self.effs:
            self.convertToDict()

    def convertToDict(self):
        """
        Converts the efficiency list to a dictionary keyed by the bins tuple for easy lookup.
        Each entry: bins_tuple -> {'bindef': {var_name: (min,max), ...}, 'efficiency': eff_obj}
        """
        output_dict = {}
        for bins, eff in self.effs.items():
            bindef = {name: eff.bins[self.var_map[name]] for name in self.var_names}
            output_dict[bins] = {
                'bindef': bindef,
                'efficiency' : eff
            }
        self.efficiency_dict = output_dict
        return self.efficiency_dict

    def __str__(self):
        header = "#" + "\t".join([f"{name}_min\t{name}_max" for name in self.var_names])
        header += "\teffData\terrData\teffMC\terrMC\taltBkg\taltSig\taltMC\taltTag"
        self._ensure_eff_dict()
        lines = []
        for entry in self.efficiency_dict.values():
            lines.append(str(entry['efficiency']))
        return header + "\n" + "\n".join(lines)

    def addEfficiency(self, eff):
        # keep self.effs (for backward compatibility) and update dictionary view
        self.effs[eff.bins] = eff
        bindef = {name: eff.bins[self.var_map[name]] for name in self.var_names}
        self.efficiency_dict[eff.bins] = {
            'bindef': bindef,
            'efficiency': eff
        }

    def combine_syst(self):
        """
        Loop over all efficiency points and calculate their combined systematics.
        """
        if self._syst_combined:
            return
        self._ensure_eff_dict()
        for entry in self.efficiency_dict.values():
            entry['efficiency'].combineSyst()
        self._syst_combined = True

    def get_1d_projection(self, plot_var, slice_vars={}, do_sf=False, do_mc=False, bin_name_dict=None):
        """
        Projects the efficiencies onto one variable, optionally slicing others.
        """
        self.combine_syst() # Ensure systematics are calculated
        if plot_var not in self.var_map:
            raise ValueError(f"Variable '{plot_var}' not found. Available: {self.var_names}")

        plot_idx = self.var_map[plot_var]
        
        # Use a dict to group projections. Key is a stable string, value stores the dict of other bins and list of effs.
        projected_effs = {}  # key_str -> {'bindef': {var_name: (min,max), ...}, 'effs': [eff,...]}
        self._ensure_eff_dict()
        for entry in self.efficiency_dict.values():
            eff = entry['efficiency']
            in_slice = all(
                var_name == plot_var or (eff.bins[self.var_map[var_name]][0] >= slice_min and eff.bins[self.var_map[var_name]][1] <= slice_max)
                for var_name, (slice_min, slice_max) in slice_vars.items()
            )
            if not in_slice: 
                continue

            # build dictionary of other variable bins
            other_bins_dict = {k: v for k, v in entry['bindef'].items() if k != plot_var}
            # stable string key for grouping (safe and readable)
            final_bin_names = bin_name_conversion(other_bins_dict, bin_name_dict) 
            key_str = ", ".join(final_bin_names) if final_bin_names else "Full Range"

            if key_str not in projected_effs:
                projected_effs[key_str] = []
            projected_effs[key_str].append(eff)

        graphs = {}
        for key_str, eff_list in projected_effs.items():
            graph_points = []
            for eff in eff_list:
                val, err = eff.effData, eff.systCombined
                if do_sf:
                    val = eff.effData / eff.effMC if eff.effMC > 0 else 1
                    err = val * math.sqrt((eff.systCombined/eff.effData)**2 + (eff.errEffMC/eff.effMC)**2) if eff.effData > 0 and eff.effMC > 0 else 0
                elif do_mc:
                    val, err = eff.effMC, eff.errEffMC

                graph_points.append({'min': eff.bins[plot_idx][0], 'max': eff.bins[plot_idx][1], 'val': val, 'err': err})
            
            # build a readable legend from the stored dict (order follows var_names)
            graphs[key_str] = graph_points
            
        return graphs

    def get_2d_histogram(self, x_var, y_var, slice_vars={}, content_type='sf'):
        """
        Creates a 2D histogram for two specified variables.
        x_var, y_var: Names of variables for X and Y axes.
        slice_vars: Dictionary to slice other dimensions.
        content_type: 'sf', 'eff_data', 'eff_mc', 'unc_abs', 'unc_rel'
                      Additionally supports per-systematic forms:
                        'unc_rel_<systName>'  -> fraction of combined syst from that component
                        'unc_abs_<systName>'  -> absolute contribution from that component (normalized to MC)
        """
        self.combine_syst()
        if x_var not in self.var_map or y_var not in self.var_map:
            raise ValueError(f"One or both variables not found. Available: {self.var_names}")

        x_idx, y_idx = self.var_map[x_var], self.var_map[y_var]
        
        # Support content types like 'unc_rel_statData' or 'unc_abs_statData'
        per_syst_idx = None
        per_syst_mode = None
        if isinstance(content_type, str):
            if content_type.startswith('unc_rel_'):
                per_syst_mode = 'rel'
                syst_name = content_type.split('unc_rel_',1)[1]
                try:
                    per_syst_idx = efficiency.getSystematicNames().index(syst_name)
                except ValueError:
                    per_syst_idx = None
            elif content_type.startswith('unc_abs_'):
                per_syst_mode = 'abs'
                syst_name = content_type.split('unc_abs_',1)[1]
                try:
                    per_syst_idx = efficiency.getSystematicNames().index(syst_name)
                except ValueError:
                    per_syst_idx = None

        # Collect bins and filter efficiencies
        x_bins, y_bins = set(), set()
        filtered_effs = []
        self._ensure_eff_dict()
        for entry in self.efficiency_dict.values():
            eff = entry['efficiency']
            in_slice = all(
                var_name in [x_var, y_var] or (eff.bins[self.var_map[var_name]][0] >= slice_min and eff.bins[self.var_map[var_name]][1] <= slice_max)
                for var_name, (slice_min, slice_max) in slice_vars.items()
            )
            if not in_slice: continue
            
            filtered_effs.append(eff)
            x_bins.add(eff.bins[x_idx][0]); x_bins.add(eff.bins[x_idx][1])
            y_bins.add(eff.bins[y_idx][0]); y_bins.add(eff.bins[y_idx][1])

        if not filtered_effs: return None

        x_bins_np = np.array(sorted(list(x_bins)))
        y_bins_np = np.array(sorted(list(y_bins)))

        h_name = f"h2_{content_type}_{x_var}_vs_{y_var}"
        h_title = f"{content_type} vs {x_var} and {y_var}"
        # If per-systematic, make the title clearer
        if per_syst_idx is not None and per_syst_mode is not None:
            syst_list = efficiency.getSystematicNames()
            sname = syst_list[per_syst_idx] if per_syst_idx < len(syst_list) else str(per_syst_idx)
            h_title = f"{per_syst_mode.upper()} syst '{sname}' vs {x_var} and {y_var}"

        h2 = rt.TH2F(h_name, h_title, len(x_bins_np) - 1, x_bins_np, len(y_bins_np) - 1, y_bins_np)

        for eff in filtered_effs:
            bin_x_center = (eff.bins[x_idx][0] + eff.bins[x_idx][1]) / 2
            bin_y_center = (eff.bins[y_idx][0] + eff.bins[y_idx][1]) / 2
            
            val, err = 0, 0
            sf = eff.effData / eff.effMC if eff.effMC > 0 else 1
            sf_err = 0
            if eff.effData > 0 and eff.effMC > 0:
                # total SF uncertainty propagated (approx)
                sf_err = sf * math.sqrt((eff.systCombined/eff.effData)**2 + (eff.errEffMC/eff.effMC)**2)

            if per_syst_idx is not None and per_syst_idx < len(eff.syst):
                # Produce per-systematic maps
                if per_syst_mode == 'rel':
                    # fraction of combined systematic contributed by this component
                    val = abs(eff.syst[per_syst_idx]) / eff.systCombined if eff.systCombined > 0 else 0
                    err = 0
                elif per_syst_mode == 'abs':
                    # absolute contribution normalized to MC (similar to old implementation)
                    val = abs(eff.syst[per_syst_idx]) / eff.effMC if eff.effMC > 0 else 0
                    err = 0
                else:
                    val, err = 0, 0
            else:
                # Standard content types
                if content_type == 'sf':
                    val, err = sf, sf_err
                elif content_type == 'eff_data':
                    val, err = eff.effData, eff.systCombined
                elif content_type == 'eff_mc':
                    val, err = eff.effMC, eff.errEffMC
                elif content_type == 'unc_abs':
                    val, err = sf_err, 0
                elif content_type == 'unc_rel':
                    val, err = (sf_err / sf) if sf > 0 else 0, 0
                else:
                    # Unknown content type -> zeros (safe fallback)
                    val, err = 0, 0
            
            bin_idx = h2.FindBin(bin_x_center, bin_y_center)
            h2.SetBinContent(bin_idx, val)
            h2.SetBinError(bin_idx, err)
            
        return h2