include "ROOT.pxi"
import math
import correctionlib  # make sure loaded
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

# Import the new C++ class
cdef extern from "HistManager.h":
    cdef struct RangeCut:
        int var_idx
        float min_val
        float max_val

    cdef struct BranchMapping:
        string pair_mass
        string run
        string probe_pt
        string probe_eta
        string probe_phi
        string probe_sc_eta
        string probe_r9
        string probe_seedGain
        string probe_et
        string tag_pt
        string tag_eta
        string tag_phi
        string tag_sc_eta
        string tag_r9
        string tag_seedGain

    cdef cppclass HistManager:
        HistManager(TChain* chain, TFile* outfile, bool isMC, string flag_selection, string base_selection, string sample_name, int max_event) except +
        void setBranchMapping(const BranchMapping& mapping)
        void addBin(string name, string title, int nbins, double min, double max, string cut_string, vector[RangeCut] range_cuts)
        void setReweightMap(vector[double] r9_bins, vector[double] eta_bins, vector[vector[double]] ratio_map)
        void configureCorrections(string json_path)
        void setBranches(vector[string] branches)
        void process()
        void finalize()

##################
# Helper functions
##################

# Check if a string can be a number
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# Default branch mapping configurations
def get_default_branch_mapping(tnp_type="electron"):
    """
    Returns default branch mapping dictionary based on TnP type.
    tnp_type: "electron" (default) or "legacy_electron" or custom
    """
    if tnp_type == "legacy_electron":
        # Legacy naming convention (el_*, tag_Ele_*)
        return {
            'pair_mass': 'pair_mass',
            'run': 'run',
            'probe_pt': 'el_pt',
            'probe_eta': 'el_eta',
            'probe_phi': 'el_phi',
            'probe_sc_eta': 'el_sc_eta',
            'probe_r9': 'el_r9',
            'probe_seedGain': 'el_seedGain',
            'probe_et': 'el_et',
            'tag_pt': 'tag_Ele_pt',
            'tag_eta': 'tag_Ele_eta',
            'tag_phi': 'tag_Ele_phi',
            'tag_sc_eta': 'tag_sc_eta',
            'tag_r9': 'tag_Ele_r9',
            'tag_seedGain': 'tag_Ele_seedGain',
        }
    else:
        # Default naming convention (prob_*, tag_*)
        return {
            'pair_mass': 'pair_mass',
            'run': 'run',
            'probe_pt': 'prob_pt',
            'probe_eta': 'prob_eta',
            'probe_phi': 'prob_phi',
            'probe_sc_eta': 'prob_sc_eta',
            'probe_r9': 'prob_r9',
            'probe_seedGain': 'prob_seedGain',
            'probe_et': 'prob_et',
            'tag_pt': 'tag_pt',
            'tag_eta': 'tag_eta',
            'tag_phi': 'tag_phi',
            'tag_sc_eta': 'tag_sc_eta',
            'tag_r9': 'tag_r9',
            'tag_seedGain': 'tag_seedGain',
        }

##################################
# To Fill Tag and Probe histograms
##################################

def makePassFailHistograms( sample, flag, bindef, var, do_correction=False, correction_json_path=None, branch_mapping=None, reweight_json=None, max_event = -1 ):
    import json

    # Declare all C variables at the beginning
    cdef TChain* tree
    cdef TFile* outfile
    cdef HistManager* manager
    cdef BranchMapping cpp_mapping
    cdef vector[RangeCut] vec_range_cuts
    cdef RangeCut temp_cut
    cdef int IDX_EL_ET = 0
    cdef int IDX_EL_PT = 1
    cdef int IDX_TAG_PT = 2
    cdef vector[string] vec_branches
    cdef vector[double] r9_bins_vec
    cdef vector[double] eta_bins_vec
    cdef vector[vector[double]] ratio_vec
    cdef vector[double] row_vec

    ###############################
    # Read in Tag and Probe Ntuples
    ###############################

    print(sample.tree)
    tree = new TChain(sample.tree)

    for p in sample.path:
        print(' adding rootfile: ', p)
        tree.Add(str.encode(p))

    if not sample.puTree is None:
        import os
        if os.path.exists(sample.puTree):
            print(' - Adding weight tree: %s from file %s ' % (sample.weight.split('.')[0], sample.puTree))
            tree.AddFriend(str.encode(sample.weight.split('.')[0]), str.encode(sample.puTree))
        else:
            print(' - Weight tree file not found or not readable, skipping AddFriend:', sample.puTree)

    #################################
    # Prepare hists, cuts and outfile
    #################################

    outfile = new TFile(str.encode(sample.histFile),'recreate')

    # Build base selection string
    base_selection = bindef.get('baseSelection', '')
    if sample.mcTruth:
        base_selection = f'{base_selection} && mcTrue ==1' if base_selection else 'mcTrue==1'
    if sample.cut is not None:
        base_selection = f'{base_selection} && {sample.cut}' if base_selection else sample.cut
    
    # Initialize C++ Manager
    sample_name = getattr(sample, 'name', 'Unknown')
    manager = new HistManager(
        tree, 
        outfile, 
        sample.isMC, 
        flag.encode('utf-8'), 
        base_selection.encode('utf-8'),
        sample_name.encode('utf-8'),
        max_event
    )

    # Configure Branch Mapping
    if branch_mapping is None:
        # Try to get from sample, otherwise use default
        tnp_type = getattr(sample, 'tnp_type', 'electron')
        branch_mapping = get_default_branch_mapping(tnp_type)
    
    cpp_mapping.pair_mass = branch_mapping.get('pair_mass', 'pair_mass').encode('utf-8')
    cpp_mapping.run = branch_mapping.get('run', 'run').encode('utf-8')
    cpp_mapping.probe_pt = branch_mapping.get('probe_pt', 'prob_pt').encode('utf-8')
    cpp_mapping.probe_eta = branch_mapping.get('probe_eta', 'prob_eta').encode('utf-8')
    cpp_mapping.probe_phi = branch_mapping.get('probe_phi', 'prob_phi').encode('utf-8')
    cpp_mapping.probe_sc_eta = branch_mapping.get('probe_sc_eta', 'prob_sc_eta').encode('utf-8')
    cpp_mapping.probe_r9 = branch_mapping.get('probe_r9', 'prob_r9').encode('utf-8')
    cpp_mapping.probe_seedGain = branch_mapping.get('probe_seedGain', 'prob_seedGain').encode('utf-8')
    cpp_mapping.probe_et = branch_mapping.get('probe_et', 'prob_et').encode('utf-8')
    cpp_mapping.tag_pt = branch_mapping.get('tag_pt', 'tag_pt').encode('utf-8')
    cpp_mapping.tag_eta = branch_mapping.get('tag_eta', 'tag_eta').encode('utf-8')
    cpp_mapping.tag_phi = branch_mapping.get('tag_phi', 'tag_phi').encode('utf-8')
    cpp_mapping.tag_sc_eta = branch_mapping.get('tag_sc_eta', 'tag_sc_eta').encode('utf-8')
    cpp_mapping.tag_r9 = branch_mapping.get('tag_r9', 'tag_r9').encode('utf-8')
    cpp_mapping.tag_seedGain = branch_mapping.get('tag_seedGain', 'tag_seedGain').encode('utf-8')
    
    manager.setBranchMapping(cpp_mapping)

    # Configure Corrections
    if do_correction:
        if correction_json_path:
            json_path_py = correction_json_path
            manager.configureCorrections(json_path_py.encode('utf-8'))
        else:
            print('Corrections required but no json file available')

    # Configure Reweighting from JSON
    if reweight_json:
        print(f"Loading reweighting factors from {reweight_json}...")
        with open(reweight_json, 'r') as f:
            rw_data = json.load(f)
        
        # Convert Python lists to C++ vectors
        for val in rw_data['r9_bins']:
            r9_bins_vec.push_back(val)
        
        for val in rw_data['eta_bins']:
            eta_bins_vec.push_back(val)
        
        for row in rw_data['ratio']:
            row_vec.clear()
            for val in row:
                row_vec.push_back(val)
            ratio_vec.push_back(row_vec)
        
        # Pass r9_bins as x-axis and eta_bins as y-axis
        manager.setReweightMap(r9_bins_vec, eta_bins_vec, ratio_vec)

    # Prepare Bins and Cuts
    cutBinList = []
    cor_cut_var = []

    for ib in range(len(bindef['bins'])):
        # Prepare Cut String
        bin_cut = bindef['bins'][ib]['cut'] if bindef['bins'][ib]['cut'] else ''
        if sample.isMC and not sample.weight is None:
            cutBin = '( %s ) * %s ' % (bin_cut, sample.weight) if bin_cut else sample.weight
            if sample.maxWeight < 999:
                cutBin = '( %s ) * (%s < %f ? %s : 1.0 )' % (bin_cut, sample.weight, sample.maxWeight, sample.weight) if bin_cut else '(%s < %f ? %s : 1.0 )' % (sample.weight, sample.maxWeight, sample.weight)
        else:
            cutBin = bin_cut if bin_cut else ''
        
        cutBinList.append(cutBin)

        # Prepare Range Cuts
        vec_range_cuts.clear()
        bin_cor_cut = bindef['bins'][ib]['cor_cut']
        for constr in bin_cor_cut:
            var_name = constr['var_name']
            var_name_str = var_name.decode('utf-8') if isinstance(var_name, bytes) else var_name
            cor_cut_var.append(var_name_str)
            
            if var_name == b'el_et' or var_name_str == 'el_et' or var_name_str == branch_mapping.get('probe_et', 'prob_et'):
                temp_cut.var_idx = IDX_EL_ET
            elif var_name == b'el_pt' or var_name_str == 'el_pt' or var_name_str == branch_mapping.get('probe_pt', 'prob_pt'):
                temp_cut.var_idx = IDX_EL_PT
            elif var_name == b'tag_Ele_pt' or var_name_str == 'tag_Ele_pt' or var_name_str == branch_mapping.get('tag_pt', 'tag_pt'):
                temp_cut.var_idx = IDX_TAG_PT
            else:
                print(f"Warning: Unknown variable in cor_cut: {var_name}")
                continue
            
            temp_cut.min_val = constr['min']
            temp_cut.max_val = constr['max']
            vec_range_cuts.push_back(temp_cut)

        # Add Bin to Manager
        manager.addBin(
            bindef['bins'][ib]['name'],
            bindef['bins'][ib]['title'],
            var['nbins'],
            var['min'],
            var['max'],
            cutBin.encode('utf-8'),
            vec_range_cuts
        )

    # Determine Branches to Enable
    replace_patterns = ['&', '|', '-', 'cos(', 'sqrt(', 'fabs(', 'abs(', '(', ')', '>', '<', '=', '!', '*', '/', '[', ']']
    branches_str = " ".join(cutBinList) + " pair_mass " + flag
    if base_selection:
        branches_str += " " + base_selection
    for p in replace_patterns:
        branches_str = branches_str.replace(p, ' ')
    
    branches_set = set([x for x in branches_str.split(" ") if x != '' and not is_number(x)])
    branches_set.update(cor_cut_var) # Add variables used in manual cuts

    print(f'the following branches will be activated: {branches_set}')
    for br in branches_set:
        vec_branches.push_back(br.encode('utf-8'))
    
    manager.setBranches(vec_branches)

    # Run Loop
    manager.process()
    
    # Finalize (Write Hists)
    manager.finalize()

    ##########
    # Clean up
    ##########

    del manager
    outfile.Close()
    del tree