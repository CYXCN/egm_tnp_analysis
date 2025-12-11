include "ROOT.pxi"
import math
#from fitUtils import *
import correctionlib  # make sure loaded

cdef extern from "TLorentzVector.h":
    cdef cppclass TLorentzVector:
        TLorentzVector()
        void SetPtEtaPhiM(double, double, double, double)
        double M()
        TLorentzVector operator+(TLorentzVector)

cdef extern from "scaleSmearing.h":
    cdef cppclass Run3Corrector:
        Run3Corrector(string json_path) except +
        double get_correction(bool is_mc, int run, float pt, float r9, float sc_eta, int seed_gain)

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

cdef void removeNegativeBins(TH1D* h):
    for i in xrange(h.GetNbinsX()):
        if (h.GetBinContent(i) < 0):
            h.SetBinContent(i, 0)

def get_json_path(year, obj_type="Ele"):
    import os
    REPO_DICT = {
        "2025": "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",
        "2024": "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",
        "2023": "Run3-23CSep23-Summer23-NanoAODv12",
        "2023BPIX": "Run3-23DSep23-Summer23BPix-NanoAODv12",
        "2022": "Run3-22CDSep23-Summer22-NanoAODv12",
        "2022EE": "Run3-22EFGSep23-Summer22EE-NanoAODv12",
    }
    base_path = "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM"
    repo = REPO_DICT.get(str(year), REPO_DICT["2024"])
    json_file = "electronSS_EtDependent.json.gz" if obj_type == "Ele" else "photonSS_EtDependent.json.gz"
    full_path = os.path.join(base_path, repo, "latest", json_file)
    if not os.path.exists(full_path):
        print(f"Warning: Correction file not found at {full_path}")
    return full_path

##################################
# To Fill Tag and Probe histograms
##################################

def makePassFailHistograms( sample, flag, bindef, var, do_correction=False, era='2024' ):

    #####################
    # C++ Initializations
    #####################

    # For tree branches
    cdef float pair_mass

    # variables needed for the correction
    cdef unsigned int run  # Data run number
    cdef float el_pt, el_eta, el_phi, el_sc_eta, el_r9, el_seedGain, el_et
    cdef float tag_Ele_pt, tag_Ele_eta, tag_Ele_phi, tag_sc_eta, tag_Ele_5x5_r9, tag_Ele_seedGain

    # Correction variables
    cdef double corr_factor_probe = 1.0
    cdef double corr_factor_tag = 1.0
    cdef double corrected_mass = 0.0

    # TLorentzVectors for mass recalculation
    cdef TLorentzVector vProbe
    cdef TLorentzVector vTag

    cdef Run3Corrector* cpp_corrector = NULL 

    # For the loop
    cdef int nbins = 0
    cdef int nevts
    cdef int frac_of_nevts
    cdef int index
    cdef int bnidx
    cdef int outcount = 0
    cdef double weight
    
    cdef TChain* tree

    cdef TTreeFormula* flag_formula
    cdef vector[TTreeFormula*] bin_formulas

    cdef vector[TH1D*] hPass
    cdef vector[TH1D*] hFail

    cdef TList formulas_list

    cdef double epass = -1.0
    cdef double efail = -1.0


    ###############################
    # Initialize Corrector
    ###############################
    # Assume 'year' is accessible via sample object or hardcoded
    if do_correction:
        current_year = "2024" 
        if hasattr(sample, 'year'):
            current_year = sample.year
            
        json_path_py = get_json_path(current_year, "Ele")
        print(f"Initializing C++ Run3 Corrector with {json_path_py}...")
        cpp_corrector = new Run3Corrector(str.encode(json_path_py))

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

    cdef TFile* outfile = new TFile(str.encode(sample.histFile),'recreate')

    cutBinList = []
    
    # Step 1: Create Histograms and Generate Cut Strings (Do NOT create Formulas yet)
    for ib in range(len(bindef['bins'])):
        hPass.push_back(new TH1D(b'%s_Pass' % bindef['bins'][ib]['name'],bindef['bins'][ib]['title'],var['nbins'],var['min'],var['max']))
        hFail.push_back(new TH1D(b'%s_Fail' % bindef['bins'][ib]['name'],bindef['bins'][ib]['title'],var['nbins'],var['min'],var['max']))
        hPass[ib].Sumw2()
        hFail[ib].Sumw2()

        cuts = bindef['bins'][ib]['cut']
        if sample.mcTruth :
            cuts = '%s && mcTrue==1' % cuts
        if not sample.cut is None :
            cuts = '%s && %s' % (cuts,sample.cut)

        if sample.isMC and not sample.weight is None:
            cutBin = '( %s ) * %s ' % (cuts, sample.weight)
            if sample.maxWeight < 999:
                cutBin = '( %s ) * (%s < %f ? %s : 1.0 )' % (cuts, sample.weight,sample.maxWeight,sample.weight)
        else:
            cutBin = '%s' % cuts

        cutBinList.append(cutBin)

    ######################################
    # Deactivate branches and set adresses
    ######################################

    # Step 2: Identify and Activate Branches (Optimized setstatus)
    tree.SetBranchStatus("*", 0) # Disable all first

    # Find out with variables are used to activate the corresponding branches
    replace_patterns = ['&', '|', '-', 'cos(', 'sqrt(', 'fabs(', 'abs(', '(', ')', '>', '<', '=', '!', '*', '/', '[', ']']
    branches = " ".join(cutBinList) + " pair_mass " + flag
    for p in replace_patterns:
        branches = branches.replace(p, ' ')
    
    # Note: with str.encode we convert a string to bytes, which is needed for C++ functions
    branches = set([str.encode(x) for x in branches.split(" ") if x != '' and not is_number(x)])

    # Add Correction specific branches if needed
    if do_correction:
        # Common
        branches.add(b"run")
        # Probe
        branches.update([b"el_pt", b"el_eta", b"el_phi", b"el_sc_eta", b"el_r9", b"el_seedGain", b"el_et"])
        # Tag
        branches.update([b"tag_Ele_pt", b"tag_Ele_eta", b"tag_Ele_phi", b"tag_sc_eta", b"tag_Ele_5x5_r9", b"tag_Ele_seedGain"])

    print(branches)
    for br in branches:
        tree.SetBranchStatus(br, 1)

    # Step 3: Set Branch Addresses (Optimized setbranch)
    # Essential to do this BEFORE creating TTreeFormula if we want the formula to see modified values
    tree.SetBranchAddress("pair_mass", <void*>&pair_mass)

    if do_correction:
        tree.SetBranchAddress("run", &run)

        # Probe Variables
        tree.SetBranchAddress("el_pt", &el_pt)
        tree.SetBranchAddress("el_eta", &el_eta)
        tree.SetBranchAddress("el_phi", &el_phi)
        tree.SetBranchAddress("el_sc_eta", &el_sc_eta)
        tree.SetBranchAddress("el_r9", &el_r9)
        if tree.GetBranch("el_seedGain"):
            tree.SetBranchAddress("el_seedGain", &el_seedGain)
        else:
            el_seedGain = 12 
        tree.SetBranchAddress("el_et", &el_et)

        # Tag Variables
        tree.SetBranchAddress("tag_Ele_pt", &tag_Ele_pt)
        tree.SetBranchAddress("tag_Ele_eta", &tag_Ele_eta)
        tree.SetBranchAddress("tag_Ele_phi", &tag_Ele_phi)
        tree.SetBranchAddress("tag_sc_eta", &tag_sc_eta)
        tree.SetBranchAddress("tag_Ele_5x5_r9", &tag_Ele_5x5_r9)
        if tree.GetBranch("tag_Ele_seedGain"):
            tree.SetBranchAddress("tag_Ele_seedGain", &tag_Ele_seedGain)
        else:
            tag_Ele_seedGain = 12

    ######################################
    # Create TTreeFormulas
    ######################################
    
    # Step 4: Create Formulas NOW. They will link to the addresses set above.
    flag_formula = new TTreeFormula('Flag_Selection', str.encode(flag), tree)
    formulas_list.Add(<TObject*>flag_formula)

    for i, cutBin in enumerate(cutBinList):
        bin_formulas.push_back(new TTreeFormula(b'%s_Selection' % bindef['bins'][i]['name'], str.encode(cutBin), tree))
        formulas_list.Add(<TObject*>bin_formulas[i])
    nbins = len(cutBinList)
        
    tree.SetNotify(<TObject*> &formulas_list)

    ################
    # Loop over Tree
    ################
    
    nevts = tree.GetEntries()
    frac_of_nevts = nevts/20

    if do_correction:
        print("Starting event loop to fill histograms with Run 3 corrections..")
    else:
        print("Starting event loop to fill histograms..")

    for index in range(nevts):
        #if index > 10000:
        #    break
        if index % frac_of_nevts == 0:
            print(f'{outcount}%, {sample.name}')
            outcount = outcount + 5

        tree.GetEntry(index)

        to_fill_mass = pair_mass

        bin_i = -1
        
        if do_correction:
            # --- Apply Corrections ---

            if run > 390000:
                run = 386946 # for 2025 using the last bin of 2024 for test
            
            # 1. Get Correction factor for Probe
            corr_factor_probe = cpp_corrector.get_correction(
                sample.isMC,
                <int>run,
                el_pt,
                el_r9,
                el_sc_eta,
                <int>el_seedGain
            )
            
            # 2. Get Correction factor for Tag
            corr_factor_tag = cpp_corrector.get_correction(
                sample.isMC,
                <int>run,
                tag_Ele_pt,
                tag_Ele_5x5_r9,
                tag_sc_eta,
                <int>tag_Ele_seedGain
            )
            #print(corr_factor_probe, corr_factor_tag)
            
            # 3. Reconstruct 4-vectors
            # SetPtEtaPhiM(pt, eta, phi, mass). Assuming electron mass ~ 0 for high energy or 0.000511

            for bnidx in range(nbins):
                weight = bin_formulas[bnidx].EvalInstance(0)
                if weight:
                    bin_i = bnidx
                    if flag_formula.EvalInstance(0):
                        hPass[bnidx].Fill(to_fill_mass, weight)
                    else:
                        hFail[bnidx].Fill(to_fill_mass, weight)
                    break
                    
            el_pt = el_pt * corr_factor_probe
            el_et = el_et * corr_factor_probe
            tag_Ele_pt = tag_Ele_pt * corr_factor_tag

            vProbe.SetPtEtaPhiM(el_pt, el_eta, el_phi, 0.000511)
            vTag.SetPtEtaPhiM(tag_Ele_pt, tag_Ele_eta, tag_Ele_phi, 0.000511)
            
            # 4. Calculate new mass
            corrected_mass = (vTag + vProbe).M()
            to_fill_mass = corrected_mass
        compare_namespace = {
            b'el_et' : el_et,
            b'el_pt' : el_pt,
            b'tag_Ele_pt' : tag_Ele_pt,
        }
        for bnidx in range(nbins):
            # -------------------------
            bin_info = bindef['bins'][bnidx]

            # (Manual Python Check)
            pass_corrected = True

            for constr in bin_info['cor_cut']:
                val_to_check = compare_namespace[constr['var_name']] 
                
                if not (constr['min'] < val_to_check < constr['max']):
                    pass_corrected = False
                    break
            if not pass_corrected:
                continue
            
            # TTreeFormula
            pass_static = True
            if bin_formulas[bnidx]:
                if bin_formulas[bnidx].EvalInstance(0) == 0:
                    pass_static = False
            
            if not pass_static:
                continue
            else:
                # Flag (Tag/Probe passing)
                weight = 1.0
                if flag_formula.EvalInstance(0):
                    hPass[bnidx].Fill(to_fill_mass, weight)
                else:
                    hFail[bnidx].Fill(to_fill_mass, weight)
                break

    #####################
    # Deal with the Hists
    #####################

    for ib in range(len(bindef['bins'])):
        removeNegativeBins(hPass[ib])
        removeNegativeBins(hFail[ib])

        hPass[ib].Write(hPass[ib].GetName())
        hFail[ib].Write(hFail[ib].GetName())

        bin1 = 1
        bin2 = hPass[ib].GetXaxis().GetNbins()
        passI = hPass[ib].IntegralAndError(bin1,bin2,epass)
        failI = hFail[ib].IntegralAndError(bin1,bin2,efail)
        eff   = 0
        e_eff = 0
        if passI > 0 :
            itot  = (passI+failI)
            eff   = passI / (passI+failI)
            e_eff = math.sqrt(passI*passI*efail*efail + failI*failI*epass*epass) / (itot*itot)
        #print(cuts)
        #print('    ==> pass: %.1f +/- %.1f ; fail : %.1f +/- %.1f : eff: %1.3f +/- %1.3f' % (passI,epass,failI,efail,eff,e_eff))

    ##########
    # Clean up
    ##########

    if cpp_corrector != NULL:
        del cpp_corrector

    outfile.Close()
    tree.Delete()
