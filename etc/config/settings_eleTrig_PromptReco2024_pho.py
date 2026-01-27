#############################################################
########## General settings
#############################################################
import sys
tnp_config_module = sys.modules.get('__tnp_config__')
if tnp_config_module and hasattr(tnp_config_module, 'used_flag'):
    used_flag = tnp_config_module.used_flag
else:
    used_flag = None
print(f"===> used_flag: {used_flag}")
DO_SEEDED = bool(used_flag and ('unseeded' not in used_flag.lower()))
print(f"===> Running on {'seeded' if DO_SEEDED else 'unseeded'} electrons")
# flag to be Tested
flags = {
    'passingCutBasedLoose94XV2'     : '(passingCutBasedLoose94XV2 == 1)',
    'passingCutBasedMedium94XV2'    : '(passingCutBasedMedium94XV2 == 1)',
    'passingCutBasedTight94XV2'     : '(passingCutBasedTight94XV2 == 1)',
    'passingMVA94XV2wp80'           : '(passingMVA94XV2wp80 == 1)',
    'passingMVA94XV2wp90'           : '(passingMVA94XV2wp90 == 1)',
    'passingCutBasedLoose122XV1'    : '(passingCutBasedLoose122XV1  == 1)',
    'passingCutBasedMedium122XV1'   : '(passingCutBasedMedium122XV1 == 1)',
    'passingCutBasedTight122XV1'    : '(passingCutBasedTight122XV1  == 1)',
    'passingMVA122XV1wp80' : '(passingMVA122XV1wp80 == 1)',
    'passingMVA122XV1wp90' : '(passingMVA122XV1wp90 == 1)',
    'passingMVARun3V1wp80' : '(passingMVARun3V1wp80 == 1)',
    'passingMVARun3V1wp90' : '(passingMVARun3V1wp90 == 1)',
    'passingHLTSeeded' : '(passhltEG30LR9Id85b90eHE12R9Id50b80eR9IdLastFilter == 1) || (passhltEG30LIso60CaloId15b35eHE12R9Id50b80eEcalIsoLastFilter == 1)',
    'passingHLTUnseeded' : '(passhltEG22R9Id85b90eHE12R9Id50b80eR9UnseededLastFilter == 1) || (passhltEG22Iso60CaloId15b35eHE12R9Id50b80eTrackIsoUnseededLastFilter == 1)',
    }
baseOutDir = '/eos/user/y/yucao/program/trigger/CMSSW_11_2_0/src/egm_tnp_analysis/output_photon/2024_EEMC2_5X5_diffR9_rwt_corMC/'

#############################################################
########## samples definition  - preparing the samples
#############################################################
### samples are defined in etc/inputs/tnpSampleDef.py
### not: you can setup another sampleDef File in inputs
import etc.inputs.tnpSampleDef as tnpSamples
tnpTreeDir = 'tnpEleTrig'

samplesDef = {
        'data'  : tnpSamples.Prompt2024_photon['data_EGamma0_2024_Run2024C_0'].clone(),
        'mcNom' : tnpSamples.Prompt2024_photon['DY_madgraph'].clone(), # using DYtoEE, indeed this is 2024 MC
        'tagSel': tnpSamples.Prompt2024_photon['DY_madgraph'].clone(),
        'mcAlt': None,
    }

samplesDef['data'].rename('data_Prompt2024')
## can add data sample easily
dict_2024 = {
    # C: 7.24, D: 7.96, E: 11.32, F: 27.76, G: 37.77, H: 5.44, I: 11.47
# "TnPTree_EGamma0_2024_Run2024C_0.root": 7.24/2,
"TnPTree_EGamma0_2024_Run2024D_0.root": 7.96/2,
"TnPTree_EGamma0_2024_Run2024E_0.root": 11.32/2,
"TnPTree_EGamma0_2024_Run2024F_0.root": 27.76/2,
"TnPTree_EGamma0_2024_Run2024G_0.root": 37.77/2,
"TnPTree_EGamma0_2024_Run2024H_0.root": 5.44/2,
"TnPTree_EGamma0_2024_Run2024I_0.root": 11.47/2,
"TnPTree_EGamma1_2024_Run2024C_1.root": 7.24/2,
"TnPTree_EGamma1_2024_Run2024D_1.root": 7.96/2,
"TnPTree_EGamma1_2024_Run2024E_1.root": 11.32/2,
"TnPTree_EGamma1_2024_Run2024F_1.root": 27.76/2,
"TnPTree_EGamma1_2024_Run2024G_1.root": 37.77/2,
"TnPTree_EGamma1_2024_Run2024H_1.root": 5.44/2,
"TnPTree_EGamma1_2024_Run2024I_1.root": 11.47/2,
}
for filename, lumi in dict_2024.items():
    # sed = 'seeded' if DO_SEEDED else 'unseeded'
    sample_name='data_' + filename.replace('TnPTree_', '').replace('.root', '')
    samplesDef['data'].add_sample( tnpSamples.Prompt2024_photon[sample_name] )


## some sample-based cuts... general cuts defined here after
## require mcTruth on MC DY samples and additional cuts
## all the samples MUST have different names (i.e. sample.name must be different for all)
## if you need to use 2 times the same sample, then rename the second one
#samplesDef['data'  ].set_cut('run >= 273726')
samplesDef['data' ].set_tnpTree(tnpTreeDir)
if not samplesDef['mcNom' ] is None: samplesDef['mcNom' ].set_tnpTree(tnpTreeDir)
if not samplesDef['mcAlt' ] is None: samplesDef['mcAlt' ].set_tnpTree(tnpTreeDir)
if not samplesDef['tagSel'] is None: samplesDef['tagSel'].set_tnpTree(tnpTreeDir)

if not samplesDef['mcNom' ] is None: samplesDef['mcNom' ].set_mcTruth()
if not samplesDef['mcAlt' ] is None: samplesDef['mcAlt' ].set_mcTruth()
if not samplesDef['tagSel'] is None: samplesDef['tagSel'].set_mcTruth()
if not samplesDef['tagSel'] is None:
    samplesDef['tagSel'].rename('mcAltSel_DY_madgraph')
    samplesDef['tagSel'].set_cut('tag_Pho_pt > 37')

## set MC weight, simple way (use tree weight) 
#weightName = 'totWeight'
#if not samplesDef['mcNom' ] is None: samplesDef['mcNom' ].set_weight(weightName)
#if not samplesDef['mcAlt' ] is None: samplesDef['mcAlt' ].set_weight(weightName)
#if not samplesDef['tagSel'] is None: samplesDef['tagSel'].set_weight(weightName)

## set MC weight, can use several pileup rw for different data taking 

puFile = '/eos/cms/store/group/phys_higgs/nonresonant_HH/PrivateProd/Yuxiang/TriggerSF/Pileup/DY_madgraph_ele_2024.pu.puTree.root'
puFile = '/eos/cms/store/group/phys_higgs/nonresonant_HH/PrivateProd/Yuxiang/TriggerSF/Pileup/DY_madgraph_EE_ele_pho2_2024.pu.puTree.root'
weightName = 'weights_2024_CDEFGHI.totWeight'
if not samplesDef['mcNom' ] is None: samplesDef['mcNom' ].set_weight(weightName)
if not samplesDef['mcAlt' ] is None: samplesDef['mcAlt' ].set_weight(weightName)
if not samplesDef['tagSel'] is None: samplesDef['tagSel'].set_weight(weightName)
if not samplesDef['mcNom' ] is None: samplesDef['mcNom' ].set_puTree(puFile)
if not samplesDef['mcAlt' ] is None: samplesDef['mcAlt' ].set_puTree(puFile)
if not samplesDef['tagSel'] is None: samplesDef['tagSel'].set_puTree(puFile)

#############################################################
########## bining definition  [can be nD bining]
#############################################################
biningDef = [
   { 'var' : 'ph_sc_eta' , 'type': 'float', 'bins': [-2.5,-2.0,-1.566,-1.4442, -0.8, 0.0, 0.8, 1.4442, 1.566, 2.0, 2.5] },
   { 'var' : 'ph_et' , 'type': 'float', 'bins': [20,35,50,100,200,500] },
   { 'var' : 'ph_r9', 'type': 'float', 'bins': [0.,0.5,0.56,0.8,0.85,0.9,1.05,2.0]},
#    { 'var' : 'ph_r9', 'type': 'float', 'bins': [0.,0.5,0.55,0.6,0.65,0.7,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.05,2.0]},
]

less_cut = [
    { 'var' : 'abs(ph_sc_eta)' , 'type': 'float', 'bins': [(0,1.446),(1.556001,2.5)] }, # the defalut left side is >=, so using 1.556001 to avoid =
    { 'var' : 'ph_et' , 'type': 'float', 'bins': [20,35,50,100,200,500] },
    { 'var' : 'ph_r9', 'type': 'float', 'bins': [0.,0.5,0.56,0.8,0.85,0.9,1.05,2.0]},
]

if DO_SEEDED:
    biningDef = {
        'EB': [
            { 'var': 'abs(ph_sc_eta)', 'type': 'float', 'bins': [(0, 1.5)] },
            { 'var': 'ph_full5x5x_r9', 'type': 'float', 'bins': [
                (0.5, 0.56),
                (0.56, 0.85),
                (0.85, 0.90),
                (0.90, 999.0)
            ]},
            { 'cor_var': 'ph_et', 'type': 'float', 'bins': [0.0, 35.0, 37.0, 40.0, 45.0, 50.0, 60.0, 70.0, 90.0, 999999.0] },
        ],
        'EE': [
            { 'var': 'abs(ph_sc_eta)', 'type': 'float', 'bins': [(1.5, 3.0)] },
            { 'var': 'ph_full5x5x_r9', 'type': 'float', 'bins': [
                (0.8, 0.85),
                (0.85, 0.90),
                (0.90, 999.0)
            ]},
            { 'cor_var': 'ph_et', 'type': 'float', 'bins': [0.0, 35.0, 37.0, 40.0, 45.0, 50.0, 60.0, 70.0, 90.0, 999999.0] },
        ]
    }
else:
    biningDef = {
        'EB': [
            { 'var': 'abs(ph_sc_eta)', 'type': 'float', 'bins': [(0, 1.5)] },
            { 'var': 'ph_full5x5x_r9', 'type': 'float', 'bins': [
                (0.5, 0.56),
                (0.56, 0.85),
                (0.85, 0.90),
                (0.90, 999.0)
            ]},
            { 'cor_var': 'ph_et', 'type': 'float', 'bins': [0.0, 25.0, 28.0, 31.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 90.0, 999999.0] },
        ],
        'EE': [
            { 'var': 'abs(ph_sc_eta)', 'type': 'float', 'bins': [(1.5, 3.0)] },
            { 'var': 'ph_full5x5x_r9', 'type': 'float', 'bins': [
                (0.8, 0.85),
                (0.85, 0.90),
                (0.90, 999.0)
            ]},
            { 'cor_var': 'ph_et', 'type': 'float', 'bins': [0.0, 25.0, 28.0, 31.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 90.0, 999999.0] },
        ]
    }

# 2024 json config: abs eta: 0.0, 1.5, 3.0,999.0; r9: 0.0,0.56,0.85,0.9,999.0
# if DO_SEEDED:
#     biningDef = [
#             { 'var': 'abs(ph_sc_eta)', 'type': 'float', 'bins': [0.0, 1.5, 3.0] },
#             { 'var': 'ph_full5x5x_r9', 'type': 'float', 'bins': [0.5, 0.56, 0.85, 0.9, 999.0]},
#             { 'cor_var': 'ph_et', 'type': 'float', 'bins': [0.0, 35.0, 37.0, 40.0, 45.0, 50.0, 60.0, 70.0, 90.0, 999999.0] },
#     ]
# else:
#     biningDef = [
#             { 'var': 'abs(ph_sc_eta)', 'type': 'float', 'bins': [0.0, 1.5, 3.0] },
#             { 'var': 'ph_full5x5x_r9', 'type': 'float', 'bins': [0.5, 0.56, 0.85, 0.9, 999.0]},
#             { 'cor_var': 'ph_et', 'type': 'float', 'bins': [0.0, 25.0, 28.0, 31.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 90.0, 999999.0] },
#     ]
        

#############################################################
########## Cuts definition for all samples
#############################################################
### cut
cutBase   = 'abs(tag_sc_eta) < 2.1 && (abs(tag_sc_eta) < 1.4442 || abs(tag_sc_eta) > 1.566) && (tag_Pho_mva122XV1 > -0.5) && (ph_mva122XV1 > -0.9) && (abs(ph_sc_eta) < 1.4442 || abs(ph_sc_eta) > 1.566)'
preselection = 'tag_Pho_hasPixelSeed == 1 && tag_Pho_electronIdx !=-1 && tag_Pho_hoe < 0.08 && ph_hoe < 0.08 && (tag_Pho_full5x5_r9 > 0.8 || tag_Pho_chIso < 20 || tag_Pho_chIso/tag_Pho_pt < 0.3 ) && (ph_full5x5x_r9 > 0.8 || ph_chIso < 20 || ph_chIso/ph_et < 0.3 )'

cutBase += ' && ' + preselection

if not DO_SEEDED:
    cutBase += ' && tag_MatchSeededLeg'

cor_cutBase = 'tag_Pho_pt > 35'
if DO_SEEDED:
    cor_cutBase += ' && ph_et > 30'
else:
    cor_cutBase += ' && ph_et > 22' # photon preselection for unseeded

branch_mapping =  {
            'pair_mass': 'pair_mass',
            'run': 'run',
            'probe_pt': 'ph_et',
            'probe_eta': 'ph_eta',
            'probe_phi': 'ph_sc_phi',
            'probe_sc_eta': 'ph_sc_eta',
            'probe_r9': 'ph_r9',
            'probe_seedGain': 'ph_seedGain',
            'probe_et': 'ph_et',
            'tag_pt': 'tag_Pho_pt',
            'tag_eta': 'tag_Pho_eta',
            'tag_phi': 'tag_Pho_phi',
            'tag_sc_eta': 'tag_sc_eta',
            'tag_r9': 'tag_Pho_r9',
            'tag_seedGain': 'tag_Pho_seedGain',
        }

# r9Eta reweighting json file
if DO_SEEDED:
    r9Eta_reweighting = '/eos/user/y/yucao/HiggsDNA_latest/project/rwt_results_lead.json'
else:
    r9Eta_reweighting = '/eos/user/y/yucao/HiggsDNA_latest/project/rwt_results_sublead.json'
    
# can add addtionnal cuts for some bins (first check bin number using tnpEGM --checkBins)
additionalCuts = { 
#    0 : 'sqrt( 2*event_met_pfmet*tag_Pho_pt*(1-cos(event_met_pfphi-tag_Pho_phi))) < 45',
#    1 : 'sqrt( 2*event_met_pfmet*tag_Pho_pt*(1-cos(event_met_pfphi-tag_Pho_phi))) < 45',
#    2 : 'sqrt( 2*event_met_pfmet*tag_Pho_pt*(1-cos(event_met_pfphi-tag_Pho_phi))) < 45',
#    3 : 'sqrt( 2*event_met_pfmet*tag_Pho_pt*(1-cos(event_met_pfphi-tag_Pho_phi))) < 45',
#    4 : 'sqrt( 2*event_met_pfmet*tag_Pho_pt*(1-cos(event_met_pfphi-tag_Pho_phi))) < 45',
#    5 : 'sqrt( 2*event_met_pfmet*tag_Pho_pt*(1-cos(event_met_pfphi-tag_Pho_phi))) < 45',
#    6 : 'sqrt( 2*event_met_pfmet*tag_Pho_pt*(1-cos(event_met_pfphi-tag_Pho_phi))) < 45',
#    7 : 'sqrt( 2*event_met_pfmet*tag_Pho_pt*(1-cos(event_met_pfphi-tag_Pho_phi))) < 45',
#    8 : 'sqrt( 2*event_met_pfmet*tag_Pho_pt*(1-cos(event_met_pfphi-tag_Pho_phi))) < 45',
#    9 : 'sqrt( 2*event_met_pfmet*tag_Pho_pt*(1-cos(event_met_pfphi-tag_Pho_phi))) < 45',
}

#### or remove any additional cut (default)
#additionalCuts = None

#############################################################
########## fitting params to tune fit by hand if necessary
#############################################################
AltSigFitUsingDSCB = True

tnpParNomFit = [
    "meanP[-0.0,-5.0,5.0]","sigmaP[0.9,0.5,5.0]",
    "meanF[-0.0,-5.0,5.0]","sigmaF[0.9,0.5,5.0]",
    "acmsP[60.,40.,100.]","betaP[0.05,0.01,0.08]","gammaP[0.1, -0.1, 0.5]","peakP[87.0,82.0,90.0]",
    "acmsF[60.,40.,100.]","betaF[0.05,0.01,0.08]","gammaF[0.1, -0.1, 0.5]","peakF[87.0,82.0,90.0]",
    ]

tnpParAltSigFit = [
    "meanP[-0.0,-5.0,5.0]","sigmaP[1,1,6.0]","alphaP[2.0,0.5,3.5]" ,'nP[1,0.1,5]',"sigmaP_2[1.5,1,6.0]","sosP[1,0.1,5.0]",
    "meanF[-0.0,-5.0,5.0]","sigmaF[2,1,15.0]","alphaF[2.0,0.5,3.5]",'nF[1,0.1,5]',"sigmaF_2[2.0,1,6.0]","sosF[1,0.1,5.0]",
    "acmsP[60.,50.,75.]","betaP[0.04,0.01,0.06]","gammaP[0.1, 0.005, 0.5]","peakP[89.0,82.0,90.0]",
    "acmsF[60.,50.,75.]","betaF[0.04,0.01,0.06]","gammaF[0.1, 0.005, 0.5]","peakF[89.0,82.0,90.0]",
    ]
tnpParAltSigFit_DSCB = [
    "meanP[-0.0,-5.0,5.0]", "sigmaP[1,1,6.0]", "sigmaP_2[1.5,1,6.0]",
    "alphaP[2.0,0.5,5.0]", "nP[2.0,0.1,10]",
    "alphaP_2[2.0,0.5,5.0]", "nP_2[2.0,0.1,10]", # New Right Tail Params
    "sosP[1,0.0,5.0]",

    "meanF[-0.0,-5.0,5.0]", "sigmaF[2,1,15.0]", "sigmaF_2[2.0,1,15.0]",
    "alphaF[2.0,0.5,5.0]", "nF[2.0,0.1,10]",
    "alphaF_2[2.0,0.5,5.0]", "nF_2[2.0,0.1,10]", # New Right Tail Params
    "sosF[1,0.0,5.0]",
    
    # Background parameters (unused in signal fit but kept for format)
    "acmsP[60.,40.,100.]", "betaP[0.04,0.01,0.08]", "gammaP[0.1, -1, 1]", "peakP[89.0,82.0,90.0]",
    "acmsF[60.,40.,100.]", "betaF[0.04,0.01,1]", "gammaF[0.1, -1, 1]", "peakF[89.0,82.0,90.0]",
]
     
tnpParAltBkgFit = [
    "meanP[-0.0,-5.0,5.0]","sigmaP[0.9,0.1,5.0]",
    "meanF[-0.0,-5.0,5.0]","sigmaF[0.9,0.1,5.0]",
    "alphaP[0.,-5.,5.]",
    "alphaF[0.,-5.,5.]",
    ]

def make_looser_ranges(parList, factor=1.0):
    newParList = []
    for par in parList:
        if par.find('[')<0:
            newParList.append(par)
            continue
        name = par.split('[')[0]
        vals = par.split('[')[1].replace(']','').split(',')
        if len(vals)!=3:
            newParList.append(par)
            continue
        mean = float(vals[0])
        low  = float(vals[1])
        high = float(vals[2])
        if low>0:
            new_low = low/factor
        else:
            new_low = low*factor
        new_high = high*factor
        newParList.append( f"{name}[{mean},{new_low},{new_high}]" )
    return newParList        

tnpParNomFitLooser = make_looser_ranges(tnpParNomFit)
tnpParNomFitLooser = [
    "meanP[-0.0,-5.0,5.0]","sigmaP[0.9,0.1,5.0]",
    "meanF[-0.0,-5.0,5.0]","sigmaF[0.9,0.1,5.0]",
    "acmsP[60.,40.,100.]","betaP[0.05,0.005, 0.2]","gammaP[0.1, -1, 1]","peakP[87.0,70.0,100.0]",
    "acmsF[60.,40.,100.]","betaF[0.05,0.005, 1]","gammaF[0.1, -1, 1]","peakF[87.0,70.0,100.0]",
    ]
tnpParNomFit = tnpParNomFitLooser
# [
#     "meanP[-0.0,-5.0,5.0]","sigmaP[0.9,0.5,5.0]",
#     "meanF[-0.0,-5.0,5.0]","sigmaF[0.9,0.5,5.0]",
#     "acmsP[60.,40.,120.]","betaP[0.05,0.01,0.16]","gammaP[0.1, -0.1, 0.5]","peakP[87.0,82.0,90.0]",
#     "acmsF[60.,40.,120.]","betaF[0.05,0.01,0.16]","gammaF[0.1, -0.1, 0.5]","peakF[87.0,82.0,90.0]",
#     ]
tnpParAltSigFitLooser = make_looser_ranges(tnpParAltSigFit)
tnpParAltBkgFitLooser = make_looser_ranges(tnpParAltBkgFit)
tnpParAltSigFit_DSCBLooser = make_looser_ranges(tnpParAltSigFit_DSCB)
