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
    'passingHLTSeeded' : '(passhltEG30LR9Id85b90eHE12R9Id50b80eR9IdLastFilter == 1) && (passhltEG30LIso60CaloId15b35eHE12R9Id50b80eEcalIsoLastFilter == 1)',
    'passingHLTUnseeded' : '(passhltEG22R9Id85b90eHE12R9Id50b80eR9UnseededLastFilter == 1) && (passhltEG22Iso60CaloId15b35eHE12R9Id50b80eTrackIsoUnseededLastFilter == 1)',
    }

baseOutDir = '/eos/user/y/yucao/program/trigger/CMSSW_11_2_0/src/egm_tnp_analysis/output3/'

#############################################################
########## samples definition  - preparing the samples
#############################################################
### samples are defined in etc/inputs/tnpSampleDef.py
### not: you can setup another sampleDef File in inputs
import etc.inputs.tnpSampleDef as tnpSamples
tnpTreeDir = 'tnpEleTrig'

samplesDef = {
        'data'  : tnpSamples.Prompt2025['data_EGamma0_2025_Run2025C_0'].clone(),
        'mcNom' : tnpSamples.Prompt2025['DY_madgraph'].clone(),
        'tagSel': tnpSamples.Prompt2025['DY_madgraph'].clone(),
        'mcAlt': None,
    }

samplesDef['data'].rename('data_Prompt2025')

## can add data sample easily
dict_2025 = {
    # C: 20.78, D: 25.29, E: 14, F: 16.11 
# "TnPTree_EGamma0_2025_Run2025C_0.root": 20.78/4,
"TnPTree_EGamma0_2025_Run2025D_0.root": 25.29/4,
"TnPTree_EGamma0_2025_Run2025E_0.root": 14.0/4,
"TnPTree_EGamma0_2025_Run2025F_0.root": 16.11/4,
"TnPTree_EGamma1_2025_Run2025C_1.root": 20.78/4,
"TnPTree_EGamma1_2025_Run2025D_1.root": 25.29/4,
"TnPTree_EGamma1_2025_Run2025E_1.root": 14.0/4,
"TnPTree_EGamma1_2025_Run2025F_1.root": 16.11/4,
"TnPTree_EGamma2_2025_Run2025C_2.root": 20.78/4,
"TnPTree_EGamma2_2025_Run2025D_2.root": 25.29/4,
"TnPTree_EGamma2_2025_Run2025E_2.root": 14.0/4,
"TnPTree_EGamma2_2025_Run2025F_2.root": 16.11/4,
"TnPTree_EGamma3_2025_Run2025C_3.root": 20.78/4,
"TnPTree_EGamma3_2025_Run2025D_3.root": 25.29/4,
"TnPTree_EGamma3_2025_Run2025E_3.root": 14.0/4,
"TnPTree_EGamma3_2025_Run2025F_3.root": 16.11/4,
}
for filename, lumi in dict_2025.items():
    # sed = 'seeded' if DO_SEEDED else 'unseeded'
    sample_name='data_' + filename.replace('TnPTree_', '').replace('.root', '')
    samplesDef['data'].add_sample( tnpSamples.Prompt2025[sample_name] )



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
    samplesDef['tagSel'].set_cut('tag_Ele_pt > 37')

## set MC weight, simple way (use tree weight) 
#weightName = 'totWeight'
#if not samplesDef['mcNom' ] is None: samplesDef['mcNom' ].set_weight(weightName)
#if not samplesDef['mcAlt' ] is None: samplesDef['mcAlt' ].set_weight(weightName)
#if not samplesDef['tagSel'] is None: samplesDef['tagSel'].set_weight(weightName)

## set MC weight, can use several pileup rw for different data taking 

puFile = '/eos/cms/store/group/phys_higgs/nonresonant_HH/PrivateProd/Yuxiang/TriggerSF/Pileup/DY_madgraph_ele.pu.puTree.root'
weightName = 'weights_2025_CDEF.totWeight'
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
   { 'var' : 'el_sc_eta' , 'type': 'float', 'bins': [-2.5,-2.0,-1.566,-1.4442, -0.8, 0.0, 0.8, 1.4442, 1.566, 2.0, 2.5] },
   { 'var' : 'el_et' , 'type': 'float', 'bins': [20,35,50,100,200,500] },
   { 'var' : 'el_r9', 'type': 'float', 'bins': [0.,0.5,0.56,0.8,0.85,0.9,1.05,2.0]},
#    { 'var' : 'el_r9', 'type': 'float', 'bins': [0.,0.5,0.55,0.6,0.65,0.7,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.05,2.0]},
]

less_cut = [
    { 'var' : 'abs(el_sc_eta)' , 'type': 'float', 'bins': [(0,1.446),(1.556001,2.5)] }, # the defalut left side is >=, so using 1.556001 to avoid =
    { 'var' : 'el_et' , 'type': 'float', 'bins': [20,35,50,100,200,500] },
    { 'var' : 'el_r9', 'type': 'float', 'bins': [0.,0.5,0.56,0.8,0.85,0.9,1.05,2.0]},
]

if DO_SEEDED:
    biningDef = {
        'EB': [
            { 'var': 'abs(el_sc_eta)', 'type': 'float', 'bins': [(0, 1.4442)] },
            { 'var': 'el_r9', 'type': 'float', 'bins': [
                (0.5, 0.56),
                (0.56, 0.85),
                (0.85, None)
            ]},
            { 'cor_var': 'el_et', 'type': 'float', 'bins': [0.0, 35.0, 37.0, 40.0, 45.0, 50.0, 60.0, 70.0, 90.0, 999999.0] },
        ],
        'EE': [
            { 'var': 'abs(el_sc_eta)', 'type': 'float', 'bins': [(1.566, 2.5)] },
            { 'var': 'el_r9', 'type': 'float', 'bins': [
                (0.8, 0.85),
                (0.85, 0.90),
                (0.90, None)
            ]},
            { 'cor_var': 'el_et', 'type': 'float', 'bins': [0.0, 35.0, 37.0, 40.0, 45.0, 50.0, 60.0, 70.0, 90.0, 999999.0] },
        ]
    }
else:
    biningDef = {
        'EB': [
            { 'var': 'abs(el_sc_eta)', 'type': 'float', 'bins': [(0, 1.4442)] },
            { 'var': 'el_r9', 'type': 'float', 'bins': [
                (0.5, 0.56),
                (0.56, 0.85),
                (0.85, None)
            ]},
            { 'cor_var': 'el_et', 'type': 'float', 'bins': [0.0, 25.0, 28.0, 31.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 90.0, 999999.0] },
        ],
        'EE': [
            { 'var': 'abs(el_sc_eta)', 'type': 'float', 'bins': [(1.566, 2.5)] },
            { 'var': 'el_r9', 'type': 'float', 'bins': [
                (0.8, 0.85),
                (0.85, 0.90),
                (0.90, None)
            ]},
            { 'cor_var': 'el_et', 'type': 'float', 'bins': [0.0, 25.0, 28.0, 31.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 90.0, 999999.0] },
        ]
    }


#############################################################
########## Cuts definition for all samples
#############################################################
### cut
cutBase   = 'abs(tag_sc_eta) < 2.1 && (abs(tag_sc_eta) < 1.4442 || abs(tag_sc_eta) > 1.566) && (tag_Ele_Iso122X > -0.5) && (el_IsoMVA_RunIIIWinter22 > -0.9)'

if not DO_SEEDED:
    cutBase += ' && tagMatchSeededLeg'

cor_cutBase = 'tag_Ele_pt > 35'
if DO_SEEDED:
    cor_cutBase += ' && el_pt > 30'
    
# can add addtionnal cuts for some bins (first check bin number using tnpEGM --checkBins)
additionalCuts = { 
#    0 : 'sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
#    1 : 'sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
#    2 : 'sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
#    3 : 'sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
#    4 : 'sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
#    5 : 'sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
#    6 : 'sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
#    7 : 'sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
#    8 : 'sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
#    9 : 'sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
}

#### or remove any additional cut (default)
#additionalCuts = None

#############################################################
########## fitting params to tune fit by hand if necessary
#############################################################
tnpParNomFit = [
    "meanP[-0.0,-5.0,5.0]","sigmaP[0.9,0.5,5.0]",
    "meanF[-0.0,-5.0,5.0]","sigmaF[0.9,0.5,5.0]",
    "acmsP[60.,40.,100.]","betaP[0.05,0.01,0.08]","gammaP[0.1, -0.1, 0.5]","peakP[87.0,82.0,90.0]",
    "acmsF[60.,40.,100.]","betaF[0.05,0.01,0.08]","gammaF[0.1, -0.1, 0.5]","peakF[87.0,82.0,90.0]",
    ]

tnpParAltSigFit = [
    "meanP[-0.0,-5.0,5.0]","sigmaP[1,0.7,6.0]","alphaP[2.0,1.2,3.5]" ,'nP[3,-5,5]',"sigmaP_2[1.5,0.5,6.0]","sosP[1,0.5,5.0]",
    "meanF[-0.0,-5.0,5.0]","sigmaF[2,0.7,15.0]","alphaF[2.0,1.2,3.5]",'nF[3,-5,5]',"sigmaF_2[2.0,0.5,6.0]","sosF[1,0.5,5.0]",
    "acmsP[60.,50.,75.]","betaP[0.04,0.01,0.06]","gammaP[0.1, 0.005, 0.5]","peakP[89.0,82.0,90.0]",
    "acmsF[60.,50.,75.]","betaF[0.04,0.01,0.06]","gammaF[0.1, 0.005, 0.5]","peakF[89.0,82.0,90.0]",
    ]
     
tnpParAltBkgFit = [
    "meanP[-0.0,-5.0,5.0]","sigmaP[0.9,0.5,5.0]",
    "meanF[-0.0,-5.0,5.0]","sigmaF[0.9,0.5,5.0]",
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
tnpParAltSigFitLooser = make_looser_ranges(tnpParAltSigFit)
tnpParAltBkgFitLooser = make_looser_ranges(tnpParAltBkgFit)