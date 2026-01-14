import ROOT as rt
rt.gROOT.LoadMacro('./libCpp/histFitter.C+')
rt.gROOT.LoadMacro('./libCpp/RooCBExGaussShape.cc+')
rt.gROOT.LoadMacro('./libCpp/RooCMSShape.cc+')
rt.gROOT.LoadMacro('./libCpp/RooTwoSigmaDSCBShape.cc+')
rt.gROOT.SetBatch(1)

from ROOT import tnpFitter

import re
import math


minPtForSwitch = 70

def ptMin( tnpBin ):
    ptmin = 1
    if tnpBin['name'].find('pt_') >= 0:
        ptmin = float(tnpBin['name'].split('pt_')[1].split('p')[0])
    elif tnpBin['name'].find('et_') >= 0:
        ptmin = float(tnpBin['name'].split('et_')[1].split('p')[0])
    return ptmin

def createWorkspaceForAltSig( sample, tnpBin, tnpWorkspaceParam, tailleft=0, useDSCB=False ):

    ### tricky: use n < 0 for high pT bin (so need to remove param and add it back)
    cbNList = ['tailLeft']
    ptmin = ptMin(tnpBin)       
    def _setting_tailright():
        print('--- Using Right tail for fitting ---')
        for par in cbNList:
            for ip in range(len(tnpWorkspaceParam)):
                x=re.compile('%s.*?' % par)
                listToRM = list(filter(x.match, tnpWorkspaceParam))
                for ir in listToRM :
                    print('**** remove', ir)
                    tnpWorkspaceParam.remove(ir)                    
            tnpWorkspaceParam.append( 'tailLeft[-1]' )
    
    if not useDSCB:
        if tailleft==0:
            if ptmin >= 35 :
                _setting_tailright()
        elif tailleft==-1:
            _setting_tailright()

    if sample.isMC:
        return tnpWorkspaceParam

    
    fileref = sample.mcRef.altSigFit
    filemc  = rt.TFile(fileref,'read')

    from ROOT import RooFit,RooFitResult
    fitresP = filemc.Get( '%s_resP' % tnpBin['name']  )
    fitresF = filemc.Get( '%s_resF' % tnpBin['name'] )

    listOfParam = ['nF','alphaF','nP','alphaP','sigmaP','sigmaF','sigmaP_2','sigmaF_2','meanGF','sigmaGF', 'sigFracF']
    if useDSCB:
        listOfParam += ['nF_2','alphaF_2','nP_2','alphaP_2']
    
    fitPar = fitresF.floatParsFinal()
    for ipar in range(len(fitPar)):
        pName = fitPar[ipar].GetName()
        print('%s[%2.3f]' % (pName,fitPar[ipar].getVal()))
        for par in listOfParam:
            if pName == par:
                x=re.compile('%s.*?' % pName)
                listToRM = list(filter(x.match, tnpWorkspaceParam))
                for ir in listToRM :
                    tnpWorkspaceParam.remove(ir)                    
                tnpWorkspaceParam.append( '%s[%2.3f]' % (pName,fitPar[ipar].getVal()) )
                              
  
    fitPar = fitresP.floatParsFinal()
    for ipar in range(len(fitPar)):
        pName = fitPar[ipar].GetName()
        print('%s[%2.3f]' % (pName,fitPar[ipar].getVal()))
        for par in listOfParam:
            if pName == par:
                x=re.compile('%s.*?' % pName)
                listToRM = list(filter(x.match, tnpWorkspaceParam))
                for ir in listToRM :
                    tnpWorkspaceParam.remove(ir)
                tnpWorkspaceParam.append( '%s[%2.3f]' % (pName,fitPar[ipar].getVal()) )

    filemc.Close()

    return tnpWorkspaceParam


#############################################################
########## nominal fitter
#############################################################
def getHistEntries(histFile, binName, category):
    """
    Get the number of entries in a histogram
    
    Args:
        histFile: Path to ROOT file
        binName: Name of the bin
        category: 'Pass' or 'Fail'
    
    Returns:
        Number of entries in the histogram, 0 if error
    """
    try:
        infile = rt.TFile(histFile, "read")
        hist = infile.Get(f'{binName}_{category}')
        entries = hist.GetEntries() if hist else 0
        infile.Close()
        return entries
    except:
        return 0

def selectZLineShape(sample, tnpBin, minEntries=500, useGenZee=False):
    """
    Intelligently select Z line shape based on histogram entries and pT threshold
    
    Priority logic:
    1. If both Pass and Fail entries < minEntries: Use ZeeGenLevel (most critical)
    2. If Pass entries < minEntries: Use ZeeGenLevel (Pass is essential)
    3. If Fail entries < minEntries AND pT > minPtForSwitch: Use Pass for Fail (high pT + low Fail stats)
    4. If Fail entries < minEntries: Use Pass for Fail (only Fail has low stats)
    5. If pT > minPtForSwitch: Use Pass for Fail (high pT fallback, but both have enough stats)
    6. Otherwise: Use respective histograms normally
    
    Args:
        sample: Sample object containing MC reference
        tnpBin: TnP bin configuration
        minEntries: Minimum entries threshold (default: 500)
    
    Returns:
        tuple: (histZLineShapeP, histZLineShapeF, fileHandle)
    """
    binName = tnpBin['name']
    ptmin = ptMin(tnpBin)
    
    # Get entries for Pass and Fail categories
    entriesP = getHistEntries(sample.mcRef.histFile, binName, 'Pass')
    entriesF = getHistEntries(sample.mcRef.histFile, binName, 'Fail')
    
    print(f'--- Z LineShape Selection ---')
    print(f'Bin: {binName}, pT_min: {ptmin} GeV')
    print(f'Pass entries: {entriesP}, Fail entries: {entriesF}, Threshold: {minEntries}')

    if useGenZee:
        print('Using ZeeGenLevel as per configuration')
        fileGenLevel = rt.TFile('etc/inputs/ZeeGenLevel.root', 'read')
        histZLineShape = fileGenLevel.Get('Mass')
        return histZLineShape, histZLineShape, fileGenLevel
    
    fileTruth = rt.TFile(sample.mcRef.histFile, 'read')
    
    # Priority 1: Both Pass and Fail have insufficient statistics -> Use ZeeGenLevel
    if entriesP < minEntries and entriesF < minEntries:
        print('Both Pass and Fail entries below threshold, switching to ZeeGenLevel')
        fileTruth.Close()
        fileGenLevel = rt.TFile('etc/inputs/ZeeGenLevel.root', 'read')
        histZLineShape = fileGenLevel.Get('Mass')
        return histZLineShape, histZLineShape, fileGenLevel
    
    # Priority 2: Pass has insufficient statistics -> Use ZeeGenLevel
    # (Pass is critical, cannot use Fail for Pass)
    elif entriesP < minEntries:
        print('Pass entries below threshold, switching to ZeeGenLevel')
        fileTruth.Close()
        fileGenLevel = rt.TFile('etc/inputs/ZeeGenLevel.root', 'read')
        histZLineShape = fileGenLevel.Get('Mass')
        return histZLineShape, histZLineShape, fileGenLevel
    
    # At this point, Pass has sufficient statistics
    histZLineShapeP = fileTruth.Get(f'{binName}_Pass')
    
    # Priority 3: Fail has insufficient statistics -> Use Pass for Fail
    if entriesF < minEntries:
        print('Fail entries below threshold, using Pass histogram for Fail')
        histZLineShapeF = fileTruth.Get(f'{binName}_Pass')
        return histZLineShapeP, histZLineShapeF, fileTruth
    
    # Priority 4: High pT bin (both have sufficient statistics) -> Use Pass for Fail
    elif ptmin > minPtForSwitch:
        print(f'High pT bin (pT > {minPtForSwitch} GeV) with sufficient statistics, using Pass histogram for Fail')
        histZLineShapeF = fileTruth.Get(f'{binName}_Pass')
        return histZLineShapeP, histZLineShapeF, fileTruth
    
    # Priority 5: Normal case -> Use respective histograms
    else:
        print('Using normal Pass and Fail histograms')
        histZLineShapeF = fileTruth.Get(f'{binName}_Fail')
        return histZLineShapeP, histZLineShapeF, fileTruth
    
def estimatedMassCut(tnpBin, TagPtOffset=35):
    ptmin = ptMin(tnpBin)
    estimatedMass = ptmin + TagPtOffset
    if estimatedMass < 60:
        estimatedMass = 60
    if estimatedMass > 80:
        estimatedMass = 80
    return estimatedMass

def histFitterNominal( sample, tnpBin, tnpWorkspaceParam, minEntries=500, useGenZee=False ):
        
    tnpWorkspaceFunc = [
        "Gaussian::sigResPass(x,meanP,sigmaP)",
        "Gaussian::sigResFail(x,meanF,sigmaF)",
        "RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)",
        "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
        ]

    tnpWorkspace = []
    tnpWorkspace.extend(tnpWorkspaceParam)
    tnpWorkspace.extend(tnpWorkspaceFunc)
    
    ## init fitter
    infile = rt.TFile( sample.histFile, "read")
    hP = infile.Get('%s_Pass' % tnpBin['name'] )
    hF = infile.Get('%s_Fail' % tnpBin['name'] )
    fitter = tnpFitter( hP, hF, tnpBin['name'] )
    infile.Close()

    ## setup
    fitter.useMinos()
    rootpath = sample.nominalFit.replace('.root', '-%s.root' % tnpBin['name'])
    rootfile = rt.TFile(rootpath,'update')
    fitter.setOutputFile( rootfile )
    
    ## generated Z LineShape
    ## Apply intelligent switching only for data (nominal fit only applies to data)
    histZLineShapeP, histZLineShapeF, fileTruth = selectZLineShape(sample, tnpBin, minEntries, useGenZee=useGenZee)
    estimatedMassCutValue = estimatedMassCut(tnpBin)
    fitter.setZLineShapes(histZLineShapeP, histZLineShapeF, True, estimatedMassCutValue)
    # check the fail numbers to decide whether to fix the bkgPassToFail ratio
    failEntries = getHistEntries(sample.histFile, tnpBin['name'], 'Fail')

    not_fixBkgPassToFail = False

    not_fixBkgPassToFail = True

    if failEntries > minEntries and not not_fixBkgPassToFail:
        print(f'---- Fixing bkg parameters to Fail for Pass fitting in Bin {tnpBin["name"]} ----')
        fitter.fixBkgPassToFail()
    fileTruth.Close()

    ### set workspace
    workspace = rt.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace )

    title = tnpBin['title'].replace(';',' - ')
    title = title.replace('probe_sc_eta','#eta_{SC}')
    title = title.replace('probe_Ele_pt','p_{T}')
    fitter.fits(sample.mcTruth,sample.isMC,title)
    rootfile.Close()



#############################################################
########## alternate signal fitter
#############################################################
def histFitterAltSig( sample, tnpBin, tnpWorkspaceParam, isaddGaus=0, useDSCB=False, doMcPreFit=True ):

    tailLeft = 0
    if doMcPreFit:
        if sample.isMC:
            print('---- Doing MC pre-fit to get signal parameters ----')
            if useDSCB:
                from . import mcPreFit_DCB
                newParam, setting = mcPreFit_DCB.perform_pre_AltSigFit( sample, tnpBin, tnpWorkspaceParam, do_plot=True )
                mcFitLow = mcPreFit_DCB.extract_fit_range_from_bin(tnpBin.get('name',''))
            else:
                from . import mcPreFit
                newParam, setting = mcPreFit.perform_pre_AltSigFit( sample, tnpBin, tnpWorkspaceParam, do_plot=True )
            tnpWorkspaceParam = newParam
            print('---- Using MC pre-fit settings ----')
            print(tnpWorkspaceParam)
            tailLeft = setting.get('tailLeft', 0)
        elif not useDSCB: # only in data and 
            import os
            pathToCheck = os.path.join(os.path.dirname(sample.mcRef.histFile), 'python_prefit', sample.mcRef.name, 'result')
            tnpBin_name = tnpBin['name']
            tailLeft = 0
            if os.path.exists(os.path.join(pathToCheck, f"{tnpBin_name}_tailLeft")):
                tailLeft = 1
            elif os.path.exists(os.path.join(pathToCheck, f"{tnpBin_name}_tailRight")):
                tailLeft = -1
            if tailLeft != 0:
                print(f'---- Using MC pre-fit settings from {pathToCheck}, tailLeft={tailLeft} ----')
    

    tnpWorkspacePar = createWorkspaceForAltSig( sample,  tnpBin, tnpWorkspaceParam, tailleft=tailLeft, useDSCB=useDSCB )

    if useDSCB:
        print('--- Using Double Crystal Ball + Exponential for fitting ---')
        tnpWorkspaceFunc = [
            "RooTwoSigmaDSCBShape::sigResPass(x,meanP,expr('sqrt(sigmaP*sigmaP+sosP*sosP)',{sigmaP,sosP}),expr('sqrt(sigmaP_2*sigmaP_2+sosP*sosP)',{sigmaP_2,sosP}),alphaP,nP,alphaP_2,nP_2)",
            "RooTwoSigmaDSCBShape::sigResFail(x,meanF,expr('sqrt(sigmaF*sigmaF+sosF*sosF)',{sigmaF,sosF}),expr('sqrt(sigmaF_2*sigmaF_2+sosF*sosF)',{sigmaF_2,sosF}),alphaF,nF,alphaF_2,nF_2)",
            "RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)",
            "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
        ]
    else:
        tnpWorkspaceFunc = [
            "tailLeft[1]",
            "RooCBExGaussShape::sigResPass(x,meanP,expr('sqrt(sigmaP*sigmaP+sosP*sosP)',{sigmaP,sosP}),alphaP,nP, expr('sqrt(sigmaP_2*sigmaP_2+sosP*sosP)',{sigmaP_2,sosP}),tailLeft)",
            "RooCBExGaussShape::sigResFail(x,meanF,expr('sqrt(sigmaF*sigmaF+sosF*sosF)',{sigmaF,sosF}),alphaF,nF, expr('sqrt(sigmaF_2*sigmaF_2+sosF*sosF)',{sigmaF_2,sosF}),tailLeft)",
            "RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)",
            "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
            ]
    if isaddGaus==1:
        tnpWorkspaceFunc += [ "Gaussian::sigGaussFail(x,meanGF,sigmaGF)", ]
        if sample.isMC:
            tnpWorkspaceFunc += [ "sigFracF[0.5,0.0,1.0]", ]

    tnpWorkspace = []
    tnpWorkspace.extend(tnpWorkspacePar)
    tnpWorkspace.extend(tnpWorkspaceFunc)

    print('-'*30)
    print(tnpWorkspace)
        
    ## init fitter
    infile = rt.TFile( sample.histFile, "read")
    hP = infile.Get('%s_Pass' % tnpBin['name'] )
    hF = infile.Get('%s_Fail' % tnpBin['name'] )
    ## for high pT change the failing spectra to passing probe to get statistics 
    ## MC only: this is to get MC parameters in data fit!
    if sample.isMC and ptMin( tnpBin ) > minPtForSwitch:     
        hF = infile.Get('%s_Pass' % tnpBin['name'] )
    if sample.isMC and useDSCB:
        fitter = tnpFitter( hP, hF, tnpBin['name'], mcFitLow)
    else:
        fitter = tnpFitter( hP, hF, tnpBin['name'] )
#    fitter.fixSigmaFtoSigmaP()
    infile.Close()

    ## setup
    rootpath = sample.altSigFit.replace('.root', '-%s.root' % tnpBin['name'])
    rootfile = rt.TFile(rootpath,'update')
    fitter.setOutputFile( rootfile )
    
    ## generated Z LineShape
    fileTruth = rt.TFile('etc/inputs/ZeeGenLevel.root','read')
    histZLineShape = fileTruth.Get('Mass')
    estimatedMassCutValue = estimatedMassCut(tnpBin)
    fitter.setZLineShapes(histZLineShape,histZLineShape, True, estimatedMassCutValue)
    fileTruth.Close()

    ### set workspace
    workspace = rt.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    print(workspace)
    fitter.setWorkspace( workspace, isaddGaus )

    title = tnpBin['title'].replace(';',' - ')
    title = title.replace('probe_sc_eta','#eta_{SC}')
    title = title.replace('probe_Ele_pt','p_{T}')
    fitter.fits(sample.mcTruth,sample.isMC,title, isaddGaus)

    rootfile.Close()



#############################################################
########## alternate background fitter
#############################################################
def histFitterAltBkg( sample, tnpBin, tnpWorkspaceParam, minEntries=500 ):

    tnpWorkspaceFunc = [
        "Gaussian::sigResPass(x,meanP,sigmaP)",
        "Gaussian::sigResFail(x,meanF,sigmaF)",
        "Exponential::bkgPass(x, alphaP)",
        "Exponential::bkgFail(x, alphaF)",
        ]

    tnpWorkspace = []
    tnpWorkspace.extend(tnpWorkspaceParam)
    tnpWorkspace.extend(tnpWorkspaceFunc)
            
    ## init fitter
    infile = rt.TFile(sample.histFile,'read')
    hP = infile.Get('%s_Pass' % tnpBin['name'] )
    hF = infile.Get('%s_Fail' % tnpBin['name'] )
    fitter = tnpFitter( hP, hF, tnpBin['name'] )
    infile.Close()

    ## setup
    rootpath = sample.altBkgFit.replace('.root', '-%s.root' % tnpBin['name'])
    rootfile = rt.TFile(rootpath,'update')
    fitter.setOutputFile( rootfile )

    ## generated Z LineShape
    ## Apply intelligent switching only for data (altBkg fit only applies to data)
    if not sample.isMC:
        histZLineShapeP, histZLineShapeF, fileTruth = selectZLineShape(sample, tnpBin, minEntries)
        estimatedMassCutValue = estimatedMassCut(tnpBin)
        fitter.setZLineShapes(histZLineShapeP, histZLineShapeF, True, estimatedMassCutValue)
        fileTruth.Close()
    else:
        # For MC, use standard approach
        fileTruth = rt.TFile(sample.mcRef.histFile,'read')
        histZLineShapeP = fileTruth.Get('%s_Pass'%tnpBin['name'])
        histZLineShapeF = fileTruth.Get('%s_Fail'%tnpBin['name'])
        if ptMin( tnpBin ) > minPtForSwitch: 
            histZLineShapeF = fileTruth.Get('%s_Pass'%tnpBin['name'])
        fitter.setZLineShapes(histZLineShapeP, histZLineShapeF)
        fileTruth.Close()

    ### set workspace
    workspace = rt.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace )

    title = tnpBin['title'].replace(';',' - ')
    title = title.replace('probe_sc_eta','#eta_{SC}')
    title = title.replace('probe_Ele_pt','p_{T}')
    fitter.fits(sample.mcTruth,sample.isMC,title)
    rootfile.Close()


