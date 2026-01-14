### python specific import
import argparse
import os
import sys
import pickle
import shutil
import multiprocessing as mp
import datetime


parser = argparse.ArgumentParser(description='tnp EGM fitter')
parser.add_argument('--checkBins'  , action='store_true'  , help = 'check  bining definition')
parser.add_argument('--createBins' , action='store_true'  , help = 'create bining definition')
parser.add_argument('--createHists', action='store_true'  , help = 'create histograms')
parser.add_argument('--sample'     , default='all'        , help = 'create histograms (per sample, expert only)')
parser.add_argument('--altSig'     , action='store_true'  , help = 'alternate signal model fit')
parser.add_argument('--addGaus'    , action='store_true'  , help = 'add gaussian to alternate signal model failing probe')
parser.add_argument('--altBkg'     , action='store_true'  , help = 'alternate background model fit')
parser.add_argument('--doFit'      , action='store_true'  , help = 'fit sample (sample should be defined in settings.py)')
parser.add_argument('--mcSig'      , action='store_true'  , help = 'fit MC nom [to init fit parama]')
parser.add_argument('--doPlot'     , action='store_true'  , help = 'plotting')
parser.add_argument('--sumUp'      , action='store_true'  , help = 'sum up efficiencies')
parser.add_argument('--iBin'       , dest = 'binNumber'   , type = str,  default='-1', help='bin number (to refit individual bin), support multiple bins separated by comma, e.g. "0,1,2" or range "0-5"')
parser.add_argument('--flag'       , default = None       , help ='WP to test')
parser.add_argument('settings'     , default = None       , help = 'setting file [mandatory]')

parser.add_argument('--doCor'      , dest = 'doCorrection', action='store_true'  , help = 'do scale and smearing correction')
parser.add_argument('--year'       , default = '2024', help = 'year config used for correction' )
parser.add_argument('--useGenZee'  , action='store_true'  , help = 'use ZeeGenLevel for Z line shape in nominal fit')


args = parser.parse_args()

import types
tnp_config_module = types.ModuleType('__tnp_config__')
tnp_config_module.used_flag = args.flag
sys.modules['__tnp_config__'] = tnp_config_module

print('===> settings %s <===' % args.settings)
importSetting = 'import %s as tnpConf' % args.settings.replace('/','.').split('.py')[0]
print(importSetting)
exec(importSetting)

### tnp library
import libPython.binUtils  as tnpBiner
import libPython.rootUtils as tnpRoot


if args.flag is None:
    print('[tnpEGM_fitter] flag is MANDATORY, this is the working point as defined in the settings.py')
    sys.exit(0)
    
if not args.flag in tnpConf.flags.keys() :
    print('[tnpEGM_fitter] flag %s not found in flags definitions' % args.flag)
    print('  --> define in settings first')
    print('  In settings I found flags: ')
    print(tnpConf.flags.keys())
    sys.exit(1)

outputDirectory = '%s/%s/' % (tnpConf.baseOutDir,args.flag)

print('===>  Output directory: ')
print(outputDirectory)


####################################################################
##### Create (check) Bins
####################################################################
if args.checkBins:
    tnpBins = tnpBiner.createBins(tnpConf.biningDef,tnpConf.cutBase,tnpConf.cor_cutBase)
    tnpBiner.tuneCuts( tnpBins, tnpConf.additionalCuts )
    for ib in range(len(tnpBins['bins'])):
        print(tnpBins['bins'][ib]['name'])
        print('  - cut: ',tnpBins['bins'][ib]['cut'])
    sys.exit(0)
    
if args.createBins:
    if os.path.exists( outputDirectory ):
            shutil.rmtree( outputDirectory )
    os.makedirs( outputDirectory )
    tnpBins = tnpBiner.createBins(tnpConf.biningDef,tnpConf.cutBase,tnpConf.cor_cutBase)
    tnpBiner.tuneCuts( tnpBins, tnpConf.additionalCuts )
    pickle.dump( tnpBins, open( '%s/bining.pkl'%(outputDirectory),'wb') )
    print('created dir: %s ' % outputDirectory)
    print('bining created successfully... ')
    print('Note than any additional call to createBins will overwrite directory %s' % outputDirectory)
    sys.exit(0)

tnpBins = pickle.load( open( '%s/bining.pkl'%(outputDirectory),'rb') )


####################################################################
##### Create Histograms
####################################################################
for s in tnpConf.samplesDef.keys():
    sample =  tnpConf.samplesDef[s]
    if sample is None: continue
    setattr( sample, 'tree'     ,'%s/fitter_tree' % tnpConf.tnpTreeDir )
    setattr( sample, 'histFile' , '%s/%s_%s.root' % ( outputDirectory , sample.name, args.flag ) )


if args.createHists:

    print(" ======== Creating Histograms (Final Precise Configuration) ========")
    import libPython.histUtils as tnpHist
    import copy

    do_correction = args.doCorrection
    correction_year = args.year
    if do_correction:
        print(f'apply scale and smearing correction using {correction_year} config')

    def tnpBins_converter(obj):
        """
        A highly specific converter based on the exact needs of histUtils.pyx.
        - Keeps all keys as strings.
        - Converts values of 'name' and 'title' keys to bytes.
        - Keeps values of 'cut' keys as strings.
        - Recursively applies this logic.
        """
        if isinstance(obj, dict):
            new_dict = {}
            for k, v in obj.items():
                if k in ['name', 'title'] and isinstance(v, str):
                    new_dict[k] = v.encode('utf-8')
                elif k in ['cut', 'baseSelection'] and isinstance(v, str):
                    new_dict[k] = v
                else:
                    new_dict[k] = tnpBins_converter(v)
            return new_dict
        elif isinstance(obj, list):
            return [tnpBins_converter(elem) for elem in obj]
        elif isinstance(obj, str):
            return obj.encode('utf-8')
        return obj

    branch_mapping = tnpConf.branch_mapping if hasattr(tnpConf, 'branch_mapping') else None
    r9Eta_reweighting = tnpConf.r9Eta_reweighting if hasattr(tnpConf, 'r9Eta_reweighting') else None
    if r9Eta_reweighting:
        print(f'apply R9-Eta reweighting using {r9Eta_reweighting} file')

    tnpBins_to_pass = copy.deepcopy(tnpBins)
    tnpBins_to_pass = tnpBins_converter(tnpBins_to_pass)
    def parallel_hists(sampleType):
        sample =  tnpConf.samplesDef[sampleType]
        if sample is None : return
        if sampleType == args.sample or args.sample == 'all' :
            print('creating histogram for sample ')
            sample.dump()
            var = { 'name' : b'pair_mass', 'nbins' : 80, 'min' : 50, 'max': 130 }
            if sample.mcTruth:
                var = { 'name' : b'pair_mass', 'nbins' : 80, 'min' : 50, 'max': 130 }
            # 1. `sample.tree` must be bytes
            if hasattr(sample, 'tree') and isinstance(getattr(sample, 'tree'), str):
                setattr(sample, 'tree', getattr(sample, 'tree').encode('utf-8'))
            # 2. `tnpBins` must be converted
            tnpHist.makePassFailHistograms( sample, tnpConf.flags[args.flag], tnpBins_to_pass, var, do_correction, correction_year, branch_mapping, r9Eta_reweighting)#, max_event = 10000000)
    
    num_samples = len(list(tnpConf.samplesDef.keys()))
    num_processes = min(num_samples, mp.cpu_count())
    
    with mp.Pool(processes=num_processes) as pool:
        pool.map(parallel_hists, tnpConf.samplesDef.keys())
    # # debug for MC
    # parallel_hists('data')

    sys.exit(0)


####################################################################
##### Actual Fitter
####################################################################
sampleToFit = tnpConf.samplesDef['data']
if sampleToFit is None:
    print('[tnpEGM_fitter, prelim checks]: sample (data or MC) not available... check your settings')
    sys.exit(1)

sampleMC = tnpConf.samplesDef['mcNom']

if sampleMC is None:
    print('[tnpEGM_fitter, prelim checks]: MC sample not available... check your settings')
    sys.exit(1)
for s in tnpConf.samplesDef.keys():
    sample =  tnpConf.samplesDef[s]
    if sample is None: continue
    setattr( sample, 'mcRef'     , sampleMC )
    setattr( sample, 'nominalFit', '%s/%s_%s.nominalFit.root' % ( outputDirectory , sample.name, args.flag ) )
    setattr( sample, 'altSigFit' , '%s/%s_%s.altSigFit.root'  % ( outputDirectory , sample.name, args.flag ) )
    setattr( sample, 'altBkgFit' , '%s/%s_%s.altBkgFit.root'  % ( outputDirectory , sample.name, args.flag ) )



### change the sample to fit is mc fit
if args.mcSig :
    sampleToFit = tnpConf.samplesDef['mcNom']

test_bin_number = len(tnpBins['bins'])
test_list = range(test_bin_number)

# Parse bin numbers
def parse_bin_numbers(bin_str, max_bin):
    """Parse bin number string, support formats like '0,1,2' or '0-5' or '-1' for all bins"""
    bin_str = bin_str.strip()
    if bin_str == '-1':
        return list(range(max_bin))
    
    bins = set()
    for part in bin_str.split(','):
        part = part.strip()
        if '-' in part and part != '-1':
            # Range format: "0-5"
            start, end = part.split('-')
            bins.update(range(int(start), int(end) + 1))
        else:
            # Single bin
            bins.add(int(part))
    
    # Filter valid bins
    return sorted([b for b in bins if 0 <= b < max_bin])

selected_bins = parse_bin_numbers(args.binNumber, test_bin_number)
print(f'Selected bins to process: {selected_bins}')

if  args.doFit:
    import glob
    print(" ======== Fitting ========")
    sampleToFit.dump()
    def parallel_fit(ib):
        # delete the previous fit file if exists
        file_name = sampleToFit.nominalFit
        if args.altSig:
            file_name = sampleToFit.altSigFit
        if args.altBkg:
            file_name = sampleToFit.altBkgFit
        files = glob.glob( file_name.replace('.root', f"-*{tnpBins['bins'][ib]['name']}.root") )
        print(f'[WARNING] removing previous fit files for bin {ib}: {files}')
        os.system(f'rm -f {files}')

        is_looser = len(selected_bins) < test_bin_number
        if is_looser: print(f'using looser fit for bin {ib}')
        
        curr_bin = tnpBins['bins'][ib]
        fit_func = None
        fit_args = [sampleToFit, curr_bin]
        fit_kwargs = {}

        # determine fit function and parameters
        if args.altSig:
            fit_func = tnpRoot.histFitterAltSig
            use_dscb = hasattr(tnpConf, 'AltSigFitUsingDSCB') and tnpConf.tnpParAltSigFit_DSCB
            fit_kwargs['useDSCB'] = use_dscb
            if use_dscb:
                params = "tnpParAltSigFit_DSCB"
            else: 
                params = "tnpParAltSigFit"
        elif args.altBkg:
            fit_func = tnpRoot.histFitterAltBkg
            params = "tnpParAltBkgFit"
        else:
            fit_func = tnpRoot.histFitterNominal
            params = "tnpParNomFit"
            fit_kwargs['useGenZee'] = args.useGenZee
        
        if is_looser:
            params = params + "Looser"
        
        if args.altSig and args.addGaus:
            params += "_addGaus"
            fit_kwargs['isaddGaus'] = True
        
        fit_args.append( getattr(tnpConf, params) )

        # perform fit
        if fit_func:
            fit_func(*fit_args, **fit_kwargs)

    timeout_seconds = 10000 if len(selected_bins) < test_bin_number else 2000
    timeout_log = os.path.join(outputDirectory, 'timeout.txt')
    
    tasks_to_run = []
    for ib in test_list:
        if ib in selected_bins:
            tasks_to_run.append(ib)
    import time
    max_procs = mp.cpu_count()
    running_procs = [] 
    
    task_idx = 0
    total_tasks = len(tasks_to_run)

    while task_idx < total_tasks or len(running_procs) > 0:
        while len(running_procs) < max_procs and task_idx < total_tasks:
            ib = tasks_to_run[task_idx]
            p = mp.Process(target=parallel_fit, args=(ib,))
            p.start()
            running_procs.append((p, ib, time.time()))
            task_idx += 1
        
        active_procs = []
        for p, ib, start_t in running_procs:
            if not p.is_alive():
                p.join()
            
            elif (time.time() - start_t) > timeout_seconds:
                try:
                    print(f"[tnpEGM_fitter] Killing bin {ib} due to timeout ({timeout_seconds}s)")
                    p.terminate()
                    p.join(5)
                except Exception as e:
                    print(f"Error terminating process: {e}")
                
                # timeout record
                ts = datetime.datetime.now().isoformat()
                ft_param = 'altSig' if args.altSig else ('altBkg' if args.altBkg else 'nominal')
                info = f"{ts}\tbin={ib}\tflag={args.flag}\tsample={sampleToFit.name}\ttimeout={timeout_seconds}s\tfitParam={ft_param}\n"
                with open(timeout_log, 'a') as fout:
                    fout.write(info)
            
            else:
                active_procs.append((p, ib, start_t))
        
        running_procs = active_procs
        if len(running_procs) > 0:
            time.sleep(0.1)

    # parallel_fit(121)
    # exit()
    if len(selected_bins) < test_bin_number:
        print(f'fitting for bins {selected_bins} done')
        sys.exit(0)
    args.doPlot = True
####################################################################
##### dumping plots
####################################################################
if  args.doPlot:
    fileName = sampleToFit.nominalFit
    fitType  = 'nominalFit'
    if args.altSig : 
        fileName = sampleToFit.altSigFit
        fitType  = 'altSigFit'
    if args.altBkg : 
        fileName = sampleToFit.altBkgFit
        fitType  = 'altBkgFit'
    
    import glob
    to_merge_file_list = glob.glob( fileName.replace('.root', '-*.root') )
    # check file size, if file size < 5KB, consider it as failed fit and remove it from merging
    valid_file_list = []
    invalid_bin_list = []
    for f in to_merge_file_list:
        if os.path.getsize(f) > 5*1024:
            valid_file_list.append(f)
        else:
            print(f"[tnpEGM_fitter] removing invalid fit file {f} (size < 5KB)")
            # get bin number from file name
            # e.g. data_EGamma0_2025_Run2025C_0_unseeded_passingHLTUnseeded.nominalFit-bin349_el_sc_eta_2p00To2p50_el_et_200p00To500p00_el_r9_1p05To2p00
            bin_str = f.split('.root')[0].split('-bin')[-1].split('_')[0]
            invalid_bin_list.append(int(bin_str))
        
    os.system('hadd -f %s %s' % (fileName, ' '.join(valid_file_list)))

    plottingDir = '%s/plots/%s/%s' % (outputDirectory,sampleToFit.name,fitType)
    if not os.path.exists( plottingDir ):
        os.makedirs( plottingDir )
    shutil.copy('etc/inputs/index.php.listPlots','%s/index.php' % plottingDir)

    for ib in range(test_bin_number):#range(len(tnpBins['bins'])):
        if ib in selected_bins:
            if ib in invalid_bin_list:
                print(f"[tnpEGM_fitter] skipping plotting for invalid bin {ib}")
                continue
            tnpRoot.histPlotter( fileName, tnpBins['bins'][ib], plottingDir )

    print(' ===> Plots saved in <=======')
#    print 'localhost/%s/' % plottingDir
    exit(0)


####################################################################
##### dumping egamma txt file 
####################################################################
if args.sumUp:
    sampleToFit.dump()
    info = {
        'data'        : sampleToFit.histFile,
        'dataNominal' : sampleToFit.nominalFit,
        'dataAltSig'  : sampleToFit.altSigFit ,
        'dataAltBkg'  : sampleToFit.altBkgFit ,
        'mcNominal'   : sampleToFit.mcRef.histFile,
        'mcAlt'       : None,
        'tagSel'      : None
        }

    #if not tnpConf.samplesDef['mcAlt' ] is None:
    #    info['mcAlt'    ] = tnpConf.samplesDef['mcAlt' ].histFile
    if not tnpConf.samplesDef['tagSel'] is None:
        info['tagSel'   ] = tnpConf.samplesDef['tagSel'].histFile

    effis = None
    out_format = 'json'
    if out_format == 'txt':
        effFileName ='%s/egammaEffi.txt' % outputDirectory 
        fOut = open(effFileName, 'w')

        if len(tnpBins['bins']) > 0:
            n_vars = len(tnpBins['bins'][0]['title'].split(';'))
            var_headers = []
            for i in range(1, n_vars + 1):
                var_headers.extend([f'var{i}_low', f'var{i}_high'])
            eff_headers = [
                'dataNom_eff', 'dataNom_err',
                'mcNom_eff', 'mcNom_err',
                'dataAltBkg_eff',
                'dataAltSig_eff',
                'mcAlt_eff',
                'tagSel_eff'
            ]
            header_line = '\t'.join(var_headers + eff_headers)
            fOut.write('# ' + header_line + '\n')
        # ====================

        for ib in range(len(tnpBins['bins'])):
            print('Summing up bin %d / %d ' % (ib, len(tnpBins['bins']) ) )
            effis = tnpRoot.getAllEffi( info, tnpBins['bins'][ib] )

            ### formatting for arbitrary N-D binning
            title_parts = tnpBins['bins'][ib]['title'].split(';')
            n_vars = len(title_parts)
            
            var_ranges = []
            for i in range(0, n_vars ):
                var_range = title_parts[i].split('<')
                if '>' in title_parts[i]:
                    var_range = title_parts[i].replace('=','').split('>')
                    var_range = [var_range[1], 'xxx', 9999]
                var_ranges.append(var_range)
            
            if ib == 0:
                print(tnpBins['bins'][ib]['title'])
                for i, var_range in enumerate(var_ranges, 1):
                    astr = '### var%d : %s' % (i, var_range[1])
                    print(astr)
                    fOut.write(astr + '\n')

            format_parts = []
            values = []
            
            for var_range in var_ranges:
                format_parts.extend(['%+8.5f', '%+8.5f'])
                values.extend([float(var_range[0]), float(var_range[2])])
            
            # add efficiency values
            # value handel: if XXX[0] == 0.5 then set to 0 as 0.5 is default for almost zero stat
            for key in ['dataNominal', 'mcNominal', 'dataAltBkg', 'dataAltSig', 'mcAlt', 'tagSel']:
                if key in effis:
                    if abs(effis[key][0] - 0.50000) < 1e-5:
                        effis[key][0] = 0.0
            format_parts.extend(['%5.5f'] * 8)
            values.extend([
                effis['dataNominal'][0], effis['dataNominal'][1],
                effis['mcNominal'][0], effis['mcNominal'][1],
                effis['dataAltBkg'][0],
                effis['dataAltSig'][0],
                effis['mcAlt'][0],
                effis['tagSel'][0],
            ])
            astr = '\t'.join(format_parts) % tuple(values)
            print(astr)
            fOut.write(astr + '\n')

        fOut.close()

        print('Effis saved in file : ',  effFileName)
        import libPython.new_EGammaID_scaleFactors as egm_sf
        egm_sf.doEGM_SFs(effFileName,sampleToFit.lumi)
        exit(0)
    elif out_format == 'json':
        effFileName ='%s/egammaEffi.json' % outputDirectory 
        import json
        effis_all = {}
        for ib in range(len(tnpBins['bins'])):
            print('Summing up bin %d / %d ' % (ib, len(tnpBins['bins']) ) )
            effis = tnpRoot.getAllEffi( info, tnpBins['bins'][ib], out_level='json' )
            effis_all[tnpBins['bins'][ib]['name']] = {
                'tittle'      : tnpBins['bins'][ib]['title'],
                'def'         : tnpBins['bins'][ib]['vars'],
                'info' : effis
            }
        with open(effFileName, 'w') as fOut:
            json.dump( effis_all, fOut, indent=4)
        print('Effis saved in file : ',  effFileName)
        import libPython.new_EGammaID_scaleFactors as egm_sf
        egm_sf.doEGM_SFs(effFileName,sampleToFit.lumi)
        exit(0)
