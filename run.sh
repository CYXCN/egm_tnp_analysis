#!/bin/bash
echo $PWD
echo "Setting environment"

# using cmssw-el7 container

source /cvmfs/cms.cern.ch/cmsset_default.sh

# for condor
# cd /eos/user/y/yucao/program/trigger/CMSSW_11_2_0/src/egm_tnp_analysis/
cmsenv

# Uncomment the following lines if you need to set up the CMSSW environment
# workdir=$PWD;
# cmssw_path="$CMSSW_BASE/src";
# cd $cmssw_path; cmsenv; cd $workdir
# export PYTHONPATH=$workdir;
# make

# --- Arguments ---
# $1: config file (e.g., etc/config/settings_pho_run2022FG.py)
# $2: ID flag
# $3: run_mode: 0=local, 1=prepare condor, 2=run on condor, 3=collect from condor
# $4: iBin: -1 for all bins, >=0 for a specific bin (used for condor)

CONFIG=$1
FLAG=$2
RUN_MODE=${3:-0} # Default to 0 (local) if not provided
IBIN=${4:--1}    # Default to -1 if not provided

echo "Config: $CONFIG"
echo "Flag: $FLAG"
echo "Run Mode: $RUN_MODE"
echo "iBin: $IBIN"

# --- Main Logic ---

# Mode 0: Full local run
if [ "$RUN_MODE" -eq 0 ]; then
    echo "Running in local mode"
    # python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --checkBins
    # python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --createBins
    # python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --createHists --doCor
    python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --doFit
    python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --doFit --mcSig --altSig
    # python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --doFit --altSig
    # python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --doFit --altBkg
    # python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --sumUp

# Mode 1: Prepare for condor submission
elif [ "$RUN_MODE" -eq 1 ]; then
    echo "Preparing for condor"
    python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --checkBins
    python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --createBins
    python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --createHists --doCor
    # You might need an additional script here to generate and submit condor jobs

# Mode 2: Run on a condor node for a specific bin
elif [ "$RUN_MODE" -eq 2 ]; then
    if [ "$IBIN" -lt 0 ]; then
        echo "Error: For condor run (mode 2), a specific bin index (>=0) must be provided as the 4th argument."
        exit 1
    fi
    echo "Running on condor for bin $IBIN"
    # python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --doFit --iBin $IBIN
    # python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --doFit --mcSig --altSig --iBin $IBIN
    python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --doFit --altSig --iBin $IBIN
    # python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --doFit --altBkg --iBin $IBIN
# Mode 3: Collect results from condor jobs
elif [ "$RUN_MODE" -eq 3 ]; then
    echo "Collecting results from condor"
    python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --doPlot
    python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --doPlot --mcSig --altSig
    python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --doPlot --altSig
    python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --doPlot --altBkg
    python3 tnpEGM_fitter.py $CONFIG --flag $FLAG --sumUp

else
    echo "Error: Invalid run mode '$RUN_MODE'. Use 0 for local, 1 for prepare, 2 for condor run, 3 for collect."
    exit 1
fi

echo "Done."