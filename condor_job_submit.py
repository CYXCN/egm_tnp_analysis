import os
import sys
import pickle
import argparse
import importlib.util

# a help script to submit condor jobs for EGM TnP analysis Fit

# submit example:
# executable             = run.sh
# Arguments = $(config_file) $(flag) 2 $(ProcId)
# getenv                 =  True
# output                  = gen.out
# error                   = gen.err
# log                     = gen.log
# output_destination = root://eosuser.cern.ch//eos/user/y/yucao/program/trigger/CMSSW_11_2_0/src/egm_tnp_analysis/condor_output/
# MY.XRDCP_CREATE_DIR = True
# x509userproxy = $ENV(X509_USER_PROXY)
# request_cpus = 1
# +JobFlavour = "tomorrow"
# queue binnumber

# read the config file to get the baseOutDir
# the workdir will be the baseOutDir/flag
# the bin pkl file will be the workdir/bining.pkl

def load_config(config_path):
    """Loads a Python config file as a module."""
    spec = importlib.util.spec_from_file_location("config", config_path)
    config_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config_module)
    return config_module

def main():
    parser = argparse.ArgumentParser(description="Submit Condor jobs for EGM TnP analysis")
    parser.add_argument("config_file", help="Path to the configuration file (e.g., etc/config/settings_pho_run2022FG.py)")
    parser.add_argument("flag", help="ID flag for the analysis")
    parser.add_argument("--dry-run", action="store_true", help="Only generate the .sub file, do not submit.")
    args = parser.parse_args()

    # --- 1. Load configuration and define paths ---
    try:
        config = load_config(args.config_file)
        base_out_dir = config.baseOutDir
    except (FileNotFoundError, AttributeError) as e:
        print(f"Error: Could not load 'baseOutDir' from config file '{args.config_file}'.")
        print(f"Details: {e}")
        sys.exit(1)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    workdir = os.path.join(base_out_dir, args.flag)
    bin_pkl_file = os.path.join(workdir, 'bining.pkl')
    
    # --- 2. Determine the number of bins/jobs ---
    try:
        with open(bin_pkl_file, 'rb') as f:
            bining = pickle.load(f)
            bins = bining.get('bins', [])
            num_bins = len(bins)
    except FileNotFoundError:
        print(f"Error: Binning file not found at '{bin_pkl_file}'.")
        print("Please run the preparation step first (run.sh <config> <flag> 1)")
        sys.exit(1)
    
    if num_bins == 0:
        print("Error: Found 0 bins. Nothing to submit.")
        sys.exit(1)

    print(f"Found {num_bins} bins. Preparing Condor submission.")

    # --- 3. Create directories and sub file content ---
    condor_log_dir = os.path.join(workdir, 'condor_logs')
    os.makedirs(condor_log_dir, exist_ok=True)

    # Note: Condor needs absolute paths for some configurations
    executable_path = os.path.join(script_dir, 'run.sh')
    config_file_path = args.config_file
    
    sub_content = f"""
executable              = {executable_path}
Arguments               = {config_file_path} {args.flag} 2 $(ProcId)
getenv                  = True

output                  = {condor_log_dir}/job_$(ProcId).out
error                   = {condor_log_dir}/job_$(ProcId).err
log                     = {condor_log_dir}/job.log

# Transfer output files to a specific directory
# The output files from the python script will be in the workdir/{args.flag}
# We will collect them later with run_mode 3
# transfer_output_files   = ""

# Use this for EOS destination
output_destination      = root://eosuser.cern.ch/{workdir}/
MY.XRDCP_CREATE_DIR     = True
MY.SingularityImage     = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/batch-team/containers/plusbatch/el7-full:latest"

request_cpus            = 1
+JobFlavour             = "tomorrow"

queue {num_bins}
"""

    sub_file_path = os.path.join(workdir, 'condor_submit.sub')
    with open(sub_file_path, 'w') as f:
        f.write(sub_content)

    print(f"Condor submission file created at: {sub_file_path}")

    # --- 4. Submit to Condor ---
    if not args.dry_run:
        submit_command = f"condor_submit  --spool {sub_file_path}"
        print(f"Executing: {submit_command}")
        os.system(submit_command)
    else:
        print("Dry run requested. Job not submitted.")

if __name__ == "__main__":
    main()