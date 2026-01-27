import json
import os

def get_json_path(year, obj_type="Ele"):
    import os
    REPO_DICT = {
        "2025": "Run3-25Prompt-Summer24-NanoAODv15",
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