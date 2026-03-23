#ifndef HISTMANAGER_H
#define HISTMANAGER_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <set>
#include <map>

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"
#include "TList.h"

// Include the corrector header
#include "scaleSmearing.h"

// Struct for manual range cuts (optimization)
struct RangeCut {
    int var_idx;
    float min_val;
    float max_val;
};

// Indices for manual cuts
const int IDX_EL_ET = 0;
const int IDX_EL_PT = 1;
const int IDX_TAG_PT = 2;

// Struct for branch name mapping
struct BranchMapping {
    std::string pair_mass;
    std::string run;
    // Probe branches
    std::string probe_pt;
    std::string probe_eta;
    std::string probe_phi;
    std::string probe_sc_eta;
    std::string probe_r9;
    std::string probe_seedGain;
    std::string probe_et;
    // Tag branches
    std::string tag_pt;
    std::string tag_eta;
    std::string tag_phi;
    std::string tag_sc_eta;
    std::string tag_r9;
    std::string tag_seedGain;

    std::string rho;
    std::string probe_pfPhoIso03;
    std::string probe_sieie;
    std::string probe_iso;
    std::string probe_rel_iso;
    std::string probe_electronVeto;
    std::string probe_mvaID;
    std::string probe_hoe;
    std::string tag_pfPhoIso03;
    std::string tag_sieie;
    std::string tag_iso;
    std::string tag_rel_iso;
    std::string tag_electronVeto;
    std::string tag_mvaID;
    std::string tag_hoe;
    
    // Default constructor with default branch names
    BranchMapping() :
        pair_mass("pair_mass"),
        run("run"),
        probe_pt("prob_pt"),
        probe_eta("prob_eta"),
        probe_phi("prob_phi"),
        probe_sc_eta("prob_sc_eta"),
        probe_r9("prob_r9"),
        probe_seedGain("prob_seedGain"),
        probe_et("prob_et"),
        tag_pt("tag_pt"),
        tag_eta("tag_eta"),
        tag_phi("tag_phi"),
        tag_sc_eta("tag_sc_eta"),
        tag_r9("tag_r9"),
        tag_seedGain("tag_seedGain"),

        rho("rho"),
        probe_pfPhoIso03("probe_pfPhoIso03"),
        probe_sieie("probe_sieie"),
        probe_iso("probe_iso"),
        probe_rel_iso("probe_rel_iso"),
        probe_electronVeto("probe_electronVeto"),
        probe_mvaID("probe_mvaID"),
        probe_hoe("probe_hoe"),
        tag_pfPhoIso03("tag_pfPhoIso03"),
        tag_sieie("tag_sieie"),
        tag_iso("tag_iso"),
        tag_rel_iso("tag_rel_iso"),
        tag_electronVeto("tag_electronVeto"),
        tag_mvaID("tag_mvaID"),
        tag_hoe("tag_hoe")
    {}
};

struct HGGSelection {

    float EA1_EB1 = 0.102056;
    float EA2_EB1 = -0.000398112;
    float EA1_EB2 = 0.0820317;
    float EA2_EB2 = -0.000286224;
    float EA1_EE1 = 0.0564915;
    float EA2_EE1 = -0.000248591;
    float EA1_EE2 = 0.0428606;
    float EA2_EE2 = -0.000171541;
    float EA1_EE3 = 0.0395282;
    float EA2_EE3 = -0.000121398;
    float EA1_EE4 = 0.0369761;
    float EA2_EE4 = -8.10369e-05;
    float EA1_EE5 = 0.0369417;
    float EA2_EE5 = -2.76885e-05;

    float max_pho_iso_EB_low_r9 = 4.0;
    float max_pho_iso_EE_low_r9 = 4.0;

    float min_full5x5_r9_EB_high_r9 = 0.85;
    float min_full5x5_r9_EE_high_r9 = 0.9;
    float min_full5x5_r9_EB_low_r9 = 0.5;
    float min_full5x5_r9_EE_low_r9 = 0.8;
    
    float max_trkSumPtHollowConeDR03_EB_low_r9 = 6.0;
    float max_trkSumPtHollowConeDR03_EE_low_r9 = 6.0;
    
    float max_sieie_EB_low_r9 = 0.015;
    float max_sieie_EE_low_r9 = 0.035;
    
    float min_pt_photon = 25.0;
    float min_mvaid = -0.9;
    float max_hovere = 0.08;
    float min_full5x5_r9 = 0.8;
    
    float max_chad_iso = 20.0;
    float max_chad_rel_iso = 0.3;
};

class HistManager {
public:
    HistManager(TChain* chain, TFile* outfile, bool isMC, std::string flag_selection, std::string base_selection = "", std::string sample_name = "", int max_events = -1) 
        : fChain(chain), fOutfile(outfile), fIsMC(isMC), fFlagSelection(flag_selection), fBaseSelection(base_selection), 
          fSampleName(sample_name), fMaxEvents(max_events), fCorrector(nullptr), fDoCorrection(false), fDoReweight(false), fApplyHGGPreselection(false), fBaseFormula(nullptr), fFlagFormula(nullptr) {
        
        // Initialize variable values vector
        fCurrentVarValues.resize(3);
    }

    ~HistManager() {
        if (fChain) {
            fChain->SetNotify(nullptr);
        }
        if (fCorrector) delete fCorrector;
        if (fBaseFormula) delete fBaseFormula;
        if (fFlagFormula) delete fFlagFormula;
        for (auto& bin : fBins) {
            if (bin.formula) delete bin.formula;
        }
    }

    // Configure branch name mapping
    void setBranchMapping(const BranchMapping& mapping) {
        fBranchMapping = mapping;
    }

    // Add a bin configuration
    void addBin(std::string name, std::string title, int nbins, double min, double max, std::string cut_string, std::vector<RangeCut> range_cuts) {
        BinDef bin;
        bin.name = name;
        bin.title = title;
        
        fOutfile->cd();
        bin.hPass = new TH1D((name + "_Pass").c_str(), title.c_str(), nbins, min, max);
        bin.hFail = new TH1D((name + "_Fail").c_str(), title.c_str(), nbins, min, max);
        bin.hPass->Sumw2();
        bin.hFail->Sumw2();

        bin.cut_string = cut_string;
        bin.range_cuts = range_cuts;
        bin.formula = nullptr;

        fBins.push_back(bin);
    }

    // Configure Run3 Corrections
    void configureCorrections(std::string json_path) {
        std::cout << "Initializing C++ Run3 Corrector with " << json_path << "..." << std::endl;
        fCorrector = new Run3Corrector(json_path);
        fDoCorrection = true;
    }

    // Set branches to be enabled (optimization)
    void setBranches(std::vector<std::string> branches) {
        fChain->SetBranchStatus("*", 0);
        for (const auto& br : branches) {
            fChain->SetBranchStatus(br.c_str(), 1);
        }
        
        // Always enable correction branches if needed
        if (fDoCorrection) {
            std::vector<std::string> corr_branches = {
                fBranchMapping.run,
                fBranchMapping.probe_pt, fBranchMapping.probe_eta, fBranchMapping.probe_phi,
                fBranchMapping.probe_sc_eta, fBranchMapping.probe_r9, fBranchMapping.probe_seedGain, fBranchMapping.probe_et,
                fBranchMapping.tag_pt, fBranchMapping.tag_eta, fBranchMapping.tag_phi,
                fBranchMapping.tag_sc_eta, fBranchMapping.tag_r9, fBranchMapping.tag_seedGain
            };
            for (const auto& br : corr_branches) {
                fChain->SetBranchStatus(br.c_str(), 1);
            }
        }

        // Enable reweight branches independently (probe_r9, probe_sc_eta are needed
        // even when correction is disabled)
        if (fDoReweight && !fDoCorrection) {
            fChain->SetBranchStatus(fBranchMapping.probe_r9.c_str(), 1);
            fChain->SetBranchStatus(fBranchMapping.probe_sc_eta.c_str(), 1);
        }
        
        // Always enable pair_mass and flag
        fChain->SetBranchStatus(fBranchMapping.pair_mass.c_str(), 1);
    }

    // Configure R9-Eta reweighting: x = R9, y = Eta
    void setReweightMap(std::vector<double> r9_bins, std::vector<double> eta_bins, std::vector<std::vector<double>> ratio_map) {
        fR9Bins = r9_bins;
        fEtaBins = eta_bins;
        fRatioMap = ratio_map;
        fDoReweight = true;
    }

    void setHGGSelection( bool apply_preselection) {
        fApplyHGGPreselection = apply_preselection;
        // cout
        std::cout << "HGG Preselection " << (fApplyHGGPreselection ? "enabled" : "disabled") << "." << std::endl;
    }

    HGGSelection HGG_Selection; // Public member to hold selection parameters
    bool passPhotonSelection(
        float eta, 
        float pfPhoIso03, 
        float sc_eta,
        float r9, 
        float sieie, 
        float iso, 
        float rel_iso, 
        float electronVeto, 
        float pt, 
        float mvaID, 
        float hoe, 
        float rho, 
        bool apply_electron_veto = false,   // corresponding to electron_veto
        bool revert_electron_veto = false
    ) const {
        
        float photon_abs_eta = std::abs(eta);
        float photon_abs_sc_eta = std::abs(sc_eta);
        bool isScEtaEB = photon_abs_sc_eta < 1.4442;
        bool isScEtaEE = photon_abs_sc_eta > 1.566 && photon_abs_sc_eta < 2.5;

        // 1. cal pass_phoIso_rho_corr_EB
        bool pass_phoIso_rho_corr_EB = false;
        if (photon_abs_eta > 0.0 && photon_abs_eta < 1.0) {
            pass_phoIso_rho_corr_EB = (pfPhoIso03 - (rho * HGG_Selection.EA1_EB1) - (rho * rho * HGG_Selection.EA2_EB1)) < HGG_Selection.max_pho_iso_EB_low_r9;
        } else if (photon_abs_eta > 1.0 && photon_abs_eta < 1.4442) {
            pass_phoIso_rho_corr_EB = (pfPhoIso03 - (rho * HGG_Selection.EA1_EB2) - (rho * rho * HGG_Selection.EA2_EB2)) < HGG_Selection.max_pho_iso_EB_low_r9;
        }

        // 2. cla pass_phoIso_rho_corr_EE
        bool pass_phoIso_rho_corr_EE = false;
        if (photon_abs_eta > 1.566 && photon_abs_eta < 2.0) {
            pass_phoIso_rho_corr_EE = (pfPhoIso03 - (rho * HGG_Selection.EA1_EE1) - (rho * rho * HGG_Selection.EA2_EE1)) < HGG_Selection.max_pho_iso_EE_low_r9;
        } else if (photon_abs_eta > 2.0 && photon_abs_eta < 2.2) {
            pass_phoIso_rho_corr_EE = (pfPhoIso03 - (rho * HGG_Selection.EA1_EE2) - (rho * rho * HGG_Selection.EA2_EE2)) < HGG_Selection.max_pho_iso_EE_low_r9;
        } else if (photon_abs_eta > 2.2 && photon_abs_eta < 2.3) {
            pass_phoIso_rho_corr_EE = (pfPhoIso03 - (rho * HGG_Selection.EA1_EE3) - (rho * rho * HGG_Selection.EA2_EE3)) < HGG_Selection.max_pho_iso_EE_low_r9;
        } else if (photon_abs_eta > 2.3 && photon_abs_eta < 2.4) {
            pass_phoIso_rho_corr_EE = (pfPhoIso03 - (rho * HGG_Selection.EA1_EE4) - (rho * rho * HGG_Selection.EA2_EE4)) < HGG_Selection.max_pho_iso_EE_low_r9;
        } else if (photon_abs_eta > 2.4 && photon_abs_eta < 2.5) {
            pass_phoIso_rho_corr_EE = (pfPhoIso03 - (rho * HGG_Selection.EA1_EE5) - (rho * rho * HGG_Selection.EA2_EE5)) < HGG_Selection.max_pho_iso_EE_low_r9;
        }

        // 3. High R9
        bool isEB_high_r9 = isScEtaEB && (r9 > HGG_Selection.min_full5x5_r9_EB_high_r9);
        bool isEE_high_r9 = isScEtaEE && (r9 > HGG_Selection.min_full5x5_r9_EE_high_r9);

        // 4. Low R9
        bool isEB_low_r9 = isScEtaEB 
                        && (r9 > HGG_Selection.min_full5x5_r9_EB_low_r9) 
                        && (r9 < HGG_Selection.min_full5x5_r9_EB_high_r9) 
                        && (iso < HGG_Selection.max_trkSumPtHollowConeDR03_EB_low_r9) 
                        && (sieie < HGG_Selection.max_sieie_EB_low_r9) 
                        && pass_phoIso_rho_corr_EB;

        bool isEE_low_r9 = isScEtaEE 
                        && (r9 > HGG_Selection.min_full5x5_r9_EE_low_r9) 
                        && (r9 < HGG_Selection.min_full5x5_r9_EE_high_r9) 
                        && (iso < HGG_Selection.max_trkSumPtHollowConeDR03_EE_low_r9) 
                        && (sieie < HGG_Selection.max_sieie_EE_low_r9) 
                        && pass_phoIso_rho_corr_EE;

        // 5. electron Veto (Electron Veto)
        bool e_veto_cut = true;
        if (apply_electron_veto) {
            e_veto_cut = (electronVeto > 0.5);
        } else if (revert_electron_veto) {
            e_veto_cut = (electronVeto < 0.5);
        }

        if (!e_veto_cut) return false;
        if (pt <= HGG_Selection.min_pt_photon) return false;
        if (!(isScEtaEB || isScEtaEE)) return false;
        if (mvaID <= HGG_Selection.min_mvaid) return false;
        if (hoe >= HGG_Selection.max_hovere) return false;

        bool pass_r9_or_iso = (r9 > HGG_Selection.min_full5x5_r9) 
                           || ((rel_iso) < HGG_Selection.max_chad_iso) // note, rel_iso‘s definition is different from NANOAOD, it is PFIsoChgQuadratic but not PFIsoChg/pt.
                           || (rel_iso / pt < HGG_Selection.max_chad_rel_iso);
        if (!pass_r9_or_iso) return false;

        bool pass_category = isEB_high_r9 || isEB_low_r9 || isEE_high_r9 || isEE_low_r9;
        if (!pass_category) return false;

        // if all pass
        if (r9 < 0.4){
            std::cout << "Photon passed Selection cut with: R9 = " << r9 << ", Eta = " << eta << ", SC Eta = " << sc_eta << ", pt = " << pt << std::endl;
        }
        return true; 
    }

    // Main event loop
    void process() {
        setupAddresses();
        setupFormulas();

        Long64_t nevts = fMaxEvents < 0 ? fChain->GetEntries() : std::min(fMaxEvents, (int)fChain->GetEntries());
        Long64_t frac_of_nevts = nevts / 20;
        if (frac_of_nevts == 0) frac_of_nevts = 1;

        std::cout << "==================================================" << std::endl;
        if (!fSampleName.empty()) {
            std::cout << "Processing sample: " << fSampleName << std::endl;
        }
        std::cout << "Total events: " << nevts << std::endl;
        std::cout << (fDoCorrection ? "Starting event loop with Run 3 corrections..." : "Starting event loop...") << std::endl;
        std::cout << "==================================================" << std::endl;

        TLorentzVector vProbe, vTag, vSum;

        Long64_t pass_base_count = 0;
        Long64_t enter_analysis_count = 0;

        bool check_correction = false;
        int max_check = 100;

        for (Long64_t index = 0; index < nevts; ++index) {
            if (index % frac_of_nevts == 0) {
                int progress = (int)(100.0 * index / nevts);
                if (!fSampleName.empty()) {
                    std::cout << "[" << fSampleName << "] " << progress << "% processed (" << index << "/" << nevts << ")" << std::endl;
                } else {
                    std::cout << progress << "% processed (" << index << "/" << nevts << ")" << std::endl;
                }
            }

            fChain->GetEntry(index);

            if (fBranchMapping.probe_pt == fBranchMapping.probe_et) {
                probe_et = probe_pt;
            }


            // Base Selection
            if (fBaseFormula && fBaseFormula->EvalInstance(0) == 0) {
                continue;
            }
            pass_base_count++;

            float current_probe_pt = probe_pt;
            float current_probe_et = probe_et;
            float current_tag_pt = tag_pt;
            double to_fill_mass = pair_mass;

            if (fDoCorrection && fCorrector) {
                unsigned int current_run = run;

                // Probe Correction
                double corr_factor_probe = fCorrector->get_correction(
                    fIsMC, current_run, probe_pt, probe_r9, probe_sc_eta, (int)probe_seedGain
                );

                // Tag Correction
                double corr_factor_tag = fCorrector->get_correction(
                    fIsMC, current_run, tag_pt, tag_r9, tag_sc_eta, (int)tag_seedGain
                );

                current_probe_pt = probe_pt * corr_factor_probe;
                current_probe_et = probe_et * corr_factor_probe;
                current_tag_pt = tag_pt * corr_factor_tag;

                // Reconstruct Mass
                vProbe.SetPtEtaPhiM(current_probe_pt, probe_eta, probe_phi, 0.000511);
                vTag.SetPtEtaPhiM(current_tag_pt, tag_eta, tag_phi, 0.000511);
                vSum = vTag + vProbe;
                to_fill_mass = vSum.M();
                if (check_correction && pass_base_count < max_check) {
                    std::cout << "Event " << index << ": Run " << current_run 
                              << ", Original Probe PT " << probe_pt << ", Corrected PT " << current_probe_pt
                              << "; Original Tag PT " << tag_pt << ", Corrected PT " << current_tag_pt
                              << "; Original Mass " << pair_mass << ", Corrected Mass " << to_fill_mass << std::endl;
                }
            }

            // do HGG preselection
            if (fApplyHGGPreselection) {
                bool probe_pass_preselection = passPhotonSelection(
                    probe_eta, 
                    probe_pfPhoIso03, 
                    probe_sc_eta,
                    probe_r9, 
                    probe_sieie, 
                    probe_iso, 
                    probe_rel_iso, 
                    probe_electronVeto, 
                    current_probe_pt, 
                    probe_mvaID, 
                    probe_hoe, 
                    rho,
                    false, // apply electron veto for preselection
                    true  // revert electron veto for probe (probe should fail electron veto) meaning probe_electronVeto should be a electron
                );
                if (nevts > 100 and !probe_pass_preselection)
                    continue;
                bool tag_pass_preselection = passPhotonSelection(
                    tag_eta, 
                    tag_pfPhoIso03, 
                    tag_sc_eta,
                    tag_r9, 
                    tag_sieie, 
                    tag_iso, 
                    tag_rel_iso, 
                    tag_electronVeto, 
                    current_tag_pt, 
                    tag_mvaID, 
                    tag_hoe, 
                    rho,
                    false // apply electron veto for preselection
                );
                if (nevts < 100) 
                    std::cout << "Event " << index << ": HGG preselection. Probe pass: " << probe_pass_preselection << ", Tag pass: " << tag_pass_preselection << std::endl;
                if (!probe_pass_preselection || !tag_pass_preselection){
                    if (nevts < 100) 
                    std::cout << "Event " << index << ": Failed HGG preselection, skipping event." << std::endl;
                    continue;
                }
            }


            double reweight_factor = 1.0;
            if (fDoReweight) {
                // Find R9 index (x-axis)
                int r9_idx = -1;
                for (size_t i = 0; i < fR9Bins.size() - 1; ++i) {
                    if (probe_r9 >= fR9Bins[i] && probe_r9 < fR9Bins[i+1]) {
                        r9_idx = (int)i;
                        break;
                    }
                }

                // Find Eta index (y-axis)
                double abs_eta = std::abs(probe_sc_eta);
                int eta_idx = -1;
                for (size_t i = 0; i < fEtaBins.size() - 1; ++i) {
                    if (abs_eta >= fEtaBins[i] && abs_eta < fEtaBins[i+1]) {
                        eta_idx = (int)i;
                        break;
                    }
                }

                if (r9_idx != -1 && eta_idx != -1) {
                    reweight_factor = fRatioMap[r9_idx][eta_idx];
                    if (reweight_factor < 0) reweight_factor = 1.0; // Handle invalid bins
                }
                if (nevts < 100) {
                    std::cout << "Event " << index << ": Probe R9 " << probe_r9 << ", SC Eta " << probe_sc_eta 
                              << "; R9 bin " << r9_idx << ", Eta bin " << eta_idx 
                              << "; Reweight factor " << reweight_factor << std::endl;
                }
            }

            // Update values for manual cuts
            fCurrentVarValues[IDX_EL_ET] = current_probe_et;
            fCurrentVarValues[IDX_EL_PT] = current_probe_pt;
            fCurrentVarValues[IDX_TAG_PT] = current_tag_pt;

            // Loop over bins
            for (auto& bin : fBins) {
                // Manual Range Check
                bool pass_manual = true;
                for (const auto& cut : bin.range_cuts) {
                    if (cut.var_idx >= 0 && cut.var_idx < 3) {
                        float val = fCurrentVarValues[cut.var_idx];
                        if (!(val >= cut.min_val && val < cut.max_val)) {
                            pass_manual = false;
                            break;
                        }
                    }
                }
                if (!pass_manual) continue;

                // Formula Check (Weight)
                double weight = bin.formula->EvalInstance(0) * reweight_factor;
                if (weight != 0) {
                    if (fFlagFormula->EvalInstance(0)) {
                        bin.hPass->Fill(to_fill_mass, weight);
                    } else {
                        bin.hFail->Fill(to_fill_mass, weight);
                    }
                    enter_analysis_count++;
                    break; 
                }
            }
        }
        std::cout << "[" << fSampleName << "] Total events passing base selection: " << pass_base_count << std::endl;
        std::cout << "[" << fSampleName << "] Total events entering analysis: " << enter_analysis_count << std::endl;
    }

    // Finalize: remove negative bins and write
    void finalize() {
        fOutfile->cd();
        for (auto& bin : fBins) {
            removeNegativeBins(bin.hPass);
            removeNegativeBins(bin.hFail);
            bin.hPass->Write();
            bin.hFail->Write();
        }
    }

private:
    TChain* fChain;
    TFile* fOutfile;
    bool fIsMC;
    std::string fFlagSelection;
    std::string fBaseSelection;
    std::string fSampleName;
    int fMaxEvents;
    
    Run3Corrector* fCorrector;
    bool fDoCorrection;
    bool fDoReweight;
    bool fApplyHGGPreselection;
    std::vector<double> fEtaBins;
    std::vector<double> fR9Bins;
    std::vector<std::vector<double>> fRatioMap;
    
    BranchMapping fBranchMapping;

    struct BinDef {
        std::string name;
        std::string title;
        TH1D* hPass;
        TH1D* hFail;
        std::string cut_string;
        TTreeFormula* formula;
        std::vector<RangeCut> range_cuts;
    };
    std::vector<BinDef> fBins;

    // Tree Variables (using generic names internally)
    float pair_mass;
    unsigned int run;
    float probe_pt, probe_eta, probe_phi, probe_sc_eta, probe_r9, probe_seedGain, probe_et;
    float tag_pt, tag_eta, tag_phi, tag_sc_eta, tag_r9, tag_seedGain;

    float rho;
    float probe_pfPhoIso03, probe_sieie, probe_iso, probe_rel_iso, probe_mvaID, probe_hoe;
    float probe_electronVeto;
    float tag_pfPhoIso03, tag_sieie, tag_iso, tag_rel_iso, tag_mvaID, tag_hoe;
    float tag_electronVeto;

    // Formulas
    TTreeFormula* fBaseFormula;
    TTreeFormula* fFlagFormula;
    TList fFormulasList;

    std::vector<double> fCurrentVarValues;

    void setupAddresses() {
        fChain->SetBranchAddress(fBranchMapping.pair_mass.c_str(), &pair_mass);
        
        if (fDoCorrection) {
            fChain->SetBranchAddress(fBranchMapping.run.c_str(), &run);
            fChain->SetBranchAddress(fBranchMapping.probe_pt.c_str(), &probe_pt);
            fChain->SetBranchAddress(fBranchMapping.probe_eta.c_str(), &probe_eta);
            fChain->SetBranchAddress(fBranchMapping.probe_phi.c_str(), &probe_phi);
            fChain->SetBranchAddress(fBranchMapping.probe_sc_eta.c_str(), &probe_sc_eta);
            fChain->SetBranchAddress(fBranchMapping.probe_r9.c_str(), &probe_r9);
            fChain->SetBranchAddress(fBranchMapping.probe_seedGain.c_str(), &probe_seedGain);
            if (fBranchMapping.probe_et != fBranchMapping.probe_pt) {
                fChain->SetBranchAddress(fBranchMapping.probe_et.c_str(), &probe_et);
            }

            fChain->SetBranchAddress(fBranchMapping.tag_pt.c_str(), &tag_pt);
            fChain->SetBranchAddress(fBranchMapping.tag_eta.c_str(), &tag_eta);
            fChain->SetBranchAddress(fBranchMapping.tag_phi.c_str(), &tag_phi);
            fChain->SetBranchAddress(fBranchMapping.tag_sc_eta.c_str(), &tag_sc_eta);
            fChain->SetBranchAddress(fBranchMapping.tag_r9.c_str(), &tag_r9);
            fChain->SetBranchAddress(fBranchMapping.tag_seedGain.c_str(), &tag_seedGain);
        } else {
            // If not doing correction, we still need these for manual cuts if they are used
            probe_pt = 0; probe_et = 0; tag_pt = 0;
            probe_r9 = 0; probe_sc_eta = 0;
            if(fChain->GetBranchStatus(fBranchMapping.probe_pt.c_str())) 
                fChain->SetBranchAddress(fBranchMapping.probe_pt.c_str(), &probe_pt);
            if(fChain->GetBranchStatus(fBranchMapping.probe_et.c_str())) 
                fChain->SetBranchAddress(fBranchMapping.probe_et.c_str(), &probe_et);
            if(fChain->GetBranchStatus(fBranchMapping.tag_pt.c_str())) 
                fChain->SetBranchAddress(fBranchMapping.tag_pt.c_str(), &tag_pt);

            // Set up reweight branches when reweighting is enabled without correction
            if (fDoReweight) {
                if(fChain->GetBranchStatus(fBranchMapping.probe_r9.c_str()))
                    fChain->SetBranchAddress(fBranchMapping.probe_r9.c_str(), &probe_r9);
                if(fChain->GetBranchStatus(fBranchMapping.probe_sc_eta.c_str()))
                    fChain->SetBranchAddress(fBranchMapping.probe_sc_eta.c_str(), &probe_sc_eta);
            }
        }
        if (fApplyHGGPreselection) {
            fChain->SetBranchAddress(fBranchMapping.rho.c_str(), &rho);
            fChain->SetBranchAddress(fBranchMapping.probe_pfPhoIso03.c_str(), &probe_pfPhoIso03);
            fChain->SetBranchAddress(fBranchMapping.probe_sieie.c_str(), &probe_sieie);
            fChain->SetBranchAddress(fBranchMapping.probe_iso.c_str(), &probe_iso);
            fChain->SetBranchAddress(fBranchMapping.probe_rel_iso.c_str(), &probe_rel_iso);
            fChain->SetBranchAddress(fBranchMapping.probe_electronVeto.c_str(), &probe_electronVeto);
            fChain->SetBranchAddress(fBranchMapping.probe_mvaID.c_str(), &probe_mvaID);
            fChain->SetBranchAddress(fBranchMapping.probe_hoe.c_str(), &probe_hoe);
            fChain->SetBranchAddress(fBranchMapping.tag_pfPhoIso03.c_str(), &tag_pfPhoIso03);
            fChain->SetBranchAddress(fBranchMapping.tag_sieie.c_str(), &tag_sieie);
            fChain->SetBranchAddress(fBranchMapping.tag_iso.c_str(), &tag_iso);
            fChain->SetBranchAddress(fBranchMapping.tag_rel_iso.c_str(), &tag_rel_iso);
            fChain->SetBranchAddress(fBranchMapping.tag_electronVeto.c_str(), &tag_electronVeto);
            fChain->SetBranchAddress(fBranchMapping.tag_mvaID.c_str(), &tag_mvaID);
            fChain->SetBranchAddress(fBranchMapping.tag_hoe.c_str(), &tag_hoe);
        }

    }

    void setupFormulas() {
        if (!fBaseSelection.empty()) {
            fBaseFormula = new TTreeFormula("Base_Selection", fBaseSelection.c_str(), fChain);
            fFormulasList.Add(fBaseFormula);
        }

        fFlagFormula = new TTreeFormula("Flag_Selection", fFlagSelection.c_str(), fChain);
        fFormulasList.Add(fFlagFormula);

        for (auto& bin : fBins) {
            std::string cut = bin.cut_string.empty() ? "1" : bin.cut_string;
            bin.formula = new TTreeFormula((bin.name + "_Selection").c_str(), cut.c_str(), fChain);
            fFormulasList.Add(bin.formula);
        }
        
        fChain->SetNotify(&fFormulasList);
    }

    void removeNegativeBins(TH1D* h) {
        for (int i = 0; i <= h->GetNbinsX() + 1; ++i) {
            if (h->GetBinContent(i) < 0) h->SetBinContent(i, 0);
        }
    }
};

#endif
