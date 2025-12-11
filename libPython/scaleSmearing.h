#ifndef SCALE_SMEARING_H
#define SCALE_SMEARING_H

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <random>
#include <cmath>
#include <variant>
#include <unordered_map>

#include <correction.h>

class Run3Corrector {
public:
    Run3Corrector(std::string json_path) {
        rng = std::mt19937(125);
        dist = std::normal_distribution<double>(0.0, 1.0);

        try {
            auto cset = correction::CorrectionSet::from_file(json_path);
            
            // Scale is CompoundCorrection
            scale_evaluator = cset->compound().at("Scale");
            // SmearAndSyst is Correction
            smear_evaluator = cset->at("SmearAndSyst");
            
            std::cout << "[C++] Loaded corrections from: " << json_path << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "[C++] Error loading corrections: " << e.what() << std::endl;
            throw;
        }
    }

    double get_correction(bool is_mc, int run, float pt, float r9, float sc_eta, int seed_gain) {
        if (!is_mc) {
            // Data: Scale correction (Compound)
            std::vector<correction::Variable::Type> inputs = {
                std::string("scale"),
                static_cast<double>(run),
                static_cast<double>(sc_eta),
                static_cast<double>(r9),
                static_cast<double>(pt),
                static_cast<double>(seed_gain)
            };
            return scale_evaluator->evaluate(inputs);
        } else {
            // MC: Smearing
            std::vector<correction::Variable::Type> smear_inputs = {
                std::string("smear"),
                static_cast<double>(pt),
                static_cast<double>(r9),
                static_cast<double>(sc_eta)
            };
            double smear_width = smear_evaluator->evaluate(smear_inputs);
            double random_number = dist(rng);
            return 1.0 + smear_width * random_number;
        }
    }

private:
    std::shared_ptr<const correction::CompoundCorrection> scale_evaluator;
    std::shared_ptr<const correction::Correction> smear_evaluator;
    std::mt19937 rng;
    std::normal_distribution<double> dist;
};

struct CutRange {
    std::string var_name;
    double min = -9999;
    double max = 9999;
};

bool passCorCut(
    const std::unordered_map<std::string, double>& vars,
    const std::vector<CutRange>& cuts
) {
    for (const auto& c : cuts) {
        auto it = vars.find(c.var_name);
        if (it == vars.end()) return false;
        double v = it->second;
        if (v <= c.min || v >= c.max) return false;
    }
    return true;
}

#endif