#ifndef RNA_ENERGY_EVALUATOR_HPP
#define RNA_ENERGY_EVALUATOR_HPP

#include <iostream>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>
#include "utilities.hpp"
#include "global_variables.hpp"

extern "C" {
    #include <ViennaRNA/fold.h>
    #include <ViennaRNA/utils/basic.h>
    #include <ViennaRNA/utils/strings.h>
    #include <ViennaRNA/model.h>
    #include <ViennaRNA/loop_energies.h>
    #include <ViennaRNA/pair_mat.h>
    #include <ViennaRNA/fold_vars.h>
}

using Matrix = std::vector<std::vector<int>>;

class RNAEnergyEvaluator {
public:
    // --- Constructor ---
    RNAEnergyEvaluator(const std::unordered_map<int, std::string> &strands);

    // --- Main DP computation ---
    void fill_single_strand(int s);
    std::unordered_map<int, std::string> strands_1b;   // with leading '$'
    std::unordered_map<int, std::string> vrna_strands; // without '$' (for ViennaRNA)


    // --- Low-level energy evaluation functions (keep all) ---
    int hairpin_energy(int i, int j, int s) const;
    int interior_loop_energy(int i, int j, int k, int l, int s) const;
    int exterior_loop_energy(int i, int j, int s) const;
    int multiloop_stem_energy(int i, int j, int s, bool is_ext = false) const;
    int stem_energy(int i, int j, int s) const;

    // --- Matrix getters ---
    const Matrix& get_F(int s) const { return F_s.at(s); }
    const Matrix& get_C(int s) const { return C_s.at(s); }
    const Matrix& get_M(int s) const { return M_s.at(s); }
    const Matrix& get_M1(int s) const { return M1_s.at(s); }
    const std::unordered_map<int, std::string>& get_strands() const {
        return strands;
    }
    const std::unordered_map<int, short*>& get_S1_map() const {
        return S1_map;
    }
        
    vrna_param_t* get_params() const { return params; }


private:
    // --- Data structures ---
    std::unordered_map<int, std::string> strands;
    std::unordered_map<int, Matrix> C_s, M_s, M1_s, F_s;
    std::unordered_map<int, vrna_fold_compound_t*> fold_compounds;
    std::unordered_map<int, short*> S1_map;
    vrna_param_t *params = nullptr;
    const int inf_energy = 1000000;
};

#endif
