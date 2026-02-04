#ifndef RNA_ENERGY_EVALUATOR_HPP
#define RNA_ENERGY_EVALUATOR_HPP

#include <string>
#include <unordered_map>
#include <vector>

#include "utilities.hpp"
#include "global_variables.hpp"

extern "C" {
  #include <ViennaRNA/fold_compound.h>
  #include <ViennaRNA/model.h>
  #include <ViennaRNA/params/basic.h>
  #include <ViennaRNA/loop_energies.h>
  #include <ViennaRNA/utils/basic.h>
}


using Matrix = std::vector<std::vector<int>>;

class RNAEnergyEvaluator {
public:
    ~RNAEnergyEvaluator();

    RNAEnergyEvaluator(const std::unordered_map<int, std::string> &strands_1b_input);

    void fill_single_strand(int s);

    // Energies
    int hairpin_energy(int i, int j, int s) const;
    int interior_loop_energy(int i, int j, int k, int l, int s) const;
    int interior_loop_energy_cross(int s, int i, int r, int j, int k, int l) const;

    int exterior_loop_energy(int i, int j, int s) const;
    int multiloop_stem_energy(int i, int j, int s, bool is_ext = false) const;
    int stem_energy(int i, int j, int s) const;

    // Getters
    const Matrix& get_F(int s)  const { return F_s.at(s); }
    const Matrix& get_C(int s)  const { return C_s.at(s); }
    const Matrix& get_M(int s)  const { return M_s.at(s); }
    const Matrix& get_M1(int s) const { return M1_s.at(s); }
    const std::unordered_map<int, std::string>& get_strands() const { return strands_1b; }


    const std::unordered_map<int, std::string>& get_strands_1b() const { return strands_1b; }
    const std::unordered_map<int, std::string>& get_vrna_strands() const { return vrna_strands; }

    const std::unordered_map<int, int>& get_len_map() const { return len_map; }
    const std::unordered_map<int, short*>& get_S1_map() const { return S1_map; }

    vrna_param_t* get_params() const { return params; }

private:
    // Input strands
    std::unordered_map<int, std::string> strands_1b;   // with '$'
    std::unordered_map<int, std::string> vrna_strands; // without '$'
    std::unordered_map<int, int>         len_map;      // clean length (no '$')

    // Vienna
    std::unordered_map<int, vrna_fold_compound_t*> fold_compounds;
    std::unordered_map<int, short*> S1_map;
    vrna_param_t *params = nullptr;

    // Single-strand DP matrices (Vienna-like)
    std::unordered_map<int, Matrix> C_s, M_s, M1_s, F_s;
};

    

#endif
