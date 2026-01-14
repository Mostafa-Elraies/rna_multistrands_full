#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <cassert>
#include "nussinov.hpp"
#include "utilities.hpp"
#include "RNAEnergyEvaluator.hpp"
extern "C" {
    #include <ViennaRNA/fold.h>
    #include <ViennaRNA/params/basic.h>
    #include <ViennaRNA/fold_compound.h>
    #include <ViennaRNA/loop_energies.h>
    #include <ViennaRNA/eval.h>
    #include <ViennaRNA/utils.h>
}

RNAEnergyEvaluator::RNAEnergyEvaluator(const std::unordered_map<int, std::string> &input_strands)
{
    strands_1b = input_strands; // keep the 1-based (with '$') copy

        vrna_md_t md;
vrna_md_set_default(&md); // Initialize model defaults

params = vrna_params(&md);
for (const auto &[s, seq_with_dollar] : strands_1b) {
    // Make a clean sequence for ViennaRNA
    const std::string clean = (!seq_with_dollar.empty() && seq_with_dollar[0] == '$')
                              ? seq_with_dollar.substr(1)
                              : seq_with_dollar;

    vrna_strands[s] = clean;  // store clean copy

    // Build fold_compound on the clean sequence
    vrna_fold_compound_t *fc = vrna_fold_compound(clean.c_str(), nullptr, VRNA_OPTION_MFE);
    fold_compounds[s] = fc;

    vrna_mfe(fc, nullptr);

    int len = static_cast<int>(clean.length()); // no '$' here
        
        Matrix C(len + 1, std::vector<int>(len + 1, inf_energy));
        Matrix M(len + 1, std::vector<int>(len + 1, inf_energy));
        Matrix M1(len + 1, std::vector<int>(len + 1, inf_energy));
        Matrix F(len + 1, std::vector<int>(len + 1, inf_energy));
        
        for (int i = 1; i <= len; ++i) {
            for (int j = i; j <= len; ++j) {
                int idx = fc->iindx[i] - j; // ViennaRNA 1D index (1-based)
                C[i][j] = fc->matrices->c[idx];
                M[i][j] = fc->matrices->fML[idx];
            }
        }
          // here you store and call the DP
          C_s[s] = std::move(C);
          M_s[s] = std::move(M);
          M1_s[s] = std::move(M1);
          F_s[s] = std::move(F);
          S1_map[s] = fc->sequence_encoding; // 1-based in Vienna

          fill_single_strand(s);   // <--- this line calls the function we’ll define next
        }                            // closes the for over strands
  }                                   // <-- ADD THIS: closes the constructor
  
  void RNAEnergyEvaluator::fill_single_strand(int s) {
      const std::string& seq = strands[s];
      auto& F  = F_s[s];
      auto& C  = C_s[s];
      auto& M  = M_s[s];
      auto& M1 = M1_s[s];
      const auto& S1 = S1_map[s];
  
      int n = seq.size() - 1; // 1-based indexing
      const int INF_ENERGY = inf_energy;
  
      // --- initialize ---
      for (int i = 1; i <= n; ++i)
          for (int j = 1; j <= n; ++j) {
              C[i][j]  = INF_ENERGY;
              M[i][j]  = INF_ENERGY;
              M1[i][j] = INF_ENERGY;
              F[i][j]  = (i >= j ? 0 : INF_ENERGY);
          }
  
      // --- main DP loops ---
      for (int span = theta + 1; span <= n; ++span) {
          for (int i = 1; i + span <= n; ++i) {
              int j = i + span;
  
              // ---- C[i][j] ----
              if (can_pair(seq[i - 1], seq[j - 1])) {
                  int best = INF_ENERGY;
  
                  best = std::min(best, hairpin_energy(i, j, s));
  
                  for (int p = i + 1; p < j - theta - 1; ++p) {
                      for (int q = p + theta + 1; q < j; ++q) {
                          if (!can_pair(seq[p - 1], seq[q - 1])) continue;
                          int inner = C[p][q];
                          if (inner >= INF_ENERGY) continue;
                          int e_int = interior_loop_energy(i, j, p, q, s);
                          best = std::min(best, inner + e_int);
                      }
                  }
  
                  for (int u = i + 1; u < j; ++u) {
                      int left  = (u - 1 >= i + 1) ? M1[i + 1][u - 1] : INF_ENERGY;
                      int right = M[u][j - 1];
                      int cand  = (left >= INF_ENERGY || right >= INF_ENERGY) ? INF_ENERGY : left + right + params->MLclosing;
                      if (cand < best) best = cand;
                  }
  
                  C[i][j] = best;
              }
  
              // ---- M[i][j] ----
              {
                  int best = INF_ENERGY;
  
                  // j unpaired
                  best = std::min(best, M[i][j - 1] + params->MLbase);
  
                  // split (M | C)
                  for (int u = i; u < j; ++u) {
                    if (M[i][u] < INF_ENERGY && C[u + 1][j] < INF_ENERGY) {
                        int type = vrna_get_ptype_md(S1_map[s][u + 1], S1_map[s][j], &params->model_details);
                        int penalty = params->MLintern[type];
                        best = std::min(best, M[i][u] + C[u + 1][j] + penalty);
                    }
                }
  
                  // start from C
                  // case 3: start directly from a C
            for (int u = i; u < j; ++u) {
                if (C[u + 1][j] < INF_ENERGY) {
                int type = vrna_get_ptype_md(S1_map[s][u + 1], S1_map[s][j], &params->model_details);
                int penalty = params->MLintern[type];
                best = std::min(best, C[u + 1][j] + penalty);
    }
}

  
                  M[i][j] = best;
              }
  
              // ---- M1[i][j] ----
              {
                  int best = INF_ENERGY;
                  best = std::min(best, M1[i][j - 1] + params->MLbase);
                  // case 2: helix at (i,j)
if (C[i][j] < INF_ENERGY) {
    int type = vrna_get_ptype_md(S1_map[s][i], S1_map[s][j], &params->model_details);
    int penalty = params->MLintern[type];
    best = std::min(best, C[i][j] + penalty);
}

                  M1[i][j] = best;
              }
  
              // ---- F[i][j] ----
              {
                  int best = INF_ENERGY;
  
                  if (i + 1 <= j)
                      best = std::min(best, F[i + 1][j]);
  
                  for (int k = i + theta + 1; k <= j; ++k) {
                      if (C[i][k] < INF_ENERGY) {
                          int cand = C[i][k];
                          if (k + 1 <= j && F[k + 1][j] < INF_ENERGY) cand += F[k + 1][j];
                          if (cand < best) best = cand;
                      }
                  }
  
                  F[i][j] = best;
              }
          }
      }
  }
  
        

  int RNAEnergyEvaluator::hairpin_energy(int i, int j, int s) const {
    int u = j - i - 1;
    int type = vrna_get_ptype_md(S1_map.at(s)[i], S1_map.at(s)[j], &params->model_details);
    return vrna_E_hairpin(u, type,
                          S1_map.at(s)[i + 1], S1_map.at(s)[j - 1],
                          vrna_strands.at(s).c_str() + (i - 1),
                          params);
}


int RNAEnergyEvaluator::interior_loop_energy(int i, int j, int k, int l, int s) const {
    int u1 = k - i - 1;
    int u2 = j - l - 1;
    int type = vrna_get_ptype_md(S1_map.at(s)[i], S1_map.at(s)[j], &params->model_details);
    int type2 = vrna_get_ptype_md(S1_map.at(s)[l], S1_map.at(s)[k], &params->model_details);
    return E_IntLoop(u1, u2, type, type2, 
        S1_map.at(s)[i + 1], S1_map.at(s)[j - 1],
        S1_map.at(s)[k - 1], S1_map.at(s)[l + 1],
        params);

}

int RNAEnergyEvaluator::exterior_loop_energy(int i, int j, int s) const {
    int type = vrna_get_ptype_md(S1_map.at(s)[i], S1_map.at(s)[j], &params->model_details);
    return E_Stem(type, S1_map.at(s)[i + 1], S1_map.at(s)[j - 1], 1, params);
}

int RNAEnergyEvaluator::multiloop_stem_energy(int i, int j, int s, bool is_ext) const {
    int type = vrna_get_ptype_md(S1_map.at(s)[i], S1_map.at(s)[j], &params->model_details);
    return E_Stem(type, S1_map.at(s)[i + 1], S1_map.at(s)[j - 1], is_ext ? 1 : 0, params);
}

int RNAEnergyEvaluator::stem_energy(int i, int j, int s) const {
    int type = vrna_get_ptype_md(S1_map.at(s)[i], S1_map.at(s)[j], &params->model_details);
    return E_Stem(type, S1_map.at(s)[i + 1], S1_map.at(s)[j - 1], 1, params);
}