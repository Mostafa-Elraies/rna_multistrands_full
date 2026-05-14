/**
 * @file RNAEnergyEvaluator.cpp
 * @brief This acts as a wrapper to interface with Vienna RNA package to extract 
 * thermodynamic parameters. 
 *
 * Here are also the single strands folding
 * 
 */


#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <cassert>
#include "utilities.hpp"
#include "RNAEnergyEvaluator.hpp"
extern "C" {
    #include <ViennaRNA/fold.h>
    #include <ViennaRNA/params/basic.h>
    #include <ViennaRNA/fold_compound.h>
    #include <ViennaRNA/loop_energies.h>
    #include <ViennaRNA/loops/all.h>
    #include <ViennaRNA/loops/multibranch.h>
    #include <ViennaRNA/loops/external.h>
  
    #include <ViennaRNA/eval.h>
    #include <ViennaRNA/utils.h>
}

// E_MLstem with no dangles (-d 0): MLintern[type] + TerminalAU for non-GC pairs (type>2).
static inline int E_MLstem_nd(int type, const vrna_param_t* p) {
    if (type == 0) return 0;
    return p->MLintern[type] + (type > 2 ? p->TerminalAU : 0);
}

RNAEnergyEvaluator::RNAEnergyEvaluator(const std::unordered_map<int, std::string> &input_strands) {
    strands_1b = input_strands; // with '$'
  
    vrna_md_t md;
    vrna_md_set_default(&md);

    // ── Conditions ──────────────────────────────────────────────
    md.salt        = 1.021;   // M  (E. coli physiological: 0.15 | default: 1.021)
    md.temperature = 37.0;   // °C (E. coli: 37 | default: 37)
    // ────────────────────────────────────────────────────────────

    params = vrna_params(&md);
    for (const auto &[s, seq] : input_strands) {
        if (seq.empty() || seq[0] != '$') {
            throw std::runtime_error("RNAEnergyEvaluator: strand " + std::to_string(s) + " is missing leading '$'");
        }
    }
    for (const auto &[s, seq_with_dollar] : strands_1b) {
  
      // Removing the $...
      const std::string clean =
          (!seq_with_dollar.empty() && seq_with_dollar[0] == '$')
          ? seq_with_dollar.substr(1)
          : seq_with_dollar;
  
      vrna_strands[s] = clean;
      len_map[s] = (int)clean.size();
  
      // build fold compound
      vrna_fold_compound_t *fc = vrna_fold_compound(clean.c_str(), &md, VRNA_OPTION_MFE);
      fold_compounds[s] = fc;
  
      // fill mfe matrices (Vienna)
      vrna_mfe(fc, nullptr);
  
      int len = (int)clean.size();
      
static constexpr int INF_INT = std::numeric_limits<int>::max() / 4;


  
        Matrix C(len + 1, std::vector<int>(len + 1, INF_INT));
        Matrix M(len + 1, std::vector<int>(len + 1, INF_INT));
        Matrix M1(len + 1, std::vector<int>(len + 1, INF_INT));
        Matrix F(len + 1, std::vector<int>(len + 1, INF_INT));

  
      C_s[s]  = std::move(C);
      M_s[s]  = std::move(M);
      M1_s[s] = std::move(M1);
      F_s[s]  = std::move(F);
  
      //  Vienna encoding pointer:
      S1_map[s] = fc->sequence_encoding;
  
      // single-strand DP for this strand
      fill_single_strand(s);
    }
  }
  RNAEnergyEvaluator::~RNAEnergyEvaluator() {
    for (auto &kv : fold_compounds) {
      if (kv.second) vrna_fold_compound_free(kv.second);
      kv.second = nullptr;
    }
    if (params) {
      free(params);
      params = nullptr;
    }
  }
      
  
  void RNAEnergyEvaluator::fill_single_strand(int s) {
    const std::string& seq = strands_1b.at(s);   // with $ at [0]
    auto& F  = F_s[s];
    auto& C  = C_s[s];
    auto& M  = M_s[s];
    auto& M1 = M1_s[s];
    const short* S1 = S1_map.at(s);
    static constexpr int INF_INT = std::numeric_limits<int>::max() / 4;

    const int n = (int)seq.size() - 1;  

    // Initialization:
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            C[i][j]  = INF_INT;
            M1[i][j] = INF_INT;
            F[i][j]  = INF_INT;

            M[i][j]  = (i > j) ? 0 : INF_INT;
        }
    }

    // -----------------------------
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (i >= j) F[i][j] = 0;
        }
    }

    // -----------------------------
    // As we did in the 6D matrix filling, we want to fill the small problems first:
    // span = j-i
    // and then we go for bigger and bigger problems. 

    // -----------------------------
    for (int span = 0; span <= n; ++span) {
        for (int i = n; i >= 1; --i) {
            int j = i + span;
            if (j < 1 || j > n) continue;

            // 1) C[i][j] : closed pair region

            if (j >= i + theta + 1 && can_pair(seq[i], seq[j])) {
                int best = INF_INT;

                // (a) hairpin
                best = std::min(best, hairpin_energy(i, j, s));

                // (b) interior loops: (i,j) encloses (p,q)
                for (int p = i + 1; p < j - theta - 1; ++p) {
                    for (int q = p + theta + 1; q < j; ++q) {
                        if (!can_pair(seq[p], seq[q])) continue;
                        int inner = C[p][q];
                        if (inner >= INF_INT) continue;

                        int e_int = interior_loop_energy(i, j, p, q, s);
                        best = std::min(best, inner + e_int);
                    }
                }

                // (c) multiloop closure: MLclosing + MLintern[type(i,j)] + M1[i+1][u] + M[u+1][j-1]
                // MLintern here is the stem cost for the closing pair (i,j) itself
                {
                    int type_ij = vrna_get_ptype_md(S1[i], S1[j], &params->model_details);
                    int mlintern_ij = (type_ij != 0) ? E_MLstem_nd(type_ij, params) : INF_INT;
                    for (int u = i + 1; u <= j - 1; ++u) {
                        int left  = M1[i + 1][u];
                        int right = M[u + 1][j - 1];

                        if (left >= INF_INT || right >= INF_INT) continue;

                        int cand = left + right + params->MLclosing + mlintern_ij;
                        best = std::min(best, cand);
                    }
                }

                C[i][j] = best;
            } else {
                C[i][j] = INF_INT;
            }

            // 2) M1[i][j] 
            {
                int best = INF_INT;

                // (a) i unpaired 
                if (i + 1 <= j && M1[i + 1][j] < INF_INT)
                    best = std::min(best, M1[i + 1][j] + params->MLbase);

                // (b) start helix C[i][u] and pay cost b, in vienne it's --> MLintern(type(i,u))
                if (C[i][j] < INF_INT) {
                    int type = vrna_get_ptype_md(S1[i], S1[j], &params->model_details);
                    int penalty = E_MLstem_nd(type, params);
                    best = std::min(best, C[i][j] + penalty);
                }

                M1[i][j] = best;
            }

            // 3) M[i][j]
            {
                // i == j: single trailing unpaired base costs MLbase (= 0 here). Must be 0
                // so that multiloop arms with one unpaired trailing base are not blocked.
                int best = (i >= j) ? 0 : INF_INT;

                if (i <= j) {
                    // (a) i unpaired
                    if (i + 1 <= j && M[i + 1][j] < INF_INT)
                        best = std::min(best, M[i + 1][j] + params->MLbase);

                    // (b) split (C | M): C[i][u] + penalty(i,u) + M[u+1][j]
                    for (int u = i; u < j; ++u) {
                        if (C[i][u] >= INF_INT || M[u + 1][j] >= INF_INT) continue;
                        int type = vrna_get_ptype_md(S1[i], S1[u], &params->model_details);
                        int penalty = E_MLstem_nd(type, params);
                        best = std::min(best, C[i][u] + penalty + M[u + 1][j]);
                    }

                    // (c) same as (b) but without the M continuation on the right — mirroring Vienna's routine
                    for (int u = i; u <= j; ++u) {
                        if (C[i][u] >= INF_INT) continue;
                        int type = vrna_get_ptype_md(S1[i], S1[u], &params->model_details);
                        int penalty = E_MLstem_nd(type, params);
                        best = std::min(best, C[i][u] + penalty);
                    }
                }

                M[i][j] = best;
            }

            // 4) F[i][j] 
            {
                int best = (i >= j) ? 0 : INF_INT;

                // (a) i unpaired
                if (i + 1 <= j && F[i + 1][j] < INF_INT)
                    best = std::min(best, F[i + 1][j]);

                // (b) i paired with k -> C[i][k] + TerminalAU(i,k) + F[k+1][j]
                // External pairs pay TerminalAU (= E_ext_stem with -d 0) for non-GC pairs.
                for (int k = i + theta + 1; k <= j; ++k) {
                    if (C[i][k] >= INF_INT) continue;
                    int type_ik = vrna_get_ptype_md(S1[i], S1[k], &params->model_details);
                    int ext_pen = (type_ik > 2) ? params->TerminalAU : 0;
                    int cand = C[i][k] + ext_pen;
                    if (k + 1 <= j && F[k + 1][j] < INF_INT) cand += F[k + 1][j];
                    best = std::min(best, cand);
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
    static constexpr int INF_INT = std::numeric_limits<int>::max() / 4;

    int u1 = k - i - 1;
    int u2 = j - l - 1;

    if (u1 < 0 || u2 < 0) return INF_INT;
    if (u1 + u2 > 30) return INF_INT;

    int type  = vrna_get_ptype_md(S1_map.at(s)[i], S1_map.at(s)[j], &params->model_details);
    int type2 = vrna_get_ptype_md(S1_map.at(s)[l], S1_map.at(s)[k], &params->model_details);

    if (type == 0 || type2 == 0) return INF_INT;

    return E_IntLoop(u1, u2, type, type2,
                     S1_map.at(s)[i + 1], S1_map.at(s)[j - 1],
                     S1_map.at(s)[k - 1], S1_map.at(s)[l + 1],
                     params);
}
int RNAEnergyEvaluator::interior_loop_energy_cross(int s, int i, int r, int j, int k, int l) const {
    // Note: s == r is valid in soupfold (two copies of the same strand type).
    const auto& S1s = S1_map;

    const int ns = len_map.at(s);
    const int nr = len_map.at(r);

    if (i < 1 || i > ns) return inf_energy;
    if (k < 1 || k > ns) return inf_energy;
    if (j < 1 || j > nr) return inf_energy;
    if (l < 1 || l > nr) return inf_energy;

    if (!(i < k && l < j)) return inf_energy;

    int u1 = k - i - 1;
    int u2 = j - l - 1;
    if (u1 < 0 || u2 < 0) return inf_energy;

    // Vienna-style interior loop size limit
    if (u1 + u2 > 30) return inf_energy;

    // Need adjacent nucleotides for E_IntLoop
    if (i + 1 > ns) return inf_energy;
    if (k - 1 < 1)  return inf_energy;
    if (j - 1 < 1)  return inf_energy;
    if (l + 1 > nr) return inf_energy;

    int type  = vrna_get_ptype_md(S1s.at(s)[i], S1s.at(r)[j], &params->model_details);
    int type2 = vrna_get_ptype_md(S1s.at(r)[l], S1s.at(s)[k], &params->model_details);

    // Reject illegal closing pairs explicitly
    if (type == 0 || type2 == 0) return inf_energy;

    int si1 = S1s.at(s)[i + 1];
    int sj1 = S1s.at(r)[j - 1];
    int sp1 = S1s.at(s)[k - 1];
    int sq1 = S1s.at(r)[l + 1];

    return E_IntLoop(u1, u2, type, type2, si1, sj1, sp1, sq1, params);
}


int RNAEnergyEvaluator::exterior_loop_energy(int i, int j, int s) const {
    int type = vrna_get_ptype_md(S1_map.at(s)[i], S1_map.at(s)[j], &params->model_details);
    int si1  = S1_map.at(s)[i + 1];
    int sj1  = S1_map.at(s)[j - 1];
    return vrna_E_exterior_stem(type, si1, sj1, params);
}
int RNAEnergyEvaluator::multiloop_stem_energy(int i, int j, int s, bool /*is_ext*/) const {
    // "Stem energy when a helix is inside a multiloop"
    int type = vrna_get_ptype_md(S1_map.at(s)[i], S1_map.at(s)[j], &params->model_details);
    return E_MLstem_nd(type, params);
}

int RNAEnergyEvaluator::stem_energy(int i, int j, int s) const {
    int type = vrna_get_ptype_md(S1_map.at(s)[i], S1_map.at(s)[j], &params->model_details);
    int si1  = S1_map.at(s)[i + 1];
    int sj1  = S1_map.at(s)[j - 1];
    return vrna_E_exterior_stem(type, si1, sj1, params);
}