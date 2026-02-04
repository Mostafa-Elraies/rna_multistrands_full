// validate_cofold.cpp
// Goal: Compare a 2-strand strand_soup folding vs ViennaRNA cofold MFE.

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <algorithm>

#define STRAND_SOUP_NO_MAIN
#include "strand_soup.cpp"

// ViennaRNA
extern "C" {
  #include <ViennaRNA/fold_compound.h>
  #include <ViennaRNA/model.h>
  #include <ViennaRNA/mfe.h>
  #include <ViennaRNA/utils.h>
}

// ----------------------------- helpers -----------------------------

// Build dot-bracket per strand and join by '&'.
// strands map contains sequences with '$' at index 0 (1-based positions in pairs).
static std::string dotbracket_from_output(
    const std::unordered_map<int, std::string>& strands,
    const output_backtrack& bt
){
  const int S = (int)bt.list_of_sequences.size();
  if (S == 0) return "";

  // Per-strand dot strings, in the REALIZED order of bt.list_of_sequences
  std::vector<std::string> db(S);
  for (int x = 1; x <= S; ++x) {
    const int seq_id = bt.list_of_sequences[x - 1];

    auto it = strands.find(seq_id);
    if (it == strands.end())
      throw std::runtime_error("dotbracket_from_output: seq_id not found in strands");

    // length without '$'
    const int n = std::max(0, (int)it->second.size() - 1);
    db[x - 1] = std::string(n, '.');
  }

  // Place parentheses
  for (const auto& p : bt.list_of_pairs) {
    if (p.size() != 4) continue;

    int x = p[0], i = p[1];
    int y = p[2], j = p[3];

    if (x < 1 || x > S || y < 1 || y > S)
      throw std::runtime_error("dotbracket_from_output: pair strand index out of range");

    if (i < 1 || i > (int)db[x - 1].size() || j < 1 || j > (int)db[y - 1].size())
      throw std::runtime_error("dotbracket_from_output: pair nucleotide index out of range");

    // detect double-use early
    if (db[x - 1][i - 1] != '.' || db[y - 1][j - 1] != '.')
      throw std::runtime_error("dotbracket_from_output: nucleotide used in >1 pair");

    if (x == y) {
      // Same strand: enforce i<j
      if (i < j) { db[x - 1][i - 1] = '('; db[y - 1][j - 1] = ')'; }
      else if (j < i) { db[x - 1][j - 1] = '('; db[y - 1][i - 1] = ')'; }
      else throw std::runtime_error("dotbracket_from_output: i==j in same-strand pair");
    } else {
      // Cross-strand: put '(' on earlier strand in realized order
      if (x < y) { db[x - 1][i - 1] = '('; db[y - 1][j - 1] = ')'; }
      else       { db[y - 1][j - 1] = '('; db[x - 1][i - 1] = ')'; }
    }
  }

  // Join with '&'
  std::string out = db[0];
  for (int t = 1; t < S; ++t) out += "&" + db[t];
  return out;
}

// Run Vienna cofold MFE on "A&B" and return (structure, mfe)
static std::pair<std::string, float> vienna_cofold_mfe(
    const std::string& A,
    const std::string& B
) {
  std::string dimer = A + "&" + B;

  vrna_md_t md;
  vrna_md_set_default(&md);

  vrna_fold_compound_t* fc = vrna_fold_compound(dimer.c_str(), &md, VRNA_OPTION_DEFAULT);
  if (!fc) throw std::runtime_error("Vienna: vrna_fold_compound returned null");

  std::string structure(dimer.size(), '.');
  float mfe = vrna_mfe(fc, structure.data());

  vrna_fold_compound_free(fc);
  return {structure, mfe};
}

// ----------------------------- main test -----------------------------

int main() {
  try {
    struct Test { std::string A, B; };
    std::vector<Test> tests = {
      {"GCGCGC", "GCGCGC"},     // strong duplex
      {"GGGCCC", "GGGCCC"},     // strong duplex
      {"AAAAAA", "AAAAAA"},     // should NOT pair (expect mostly dots)
      {"CCCCCC", "CCCCCC"},     // should NOT pair
      {"AAAAAA", "CCCCCC"},     // should NOT pair (A-C not canonical)
      {"GCAAAA", "UUUGCU"},     // mixed
    };

    for (const auto& T : tests) {
      const std::string& A = T.A;
      const std::string& B = T.B;

      std::cout << "\n============================\n";
      std::cout << "A = " << A << "\nB = " << B << "\n";

      // Build 2-strand soup (with '$')
      std::unordered_map<int, std::string> strands;
      strands[1] = with_dollar(A);
      strands[2] = with_dollar(B);

      // Single-strand evaluator
      RNAEnergyEvaluator evaluator(strands);
      vrna_param_t* params = evaluator.get_params();

      // Allocate multi-strand DP matrices for m_start=2
      int m_start = 2;
      int m_size = m_start - 1;          // => 1
      int s_size = (int)strands.size();  // 2
      int r_size = (int)strands.size();  // 2
      int i_size = (int)A.size();
      int j_size = (int)B.size();
      int c_size = 2;

      Matrix5D C(m_size, s_size, i_size, r_size, j_size);
      Matrix5D M(m_size, s_size, i_size, r_size, j_size);
      Matrix5D M1_multi(m_size, s_size, i_size, r_size, j_size);
      Matrix6D F(m_size, s_size, i_size, r_size, j_size, c_size);

      // Fill DP
      MainAuxiliaryMatrix(strands, C, M, M1_multi, F, evaluator, params);

      // Backtrack
      auto structs = full_backtrack(strands, C, M, M1_multi, F, evaluator, params, theta);

      std::string my_db = "(no structure)";
      if (!structs.empty()) {
        my_db = dotbracket_from_output(strands, structs[0]);
      }

      // Vienna cofold
      auto [vrna_db, vrna_mfe] = vienna_cofold_mfe(A, B);

      // Print
      std::cout << "\n--- My model (2-strand) ---\n";
      std::cout << "structure: " << my_db << "\n";

      std::cout << "\n--- ViennaRNA cofold ---\n";
      std::cout << "structure: " << vrna_db << "\n";
      std::cout << "MFE: " << vrna_mfe << " kcal/mol\n";
    }

    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Fatal: " << e.what() << "\n";
    return 1;
  }
}
