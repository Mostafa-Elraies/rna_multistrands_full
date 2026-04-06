/**
 * @file strand_soup.cpp
 * @brief Implements the turner enrgy model of the strand_soup algorithm for RNA secondary structure prediction.
 *
 * This file contains the implementation of the dynamic programming approach 
 * for RNA structure prediction, including the matrix filling and the backtracking.
 * 
 * 
 * COMPILATION of this file for testing:
 * Go to the strand_soup.ccp directory where the Makefile should also be and run the following command in the terminal:
 *      make clean
 *      make strand_soup
 *      
 * 
 * After compilation, run the program with:
 *      ./strand_soup.exe
 *
 * Dependencies:
 * "global_variables.hpp" for global variables.
 * "RNAEnergyEvaluator.hpp" for the single strand folding and extracting energy parameters from Vienna RNA package.
 * "utilities.hpp" for helper functions.
 * 
 * The linking is done by the Makefile
 *
 * 
 * @date 2026-02
 */


 #include <iostream>
 #include <fstream>
 #include <vector>
 #include <string>
 #include <cstdlib>
 #include <ctime>
 #include <limits>
 #include <chrono> 
 #include <iomanip> 
 #include <sstream>
 #include <cmath>
 #include <algorithm>
 #include <unordered_set>
 #include "RNAEnergyEvaluator.hpp"
 #include "utilities.hpp"
 
 extern "C" {
     #include <ViennaRNA/loop_energies.h>
     #include <ViennaRNA/utils.h>
     #include <ViennaRNA/model.h>
   }
 // just for the colors in the terminal
 #define RESET   "\033[0m"
 #define RED     "\033[31m"
 #define GREEN   "\033[32m"

// E_MLstem with no dangles (-d 0): MLintern[type] + TerminalAU for non-GC pairs (type>2).
// Matches ViennaRNA's E_MLstem(type, -1, -1, params).
static inline float E_MLstem_nd(int type, const vrna_param_t* p) {
    if (type == 0) return 0.0f;
    return (float)(p->MLintern[type] + (type > 2 ? p->TerminalAU : 0));
}
 
 
 
 //======================================    The Classes    ==============================================//
 
 
 
 /**
  * @class Matrix6D
  * @brief A class representing a 6D matrix used in the strand soup algorithm.
  * 
  * This class handles a 6-dimensional matrix with specific sizes for each dimension.
  * The matrix corresponds to the full energy matrix (free) used in the RNA folding calculations.
  * It provides methods for accessing and modifying the matrix elements.
  * There are 6 dimensions:
  * - m: Number of dimensions needed in total for m. If we start with m_start strands, we need m = m_start-1 dimensions (easy to see with m_start = 2. once s and r are selected we only have the case m=0,so this dim is of size 2-1 = 1
  * - s: Strand type count (1-based)
  * - i: Max Nucleotide count in one of the strands (1-based)
  * - r: Strand type count (1-based)
  * - j: Max Nucleotide count in one of the strands (1-based)
  * - c: Connectivity bit (0 or 1)
  * 
  */
 
 
 
 
 class Matrix6D{
     public:
     /**
      * @brief Constructor that initializes a 6D matrix with specified sizes.
      * 
      * @param m_size Size of the first dimension. Strand count in the soup -2 (the strands s and r are not counted)
      * @param s_size Size of the second dimension : Strand type count (1-based)
      * @param i_size Size of the third dimension : Nucleotide count in the first strand (1-based)
      * @param r_size Size of the fourth dimension : Strand type count (1-based)
      * @param j_size Size of the fifth dimension : Nucleotide count in the second strand (1-based)
      * @param c_size Size of the sixth dimension = 2  : Connectivity bit (0 or 1).
      */
     Matrix6D(int m_size, int s_size, int i_size, int r_size, int j_size, int c_size)
         : m_size(m_size), s_size(s_size), i_size(i_size), r_size(r_size), j_size(j_size), c_size(c_size) {
         data = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>>>(
             m_size, std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>>(
             s_size, std::vector<std::vector<std::vector<std::vector<float>>>>(
             i_size + 2, std::vector<std::vector<std::vector<float>>>( //! we add one column before and after for border effect
             r_size, std::vector<std::vector<float>>(
             j_size + 2, std::vector<float>( //! we add one column before and after for border effect
             c_size, 0))))));
     }
 
     /**
      * @brief Accessor method for the matrix elements.
      * 
      * @param m Index of the first dimension
      * @param s Index of the second dimension
      * @param i Index of the third dimension
      * @param r Index of the fourth dimension
      * @param j Index of the fifth dimension
      * @param c Index of the sixth dimension
      * @return Reference to the element in the matrix
      * 
      * Throws an out_of_range exception if indices are out of bounds.
      */
 
 float& operator()(int m, int s, int i, int r, int j, int c) {
     if (s <= 0 || s > s_size) {
         std::ostringstream oss;
         oss << "Index s out of range: s=" << s << " s_size=" << s_size;
         throw std::out_of_range(oss.str());
     }
     if (r <= 0 || r > r_size) {
         std::ostringstream oss;
         oss << "Index r out of range: r=" << r << " r_size=" << r_size;
         throw std::out_of_range(oss.str());
     }
     if (i < 0 || i >= i_size + 2) {
         std::ostringstream oss;
         oss << "Index i out of range: i=" << i
             << " allowed=[0.." << (i_size + 1) << "]"
             << " (i_size=" << i_size << ")"
             << " at F(m=" << m << ", s=" << s << ", i=" << i
             << ", r=" << r << ", j=" << j << ", c=" << c << ")";
         throw std::out_of_range(oss.str());
     }
     if (j < 0 || j >= j_size + 2) {
         std::ostringstream oss;
         oss << "Index j out of range: j=" << j
             << " allowed=[0.." << (j_size + 1) << "]"
             << " (j_size=" << j_size << ")"
             << " at F(m=" << m << ", s=" << s << ", i=" << i
             << ", r=" << r << ", j=" << j << ", c=" << c << ")";
         throw std::out_of_range(oss.str());
     }
     if (m < 0 || m >= m_size) {
         std::ostringstream oss;
         oss << "Index m out of range: m=" << m << " m_size=" << m_size;
         throw std::out_of_range(oss.str());
     }
     if (c < 0 || c >= c_size) {
         std::ostringstream oss;
         oss << "Index c out of range: c=" << c << " c_size=" << c_size;
         throw std::out_of_range(oss.str());
     }
 
     return data[m][s-1][i][r-1][j][c];
 }
 
 
     int get_m_size() const { return m_size; }
     int get_s_size() const { return s_size; }
     int get_i_size() const { return i_size; } // i_size is the real number of bases in the strand (It's the size of the 3rd dimension -2)
     int get_r_size() const { return r_size; }
     int get_j_size() const { return j_size; } // j_size is the real number of bases in the strand (It's the size of the 5th dimension -2)
     int get_c_size() const { return c_size; }
 
 private:
     int m_size, s_size, i_size, r_size, j_size, c_size;
     std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>>> data;
 };
 
 
 /**
  * @class Matrix5D
  * @brief A class representing a 5D matrix used in the strand soup algorithm.
  * 
  * This class handles a 5-dimensional matrix with specific sizes for each dimension.
  * The matrix corresponds to the C (close) and the M (multi) matrices used in the full model RNA prediction.
  * It provides methods for accessing and modifying the matrix elements.
  * There are 5 dimensions:
  * - m: Number of dimensions needed in total for m. If we start with m_start strands, we need m = m_start-1 dimensions (easy to see with m_start = 2. once s and r are selected we only have the case m=0,so this dim is of size 2-1 = 1
  * - s: Strand type count (1-based)
  * - i: Max Nucleotide count in one of the strands (1-based)
  * - r: Strand type count (1-based)
  * - j: Max Nucleotide count in one of the strands (1-based)
  */
 class Matrix5D{
     public:
     /**
      * @brief Constructor that initializes a 5D matrix with specified sizes.
      * 
      * @param m_size Size of the first dimension. Strand count in the soup -2 (the strands s and r are not counted)
      * @param s_size Size of the second dimension : Strand type count (1-based)
      * @param i_size Size of the third dimension : Nucleotide count in the first strand (1-based)
      * @param r_size Size of the fourth dimension : Strand type count (1-based)
      * @param j_size Size of the fifth dimension : Nucleotide count in the second strand (1-based)
      */
     Matrix5D(int m_size, int s_size, int i_size, int r_size, int j_size)
         : m_size(m_size), s_size(s_size), i_size(i_size), r_size(r_size), j_size(j_size) {
         data = std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>>(
             m_size, std::vector<std::vector<std::vector<std::vector<float>>>>(
             s_size, std::vector<std::vector<std::vector<float>>>(
             i_size + 2, std::vector<std::vector<float>>( //! we add one column before and after for border effect
             r_size, std::vector<float>(
             j_size + 2,  //! we add one column before and after for border effect
             0)))));
     }
 
     /**
      * @brief Accessor method for the matrix elements.
      * 
      * @param m Index of the first dimension
      * @param s Index of the second dimension
      * @param i Index of the third dimension
      * @param r Index of the fourth dimension
      * @param j Index of the fifth dimension
      * @return Reference to the element in the matrix
      * 
      * Throws an out_of_range exception if indices are out of bounds.
      */
     float& operator()(int m, int s, int i, int r, int j) {
         if (s <= 0 || s > s_size) {
             throw std::out_of_range("Index s is out of range");
         }
         if (r <= 0 || r > r_size) {
             throw std::out_of_range("Index r is out of range");
         }
         if (i < 0 || i >= i_size + 2) {
             throw std::out_of_range("Index i is out of range");
         }
         if (j < 0 || j >= j_size + 2) {
             throw std::out_of_range("Index j is out of range");
         }
         return data[m][s-1][i][r-1][j]; // 1-based indexing for strands and bases
     }
 
     int get_m_size() const { return m_size; }
     int get_s_size() const { return s_size; }
     int get_i_size() const { return i_size; } // i_size is the real number of bases in the strand (It's the size of the 3rd dimension -2)
     int get_r_size() const { return r_size; }
     int get_j_size() const { return j_size; } // j_size is the real number of bases in the strand (It's the size of the 5th dimension -2)
    
 private:
     int m_size, s_size, i_size, r_size, j_size;
     std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> data;
 };
 
 
 
 
 
 /**
  * @class output_backtrack
  * @brief Class to manage the backtracking output from RNA folding calculations of the strand_soup algorithm.
  * 
  * This class stores the sequences and pairs of nucleotides that are part of the optimal secondary structure.
  * It provides methods to add sequences and pairs, shift indices, merge outputs from subproblems, and print the results.
  * 
  */
 class output_backtrack{
 
     public:
     std::vector<int> list_of_sequences;
     std::vector<std::vector<int>> list_of_pairs;
 
     /**
      * @brief Add a sequence to the end of the list of sequences.
      * 
      * @param sequence The sequence to be added.
      */
     void add_sequence(int sequence){
         list_of_sequences.push_back(sequence);
     }
 
     /**
      * @brief Add a sequence to the front of the list  of sequences.
      * 
      * @param sequence The sequence to be added at the front.
      */
     void add_sequence_front(int sequence){
         list_of_sequences.insert(list_of_sequences.begin(),sequence);
     }
 
     /**
      * @brief Add a pair of nucleotides to the output.
      * 
      * @param x Index of the first strand
      * @param i Index of the first nucleotide in strand x
      * @param y Index of the second strand
      * @param j Index of the second nucleotide in strand y
      */
     void add_pair(int x, int i, int y, int j) {
         // std::cout << "Adding pair: (" << x << "," << i << "," << y << "," << j << ")" << std::endl;
         list_of_pairs.push_back({x, i, y, j});
     }
 
     /**
      * @brief Shift the indices of the pairs by a given value.
      * 
      * This function is usefull when backtracking the output of a subproblem.
      * 
      * @param shift The value by which to shift the indices.
      */
     void shift(int shift){
         for (int i=0; i<int(list_of_pairs.size()); i++){
             list_of_pairs[i][0] += shift;
             list_of_pairs[i][2] += shift;
         }
     }
 
     /**
      * @brief Merge the current output with another output.
      * 
      * @param sub_output The output of the sub-problem to be merged.
      */
     void merge(const output_backtrack& sub_output){
         for (const auto& seq : sub_output.list_of_sequences){
             list_of_sequences.push_back(seq);
         }
         for (const auto& pair : sub_output.list_of_pairs){
             list_of_pairs.push_back({pair[0], pair[1], pair[2], pair[3]});
         }}
 
     /**
      * @brief Print the sequences and pairs in the output.
      */
     void print() const{
         std::cout << "Ordered list of sequences: ";
         for (const auto& seq : list_of_sequences){
             std::cout << seq << " ";
         }
         std::cout << std::endl;
 
         std::cout << "Pairs:  ";
         for (const auto& pair : list_of_pairs){
             std::cout << "(" << pair[0] << "," << pair[1] << "," << pair[2] << "," << pair[3] << ") ";
         }
         std::cout << std::endl;
     }
 
     void print_dot_bracket(const std::unordered_map<int, std::string>& strands) const {
         int n = (int)list_of_sequences.size();
         // Build per-strand sequences (strip leading '$')
         std::vector<std::string> seqs;
         std::vector<int> lens;
         for (int seq_id : list_of_sequences) {
             const std::string& raw = strands.at(seq_id);
             std::string clean = (raw[0] == '$') ? raw.substr(1) : raw;
             seqs.push_back(clean);
             lens.push_back((int)clean.size());
         }
         // pair_of[slot][pos] = {partner_slot, partner_pos}, 1-based; {0,0} = unpaired
         std::vector<std::vector<std::pair<int,int>>> pair_of(n + 1);
         for (int k = 1; k <= n; ++k)
             pair_of[k].assign(lens[k-1] + 1, {0, 0});
         for (const auto& p : list_of_pairs) {
             int sx = p[0], px = p[1], sy = p[2], py = p[3];
             if (sx >= 1 && sx <= n && px >= 1 && px <= lens[sx-1] &&
                 sy >= 1 && sy <= n && py >= 1 && py <= lens[sy-1]) {
                 pair_of[sx][px] = {sy, py};
                 pair_of[sy][py] = {sx, px};
             }
         }
         std::string seq_str, db_str;
         for (int k = 1; k <= n; ++k) {
             if (k > 1) { seq_str += '&'; db_str += '&'; }
             for (int p = 1; p <= lens[k-1]; ++p) {
                 seq_str += seqs[k-1][p-1];
                 auto [ps, pp] = pair_of[k][p];
                 if (ps == 0) {
                     db_str += '.';
                 } else if (ps > k || (ps == k && pp > p)) {
                     db_str += '(';
                 } else {
                     db_str += ')';
                 }
             }
         }
         std::cout << seq_str << '\n' << db_str << '\n';
     }
 };
 
 
 
 
 
 //======================================    Energy matrix functions    ==============================================//
 /**
  * @brief Handles the free structure case 'F' with the allonce of disconnect of outermost fragments.
  * 
  * This function computes the minimum energy for the bubble case in the RNA strand soup problem. It considers the 4 cases and returns the minimum energy.
  * 
  * @param m The number of strands remaining in the soup.
  * @param s The index of the starting strand. (1-based)
  * @param i The index of the starting nucleotide in strand s. (1-based)
  * @param r The index of the ending strand. (1-based)
  * @param j The index of the ending nucleotide in strand r. (1-based)
  * @param c The connectivity bit.
  * @param strands A dictionary of strands (index, sequence).
  * @param C An energy matrix for closed structure.
  * @param M  An energy matrix for multiple loop structure.
  * @param F An energy matrix for free structure.
  * @param evaluator Energy evaluator with ViennaRNA functions.
  * @return float The minimum energy.
  */
 
  float GeneralCaseMinimization (int m, int s, int i, int r, int j, int c,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M, Matrix6D& F,
     const RNAEnergyEvaluator& evaluator)
 {
     // s == r is valid: two copies of the same strand type can interact (soupfold).
     float min_value = inf_energy;
     const int len_s = int(strands.at(s).size()) - 1;
 
     if (i > len_s) return inf_energy;  // We mustn't do the general case when s is empty.
     if (j < 1)     return inf_energy;  // We mustn't do the general case when r is empty.
 // those are handeled in the border case function, main auxiliary.
 
     // --- Case 1 : i unpaired ---
     // Caller ensures i <= |s|-1 and j >= 1 in general case; borders handled outside.
     if (i + 1 <= len_s + 1)  
         min_value = F(m, s, i + 1, r, j, c);
     else
         min_value = inf_energy;
     
     const auto& S1s_F = evaluator.get_S1_map();
     const vrna_param_t* params_F = evaluator.get_params();
     vrna_md_t md_F = params_F->model_details;

     // --- Case 2 : i pairs within strand the first strand 's' (i,k) ---
     for (int k = i + theta + 1; k <= int(strands.at(s).length()) - 1; ++k) {
         if (can_pair(strands.at(s)[i], strands.at(s)[k])) {
             int type_ik = vrna_get_ptype_md(S1s_F.at(s)[i], S1s_F.at(s)[k], &md_F);
             float ext_pen = (type_ik > 2) ? (float)params_F->TerminalAU : 0.0f;
             float value = evaluator.get_C(s)[i][k] + ext_pen + F(m, s, k + 1, r, j, c);
             if (value < min_value) min_value = value;
         }
     }

     // --- Case 3 : i pairs with new strand t at position k ---
     if (m >= 1) {
         for (int t = 1; t <= M.get_s_size(); ++t) {
             // Any strand type is allowed as intermediate (soupfold: types can repeat).
             for (int m1 = 0; m1 < m; ++m1) {
                 int m2 = m - m1 - 1;
                 for (int k = 1; k <= int(strands.at(t).length()) - 1; ++k) {
                     if (can_pair(strands.at(s)[i], strands.at(t)[k])) {
                         int type_ik = vrna_get_ptype_md(S1s_F.at(s)[i], S1s_F.at(t)[k], &md_F);
                         float ext_pen = (type_ik > 2) ? (float)params_F->TerminalAU : 0.0f;
                         float value = C(m1, s, i, t, k) + ext_pen + F(m2, t, k + 1, r, j, c);
                         if (value < min_value) min_value = value;
                     }
                 }
             }
         }
     }

     // --- Case 4 : i pairs with the last strand r at k (only if connected) ---
     if (c == 1) {
         for (int k = 1; k <= j; ++k) {
             if (!can_pair(strands.at(s)[i], strands.at(r)[k])) continue;
             int type_ik = vrna_get_ptype_md(S1s_F.at(s)[i], S1s_F.at(r)[k], &md_F);
             float ext_pen = (type_ik > 2) ? (float)params_F->TerminalAU : 0.0f;
             if (k == int(strands.at(r).length()) - 1) {
                 float value = C(m, s, i, r, k) + ext_pen;
                 if (value < min_value) min_value = value;
             } else {
                 float value = C(m, s, i, r, k) + ext_pen + evaluator.get_F(r)[k + 1][j];
                 if (value < min_value) min_value = value;
             }
         }
     }
 
     return min_value;
 }
 
 /**
  *  @brief Computes the 5D matrix for the M1_multi loop. 
  * 
  * This function computes the minimum energy for the first multi loop in the RNA strand soup problem. It considers the 2 cases and fills the matrix. 
  * 
  * @param m The number of strands remaining in the soup.
  * @param s The index of the starting strand. (1-based)
  * @param i The index of the starting nucleotide in strand s. (1-based)
  * @param r The index of the ending strand. (1-based)
  * @param j The index of the ending nucleotide in strand r. (1-based)
  * @param strands A dictionary of strands (index, sequence).
  * @param C An energy matrix for closed structure.
  * @param M1_multi  An energy matrix for first multi loop used to decompose mulitple loops with M.
  * @param evaluator Energy evaluator with ViennaRNA functions.
  * @param params Parameters from Vienna RNA.
 
  */
 
 
 void M1MultiMatrixMinimization(int m, int s, int i, int r, int j,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M1_multi,
     RNAEnergyEvaluator& evaluator,
     vrna_param_t* params)
 {
     if (i == int(strands.at(s).length()) || j == 0) {
         M1_multi(m,s,i,r,j) = inf_energy;
         return;
     }
 
     float best = inf_energy;
 
     // Case 1: unpaired on left boundary 
     if (i < int(strands.at(s).length()) - 1) {
         best = std::min(best, M1_multi(m,s,i+1,r,j) + params->MLbase);
     }
 
     // Case 2: start with a helix at (i,j)
     const auto& S1s = evaluator.get_S1_map();
     if (can_pair(strands.at(s)[i], strands.at(r)[j])) {
         int type = vrna_get_ptype_md(S1s.at(s)[i], S1s.at(r)[j], &params->model_details);
         if (type != 0) {
             best = std::min(best, C(m,s,i,r,j) + E_MLstem_nd(type, params));
         }
     }
 
     M1_multi(m,s,i,r,j) = best;
 }
 
 /**
  * @brief Computes the 5D matrix for the  closed structure case of the RNA strand soup problem.
  * 
  * This function computes the minimum energy for a bubble case in the RNA strand soup problem. 
  * It considers the 5 cases and fills the matrix enteries with the minimum energy.
  * 
  * @param m The number of strands remaining in the soup.
  * @param s The index of the starting strand. (1-based)
  * @param i The index of the starting nucleotide in strand s. (1-based)
  * @param r The index of the ending strand. (1-based)
  * @param j The index of the ending nucleotide in strand r. (1-based)
  * @param strands A dictionary of strands (index, sequence).
  * @param C An energy matrix for closed structure.
  * @param M1_multi  An energy matrix for first multi loop used to decompose mulitple loops with M.
  * @param F The energy matrix.
  * @param M The multi loop case  matrix.
  * @param evaluator Energy evaluator with ViennaRNA functions.
  * @param params Parameters from Vienna RNA.
  */
 
 
  void ClosedCaseMinimization(int m, int s, int i, int r, int j, 
     const std::unordered_map<int, std::string>& strands, 
     Matrix5D& C, Matrix5D& M1_multi, Matrix6D& F, Matrix5D& M,
     RNAEnergyEvaluator& evaluator, vrna_param_t* params)
 {
     const int len_s = int(strands.at(s).size()) - 1;
     const int len_r = int(strands.at(r).size()) - 1;
 
     if (i < 1 || i > len_s || j < 1 || j > len_r) {
     C(m,s,i,r,j) = inf_energy;
     return;
     }
 
      // s == r is valid: two copies of the same strand type can form a closed pair (soupfold).
     if (!can_pair(strands.at(s)[i], strands.at(r)[j])) {
         C(m, s, i, r, j) = inf_energy;
         return;
     }
     float min_value = inf_energy;

     // MLintern for the closing pair (i,j): ViennaRNA's E_MLstem includes MLintern[rtype[type]]
     // for the closing pair of every multiloop.  We use type_ij directly; for canonical pairs
     // (CG/GC/AU/UA/GU/UG) MLintern[type] == MLintern[rtype[type]], so this is equivalent.
     const auto& S1s_c = evaluator.get_S1_map();
     int type_ij = vrna_get_ptype_md(S1s_c.at(s)[i], S1s_c.at(r)[j], &params->model_details);
     float mlintern_closing = (type_ij != 0) ? E_MLstem_nd(type_ij, params) : 0.0f;

     // --- Case 1 : interior loop (i,j) closed by (k,l) ---
     // Require i < k and l < j to form a proper interior loop.
     // note: we imposed a limit on l-k = 30.
     for (int k = i + 1; k <= int(strands.at(s).length()) - 1; ++k) {
         for (int l = 1; l <= j - 1; ++l) {
             if (k <= i || l >= j) continue;
             if (!can_pair(strands.at(s)[k], strands.at(r)[l])) continue;
             float loop_energy = evaluator.interior_loop_energy_cross(s, i, r, j, k, l);
             float value = C(m, s, k, r, l) + loop_energy;
             if (value < min_value) min_value = value;
         }
     }

     // --- Case 2 : leftmost stem within s (i+1 pairs with k), plus multiloop split ---
     // Note: No can_pair filter on s[i+1]: M1 already handles unpaired leading positions.
     for (int k = i + theta + 1; k <= int(strands.at(s).length()) - 1; ++k) {
         float multiloop_energy = params->MLclosing + mlintern_closing + evaluator.get_M1(s)[i + 1][k];
         float value = multiloop_energy + M(m, s, k + 1, r, j - 1);
         if (value < min_value) min_value = value;
     }

     // --- Case 3 : leftmost stem with new strand t ---
 if (m >= 1 && i + 1 <= len_s && j >= 2) { // j>=2 so j-1>=1
     for (int t = 1; t <= F.get_s_size(); ++t) {
         // Any strand type allowed as intermediate (soupfold: types can repeat).

         for (int m1 = 0; m1 < m; ++m1) {
             int m2 = m - m1 - 1;
             for (int k = 1; k <= int(strands.at(t).length()) - 1; ++k) {
                 float value = params->MLclosing
                             + mlintern_closing
                             + M1_multi(m1, s, i + 1, t, k)
                             + M(m2, t, k + 1, r, j - 1);

                 if (value < min_value) min_value = value;
             }
         }
     }
 }


     // --- Case 4 : leftmost stem with rightmost strand r at k ---
    if (i + 1 <= len_s)
     for (int k = 1; k <= j - 1; ++k) {


         if (k == int(strands.at(r).length()) - 1) {
             float value = params->MLclosing
                         + mlintern_closing
                         + M1_multi(m, s, i + 1, r, k);
             if (value < min_value) min_value = value;
         } else {
             float value = params->MLclosing
                         + mlintern_closing
                         + M1_multi(m, s, i + 1, r, k)
                         + evaluator.get_M(r)[k + 1][j - 1];
             if (value < min_value) min_value = value;
         }
     }
 
 // --- Case 5: disconnect (external) ---
 // The pair (i,j) is in an external loop context here, so TerminalAU applies for non-GC pairs.
 /* extloop gives only TerminalAU, see from the documentation here:
 
 int E_Stem(int type, int si1, int sj1, int extLoop, vrna_param_t *P)
#include <ViennaRNA/eval/exterior.h> Compute the energy contribution of a stem branching off a
loop-region.
This function computes the energy contribution of a stem that branches off a loop region. This can be
the case in multiloops, when a stem branching off increases the degree of the loop but also immediately
interior base pairs of an exterior loop contribute free energy. To switch the behavior of the function
according to the evaluation of a multiloop- or exterior-loop-stem, you pass the flag ‘extLoop’. The
returned energy contribution consists of a TerminalAU penalty if the pair type is greater than 2, dangling
end contributions of mismatching nucleotides adjacent to the stem if only one of the si1, sj1 parameters
is greater than 0 and mismatch energies if both mismatching nucleotides are positive values. Thus, to
avoid incorporating dangling end or mismatch energies just pass a negative number, e.g. -1 to the
mismatch argument.
 

 */
 if (i + 1 <= F.get_i_size() + 1 && j - 1 >= 0) {
     // TerminalAU here is for (i,j) as the closing boundary of the internal external-loop
     // context (c=0 inside). This is a soup-model-specific case with no Zuker equivalent:
     // the pair (i,j) appears at two interfaces — outer (paid by F caller) and inner (here).
     float term_au = (type_ij > 2) ? (float)params->TerminalAU : 0.0f;
     min_value = std::min(min_value, term_au + F(m, s, i + 1, r, j - 1, 0));
 }
 
 
     C(m, s, i, r, j) = min_value;
 }
 
 
 /**
  * @brief Handles the multi loop structure case of the RNA strand soup problem.
  * 
  * This function computes the minimum energy for the multiple loop case in the RNA strand soup problem. 
  * It considers the 4 cases and returns the minimum energy.
  * @param m The number of strands remaining in the soup.
  * @param s The index of the starting strand. (1-based)
  * @param i The index of the starting nucleotide in strand s. (1-based)
  * @param r The index of the ending strand. (1-based)
  * @param j The index of the ending nucleotide in strand r. (1-based)
  * @param strands A dictionary of strands (index, sequence).
  * @param C An energy matrix for closed structure.
  * @param M1_multi  An energy matrix for first multi loop used to decompose mulitple loops with M.
  * @param M The multi loop case  matrix.
  * @param F The full energy matrix.
  * @param evaluator Energy evaluator with ViennaRNA functions.
  * @param params Parameters from Vienna RNA. 
  * @return float The minimum energy.
  */
 
 
  float MultipleCaseMinimization(int m, int s, int i, int r, int j,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& /*M1_multi*/, Matrix5D& M, Matrix6D& F,
     RNAEnergyEvaluator& evaluator, vrna_param_t *params)
 {
    // s == r is valid: two copies of the same strand type can be in a multiloop (soupfold).
 
     float min_value = inf_energy;
 
     const auto& S1s  = evaluator.get_S1_map();
     const int len_s  = int(strands.at(s).size()) - 1;
     const int len_r  = int(strands.at(r).size()) - 1;
 
     // --- Case 1 : i unpaired ---
     if (i + 1 <= len_s) {
         min_value = std::min(
             min_value,
             M(m, s, i + 1, r, j) + params->MLbase
         );
     }
 
     // --- Case 2 : i pairs within s ---
     for (int k = i + theta + 1; k <= len_s; ++k) {
         if (!can_pair(strands.at(s)[i], strands.at(s)[k])) continue;
 
         int type = vrna_get_ptype_md(
             S1s.at(s)[i], S1s.at(s)[k], &params->model_details
         );
         if (type == 0) continue;
 
         min_value = std::min(
             min_value,
             evaluator.get_C(s)[i][k]
             + E_MLstem_nd(type, params)
             + M(m, s, k + 1, r, j)
         );
     }
 
     // --- Case 3 : i pairs with new strand t ---
     if (m >= 1) {
         for (int t = 1; t <= F.get_s_size(); ++t) {
             const int len_t = int(strands.at(t).size()) - 1;
             // Any strand type allowed as intermediate (soupfold: types can repeat).
 
 
             for (int m1 = 0; m1 < m; ++m1) {
                 int m2 = m - m1 - 1;
 
                 for (int k = 1; k <= len_t; ++k) {
                     if (!can_pair(strands.at(s)[i], strands.at(t)[k])) continue;
 
                     int type = vrna_get_ptype_md(
                         S1s.at(s)[i], S1s.at(t)[k], &params->model_details
                     );
                     if (type == 0) continue;
                     min_value = std::min(
                         min_value,
                         C(m1, s, i, t, k)
                         + E_MLstem_nd(type, params)
                         + M(m2, t, k + 1, r, j)
                     );
                 }
             }
         }
     }
 
     // --- Case 4 : i pairs with last strand r at k ---
     // kmax = j (inclusive): j is the last available position on r, not exclusive
     int kmax = std::min(j, len_r);
     for (int k = 1; k <= kmax; ++k) {
         if (!can_pair(strands.at(s)[i], strands.at(r)[k])) continue;

         int type = vrna_get_ptype_md(
             S1s.at(s)[i], S1s.at(r)[k], &params->model_details
         );
         if (type == 0) continue;
         if (k == j) {
             min_value = std::min(
                 min_value,
                 C(m, s, i, r, k) + E_MLstem_nd(type, params)
             );
         } else {
             min_value = std::min(
                 min_value,
                 C(m, s, i, r, k)
                 + E_MLstem_nd(type, params)
                 + evaluator.get_M(r)[k + 1][j]
             );
         }
     }
 
     return min_value;
 }
 
 
 
 /**
  * @brief This is similar to the square case of the full energy matrix, where we check the border case and actually fill the 5D matrix.
  * This function checks for j = 0 and i = |s|+1, if not it redirects to the 4 cases computed by the precedent function. 
  * This function and similar 5D functions will all be called withing the MainAuxiliaryMatrix.
  * @param m The number of strands remaining in the soup.
  * @param s The index of the starting strand. (1-based)
  * @param i The index of the starting nucleotide in strand s. (1-based)
  * @param r The index of the ending strand. (1-based)
  * @param j The index of the ending nucleotide in strand r. (1-based)
  * @param strands A dictionary of strands (index, sequence).
  * @param C An energy matrix for closed structure.
  * @param M1_multi  An energy matrix for first multi loop used to decompose mulitple loops with M.
  * @param M The multi loop case  matrix.
  * @param F The full energy matrix.
  * @param evaluator Energy evaluator with ViennaRNA functions.
  * @param params Parameters from Vienna RNA. 
  */
 
 void MMatrixMinimization(int m, int s, int i, int r, int j,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M1_multi, Matrix5D& M, Matrix6D& F,
     RNAEnergyEvaluator& evaluator,
     vrna_param_t* params)
 {
     if (i == int(strands.at(s).length()) || j == 0) {
         M(m, s, i, r, j) = inf_energy;
     } else {
         M(m, s, i, r, j) =
             MultipleCaseMinimization(m, s, i, r, j, strands, C, M1_multi, M, F, evaluator, params);
     }
 }
 
 
 /**
  * @brief Computes the 6D matrix for the RNA strand soup problem. This corresponds to the square case.
  * Within it the C,M,M1 matrices are filled too before filling the F full energy matrix.
  * 
  * @param strands The dictionary of strands (index, sequence).
  * @param C The closed matrix.
  * @param M The multiloop matrix.
  * @param M1_multi  An energy matrix for first multi loop used to decompose mulitple loops with M.
  * @param F The full matrix (6D).
  * @param evaluator Energy evaluator with ViennaRNA functions.
  * @param params ViennaRNA thermodynamic parameter set.
  */
 void MainAuxiliaryMatrix(std::unordered_map<int, std::string> strands,
     Matrix5D& C, Matrix5D& M,Matrix5D& M1_multi, Matrix6D& F,
     RNAEnergyEvaluator& evaluator,
     vrna_param_t* params) {
         
 
     for (int m = 0; m < F.get_m_size(); m++) {
         for (int s = 1; s <= F.get_s_size(); s++) {
             for (int i = F.get_i_size() +1 ; i >= 0; --i) {
                 static auto t_last = std::chrono::steady_clock::now();
 
                 if (m == 0 && s == 1) {
                     if (i % 5 == 0) { // every 5 rows it prints the time to make sure the code is running in case we enter a long sequence. 
                         auto now = std::chrono::steady_clock::now();
                         auto dt = std::chrono::duration_cast<std::chrono::seconds>(now - t_last).count();
                         t_last = now;
                 
                         std::cerr << "PROGRESS m=" << m
                                   << " s=" << s
                                   << " i=" << i << "/" << (F.get_i_size()+1)
                                   << "  dt=" << dt << "s"
                                   << std::endl;
                     }
                 }
                 
                   
                 for (int r = 1; r <= F.get_r_size(); r++) {
                     const int len_s = int(strands.at(s).size()) - 1;
                      const int len_r = int(strands.at(r).size()) - 1;
                    // s == r is valid in soupfold: two copies of same type can interact.
                 
                     for (int j = 0; j <= F.get_j_size() +1; ++j) {
 
                           
                           
                           bool s_ok_for_closed = (i >= 1 && i <= len_s);
                           bool r_ok_for_closed = (j >= 1 && j <= len_r);
                           
                           if (s_ok_for_closed && r_ok_for_closed) {
                               ClosedCaseMinimization(m,s,i,r,j,strands,C,M1_multi,F,M,evaluator,params);
                               M1MultiMatrixMinimization(m,s,i,r,j,strands,C,M1_multi,evaluator,params);
                               MMatrixMinimization(m,s,i,r,j,strands,C,M1_multi,M,F,evaluator,params);
                           } else {
                               C(m,s,i,r,j) = inf_energy;
                               M(m,s,i,r,j) = inf_energy;
                               M1_multi(m,s,i,r,j) = inf_energy;  
                           }
                           
   
 
                         for (int c = 0; c <= 1; c++) {
                             if (i > len_s) {
                                 // === CASE: strand s is empty ===
                                 if (c == 1) {
                                     F(m, s, i, r, j, c) = inf_energy;
                                     continue;
                                 } else {
                                     if (m == 0) {
                                         if (j == 0) {
                                             F(m, s, i, r, j, c) = 0;
                                             continue;
                                         } else {
                                             F(m, s, i, r, j, c) = evaluator.get_F(r)[1][j];
                                             continue;
                                         }
                                     } else {
                                         // Minimum over all t (different sequences type): M(m-1, t, 1, r, j, 1)
                                         float min_value = inf_energy;
                                         for (int t = 1; t <= M.get_s_size(); t++) {
                                            // Any strand type allowed (soupfold).
 
                                             float val = F(m - 1, t, 1, r, j, 1);
                                             if (val < min_value) min_value = val;
                                         }
                                         F(m, s, i, r, j, c) = min_value;
                                         continue;
                                     }
                                 }
                             } else {
                                 // === CASE: strand r is empty ===
                                 if (j < 1) {
                                     if (c == 1) {
                                         F(m, s, i, r, j, c) = inf_energy;
                                         continue;
                                     } else {
                                         if (m == 0) {
                                             if (i == 0) {
                                                 F(m, s, i, r, j, c) = 0;
                                                 continue;
                                             } else {
                                                 F(m, s, i, r, j, c) =
                                                     evaluator.get_F(s)[i][strands.at(s).length() - 1];
                                                     continue;
                                             }
                                         } else {
                                             // Minimum over all t: M(m-1, s, i, t, |t|-1, 1)
                                             float min_value = inf_energy;
                                             for (int t = 1; t <= M.get_s_size(); t++) {
                                                 int len_t = int(strands.at(t).size()) - 1;
                                                // Any strand type allowed (soupfold).
 
                                                     float val = F(m - 1, s, i, t, len_t, 1);
                                                 if (val < min_value) min_value = val;
                                             }
                                             F(m, s, i, r, j, c) = min_value;
                                             continue;
                                         }
                                     }
                                 } else {
                                     // === GENERAL CASE (bubble case) ===
                                     if (j > len_r) {
                                         F(m, s, i, r, j, c) = inf_energy;
                                         continue;
                                     }
                                     F(m, s, i, r, j, c) =
                                         GeneralCaseMinimization(m, s, i, r, j, c,
                                             strands, C, M, F, evaluator);
                                 }
                             }
                         } 
                     } 
                 } 
             } 
         } 
     } 
 }
 //======================================    Backtrack functions    ==============================================//
 
 
 
 /**
  * @brief Finds all the starting point in the energy matrix, which corresponds to the minimum energy (excluding the border effect).
  * 
  * This function searches through the filled energy matrix (F) and identifies the minimum energy value. It returns all the coordinates (m,s,i,r,j,c) of the minimum energy, excluding the borders.
  * 
  * @param F The 6D energy matrix that stores the computed energy values for all possible configurations.
  * 
  * @return std::vector<std::vector<int>> A vector containing all the starting points (m,s,i,r,j,c) where the minimum energy is located.
  */
 std::vector<std::vector<int>>
 Find_all_start_backtrack(Matrix6D& F, int /*n_strands*/, int c_target = 1) {
     std::vector<std::vector<int>> starting_points;
     float min_value = inf_energy;
 
     // m_required derived from matrix size: F.get_m_size() = m_start - 1, so m_required = m_start - 2.
     int m_required = F.get_m_size() - 1;
     if (m_required < 0) {
         throw std::runtime_error("Find_all_start_backtrack: need at least 2 strands (m_size >= 1)");
     }
 
     int i_start = 1;
     int j_max   = F.get_j_size();  
 
     for (int s = 1; s <= F.get_s_size(); ++s) {
         for (int r = 1; r <= F.get_r_size(); ++r) {
             // s == r is valid: two copies of same type as first/last strand (soupfold).
 
             float val = F(m_required, s, i_start, r, j_max, c_target);
             if (val == inf_energy) continue;
 
             if (val < min_value) {
                 min_value = val;
                 starting_points.clear();
                 starting_points.push_back({m_required, s, i_start, r, j_max, c_target});
             } else if (val == min_value) {
                 starting_points.push_back({m_required, s, i_start, r, j_max, c_target});
             }
         }
     }
 
     if (min_value == inf_energy) {
         throw std::runtime_error("No valid starting point found at required m (all strands) and c=1");
     }
 
     return starting_points;
 }
 
 
 /**
  * ============================================================
  *  Single-strand thermodynamic backtracking functions
  *  (Turner/Vienna model version)
  * ============================================================
  *  These reconstruct the optimal secondary structure of a single RNA using the
  *  thermodynamic matrices: F, C, M, and M1.
  *  Each function returns an output_backtrack object so that the
  *  same data structure can be merged into multistrand backtracking.
  * ============================================================
  */
 
  output_backtrack backtrack_Fs(int s, int s_slot, int i, int j, RNAEnergyEvaluator& evaluator, int theta);
  output_backtrack backtrack_Cs(int s, int s_slot, int i, int j, RNAEnergyEvaluator& evaluator, int theta);
  output_backtrack backtrack_Ms(int s, int s_slot, int i, int j, RNAEnergyEvaluator& evaluator, int theta);
  output_backtrack backtrack_M1s(int s, int s_slot, int i, int j, RNAEnergyEvaluator& evaluator, int theta);
  
  /*
  -------------------------------------------------------------
   *  F-matrix (external / free region)
 *------------------------------------------------------------
   */
  output_backtrack backtrack_Fs(int s, int s_slot, int i, int j, RNAEnergyEvaluator& evaluator, int theta) {
     output_backtrack out;
     if (i >= j) return out;
 
     const auto& F = evaluator.get_F(s);
     const auto& C = evaluator.get_C(s);
     const std::string& seq = evaluator.get_strands().at(s);
 
     int best = F[i][j];
 
     // Case 1: i unpaired
     if (i + 1 <= j && best == F[i + 1][j]) {
         return backtrack_Fs(s, s_slot, i + 1, j, evaluator, theta);
     }
 
     // Case 2: i paired with k
     const short* S1_Fs = evaluator.get_S1_map().at(s);
     vrna_md_t md_Fs = evaluator.get_params()->model_details;
     for (int k = i + theta + 1; k <= j; ++k) {
         if (!can_pair(seq[i], seq[k])) continue;
         int type_ik = vrna_get_ptype_md(S1_Fs[i], S1_Fs[k], &md_Fs);
         int ext_pen = (type_ik > 2) ? evaluator.get_params()->TerminalAU : 0;
         int cand = C[i][k] + ext_pen;
         if (k + 1 <= j) cand += F[k + 1][j];
 
         if (best == cand) {
             output_backtrack left  = backtrack_Cs(s, s_slot, i, k, evaluator, theta);
             output_backtrack right = backtrack_Fs(s, s_slot, k + 1, j, evaluator, theta);
             left.merge(right);
             return left;
         }
     }
     std::cerr << "BT FAIL in backtrack_X: "
     << " s="<<s<<" i="<<i<<" j="<<j
     << " best="<<best << "\n";
     return out;
 }  
  
  /*
  -------------------------------------------------------------
   *  C-matrix (closed pair region)
 *------------------------------------------------------------
   */
  output_backtrack backtrack_Cs(int s, int s_slot, int i, int j, RNAEnergyEvaluator& evaluator, int theta) {
      output_backtrack out;
      const auto& C = evaluator.get_C(s);
      const auto& M = evaluator.get_M(s);
      const auto& M1 = evaluator.get_M1(s);
      const std::string& seq = evaluator.get_strands().at(s);
      const vrna_param_t* params = evaluator.get_params();
  
      int best = C[i][j];
      // Store actual strand ID s (not a slot number).
      // full_backtrack will convert IDs -> slots after backtracking completes.
      out.add_pair(s_slot, i, s_slot, j);
  
      // --- Hairpin ---
      if (best == evaluator.hairpin_energy(i, j, s)) {
          return out;
      }
  
      // --- Internal loops ---
      for (int p = i + 1; p < j - theta - 1; ++p) {
          for (int q = p + theta + 1; q < j; ++q) {
              if (!can_pair(seq[p], seq[q ])) continue;
              int e_int = evaluator.interior_loop_energy(i, j, p, q, s);
              if (best == C[p][q] + e_int) {
                  output_backtrack inner = backtrack_Cs(s, s_slot, p, q, evaluator, theta);
                  out.merge(inner);
                  return out;
              }
          }
      }
  
      // --- Multiloop ---
      {
          const short* S1_bt = evaluator.get_S1_map().at(s);
          vrna_md_t md_copy = params->model_details;
          int type_ij = vrna_get_ptype_md(S1_bt[i], S1_bt[j], &md_copy);
          int mlintern_ij = (type_ij != 0) ? E_MLstem_nd(type_ij, params) : 0;
          for (int u = i + 1; u < j; ++u) {
              int left  = M1[i + 1][u];
              int right = M[u + 1][j - 1];

              if (left >= inf_energy || right >= inf_energy) continue;

              int cand = left + right + params->MLclosing + mlintern_ij;

              if (best == cand) {
                  output_backtrack L = backtrack_M1s(s, s_slot, i + 1, u, evaluator, theta);
                  output_backtrack R = backtrack_Ms(s, s_slot, u + 1, j - 1, evaluator, theta);
                  out.merge(L);
                  out.merge(R);
                  return out;
              }
          }
      } // end multiloop block

     std::cerr << "BT FAIL in backtrack_X: "
           << " s="<<s<<" i="<<i<<" j="<<j
           << " best="<<best << "\n";
 
  
      return out; 
  }
  
  /*
  -------------------------------------------------------------
   *  M-matrix (multiloop segment)
   *------------------------------------------------------------
   */
  output_backtrack backtrack_Ms(int s, int s_slot, int i, int j,
     RNAEnergyEvaluator& evaluator, int theta)
 {
 output_backtrack out;
 if (i >= j) return out;  // i==j: single trailing unpaired base, M[i][i]=0, nothing to record
 
 const auto& M  = evaluator.get_M(s);
 const auto& C  = evaluator.get_C(s);
 const auto& S1 = evaluator.get_S1_map().at(s);
 vrna_param_t* params = evaluator.get_params();
 
 int best = M[i][j];
 
 
 // Case 1: i unpaired
 if (i + 1 <= j && best == M[i + 1][j] + params->MLbase) {
 return backtrack_Ms(s, s_slot, i + 1, j, evaluator, theta);
 }
 
 // Case 2: split (C | M)
 for (int u = i; u < j; ++u) {
 if (C[i][u] >= inf_energy || M[u + 1][j] >= inf_energy) continue;
 if (!can_pair(evaluator.get_strands().at(s)[i], evaluator.get_strands().at(s)[u])) continue;
 
 int type = vrna_get_ptype_md(S1[i], S1[u], &params->model_details);
 if (type == 0) continue;
 
 int penalty = E_MLstem_nd(type, params);
 
 if (best == C[i][u] + M[u + 1][j] + penalty) {
 output_backtrack left  = backtrack_Cs(s, s_slot, i, u, evaluator, theta);
 output_backtrack right = backtrack_Ms(s, s_slot, u + 1, j, evaluator, theta);
 left.merge(right);
 return left;
 }
 }
 
 // Case 3: start directly from C[i][u]
 for (int u = i; u <= j; ++u) {
 if (C[i][u] >= inf_energy) continue;
 if (!can_pair(evaluator.get_strands().at(s)[i], evaluator.get_strands().at(s)[u])) continue;
 
 int type = vrna_get_ptype_md(S1[i], S1[u], &params->model_details);
 if (type == 0) continue;
 
 int penalty = E_MLstem_nd(type, params);
 
 if (best == C[i][u] + penalty) {
 return backtrack_Cs(s, s_slot, i, u, evaluator, theta);
 }
 }
 std::cerr << "BT FAIL in backtrack_X: "
           << " s="<<s<<" i="<<i<<" j="<<j
           << " best="<<best << "\n";
 return out;
 }
 
  
  /*
  -------------------------------------------------------------
   *  M1-matrix (multiloop entry segment)
   *------------------------------------------------------------
   */
  output_backtrack backtrack_M1s(int s, int s_slot, int i, int j, RNAEnergyEvaluator& evaluator, int theta) {
      output_backtrack out;
      const auto& M1 = evaluator.get_M1(s);
      const auto& C = evaluator.get_C(s);
      const auto& S1 = evaluator.get_S1_map().at(s);
      vrna_param_t* params = evaluator.get_params();
  
      int best = M1[i][j];
      if (i > j) return out;
         
 
  
      // Case 1: i unpaired
      if (i + 1 <= j && best == M1[i+1][j] + params->MLbase) {
          return backtrack_M1s(s, s_slot, i+1, j , evaluator, theta);
      }
  
      // Case 2: helix at (i,j)
      const std::string& seq = evaluator.get_strands().at(s);
 if (can_pair(seq[i], seq[j])) {
     int type = vrna_get_ptype_md(S1[i], S1[j], &params->model_details);
     if (type != 0) {
         int penalty = E_MLstem_nd(type, params);
      if (best == C[i][j] + penalty) {
          output_backtrack inner = backtrack_Cs(s, s_slot, i, j, evaluator, theta);
          out.merge(inner);
          return out;}
         }
      }
      std::cerr << "BT FAIL in backtrack_X: "
           << " s="<<s<<" i="<<i<<" j="<<j
           << " best="<<best << "\n";
      return out;}
 
 /**
  * ============================================================
  
  * ============================================================
  */
 output_backtrack backtrack_F_multi_square(
     int m, int s, int s_slot, int i, int r, int j, int c,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
     RNAEnergyEvaluator& evaluator, vrna_param_t* params, int theta
 );
 
 output_backtrack backtrack_F_multi_bubble(
     int m, int s, int s_slot, int i, int r, int j, int c,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
     RNAEnergyEvaluator& evaluator, vrna_param_t* params, int theta
 );
 
 output_backtrack backtrack_C_multi(
     int m, int s, int s_slot, int i, int r, int j,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
     RNAEnergyEvaluator& evaluator,
     vrna_param_t* params,
     int theta);
 
 output_backtrack backtrack_M_multi(
     int m, int s, int s_slot, int i, int r, int j,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
     RNAEnergyEvaluator& evaluator,
     vrna_param_t* params,
     int theta
 
 );
 
 output_backtrack backtrack_M1_multi(
     int m, int s, int s_slot, int i, int r, int j,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
     RNAEnergyEvaluator& evaluator,
     vrna_param_t* params,
     int theta);
 
 
 // This corresponds to the GeneralCaseMinimization function (bubble case).
 output_backtrack backtrack_F_multi_bubble(
     int m, int s, int s_slot, int i, int r, int j, int c,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
     RNAEnergyEvaluator& evaluator, vrna_param_t* params, int theta
 ) {
     output_backtrack out;
     const int r_slot = s_slot + m + 1;
 
     const int len_s = int(strands.at(s).size()) - 1;
     const int len_r = int(strands.at(r).size()) - 1;
 
     float best = F(m, s, i, r, j, c);
     if (best >= inf_energy) return out;

     const auto& S1s_bt = evaluator.get_S1_map();
     vrna_md_t md_bt = params->model_details;

     // --- Case 1 : i unpaired ---
     if (F(m, s, i + 1, r, j, c) == best) {
         return backtrack_F_multi_square(m, s, s_slot, i + 1, r, j, c,
                                         strands, C, M, M1_multi, F,
                                         evaluator, params, theta);
     }

     // --- Case 2 : i pairs within s (single strand C_s) ---
     {
         const std::string& seq_s = strands.at(s);
         for (int k = i + theta + 1; k <= len_s; ++k) {
             if (!can_pair(seq_s[i], seq_s[k])) continue;
             int type_ik = vrna_get_ptype_md(S1s_bt.at(s)[i], S1s_bt.at(s)[k], &md_bt);
             float ext_pen = (type_ik > 2) ? (float)params->TerminalAU : 0.0f;
             float cand = evaluator.get_C(s)[i][k] + ext_pen + F(m, s, k + 1, r, j, c);
             if (cand == best) {
                 output_backtrack left = backtrack_Cs(s, s_slot, i, k, evaluator, theta);
                 output_backtrack right =
                     backtrack_F_multi_square(m, s, s_slot, k + 1, r, j, c,
                                              strands, C, M, M1_multi, F,
                                              evaluator, params, theta);
                 left.merge(right);
                 return left;
             }
         }
     }

     // --- Case 3 : i pairs with new strand t at position k ---
     if (m >= 1) {
         const std::string& seq_s = strands.at(s);

         for (int t = 1; t <= F.get_s_size(); ++t) {
             const std::string& seq_t = strands.at(t);
             int len_t = int(seq_t.size()) - 1;

             for (int m1 = 0; m1 < m; ++m1) {
                 int m2 = m - m1 - 1;

                 for (int k = 1; k <= len_t; ++k) {
                     if (!can_pair(seq_s[i], seq_t[k])) continue;
                     int type_ik = vrna_get_ptype_md(S1s_bt.at(s)[i], S1s_bt.at(t)[k], &md_bt);
                     float ext_pen = (type_ik > 2) ? (float)params->TerminalAU : 0.0f;
                     float cand = C(m1, s, i, t, k) + ext_pen + F(m2, t, k + 1, r, j, c);
                     if (cand == best) {
                         const int t_slot = s_slot + m1 + 1;
                         output_backtrack left =
                             backtrack_C_multi(m1, s, s_slot, i, t, k,
                                               strands, C, M, M1_multi, F,
                                               evaluator, params, theta);
                         left.add_sequence(t);
                         output_backtrack right =
                             backtrack_F_multi_square(m2, t, t_slot, k + 1, r, j, c,
                                                      strands, C, M, M1_multi, F,
                                                      evaluator, params, theta);
                         left.merge(right);
                         return left;
                     }
                 }
             }
         }
     }

     // --- Case 4 : i pairs with r at k ---
     if (c == 1) {
         const std::string& seq_s = strands.at(s);
         const std::string& seq_r = strands.at(r);

         for (int k = 1; k <= j; ++k) {
             if (!can_pair(seq_s[i], seq_r[k])) continue;
             int type_ik = vrna_get_ptype_md(S1s_bt.at(s)[i], S1s_bt.at(r)[k], &md_bt);
             float ext_pen = (type_ik > 2) ? (float)params->TerminalAU : 0.0f;
             if (k == len_r) {
                 float cand = C(m, s, i, r, k) + ext_pen;
                 if (cand == best) {
                     return backtrack_C_multi(m, s, s_slot, i, r, k,
                                              strands, C, M, M1_multi, F,
                                              evaluator, params, theta);
                 }
             } else {
                 float cand = C(m, s, i, r, k) + ext_pen + evaluator.get_F(r)[k + 1][j];
                 if (cand == best) {
                     output_backtrack left =
                         backtrack_C_multi(m, s, s_slot, i, r, k,
                                           strands, C, M, M1_multi, F,
                                           evaluator, params, theta);
                     output_backtrack right =
                         backtrack_Fs(r, r_slot, k + 1, j, evaluator, theta);
                     left.merge(right);
                     return left;
                 }
             }
         }
     }
 
     std::cerr << "BT FAIL in backtrack_X: "
     << "m="<<m<<" s="<<s<<" i="<<i<<" r="<<r<<" j="<<j
     << " best="<<best << "\n";
     return output_backtrack();
 }
 // This corresponds to the MainAuxiliaryMatrix function (square case).
 
 output_backtrack backtrack_F_multi_square(
     int m, int s, int s_slot, int i, int r, int j, int c,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
     RNAEnergyEvaluator& evaluator, vrna_param_t* params, int theta
 ) {
     output_backtrack out;
     const int r_slot = s_slot + m + 1;
 
     const int len_s = int(strands.at(s).size()) - 1; // '$' at 0
 
     float best = F(m, s, i, r, j, c);
     if (best >= inf_energy) return out;
 
     // =========================
     // CASE: s empty (i > |s|)
     // =========================
     if (i > len_s) {
         if (c == 1) {
             return output_backtrack(); 
         } else {
             if (m == 0) {
                 // F = single-strand folding on r from 1..j
                 if (j <= 0) return output_backtrack();
                 return backtrack_Fs(r, r_slot, 1, j, evaluator, theta);
             } else {
                 // F(m,s,|s|+1,r,j,0) = min_t F(m-1,t,1,r,j,1)
                 for (int t = 1; t <= F.get_s_size(); ++t) {
                    // Any strand type allowed as intermediate (soupfold).
                     if (F(m, s, i, r, j, c) == F(m - 1, t, 1, r, j, 1)) {
                         output_backtrack sub =
                             backtrack_F_multi_square(m - 1, t, s_slot + 1, 1, r, j, 1,
                                                      strands, C, M, M1_multi, F,
                                                      evaluator, params, theta);
                         sub.add_sequence_front(t);  // t is the new first strand; prepend
                         return sub;
                     }
                 }
             }
         }
 
         std::cerr << "BT FAIL in backtrack_X: "
         << "m="<<m<<" s="<<s<<" i="<<i<<" r="<<r<<" j="<<j
         << " best="<<best << "\n";
         return output_backtrack();
     }
 
     // =========================
     // CASE: r empty (j < 1)
     // =========================
     if (j < 1) {
         if (c == 1) {
             return output_backtrack();
         } else {
             if (m == 0) {
                 // F = single-strand folding on s from i..|s|
                 if (i <= 0) return output_backtrack();
                 return backtrack_Fs(s, s_slot, i, len_s, evaluator, theta);
             } else {
                 // F(m,s,i,r,0,0) = min_t F(m-1,s,i,t,|t|,1)
                 for (int t = 1; t <= F.get_s_size(); ++t) {
                    // Any strand type allowed as intermediate (soupfold).
                     int len_t = int(strands.at(t).size()) - 1;
                     if (F(m, s, i, r, j, c) == F(m - 1, s, i, t, len_t, 1)) {
                         output_backtrack sub =
                             backtrack_F_multi_square(m - 1, s, s_slot, i, t, len_t, 1,
                                                      strands, C, M, M1_multi, F,
                                                      evaluator, params, theta);
                         sub.add_sequence(t);
                         return sub;
                     }
                 }
             }
         }
 
         std::cerr << "BT FAIL in backtrack_X: "
         << "m="<<m<<" s="<<s<<" i="<<i<<" r="<<r<<" j="<<j
         << " best="<<best << "\n";
         return output_backtrack();
     }
 
     // =========================
     // Otherwise general case minimization
     // =========================
     return backtrack_F_multi_bubble(m, s, s_slot, i, r, j, c,
                                     strands, C, M, M1_multi, F,
                                     evaluator, params, theta);
 }
 
 
 
 output_backtrack backtrack_C_multi(
     int m, int s, int s_slot, int i, int r, int j,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
     RNAEnergyEvaluator& evaluator,
     vrna_param_t* params,
     int theta)
  {
     if (!can_pair(strands.at(s)[i], strands.at(r)[j])) {
         std::cerr << "Illegal outer pair in backtrack_C_multi\n";
         return output_backtrack();
     }
     output_backtrack out;
     const int r_slot = s_slot + m + 1;
 
     float best = C(m, s, i, r, j);
     if (best >= inf_energy) return out;
 
     const std::string& seq_s = strands.at(s);
     const std::string& seq_r = strands.at(r);
     const auto& S1s = evaluator.get_S1_map();
 
     int len_s = int(seq_s.length()) - 1;
     int len_r = int(seq_r.length()) - 1;

     // MLintern for closing pair (i,j): matches ViennaRNA's E_MLstem for multiloop closing.
     int type_ij_bt = vrna_get_ptype_md(S1s.at(s)[i], S1s.at(r)[j], &params->model_details);
     float mlintern_closing_bt = (type_ij_bt != 0) ? E_MLstem_nd(type_ij_bt, params) : 0.0f;

     out.add_pair(s_slot, i, r_slot, j);

     // --- Case 1 : interior loop, closed by (k,l) ---
     for (int k = i + 1; k <= len_s; ++k) {
         for (int l = 1; l <= j - 1; ++l) {
             if (k <= i || l >= j) continue;
             if (!can_pair(seq_s[k], seq_r[l])) continue;
             float loop_energy = evaluator.interior_loop_energy_cross(s, i, r, j, k, l);
             float cand = C(m, s, k, r, l) + loop_energy;
             if (cand == best) {
                 output_backtrack inner =
                     backtrack_C_multi(m, s, s_slot, k, r, l,
                                       strands, C, M, M1_multi, F, evaluator, params, theta);
                 out.merge(inner);
                 return out;
             }
         }
     }
 
     // --- Case 2 : leftmost stem within s ---
     for (int k = i + theta + 1; k <= len_s; ++k) {
         float ml_energy = params->MLclosing + mlintern_closing_bt + evaluator.get_M1(s)[i + 1][k];
         float cand = ml_energy + M(m, s, k + 1, r, j - 1);
         if (cand == best) {
             // left: M1 on single strand s
             output_backtrack left = backtrack_M1s(s, s_slot, i + 1, k, evaluator, theta);
             // right: multi-strand M
             output_backtrack right =
                 backtrack_M_multi(m, s, s_slot, k + 1, r, j - 1,
                                   strands, C, M, M1_multi, F, evaluator, params, theta);
             out.merge(left);
             out.merge(right);
             return out;
         }
     }
 
     // --- Case 3 : leftmost stem with new strand t ---
     if (m >= 1) {
         for (int t = 1; t <= F.get_s_size(); ++t) {
            // Any strand type allowed as intermediate (soupfold).
             const std::string& seq_t = strands.at(t);
             int len_t = int(seq_t.length()) - 1;
             for (int m1 = 0; m1 < m; ++m1) {
                 int m2 = m - m1 - 1;
                 for (int k = 1; k <= len_t; ++k) {
 
 
                     float cand =
                         params->MLclosing
                       + mlintern_closing_bt
                       + M1_multi(m1, s, i + 1, t, k)
                       + M(m2, t, k + 1, r, j - 1);
 
                       if (cand == best) {
                        const int t_slot = s_slot + m1 + 1;
                         output_backtrack left =
                             backtrack_M1_multi(m1, s, s_slot, i + 1, t, k,
                                                strands, C, M, M1_multi, F,
                                                evaluator, params, theta);
                     
                         left.add_sequence(t);
                     
                         output_backtrack right =
                             backtrack_M_multi(m2, t, t_slot, k + 1, r, j - 1,
                                               strands, C, M, M1_multi, F,
                                               evaluator, params, theta);
                     
                     
                         out.merge(left);
                         out.merge(right);
                         return out;
                     }
                                    
                 }
             }
         }
     }
 
     // --- Case 4 : leftmost stem with rightmost strand r at k ---
     for (int k = 1; k <= j - 1; ++k) {
 
         if (k == len_r) {
             float cand =
                 params->MLclosing
               + mlintern_closing_bt
               + M1_multi(m, s, i + 1, r, k);
             if (cand == best) {
                 output_backtrack left =
                     backtrack_M1_multi(m, s, s_slot, i + 1, r, k,
                                        strands, C, M, M1_multi, F, evaluator, params, theta);
                 out.merge(left);
                 return out;
             }
         } else {
             float cand =
                 params->MLclosing
               + mlintern_closing_bt
               + M1_multi(m, s, i + 1, r, k)
              
               + evaluator.get_M(r)[k + 1][j - 1];
 
             if (cand == best) {
                 output_backtrack left =
                     backtrack_M1_multi(m, s, s_slot, i + 1, r, k,
                                        strands, C, M, M1_multi, F, evaluator, params, theta);
                 output_backtrack right =
                     backtrack_Ms(r, r_slot, k + 1, j - 1, evaluator, theta);
                 out.merge(left);
                 out.merge(right);
                 return out;
             }
         }
     }
 
     // --- Case 5: external disconnect ---
     {
         float term_au_c5 = (type_ij_bt > 2) ? (float)params->TerminalAU : 0.0f;
         if (term_au_c5 + F(m, s, i + 1, r, j - 1, 0) == best) {
             output_backtrack inner =
                 backtrack_F_multi_square(m, s, s_slot, i + 1, r, j - 1, 0,
                                         strands, C, M, M1_multi, F, evaluator, params, theta);
             out.merge(inner);
             return out;
         }
     }

     std::cerr << "BT FAIL in backtrack_X: "
     << "m="<<m<<" s="<<s<<" i="<<i<<" r="<<r<<" j="<<j
     << " best="<<best << "\n";
     return out;
 }
 
 output_backtrack backtrack_M_multi(
     int m, int s, int s_slot, int i, int r, int j,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
     RNAEnergyEvaluator& evaluator,
     vrna_param_t* params,
     int theta) {
     output_backtrack out;
     const int r_slot = s_slot + m + 1;
 
     float best = M(m, s, i, r, j);
     if (best >= inf_energy) return out;
 
     const std::string& seq_s = strands.at(s);
     const std::string& seq_r = strands.at(r);
     const auto& S1s = evaluator.get_S1_map();
 
     int len_s = int(seq_s.length()) - 1;
     int len_r = int(seq_r.length()) - 1;
     // --- Case 1 : i unpaired ---
     if (M(m, s, i + 1, r, j) + params->MLbase == best) {
         return backtrack_M_multi(m, s, s_slot, i + 1, r, j,
                                  strands, C, M, M1_multi, F, evaluator, params, theta);
     }
 
     // --- Case 2 : i pairs within s ---
     for (int k = i + theta + 1; k <= len_s; ++k) {
         if (!can_pair(seq_s[i], seq_s[k])) continue;
 
         int type = vrna_get_ptype_md(
             S1s.at(s)[i], S1s.at(s)[k], &params->model_details);
             if (type == 0) continue;
         float cand = evaluator.get_C(s)[i][k] + E_MLstem_nd(type, params)
                    + M(m, s, k + 1, r, j);
         if (cand == best) {
             output_backtrack left = backtrack_Cs(s, s_slot, i, k, evaluator, theta);
             output_backtrack right =
                 backtrack_M_multi(m, s, s_slot, k + 1, r, j,
                                   strands, C, M, M1_multi, F, evaluator, params, theta);
             left.merge(right);
             return left;
         }
     }
 
     // --- Case 3 : i pairs with new strand t ---
     if (m >= 1) {
         for (int t = 1; t <= F.get_s_size(); ++t) {
            // Any strand type allowed as intermediate (soupfold).
             const std::string& seq_t = strands.at(t);
             int len_t = int(seq_t.length()) - 1;
             for (int m1 = 0; m1 < m; ++m1) {
                 int m2 = m - m1 - 1;
                 for (int k = 1; k <= len_t; ++k) {
                     if (!can_pair(seq_s[i], seq_t[k])) continue;
 
                     int type = vrna_get_ptype_md(
                         S1s.at(s)[i], S1s.at(t)[k], &params->model_details);
                         if (type == 0) continue;
 
                     float cand = C(m1, s, i, t, k)
                                + E_MLstem_nd(type, params)
                                + M(m2, t, k + 1, r, j);
                                if (cand == best) {
                                 const int t_slot = s_slot + m1 + 1;
                                 output_backtrack left =
                                     backtrack_C_multi(m1, s, s_slot, i, t, k,
                                                       strands, C, M, M1_multi, F,
                                                       evaluator, params, theta);
                             
                                 left.add_sequence(t);
                             
                                 output_backtrack right =
                                     backtrack_M_multi(m2, t, t_slot, k + 1, r, j,
                                                       strands, C, M, M1_multi, F,
                                                       evaluator, params, theta);
                             
                             
                                 left.merge(right);
                                 return left;
                             }
                             
                 }
             }
         }
     }
 
     // --- Case 4 : i pairs with r at k ---
     for (int k = 1; k <= j - 1; ++k) {
         if (!can_pair(seq_s[i], seq_r[k])) continue;
 
         int type = vrna_get_ptype_md(
             S1s.at(s)[i], S1s.at(r)[k], &params->model_details);
 
         if (k == len_r) {
             float cand = C(m, s, i, r, k) + E_MLstem_nd(type, params);
             if (cand == best) {
                 output_backtrack left =
                     backtrack_C_multi(m, s, s_slot, i, r, k,
                                       strands, C, M, M1_multi, F, evaluator, params, theta);
                 return left;
             }
         } else {
             float cand =
             C(m, s, i, r, k) + E_MLstem_nd(type, params)
               + evaluator.get_M(r)[k + 1][j];
 
             if (cand == best) {
                 output_backtrack left =
                     backtrack_C_multi(m, s, s_slot, i , r, k,
                                        strands, C, M, M1_multi, F, evaluator, params, theta);
                 output_backtrack right =
                     backtrack_Ms(r, r_slot, k + 1, j, evaluator, theta);
                 left.merge(right);
                 return left;
             }
         }
     }
 
     std::cerr << "BT FAIL in backtrack_X: "
           << "m="<<m<<" s="<<s<<" i="<<i<<" r="<<r<<" j="<<j
           << " best="<<best << "\n";
 
     return out;
 }
 
 output_backtrack backtrack_M1_multi(
     int m, int s, int s_slot, int i, int r, int j,
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
     RNAEnergyEvaluator& evaluator,
     vrna_param_t* params,
     int theta)
 {
     output_backtrack out;
 
     float best = M1_multi(m, s, i, r, j);
     if (best >= inf_energy) return out;
 
     const auto& S1s = evaluator.get_S1_map();
     int len_s = int(strands.at(s).length()) - 1;
 
     // --- Case 1: i unpaired ---
     if (i < len_s) {
         float cand = M1_multi(m, s, i + 1, r, j) + params->MLbase;
         if (cand == best) {
             return backtrack_M1_multi(m, s, s_slot, i + 1, r, j,
                                       strands, C, M, M1_multi, F,
                                       evaluator, params, theta);
         }
     }
 
     // --- Case 2: helix at (i,j) ---
     if (can_pair(strands.at(s)[i], strands.at(r)[j])) {
         int type = vrna_get_ptype_md(
             S1s.at(s)[i], S1s.at(r)[j], &params->model_details);
 
         if (type != 0) {
             float cand2 = C(m, s, i, r, j) + E_MLstem_nd(type, params);
 
             if (cand2 == best) {
                 output_backtrack inner =
                     backtrack_C_multi(m, s, s_slot, i, r, j,
                                       strands, C, M, M1_multi, F,
                                       evaluator, params, theta);
                 out.merge(inner);
                 return out;
             }
         }
     }
 
     std::cerr << "BT FAIL in backtrack_M1_multi: "
               << "m=" << m << " s=" << s << " i=" << i
               << " r=" << r << " j=" << j
               << " best=" << best << "\n";
     return out;
 }
 
 std::vector<output_backtrack> full_backtrack(
     const std::unordered_map<int, std::string>& strands,
     Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
     RNAEnergyEvaluator& evaluator,
     vrna_param_t* params,
     int theta
 ) {
     std::vector<output_backtrack> secondary_structures;
 
     const int n_strands = (int)strands.size();
 
     std::vector<std::vector<int>> starting_points;
     try {
         starting_points = Find_all_start_backtrack(F, n_strands, 1); // c_target=1 for a connected structure
     } catch (const std::runtime_error& e) {
         std::cerr << "Warning: " << e.what() << "\n";
         return secondary_structures;
     }
 
     for (auto& sp : starting_points) {
         int m = sp[0], s = sp[1], i = sp[2], r = sp[3], j = sp[4], c = sp[5];
 
         std::cout << "Starting point: F("
                   << m << "," << s << "," << i << ","
                   << r << "," << j << "," << c << ") = "
                   << F(m, s, i, r, j, c) << "\n";
 
         output_backtrack secondary =
             backtrack_F_multi_square(m, s, 1, i, r, j, c,
                                      strands, C, M, M1_multi, F,
                                      evaluator, params, theta);
 
          // -----------------------------------------------------------------
         // Build the definitive ordered strand list: slot 1 = s, slot m+2 = r,
         // slots 2..m+1 = intermediate strands in the order they were encountered.
         //
         // Pairs already store SLOT NUMBERS (set by backtracking functions).
         // bt_sequences contains intermediate TYPE IDs in order (via add_sequence calls).
         std::vector<int> bt_sequences = secondary.list_of_sequences;

         // Build ordered list: s, intermediates..., r
         secondary.list_of_sequences.clear();
         secondary.list_of_sequences.push_back(s);
         for (int x : bt_sequences)
             secondary.list_of_sequences.push_back(x);
         secondary.list_of_sequences.push_back(r);

         if ((int)secondary.list_of_sequences.size() != m + 2) {
             std::cerr << "WARNING: expected " << m+2 << " strands in list, got "
                       << secondary.list_of_sequences.size()
                       << " (m=" << m << " s=" << s << " r=" << r << ")\n";
         }
         // Pairs already have slot numbers — no id_to_slot conversion needed.
       
         float E_raw = F(m, s, i, r, j, c);
 
         const int m_used = m + 2;   
 
         float Kassoc = params->DuplexInit;
         float E_assoc = (m_used > 1) ? (m_used - 1) * Kassoc : 0.0f;
         float E_total = E_raw + E_assoc;
 
         std::cout << "m_used=" << m_used
                   << "  E_raw="   << (E_raw / 100.0f)   << " kcal/mol"
                   << "  E_assoc=" << (E_assoc / 100.0f) << " kcal/mol"
                   << "  E_total=" << (E_total / 100.0f) << " kcal/mol\n";
 
         secondary.print();
         secondary.print_dot_bracket(strands);
         secondary_structures.push_back(secondary);
     }
 
     return secondary_structures;
 }
 
 //======================================    Application on Triplet Repeats    ==============================================//
 
 
 
 
 
 /**
  * @brief Computes the repartition of pairs (internal/homogeneous/heterogeneous) in the secondary structure for all triplet repeats.
  * 
  * @param number_of_repeats 
  * @param m_start 
  * 
  * STEP 1 : This function generates all possible triplet combinations.
  * STEP 2,3 : Then for each combination, it generates the corresponding RNA strands and computes the single strand folding.
  * STEP 4 : It then fills the 6D matrix of the strand_soup problem
  * STEP 5 : Finds all the strating points in the energy matrix and backtracks to find the secondary structures.
  * STEP 6 : Computes the repartition of pairs (internal/homogeneous/heterogeneous) in the secondary structures.
  * STEP 7 : Writes the repartition to CSV files.
  */
 
 
 /**
  * @brief Computes the repartition of pairs (internal/homogeneous/heterogeneous)
  *        in the secondary structure for all triplet repeats, using the
  *        thermodynamic Vienna-based strand soup (C/M/M1/F) algorithm.
  *
  * @param number_of_repeats Number of triplet repetitions in each strand.
  * @param m_start           Total number of strands in the soup.
  *
  * STEP 1 : Generate all 64 RNA triplets.
  * STEP 2 : For each pair (triplet1, triplet2), generate two RNA strands.
  * STEP 3 : Build RNAEnergyEvaluator (single-strand F/C/M/M1) on these strands.
  * STEP 4 : Build multi-strand matrices C, M, M1_multi, F and fill them with MainAuxiliaryMatrix.
  * STEP 5 : Backtrack via full_backtrack to get secondary_structures.
  * STEP 6 : Compute repartition of pairs (internal / homogeneous / heterogeneous).
  * STEP 7 : Write the repartition to CSV files.
  */
 
  inline std::string with_dollar(const std::string& s) { return "$" + s; }
 
  void compute_all_triplet_repartition(int number_of_repeats, int m_start) {
      // STEP 1: Generate all RNA Triplets
      std::vector<std::string> triplets;
      std::string bases = "AUCG";
      for (char b1 : bases) {
          for (char b2 : bases) {
              for (char b3 : bases) {
                  triplets.push_back(std::string(1, b1) + b2 + b3);
              }
          }
      }
  
      std::string internal_filename =
          "../../results/internal_n=" + std::to_string(number_of_repeats * 3) +
          "_m=" + std::to_string(m_start) + ".csv";
  
      std::string homogeneous_filename =
          "../../results/homogeneous_n=" + std::to_string(number_of_repeats * 3) +
          "_m=" + std::to_string(m_start) + ".csv";
  
      std::string heterogeneous_filename =
          "../../results/heterogeneous_n=" + std::to_string(number_of_repeats * 3) +
          "_m=" + std::to_string(m_start) + ".csv";
  
      std::ofstream internal_file(internal_filename);
      std::ofstream homogeneous_file(homogeneous_filename);
      std::ofstream heterogeneous_file(heterogeneous_filename);
  
      if (!internal_file.is_open() || !homogeneous_file.is_open() || !heterogeneous_file.is_open()) {
          std::cerr << "Error: Unable to open one or more files for writing." << std::endl;
          return;
      }
  
      // CSV header row: first column empty, then the 64 triplets as columns
      internal_file << ",";
      homogeneous_file << ",";
      heterogeneous_file << ",";
      for (size_t i = 0; i < triplets.size(); ++i) {
          internal_file   << triplets[i];
          homogeneous_file << triplets[i];
          heterogeneous_file << triplets[i];
          if (i != triplets.size() - 1) {
              internal_file   << ",";
              homogeneous_file << ",";
              heterogeneous_file << ",";
          }
      }
      internal_file   << "\n";
      homogeneous_file << "\n";
      heterogeneous_file << "\n";
  
      int total_combinations = int(triplets.size() * triplets.size());
      int current_combination = 0;
      auto start_time = std::chrono::high_resolution_clock::now();
  
      // STEP 2: Iterate over each triplet pair
      for (const auto& triplet1 : triplets) {
          std::string strand1 = generate_triplet_repeat(triplet1, number_of_repeats); // already has '$'
  
          internal_file   << triplet1 << ",";
          homogeneous_file << triplet1 << ",";
          heterogeneous_file << triplet1 << ",";
  
          for (const auto& triplet2 : triplets) {
              std::cout << "\n========== Triplet : " << triplet1
                        << " and " << triplet2 << " ==========" << std::endl;
  
              std::string strand2 = generate_triplet_repeat(triplet2, number_of_repeats); // already has '$'
  
              // Create strands map (index -> sequence with $)
              std::unordered_map<int, std::string> strands;
              add_strand_if_unique(strands, strand1);
              add_strand_if_unique(strands, strand2);
              // NOTE: add_strand_if_unique is assumed not to add '$' itself
  
              // STEP 3: Build RNAEnergyEvaluator (computes single-strand F/C/M/M1)
              RNAEnergyEvaluator evaluator(strands);
              vrna_param_t* params = evaluator.get_params();
  
              // STEP 4: Compute multi-strand matrices
              int m_size = m_start - 1;                
              int s_size = int(strands.size());        // number of distinct strand types
              int i_size = int(strand1.length()) - 1;  // real #bases without '$'
              int r_size = int(strands.size());
              int j_size = int(strand2.length()) - 1;
              int c_size = 2; // connectivity 0 or 1
  
              Matrix5D C(m_size, s_size, i_size, r_size, j_size);
              Matrix5D M(m_size, s_size, i_size, r_size, j_size);
              Matrix5D M1_multi(m_size, s_size, i_size, r_size, j_size);
              Matrix6D F(m_size, s_size, i_size, r_size, j_size, c_size);
  
              MainAuxiliaryMatrix(strands, C, M, M1_multi, F, evaluator, params);
  
              // STEP 5: Backtrack secondary structure(s)
              std::vector<output_backtrack> secondary_structures =
                  full_backtrack(strands, C, M, M1_multi, F, evaluator, params, theta);
  
              // STEP 6: Compute repartition of pairs
              int internal = 0;
              int homogeneous = 0;
              int heterogeneous = 0;
  
              for (const auto& secondary_structure : secondary_structures) {
                  for (const auto& pair : secondary_structure.list_of_pairs) {
                      int strand_idx1 = pair[0];
                      int strand_idx2 = pair[2];
  
                      if (strand_idx1 == strand_idx2) {
                          internal++;
                      } else {
                          int seq_id1 = secondary_structure.list_of_sequences[strand_idx1 - 1];
                          int seq_id2 = secondary_structure.list_of_sequences[strand_idx2 - 1];
  
                          if (seq_id1 == seq_id2) {
                              homogeneous++;
                          } else {
                              heterogeneous++;
                          }
                      }
                  }
              }
  
              int total = internal + homogeneous + heterogeneous;
              if (total > 0) {
                  internal_file   << std::fixed << std::setprecision(3)
                                  << static_cast<float>(internal) / total;
                  homogeneous_file << std::fixed << std::setprecision(3)
                                  << static_cast<float>(homogeneous) / total;
                  heterogeneous_file << std::fixed << std::setprecision(3)
                                    << static_cast<float>(heterogeneous) / total;
  
                  std::cout << "Internal : "     << static_cast<float>(internal) / total << std::endl;
                  std::cout << "Homogeneous : "  << static_cast<float>(homogeneous) / total << std::endl;
                  std::cout << "Heterogeneous : "<< static_cast<float>(heterogeneous) / total << std::endl;
              } else {
                  internal_file   << "0";
                  homogeneous_file << "0";
                  heterogeneous_file << "0";
              }
  
              if (triplet2 != triplets.back()) {
                  internal_file   << ",";
                  homogeneous_file << ",";
                  heterogeneous_file << ",";
              }
  
              current_combination++;
              float progress = (float(current_combination) / float(total_combinations)) * 100.0f;
  
              std::cout << "\rProgress: [";
              int bar_width = 50;
              int pos = int(bar_width * progress / 100.0f);
              for (int i = 0; i < bar_width; ++i) {
                  if (i < pos) std::cout << "=";
                  else if (i == pos) std::cout << ">";
                  else std::cout << " ";
              }
              std::cout << "] " << std::fixed << std::setprecision(2)
                        << progress << "%";
              std::cout.flush();
          }
  
          // End the current row in each file
          internal_file   << "\n";
          homogeneous_file << "\n";
          heterogeneous_file << "\n";
      }
  
      // End timing
      auto end_time = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed_time = end_time - start_time;
  
      std::cout << "\nComputation completed in " << elapsed_time.count() << " seconds." << std::endl;
      std::cout << "Results saved to " << internal_filename
                << ", " << homogeneous_filename
                << ", " << heterogeneous_filename << std::endl;
  }
 
  //=================================================================//
 
 void run_homogeneous_soup(const std::string& triplet, int repeats, int m_start) {
     const int seq_len = repeats * 3;
 
     std::cout << "\n============================\n";
     std::cout << "Homogeneous soup: (" << triplet << ")_" << repeats
               << "   m_start=" << m_start
               << "   length=" << seq_len << "\n";
     std::cout << "============================\n";
 
     // Build distinct strands (even if sequences identical)
     std::unordered_map<int, std::string> strands;
     std::string strand = generate_triplet_repeat(triplet, repeats);
    // Use one entry per strand type; m_start controls total strands in complex.
    strands[1] = strand;
 
 
 
     // Single-strand evaluator
     auto t0 = std::chrono::high_resolution_clock::now();
     RNAEnergyEvaluator evaluator(strands);
     vrna_param_t* params = evaluator.get_params();
     auto t1 = std::chrono::high_resolution_clock::now();
 
     // Multi-strand matrices
     int m_size = m_start - 1;
     int s_size = int(strands.size());
     int i_size = seq_len;
     int r_size = int(strands.size());
     int j_size = seq_len;
     int c_size = 2;
 
     Matrix5D C(m_size, s_size, i_size, r_size, j_size);
     Matrix5D M(m_size, s_size, i_size, r_size, j_size);
     Matrix5D M1_multi(m_size, s_size, i_size, r_size, j_size);
     Matrix6D F(m_size, s_size, i_size, r_size, j_size, c_size);
     auto t2 = std::chrono::high_resolution_clock::now();
     try {
         MainAuxiliaryMatrix(strands, C, M, M1_multi, F, evaluator, params);
     } catch (const std::out_of_range& e) {
         std::cerr << "\nOUT_OF_RANGE during MainAuxiliaryMatrix: " << e.what() << "\n";
         throw;
     }
     auto t3 = std::chrono::high_resolution_clock::now();
     
 
     // Backtrack
     auto t4 = std::chrono::high_resolution_clock::now();
     std::vector<output_backtrack> structures =
         full_backtrack(strands, C, M, M1_multi, F, evaluator, params, theta);
     auto t5 = std::chrono::high_resolution_clock::now();
 
     // Report timings
     std::chrono::duration<double> dt_eval = t1 - t0;
     std::chrono::duration<double> dt_alloc = t2 - t1;
     std::chrono::duration<double> dt_dp   = t3 - t2;
     std::chrono::duration<double> dt_bt   = t5 - t4;
 
     std::cout << "\nTimings:\n";
     std::cout << "  evaluator (single strand): " << dt_eval.count() << " s\n";
     std::cout << "  matrix alloc:              " << dt_alloc.count() << " s\n";
     std::cout << "  MainAuxiliaryMatrix DP:    " << dt_dp.count() << " s\n";
     std::cout << "  backtrack:                 " << dt_bt.count() << " s\n";
     std::cout << "  structures found:          " << structures.size() << "\n";
 
     // pair repartition for the returned structure(s)
     int internal = 0, homogeneous = 0, heterogeneous = 0;
 
     for (const auto& sec : structures) {
         for (const auto& pair : sec.list_of_pairs) {
             int strand_idx1 = pair[0];
             int strand_idx2 = pair[2];
 
             if (strand_idx1 == strand_idx2) {
                 internal++;
             } else {
                 int seq_id1 = sec.list_of_sequences[strand_idx1 - 1];
                 int seq_id2 = sec.list_of_sequences[strand_idx2 - 1];
 
                 if (seq_id1 == seq_id2) homogeneous++;
                 else heterogeneous++;
             }
         }
     }
 
     int total = internal + homogeneous + heterogeneous;
     if (total > 0) {
         std::cout << "\nPair repartition:\n";
         std::cout << "  internal      = " << (double)internal / total << "\n";
         std::cout << "  homogeneous   = " << (double)homogeneous / total << "\n";
         std::cout << "  heterogeneous = " << (double)heterogeneous / total << "\n";
     } else {
         std::cout << "\nNo pairs found.\n";
     }
 }
 
 void run_random_homogeneous_soup(int length, int m_start, uint64_t seed) {
     std::cout << "\n============================\n";
     std::cout << "Random homogeneous soup"
               << "   length=" << length
               << "   m_start=" << m_start
               << "   seed=" << seed << "\n";
     std::cout << "============================\n";
 
     seed_random(seed);
     std::string base = generate_random_sequence(length);
 
     std::unordered_map<int, std::string> strands;
    // Use one entry per strand type; m_start controls total strands in complex.
    strands[1] = "$" + base;  // add $ for 1-based indexing
 
         RNAEnergyEvaluator evaluator(strands);
         vrna_param_t* params = evaluator.get_params();
     
         int m_size = m_start - 1;
         int s_size = (int)strands.size();
         int i_size = length;
         int r_size = (int)strands.size();
         int j_size = length;
         int c_size = 2;
     
         Matrix5D C(m_size, s_size, i_size, r_size, j_size);
         Matrix5D M(m_size, s_size, i_size, r_size, j_size);
         Matrix5D M1_multi(m_size, s_size, i_size, r_size, j_size);
         Matrix6D F(m_size, s_size, i_size, r_size, j_size, c_size);
       
     
         MainAuxiliaryMatrix(strands, C, M, M1_multi, F, evaluator, params);
     
         auto structures = full_backtrack(strands, C, M, M1_multi, F, evaluator, params, theta);
         std::cout << "structures found: " << structures.size() << "\n";
     }
 
 //=======================================//
 
 void run_homogeneous(const std::string& motif, int motif_repeats, int n_strands) {
 
     std::cout << "\n============================\n";
     std::cout << "Homogeneous soup: (" << motif << ")_" << motif_repeats
               << "   n_strands=" << n_strands;
     std::cout << "\n============================\n";
 
    // Use one entry per strand type; n_strands controls total strands in complex.
     std::unordered_map<int, std::string> strands;
     std::string strand = generate_triplet_repeat(motif, motif_repeats);
 
    strands[1] = strand;
 
    
     // Single-strand evaluator
     auto t0 = std::chrono::high_resolution_clock::now();
     RNAEnergyEvaluator evaluator(strands);
     vrna_param_t* params = evaluator.get_params();
     auto t1 = std::chrono::high_resolution_clock::now();
 
     // --- Dimensions ---
     // m counts "extra strands besides the two endpoints (s and r)"
     // So for N strands total, we need m_max = N - 2.
     int m_max  = n_strands - 2;
     if (m_max < 0) {
         throw std::runtime_error("Need at least 2 strands for a multi-strand structure.");
     }
 
     int m_size = m_max + 1;              
    int s_size = int(strands.size());  // number of distinct strand types
    int r_size = int(strands.size());  // number of distinct strand types
 
     int i_size = int(strands.at(1).size()) - 1;   
     int j_size = int(strands.at(1).size()) - 1;
 
     int c_size = 2;
 
     Matrix5D C(m_size, s_size, i_size, r_size, j_size);
     Matrix5D M(m_size, s_size, i_size, r_size, j_size);
     Matrix5D M1_multi(m_size, s_size, i_size, r_size, j_size);
     Matrix6D F(m_size, s_size, i_size, r_size, j_size, c_size);
 
     auto t2 = std::chrono::high_resolution_clock::now();
     MainAuxiliaryMatrix(strands, C, M, M1_multi, F, evaluator, params);
     auto t3 = std::chrono::high_resolution_clock::now();
 
     // Backtrack
     auto t4 = std::chrono::high_resolution_clock::now();
     std::vector<output_backtrack> structures =
         full_backtrack(strands, C, M, M1_multi, F, evaluator, params, theta);
     auto t5 = std::chrono::high_resolution_clock::now();
 
     // Report timings
     std::chrono::duration<double> dt_eval  = t1 - t0;
     std::chrono::duration<double> dt_alloc = t2 - t1;
     std::chrono::duration<double> dt_dp    = t3 - t2;
     std::chrono::duration<double> dt_bt    = t5 - t4;
 
     std::cout << "\nTimings:\n";
     std::cout << "  evaluator (single strand): " << dt_eval.count()  << " s\n";
     std::cout << "  matrix alloc:              " << dt_alloc.count() << " s\n";
     std::cout << "  MainAuxiliaryMatrix DP:    " << dt_dp.count()    << " s\n";
     std::cout << "  backtrack:                 " << dt_bt.count()    << " s\n";
     std::cout << "  structures found:          " << structures.size() << "\n";
 }
 
 //=============================================================================================//
 
 void run_heterogeneous_triplets(const std::vector<std::string>& triplets,
     int repeats) {
 std::cout << "\n============================\n";
 std::cout << "Heterogeneous triplet soup   repeats=" << repeats << "\n";
 std::cout << "============================\n";
 
 // Build strands: one strand per triplet motif
 std::unordered_map<int, std::string> strands;
 for (int id = 1; id <= (int)triplets.size(); ++id) {
 strands[id] = generate_triplet_repeat(triplets[id - 1], repeats);
 }
 
 std::cout << "Generated strands:\n";
 for (const auto& pair : strands) {
 std::cout << "  Index " << pair.first
 << " : " << pair.second.substr(1) << "\n"; // remove '$' in print
 }
 
 // Single-strand evaluator
 auto t0 = std::chrono::high_resolution_clock::now();
 RNAEnergyEvaluator evaluator(strands);
 vrna_param_t* params = evaluator.get_params();
 auto t1 = std::chrono::high_resolution_clock::now();
 
 // Multi-strand matrices
 int n_strands = (int)strands.size();
 int m_size = n_strands - 1;
 int s_size = n_strands;
 int r_size = n_strands;
 int i_size = (int)strands.at(1).length() - 1; 
 int j_size = (int)strands.at(1).length() - 1; 
 int c_size = 2;
 
 Matrix5D C(m_size, s_size, i_size, r_size, j_size);
 Matrix5D M(m_size, s_size, i_size, r_size, j_size);
 Matrix5D M1_multi(m_size, s_size, i_size, r_size, j_size);
 Matrix6D F(m_size, s_size, i_size, r_size, j_size, c_size);
 auto t2 = std::chrono::high_resolution_clock::now();
 
 try {
 MainAuxiliaryMatrix(strands, C, M, M1_multi, F, evaluator, params);
 } catch (const std::out_of_range& e) {
 std::cerr << "\nOUT_OF_RANGE during MainAuxiliaryMatrix: "
 << e.what() << "\n";
 throw;
 }
 
 auto t3 = std::chrono::high_resolution_clock::now();
 
 // Backtrack
 auto t4 = std::chrono::high_resolution_clock::now();
 std::vector<output_backtrack> structures =
 full_backtrack(strands, C, M, M1_multi, F, evaluator, params, theta);
 auto t5 = std::chrono::high_resolution_clock::now();
 
 // Timings
 std::chrono::duration<double> dt_eval = t1 - t0;
 std::chrono::duration<double> dt_alloc = t2 - t1;
 std::chrono::duration<double> dt_dp   = t3 - t2;
 std::chrono::duration<double> dt_bt   = t5 - t4;
 
 std::cout << "\nTimings:\n";
 std::cout << "  evaluator (single strand): " << dt_eval.count() << " s\n";
 std::cout << "  matrix alloc:              " << dt_alloc.count() << " s\n";
 std::cout << "  MainAuxiliaryMatrix DP:    " << dt_dp.count() << " s\n";
 std::cout << "  backtrack:                 " << dt_bt.count() << " s\n";
 std::cout << "  structures found:          " << structures.size() << "\n";
 }
 //=============================================================================================//

// k distinct strand types, m_total strands in the complex (m_total can be > k)
// CLI: ./strand_soup.exe --heterogeneous <m_total> <seq1> <seq2> ...
void run_heterogeneous_general(const std::vector<std::string>& sequences, int m_total) {

    if (sequences.empty())
        throw std::runtime_error("run_heterogeneous_general: no sequences provided.");
    if (m_total < 2)
        throw std::runtime_error("run_heterogeneous_general: need at least 2 strands.");

    int k = (int)sequences.size();   // distinct types

    // assign IDs 1..k, prepend '$'
    std::unordered_map<int, std::string> strands;
    for (int id = 1; id <= k; ++id)
        strands[id] = "$" + sequences[id - 1];

    std::cout << "\n============================\n";
    std::cout << "Heterogeneous general soup\n";
    std::cout << "  distinct types (k) = " << k      << "\n";
    std::cout << "  total strands  (m) = " << m_total << "\n";
    std::cout << "  strands:\n";
    for (int id = 1; id <= k; ++id)
        std::cout << "    [" << id << "] " << sequences[id - 1] << "  (" << sequences[id-1].size() << " nt)\n";
    std::cout << "============================\n";

    // i/j dimensions need to cover the longest strand
    int max_len = 0;
    for (const auto& seq : sequences)
        max_len = std::max(max_len, (int)seq.size());

    // same dimension convention as run_homogeneous:
    //   m_size = m_total - 1  (middle strands)
    //   s_size = r_size = k   (loop over distinct types)
    //   i_size = j_size = max_len
    int m_size = m_total - 1;
    int s_size = k;
    int r_size = k;
    int i_size = max_len;
    int j_size = max_len;
    int c_size = 2;

    auto t0 = std::chrono::high_resolution_clock::now();
    RNAEnergyEvaluator evaluator(strands);
    vrna_param_t* params = evaluator.get_params();
    auto t1 = std::chrono::high_resolution_clock::now();

    Matrix5D C       (m_size, s_size, i_size, r_size, j_size);
    Matrix5D M       (m_size, s_size, i_size, r_size, j_size);
    Matrix5D M1_multi(m_size, s_size, i_size, r_size, j_size);
    Matrix6D F       (m_size, s_size, i_size, r_size, j_size, c_size);
    auto t2 = std::chrono::high_resolution_clock::now();

    MainAuxiliaryMatrix(strands, C, M, M1_multi, F, evaluator, params);
    auto t3 = std::chrono::high_resolution_clock::now();

    std::vector<output_backtrack> structures =
        full_backtrack(strands, C, M, M1_multi, F, evaluator, params, theta);
    auto t4 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> dt_eval  = t1 - t0;
    std::chrono::duration<double> dt_alloc = t2 - t1;
    std::chrono::duration<double> dt_dp    = t3 - t2;
    std::chrono::duration<double> dt_bt    = t4 - t3;

    std::cout << "\nTimings:\n";
    std::cout << "  evaluator (single strand): " << dt_eval.count()  << " s\n";
    std::cout << "  matrix alloc:              " << dt_alloc.count() << " s\n";
    std::cout << "  MainAuxiliaryMatrix DP:    " << dt_dp.count()    << " s\n";
    std::cout << "  backtrack:                 " << dt_bt.count()    << " s\n";
    std::cout << "  structures found:          " << structures.size() << "\n";
}

 //=============================================================================================//

 int main(int argc, char* argv[]) {

     // Test mode: ./strand_soup.exe --test-homogeneous MOTIF REPEATS N_STRANDS
     if (argc == 5 && std::string(argv[1]) == "--test-homogeneous") {
         run_homogeneous(argv[2], std::atoi(argv[3]), std::atoi(argv[4]));
         return 0;
     }

     // --heterogeneous <m_total> <seq1> <seq2> ...
     if (argc >= 4 && std::string(argv[1]) == "--heterogeneous") {
         int m_total = std::atoi(argv[2]);
         std::vector<std::string> seqs;
         for (int a = 3; a < argc; ++a)
             seqs.push_back(argv[a]);
         run_heterogeneous_general(seqs, m_total);
         return 0;
     }
 
     std::cout << "\n==================== Strand Soup ====================\n";
 
     //int m_start = 5;
    // int repeats = 31;
 
     auto start_time = std::chrono::high_resolution_clock::now();
     //run_random_homogeneous_soup(30,  10, 1112);
 
     //run_random_homogeneous_soup(60, m_start, 1111);
     //run_homogeneous("AGGUACCUAAUUGCCUAGAAAACAUGAGGAUCACCCAUG", 1, m_start);
 
    run_homogeneous_soup("CGG", 16, 5);
    //run_homogeneous_soup("CAU", repeats, m_start);
 
    //run_heterogeneous_triplets({"CAG", "GAU", "UAG","CCG"}, 10);
    //run_heterogeneous_general ({"CAGCAGCAGCAG", "CGGCGGCGG", "UUUUUCGAGCAGCA","CUGCUGCUGCUGCUCCUG"}, 12);
     auto end_time = std::chrono::high_resolution_clock::now();
     std::chrono::duration<double> elapsed = end_time - start_time;
     std::cout << "\nTOTAL wall time: " << elapsed.count() << " s\n";
     return 0;
 }