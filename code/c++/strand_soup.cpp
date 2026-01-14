/**
 * @file strand_soup.cpp
 * @brief Implements the strand_soup algorithm for RNA secondary structure prediction.
 *
 * This file contains the implementation of the dynamic programming approach 
 * for RNA structure prediction, including the matrix filling and the backtracking.
 * 
 * 
 * 
 * COMPILATION of this file for testing:
 * Go to the strand_soup.ccp directory where the Makefile should also be and run the following command in the terminal:
 *     make strand_soup
 * Execution:
 * After compilation, run the program with:
 *     ./strand_soup.exe
 *
 * Dependencies:
 * "global_variables.hpp" for global variables.
 * "nussinov.hpp" for the Nussinov algorithm.
 * "utilities.hpp" for helper functions.
 * 
 * The linking is done by the Makefile
 *
 * 
 * @date 2025-03
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

#include "RNAEnergyEvaluator.hpp"
#include "nussinov.hpp"
#include "utilities.hpp"

// just for the colors in the terminal
#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"



//======================================    The Classes    ==============================================//



/**
 * @class Matrix6D
 * @brief A class representing a 6D matrix used in the strand soup algorithm.
 * 
 * This class handles a 6-dimensional matrix with specific sizes for each dimension.
 * The matrix corresponds to the energy matrix used in the RNA folding calculations.
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
extern "C" {
    #include <ViennaRNA/utils.h>
    #include <ViennaRNA/model.h>
}


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
        if (m < 0 || m >= m_size)
            throw std::out_of_range("Index m is out of range");
        if (c < 0 || c >= c_size)
            throw std::out_of_range("Index c is out of range");
        return data[m][s-1][i][r-1][j][c]; // 1-based indexing for strands and bases
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
 * The matrix corresponds to the C and the M matrices used in the full model RNA prediction.
 * It provides methods for accessing and modifying the matrix elements.
 * There are 5 dimensions:
 * - m: Number of dimensions needed in total for m. If we start with m_start strands, we need m = m_start-1 dimensions (easy to see with m_start = 2. once s and r are selected we only have the case m=0,so this dim is of size 2-1 = 1
 * - s: Strand type count (1-based)
 * - i: Max Nucleotide count in one of the strands (1-based)
 * - r: Strand type count (1-based)
 * - j: Max Nucleotide count in one of the strands (1-based)
 * 
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
     * @param c Index of the sixth dimension
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
};







//======================================    Energy matrix functions    ==============================================//




/**
 * @brief Handles the bubble case of the RNA strand soup problem.
 * 
 * This function computes the minimum energy for a bubble case in the RNA strand soup problem. It considers the 4 cases and returns the minimum energy.
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
 * @param params ViennaRNA thermodynamic parameter set.
 * @return float The minimum energy.
 */

 float GeneralCaseMinimization (int m, int s, int i, int r, int j, int c, 
    const std::unordered_map<int, std::string>& strands,
    Matrix5D& C, Matrix5D& M, Matrix6D& F,
    const RNAEnergyEvaluator& evaluator)
{
    float min_value = inf_energy;
    auto* params = evaluator.get_params(); // (not used here yet, but fine)

    // --- Case 1 : i unpaired ---
    // Caller ensures i <= |s|-1 and j >= 1 in general case; borders handled outside.
    min_value = F(m, s, i + 1, r, j, c);

    // --- Case 2 : i pairs within strand s (i,k) ---
    for (int k = i + theta + 1; k <= int(strands.at(s).length()) - 1; ++k) {
        if (can_pair(strands.at(s)[i], strands.at(s)[k])) {
            float value = evaluator.get_C(s)[i][k] + F(m, s, k + 1, r, j, c);
            if (value < min_value) min_value = value;
        }
    }

    // --- Case 3 : i pairs with new strand t at position k ---
    if (m >= 1) {
        for (int t = 1; t <= M.get_s_size(); ++t) {
            for (int m1 = 0; m1 < m; ++m1) {
                int m2 = m - m1 - 1;
                for (int k = 1; k <= int(strands.at(t).length()) - 1; ++k) {
                    if (can_pair(strands.at(s)[i], strands.at(t)[k])) {
                        float value = C(m1, s, i, t, k) + F(m2, t, k + 1, r, j, c);
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

            if (k == int(strands.at(r).length()) - 1) {
                float value = C(m, s, i, r, k);
                if (value < min_value) min_value = value;
            } else {
                float value = C(m, s, i, r, k) + evaluator.get_F(r)[k + 1][j];
                if (value < min_value) min_value = value;
            }
        }
    }

    return min_value;
}


float M1Minimization (int m, int s, int i, int r, int j, 
    const std::unordered_map<int, std::string>& strands, 
    Matrix5D& C, Matrix5D& M1_multi, 
    vrna_param_t *params, const RNAEnergyEvaluator& evaluator) 
{
    float min_value = inf_energy;

    // --- Case 1: extend unpaired ---
    if (i < int(strands.at(s).length()) - 1) {
        float cand = M1_multi(m, s, i + 1, r, j) + params->MLbase;   
        if (cand < min_value) min_value = cand;
    }

    // --- Case 2: helix at (i, j) ---
    vrna_md_t md;
    vrna_md_set_default(&md);

    // Use encoded S1 maps for ptype (cross-strand allowed: use each strand's S1)
    // NOTE: M1_multi here is cross-strand variant, so we don't have evaluator S1 for 'r' directly.
    // If you keep only evaluator's S1 for single strands, you can't get S1[r] here.
    // Minimal fix: get pair type using characters via helper (fallback):
    // Better: pass S1 maps for both strands; but for now, compute via evaluator helper if available.
    // If you stored S1 for all strands in evaluator, do:
    const auto& S1s = evaluator.get_S1_map();
    int type = vrna_get_ptype_md(S1s.at(s)[i], S1s.at(r)[j], &params->model_details);

    float cand2 = C(m, s, i, r, j) + params->MLintern[type];
    if (cand2 < min_value) min_value = cand2;

    return min_value;
}


float Ms1Minimization (int s, int i, int j,
    const std::unordered_map<int, std::string>& strands,
    RNAEnergyEvaluator& evaluator, vrna_param_t *params)
{
    float min_value = inf_energy;

    if (i >= j) return inf_energy;

    // extend unpaired
    if (i < int(strands.at(s).length()) - 1) {
        float cand = evaluator.get_M1(s)[i + 1][j] + params->MLbase;
        if (cand < min_value) min_value = cand;
    }

    // helix at (i,j)
    vrna_md_t md;
    vrna_md_set_default(&md);
    const auto& S1s = evaluator.get_S1_map();
    int type = vrna_get_ptype_md(S1s.at(s)[i], S1s.at(s)[j], &params->model_details);

    float cand2 = evaluator.get_C(s)[i][j] + params->MLintern[type];
    if (cand2 < min_value) min_value = cand2;

    return min_value;
}


/**
 * @brief Handles the closed structure case of the RNA strand soup problem.
 * 
 * This function computes the minimum energy for a bubble case in the RNA strand soup problem. It considers the 4 cases and returns the minimum energy.
 * 
 * @param m The number of strands remaining in the soup.
 * @param s The index of the starting strand. (1-based)
 * @param i The index of the starting nucleotide in strand s. (1-based)
 * @param r The index of the ending strand. (1-based)
 * @param j The index of the ending nucleotide in strand r. (1-based)
 * @param strands A dictionary of strands (index, sequence).
 * @param F The energy matrix.
 * @param nussinov_matrices A dictionary of Nussinov matrices.
 * @return float The minimum energy.
 */


 void ClosedCaseMinimization(int m, int s, int i, int r, int j, 
    const std::unordered_map<int, std::string>& strands, 
    Matrix5D& C, Matrix5D& M1_multi, Matrix6D& F, Matrix5D& M,
    RNAEnergyEvaluator& evaluator, vrna_param_t* params)
{
    float min_value = inf_energy;
    const auto& S1s = evaluator.get_S1_map();

    // --- Case 1 : interior loop (i,j) closed by (k,l) ---
    // Require i < k and l < j to form a proper interior loop.
    for (int k = i + 1; k <= int(strands.at(s).length()) - 1; ++k) {
        for (int l = 1; l <= j - 1; ++l) {
            if (k <= i || l >= j) continue;
            if (!can_pair(strands.at(s)[k], strands.at(r)[l])) continue;
            float loop_energy = evaluator.interior_loop_energy(i, j, k, l, s);
            float value = C(m, s, k, r, l) + loop_energy;
            if (value < min_value) min_value = value;
        }
    }

    // --- Case 2 : leftmost stem within s (i+1 pairs with k), plus multiloop split ---
    for (int k = i + theta + 1; k <= int(strands.at(s).length()) - 1; ++k) {
        if (!can_pair(strands.at(s)[i + 1], strands.at(s)[k])) continue;
        float multiloop_energy = params->MLclosing + evaluator.get_M1(s)[i + 1][k];
        float value = multiloop_energy + M(m, s, k + 1, r, j - 1);
        if (value < min_value) min_value = value;
    }

    // --- Case 3 : leftmost stem with new strand t ---
    if (m >= 1) {
        for (int t = 1; t <= F.get_s_size(); ++t) {
            for (int m1 = 0; m1 < m; ++m1) {
                int m2 = m - m1 - 1;
                for (int k = 1; k <= int(strands.at(t).length()) - 1; ++k) {
                    if (!can_pair(strands.at(s)[i + 1], strands.at(t)[k])) continue;

                    int type = vrna_get_ptype_md(S1s.at(s)[i + 1], S1s.at(t)[k], &params->model_details);

                    float value = params->MLclosing
                                + M1Minimization(m1, s, i + 1, t, k, strands, C, M1_multi, params, evaluator)
                                + params->MLintern[type]
                                + M(m2, t, k + 1, r, j - 1);
                    if (value < min_value) min_value = value;
                }
            }
        }
    }

    // --- Case 4 : leftmost stem with rightmost strand r at k ---
    for (int k = 1; k <= j - 1; ++k) {
        if (!can_pair(strands.at(s)[i + 1], strands.at(r)[k])) continue;

        int type = vrna_get_ptype_md(S1s.at(s)[i + 1], S1s.at(r)[k], &params->model_details);

        if (k == int(strands.at(r).length()) - 1) {
            float value = params->MLclosing
                        + M1Minimization(m, s, i + 1, r, k, strands, C, M1_multi, params, evaluator)
                        + params->MLintern[type];
            if (value < min_value) min_value = value;
        } else {
            float value = params->MLclosing
                        + M1Minimization(m, s, i + 1, r, k, strands, C, M1_multi, params, evaluator)
                        + params->MLintern[type]
                        + evaluator.get_M1(r)[k + 1][j - 1];
            if (value < min_value) min_value = value;
        }
    }

    // --- Case 5: disconnect (external) ---
    if (F(m, s, i + 1, r, j - 1, 0) < min_value) {
        min_value = F(m, s, i + 1, r, j - 1, 0);
    }

    C(m, s, i, r, j) = min_value;
}


/**
 * @brief Handles the closed structure case of the RNA strand soup problem.
 * 
 * This function computes the minimum energy for a bubble case in the RNA strand soup problem. It considers the 4 cases and returns the minimum energy.
 * 
 * @param m The number of strands remaining in the soup.
 * @param s The index of the starting strand. (1-based)
 * @param i The index of the starting nucleotide in strand s. (1-based)
 * @param r The index of the ending strand. (1-based)
 * @param j The index of the ending nucleotide in strand r. (1-based)
 * @param strands A dictionary of strands (index, sequence).
 * @param nussinov_matrices A dictionary of Nussinov matrices.
 * @return float The minimum energy.
 */

 float MultipleCaseMinimization(int m, int s, int i, int r, int j,
    const std::unordered_map<int, std::string>& strands,
    Matrix5D& C, Matrix5D& M1_multi, Matrix5D& M, Matrix6D& F,
    RNAEnergyEvaluator& evaluator, vrna_param_t *params)
{
    float min_value = inf_energy;
    const auto& S1s = evaluator.get_S1_map();

    // --- Case 1 : i unpaired ---
    min_value = M(m, s, i + 1, r, j) + params->MLbase;

    // --- Case 2 : i pairs within s ---
    for (int k = i + theta + 1; k <= int(strands.at(s).length()) - 1; ++k) {
        if (!can_pair(strands.at(s)[i], strands.at(s)[k])) continue;

        int type = vrna_get_ptype_md(S1s.at(s)[i], S1s.at(s)[k], &params->model_details);
        float value = evaluator.get_C(s)[i][k] + params->MLintern[type] + M(m, s, k + 1, r, j);
        if (value < min_value) min_value = value;
    }

    // --- Case 3 : i pairs with new strand t ---
    if (m >= 1) {
        for (int t = 1; t <= F.get_s_size(); ++t) {
            for (int m1 = 0; m1 < m; ++m1) {
                int m2 = m - m1 - 1;
                for (int k = 1; k <= int(strands.at(t).length()) - 1; ++k) {
                    if (!can_pair(strands.at(s)[i], strands.at(t)[k])) continue;

                    int type = vrna_get_ptype_md(S1s.at(s)[i], S1s.at(t)[k], &params->model_details);
                    float value = C(m1, s, i, t, k) + params->MLintern[type] + M(m2, t, k + 1, r, j);
                    if (value < min_value) min_value = value;
                }
            }
        }
    }

    // --- Case 4 : i pairs with last strand r at k ---
    for (int k = 1; k <= j - 1; ++k) {
        if (!can_pair(strands.at(s)[i], strands.at(r)[k])) continue;

        int type = vrna_get_ptype_md(S1s.at(s)[i], S1s.at(r)[k], &params->model_details);

        if (k == int(strands.at(r).length()) - 1) {
            float value = C(m, s, i, r, k) + params->MLintern[type];
            if (value < min_value) min_value = value;
        } else {
            float value = params->MLclosing
                        + M1Minimization(m, s, i + 1, r, k, strands, C, M1_multi, params, evaluator)
                        + params->MLintern[type]
                        + evaluator.get_M1(r)[k + 1][j];
            if (value < min_value) min_value = value;
        }
    }

    return min_value;
}


/// @brief this is the case of multiloops
/// @param m 
/// @param s 
/// @param i 
/// @param r 
/// @param j 
/// @param strands 
/// @param C 
/// @param M 
/// @param nussinov_matrices 

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
 * @brief Computes the 6D matrix for the RNA strand soup problem
 * 
 * @param strands The dictionary of strands (index, sequence).
 * @param C The closed matrix.
 * @param M The multiloop matrix.
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
            for (int i = F.get_i_size() + 1; i >= 0; i--) {
                for (int r = 1; r <= F.get_r_size(); r++) {
                    for (int j = 0; j <= F.get_j_size() + 1; j++) {

                          // ---- NEW: only compute C and M when both strands are non-empty ----
                          bool s_nonempty = (i >= 1 && i <= int(strands.at(s).length()) - 1);
                          bool r_nonempty = (j >= 1 && j <= int(strands.at(r).length()) - 1);
                      
                          if (s_nonempty && r_nonempty) {
                              ClosedCaseMinimization(
                                  m, s, i, r, j,
                                  strands,
                                  C, M1_multi, F, M,
                                  evaluator, params
                              );
                              MMatrixMinimization(
                                  m, s, i, r, j,
                                  strands,
                                  C, M1_multi, M, F,
                                  evaluator, params
                              );
                          } else {
                              C(m, s, i, r, j) = inf_energy;
                              M(m, s, i, r, j) = inf_energy;
                          }
  

                        for (int c = 0; c <= 1; c++) {
                            if (i > int(strands.at(s).length() - 1)) {
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
                                            F(m, s, i, r, j, c) = evaluator.get_M1(r)[1][j];
                                            continue;
                                        }
                                    } else {
                                        // Minimum over all t: M(m-1, t, 1, r, j, 1)
                                        float min_value = inf_energy;
                                        for (int t = 1; t <= M.get_s_size(); t++) {
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
                                                    evaluator.get_M1(s)[i][strands.at(s).length() - 1];
                                                    continue;
                                            }
                                        } else {
                                            // Minimum over all t: M(m-1, s, i, t, |t|-1, 1)
                                            float min_value = inf_energy;
                                            for (int t = 1; t <= M.get_s_size(); t++) {
                                                float val = F(m - 1, s, i, t,
                                                    strands[t].length() - 1, 1);
                                                if (val < min_value) min_value = val;
                                            }
                                            F(m, s, i, r, j, c) = min_value;
                                            continue;
                                        }
                                    }
                                } else {
                                    // === GENERAL CASE ===
                                    F(m, s, i, r, j, c) =
                                        GeneralCaseMinimization(m, s, i, r, j, c,
                                            strands, C, M, F, evaluator);
                                }
                            }
                        } // end for c
                    } // end for j
                } // end for r
            } // end for i
        } // end for s
    } // end for m
}//======================================    Backtrack functions    ==============================================//



/**
 * @brief Finds all the starting point in the energy matrix, which corresponds to the minimum energy (excluding the border effect).
 * 
 * This function searches through the filled energy matrix (M) and identifies the minimum energy value. It returns all the coordinates (m,s,i,r,j,c) of the minimum energy, excluding the borders.
 * 
 * @param M The 6D energy matrix that stores the computed energy values for all possible configurations.
 * 
 * @return std::vector<std::vector<int>> A vector containing all the starting points (m,s,i,r,j,c) where the minimum energy is located.
 */
std::vector<std::vector<int>> Find_all_start_backtrack(Matrix6D& F) {
    std::vector<std::vector<int>> starting_points;
    float min_value = inf_energy;

    // We assume final optimum is at m = F.get_m_size()-1, i = 1, j = max, c = 1
    int m_last = F.get_m_size() - 1;
    int i_start = 1;
    int j_max   = F.get_j_size();

    for (int s = 1; s <= F.get_s_size(); ++s) {
        for (int r = 1; r <= F.get_r_size(); ++r) {
            float val = F(m_last, s, i_start, r, j_max, 1);
            if (val < min_value) {
                min_value = val;
                starting_points.clear();
            }
            if (val == min_value) {
                starting_points.push_back({m_last, s, i_start, r, j_max, 1});
            }
        }
    }

    if (min_value == inf_energy) {
        throw std::runtime_error(
            "No valid starting point found by Find_all_start_backtrack in the given energy matrix"
        );
    }
    return starting_points;
}

/**
 * ============================================================
 *  Single-strand thermodynamic backtracking functions
 *  (Turner/Vienna model version)
 * ============================================================
 *  These replace the old Nussinov backtrack and reconstruct the
 *  optimal secondary structure of a single RNA using the
 *  thermodynamic matrices: F, C, M, and M1.
 *  Each function returns an output_backtrack object so that the
 *  same data structure can be merged into multistrand backtracking.
 * ============================================================
 */

 output_backtrack backtrack_Fs(int s, int i, int j, RNAEnergyEvaluator& evaluator, int theta);
 output_backtrack backtrack_Cs(int s, int i, int j, RNAEnergyEvaluator& evaluator, int theta);
 output_backtrack backtrack_Ms(int s, int i, int j, RNAEnergyEvaluator& evaluator, int theta);
 output_backtrack backtrack_M1s(int s, int i, int j, RNAEnergyEvaluator& evaluator, int theta);
 
 /*-------------------------------------------------------------
  *  F-matrix (external / free region)
  *------------------------------------------------------------*/
 output_backtrack backtrack_Fs(int s, int i, int j, RNAEnergyEvaluator& evaluator, int theta) {
    output_backtrack out;
    if (i >= j) return out;

    const auto& F = evaluator.get_F(s);
    const auto& C = evaluator.get_C(s);
    const std::string& seq = evaluator.get_strands().at(s);

    int best = F[i][j];

    // Case 1: i unpaired
    if (i + 1 <= j && best == F[i + 1][j]) {
        return backtrack_Fs(s, i + 1, j, evaluator, theta);
    }

    // Case 2: i paired with k
    for (int k = i + theta + 1; k <= j; ++k) {
        if (!can_pair(seq[i - 1], seq[k - 1])) continue;

        int cand = C[i][k];
        if (k + 1 <= j) cand += F[k + 1][j];

        if (best == cand) {
            output_backtrack left  = backtrack_Cs(s, i, k, evaluator, theta);
            output_backtrack right = backtrack_Fs(s, k + 1, j, evaluator, theta);
            left.merge(right);
            return left;
        }
    }

    return out;
}  // <--- only one closing brace here

 
 /*-------------------------------------------------------------
  *  C-matrix (closed pair region)
  *------------------------------------------------------------*/
 output_backtrack backtrack_Cs(int s, int i, int j, RNAEnergyEvaluator& evaluator, int theta) {
     output_backtrack out;
     const auto& C = evaluator.get_C(s);
     const auto& M = evaluator.get_M(s);
     const auto& M1 = evaluator.get_M1(s);
     const std::string& seq = evaluator.get_strands().at(s);
     const vrna_param_t* params = evaluator.get_params();
 
     int best = C[i][j];
     out.add_pair(1, i, 1, j); // record (i,j)
 
     // --- Hairpin ---
     if (best == evaluator.hairpin_energy(i, j, s)) {
         return out;
     }
 
     // --- Internal loops ---
     for (int p = i + 1; p < j - theta - 1; ++p) {
         for (int q = p + theta + 1; q < j; ++q) {
             if (!can_pair(seq[p - 1], seq[q - 1])) continue;
             int e_int = evaluator.interior_loop_energy(i, j, p, q, s);
             if (best == C[p][q] + e_int) {
                 output_backtrack inner = backtrack_Cs(s, p, q, evaluator, theta);
                 out.merge(inner);
                 return out;
             }
         }
     }
 
     // --- Multiloop closure ---
     for (int u = i + 1; u < j; ++u) {
         if (best == M1[i + 1][u - 1] + M[u][j - 1] + params->MLclosing) {
             output_backtrack left = backtrack_M1s(s, i + 1, u - 1, evaluator, theta);
             output_backtrack right = backtrack_Ms(s, u, j - 1, evaluator, theta);
             out.merge(left);
             out.merge(right);
             return out;
         }
     }
 
     return out; // fallback
 }
 
 /*-------------------------------------------------------------
  *  M-matrix (multiloop segment)
  *------------------------------------------------------------*/
 output_backtrack backtrack_Ms(int s, int i, int j, RNAEnergyEvaluator& evaluator, int theta) {
     output_backtrack out;
     const auto& M = evaluator.get_M(s);
     const auto& C = evaluator.get_C(s);
     const auto& S1 = evaluator.get_S1_map().at(s);
     vrna_param_t* params = evaluator.get_params();

 
     int best = M[i][j];
 
     // Case 1: j unpaired
     if (best == M[i][j - 1] + params->MLbase) {
         return backtrack_Ms(s, i, j - 1, evaluator, theta);
     }
 
     // Case 2: segmentation (M | C)
     for (int u = i; u < j; ++u) {
         int type = vrna_get_ptype_md(S1[u + 1], S1[j], &params->model_details);
         int penalty = params->MLintern[type];
         if (best == M[i][u] + C[u + 1][j] + penalty) {
             output_backtrack left = backtrack_Ms(s, i, u, evaluator, theta);
             output_backtrack right = backtrack_Cs(s, u + 1, j, evaluator, theta);
             left.merge(right);
             return left;
         }
     }
 
     // Case 3: start directly from a C
     for (int u = i; u < j; ++u) {
         int type = vrna_get_ptype_md(S1[u + 1], S1[j], &params->model_details);
         int penalty = params->MLintern[type];
         if (best == C[u + 1][j] + penalty) {
             return backtrack_Cs(s, u + 1, j, evaluator, theta);
         }
     }
 
     return out;
 }
 
 /*-------------------------------------------------------------
  *  M1-matrix (multiloop entry segment)
  *------------------------------------------------------------*/
 output_backtrack backtrack_M1s(int s, int i, int j, RNAEnergyEvaluator& evaluator, int theta) {
     output_backtrack out;
     const auto& M1 = evaluator.get_M1(s);
     const auto& C = evaluator.get_C(s);
     const auto& S1 = evaluator.get_S1_map().at(s);
     vrna_param_t* params = evaluator.get_params();
 
     int best = M1[i][j];
 
     // Case 1: j unpaired
     if (best == M1[i][j - 1] + params->MLbase) {
         return backtrack_M1s(s, i, j - 1, evaluator, theta);
     }
 
     // Case 2: helix at (i,j)
     int type = vrna_get_ptype_md(S1[i], S1[j], &params->model_details);
     int penalty = params->MLintern[type];
     if (best == C[i][j] + penalty) {
         output_backtrack inner = backtrack_Cs(s, i, j, evaluator, theta);
         out.merge(inner);
         return out;
     }
 
     return out;}

/**
 * ============================================================
 
 * ============================================================
 */

 output_backtrack backtrack_F_multi(
    int m, int s, int i, int r, int j, int c,
    const std::unordered_map<int, std::string>& strands,
    Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
    RNAEnergyEvaluator& evaluator,
    vrna_param_t* params,
    int theta
);

output_backtrack backtrack_C_multi(
    int m, int s, int i, int r, int j,
    const std::unordered_map<int, std::string>& strands,
    Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
    RNAEnergyEvaluator& evaluator,
    vrna_param_t* params,
    int theta
);

output_backtrack backtrack_M_multi(
    int m, int s, int i, int r, int j,
    const std::unordered_map<int, std::string>& strands,
    Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
    RNAEnergyEvaluator& evaluator,
    vrna_param_t* params,
    int theta
);

output_backtrack backtrack_M1_multi(
    int m, int s, int i, int r, int j,
    const std::unordered_map<int, std::string>& strands,
    Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
    RNAEnergyEvaluator& evaluator,
    vrna_param_t* params,
    int theta
);
output_backtrack backtrack_F_multi(
    int m, int s, int i, int r, int j, int c,
    const std::unordered_map<int, std::string>& strands,
    Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
    RNAEnergyEvaluator& evaluator,
    vrna_param_t* params,
    int theta
) {
    output_backtrack out;

    const int len_s = int(strands.at(s).length()) - 1; // because of '$'
    const int len_r = int(strands.at(r).length()) - 1;

    float best = F(m, s, i, r, j, c);

    // Border sanity: if we ever end up at an impossible cell, just return empty.
    if (best >= inf_energy) return out;

    bool s_nonempty = (i <= len_s);
    bool r_nonempty = (j >= 1);

    // ============================================================
    // CASE A: s empty (i > len_s) → mirror DP "s empty" block
    // ============================================================
    if (!s_nonempty) {
        if (c == 1) {
            // F was set to +inf here; nothing to backtrack
            return out;
        } else {
            if (m == 0) {
                if (j == 0) {
                    // F = 0; nothing to do
                    return out;
                } else {
                    // F(m,s,i,r,j,c) = evaluator.get_F(r)[1][j]
                    // => single-strand folding on strand r
                    output_backtrack single = backtrack_M1s(r, 1, j, evaluator, theta);
                    // pairs in 'single' use strand index 1;
                    // they correspond to strand r, which will be added in full_backtrack.
                    return single;
                }
            } else {
                // min_t F(m-1, t, 1, r, j, 1)
                for (int t = 1; t <= F.get_s_size(); ++t) {
                    float cand = F(m - 1, t, 1, r, j, 1);
                    if (cand == best) {
                        output_backtrack sub =
                            backtrack_F_multi(m - 1, t, 1, r, j, 1,
                                              strands, C, M, M1_multi, F, evaluator, params, theta);
                        // t is a new sequence in the middle block
                        sub.add_sequence(t);
                        sub.shift(1); // shift strand indices because in the final structure,
                                      // there will be an extra leftmost starting strand
                        return sub;
                    }
                }
            }
        }
        std::cerr << "ERROR: backtrack_F_multi reached inconsistent 's empty' branch\n";
        return out;
    }

    // ============================================================
    // CASE B: r empty (j < 1) → mirror DP "r empty" block
    // ============================================================
    if (!r_nonempty) {
        if (c == 1) {
            return out;
        } else {
            if (m == 0) {
                if (i == 0) {
                    return out;
                } else {
                    // F = single-strand folding on s from i..len_s
                    output_backtrack single = backtrack_M1s(s, i, len_s, evaluator, theta);

                    return single;
                }
            } else {
                // min_t F(m-1, s, i, t, |t|-1, 1)
                for (int t = 1; t <= F.get_s_size(); ++t) {
                    int len_t = int(strands.at(t).length()) - 1;
                    float cand = F(m - 1, s, i, t, len_t, 1);
                    if (cand == best) {
                        output_backtrack sub =
                            backtrack_F_multi(m - 1, s, i, t, len_t, 1,
                                              strands, C, M, M1_multi, F, evaluator, params, theta);
                        sub.add_sequence(t); // t appended on the right side of sequences
                        return sub;
                    }
                }
            }
        }
        std::cerr << "ERROR: backtrack_F_multi reached inconsistent 'r empty' branch\n";
        return out;
    }

    // ============================================================
    // CASE C: general case (both s and r non-empty)
    // Mirror GeneralCaseMinimization
    // ============================================================

    // --- Case 1 : i unpaired ---
    if (F(m, s, i + 1, r, j, c) == best) {
        return backtrack_F_multi(m, s, i + 1, r, j, c,
                                 strands, C, M, M1_multi, F, evaluator, params, theta);
    }

    // --- Case 2 : i pairs within s (i,k), closed by single-strand C(s) ---
    {
        const std::string& seq_s = strands.at(s);
        for (int k = i + theta + 1; k <= len_s; ++k) {
            if (!can_pair(seq_s[i], seq_s[k])) continue;  // i,k are 1-based and seq has '$' at 0
            float cand = evaluator.get_C(s)[i][k] + F(m, s, k + 1, r, j, c);
            if (cand == best) {
                // closed structure on single strand s between i..k
                output_backtrack left = backtrack_Cs(s, i, k, evaluator, theta);
                // rest from k+1.. on s,r
                output_backtrack right = backtrack_F_multi(
                    m, s, k + 1, r, j, c, strands, C, M, M1_multi, F, evaluator, params, theta);
                left.merge(right);
                return left;
            }
        }
    }

    // --- Case 3 : i pairs with new strand t at k (uses multi-strand C) ---
    if (m >= 1) {
        const std::string& seq_s = strands.at(s);
        for (int t = 1; t <= M.get_s_size(); ++t) {
            const std::string& seq_t = strands.at(t);
            int len_t = int(seq_t.length()) - 1;
            for (int m1 = 0; m1 < m; ++m1) {
                int m2 = m - m1 - 1;
                for (int k = 1; k <= len_t; ++k) {
                    if (!can_pair(seq_s[i], seq_t[k])) continue;
                    float cand = C(m1, s, i, t, k) + F(m2, t, k + 1, r, j, c);
                    if (cand == best) {
                        // multi-strand closed region between i (s) and k (t)
                        output_backtrack left =
                            backtrack_C_multi(m1, s, i, t, k,
                                              strands, C, M, M1_multi, F, evaluator, params, theta);
                        // F on remainder starting from t,k+1
                        output_backtrack right =
                            backtrack_F_multi(m2, t, k + 1, r, j, c,
                                              strands, C, M, M1_multi, F, evaluator, params, theta);

                        // Strand index mapping (like old code):
                        // s is strand 1; added 't' lies somewhere in the soup and will be
                        // accounted for via list_of_sequences and shift logic.  
                        // You may need to adjust shifts as you did in bubble_backtrack.
                        left.merge(right);
                        return left;
                    }
                }
            }
        }
    }

    // --- Case 4 : i pairs with r at k (multi-strand C(m,s,i,r,k))
    if (c == 1) {
        const std::string& seq_s = strands.at(s);
        const std::string& seq_r = strands.at(r);
        for (int k = 1; k <= j; ++k) {
            if (!can_pair(seq_s[i], seq_r[k])) continue;

            if (k == len_r) {
                float cand = C(m, s, i, r, k);
                if (cand == best) {
                    // closed region C(m,s,i,r,k)
                    output_backtrack left =
                        backtrack_C_multi(m, s, i, r, k,
                                          strands, C, M, M1_multi, F, evaluator, params, theta);
                    // no extra F on the right
                    return left;
                }
            } else {
                float cand = C(m, s, i, r, k) + evaluator.get_F(r)[k + 1][j];
                if (cand == best) {
                    output_backtrack left =
                        backtrack_C_multi(m, s, i, r, k,
                                          strands, C, M, M1_multi, F, evaluator, params, theta);
                    output_backtrack right =
                        backtrack_Fs(r, k + 1, j, evaluator, theta);
                    // right has strand index 1; in multi-strand setting r is not 1.
                    // You will likely want to shift right by (1 + m) as in original code.
                    right.shift(1); // minimal consistent shift; tune if needed
                    left.merge(right);
                    return left;
                }
            }
        }
    }

    std::cerr << "ERROR: backtrack_F_multi could not match any GeneralCase branch\n";
    return out;
}
output_backtrack backtrack_C_multi(
    int m, int s, int i, int r, int j,
    const std::unordered_map<int, std::string>& strands,
    Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
    RNAEnergyEvaluator& evaluator,
    vrna_param_t* params,
    int theta
) {
    output_backtrack out;

    float best = C(m, s, i, r, j);
    if (best >= inf_energy) return out;

    const std::string& seq_s = strands.at(s);
    const std::string& seq_r = strands.at(r);
    const auto& S1s = evaluator.get_S1_map();

    int len_s = int(seq_s.length()) - 1;
    int len_r = int(seq_r.length()) - 1;

    // First, record the closing pair (i,j) between s and r
    // In the "abstract" strand indexing, s is strand 1, r is 1+m+1
    out.add_pair(1, i, 1 + m + 1, j);

    // --- Case 1 : interior loop, closed by (k,l) ---
    for (int k = i + 1; k <= len_s; ++k) {
        for (int l = 1; l <= j - 1; ++l) {
            if (k <= i || l >= j) continue;
            if (!can_pair(seq_s[k], seq_r[l])) continue;
            float loop_energy = evaluator.interior_loop_energy(i, j, k, l, s);
            float cand = C(m, s, k, r, l) + loop_energy;
            if (cand == best) {
                output_backtrack inner =
                    backtrack_C_multi(m, s, k, r, l,
                                      strands, C, M, M1_multi, F, evaluator, params, theta);
                out.merge(inner);
                return out;
            }
        }
    }

    // --- Case 2 : leftmost stem within s (i+1 pairs k), multiloop split ---
    for (int k = i + theta + 1; k <= len_s; ++k) {
        if (!can_pair(seq_s[i + 1], seq_s[k])) continue;
        float ml_energy = params->MLclosing + evaluator.get_M1(s)[i + 1][k];
        float cand = ml_energy + M(m, s, k + 1, r, j - 1);
        if (cand == best) {
            // left: M1 on single strand s
            output_backtrack left = backtrack_M1s(s, i + 1, k, evaluator, theta);
            // right: multi-strand M
            output_backtrack right =
                backtrack_M_multi(m, s, k + 1, r, j - 1,
                                  strands, C, M, M1_multi, F, evaluator, params, theta);
            out.merge(left);
            out.merge(right);
            return out;
        }
    }

    // --- Case 3 : leftmost stem with new strand t ---
    if (m >= 1) {
        for (int t = 1; t <= F.get_s_size(); ++t) {
            const std::string& seq_t = strands.at(t);
            int len_t = int(seq_t.length()) - 1;
            for (int m1 = 0; m1 < m; ++m1) {
                int m2 = m - m1 - 1;
                for (int k = 1; k <= len_t; ++k) {
                    if (!can_pair(seq_s[i + 1], seq_t[k])) continue;

                    int type = vrna_get_ptype_md(
                        S1s.at(s)[i + 1], S1s.at(t)[k], &params->model_details);

                    float cand =
                        params->MLclosing
                      + M1Minimization(m1, s, i + 1, t, k, strands, C, M1_multi, params, evaluator)
                      + params->MLintern[type]
                      + M(m2, t, k + 1, r, j - 1);

                      if (cand == best) {
                        output_backtrack left =
                            backtrack_M1_multi(m1, s, i + 1, t, k,
                                               strands, C, M, M1_multi, F, evaluator, params, theta);
                        output_backtrack right =
                            backtrack_M_multi(m2, t, k + 1, r, j - 1,
                                              strands, C, M, M1_multi, F, evaluator, params, theta);
                    
                        // If you need shifts, apply them here before merging
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
        if (!can_pair(seq_s[i + 1], seq_r[k])) continue;

        int type = vrna_get_ptype_md(
            S1s.at(s)[i + 1], S1s.at(r)[k], &params->model_details);

        if (k == len_r) {
            float cand =
                params->MLclosing
              + M1Minimization(m, s, i + 1, r, k, strands, C, M1_multi, params, evaluator)
              + params->MLintern[type];
            if (cand == best) {
                output_backtrack left =
                    backtrack_M1_multi(m, s, i + 1, r, k,
                                       strands, C, M, M1_multi, F, evaluator, params, theta);
                out.merge(left);
                return out;
            }
        } else {
            float cand =
                params->MLclosing
              + M1Minimization(m, s, i + 1, r, k, strands, C, M1_multi, params, evaluator)
              + params->MLintern[type]
              + evaluator.get_M1(r)[k + 1][j - 1];

            if (cand == best) {
                output_backtrack left =
                    backtrack_M1_multi(m, s, i + 1, r, k,
                                       strands, C, M, M1_multi, F, evaluator, params, theta);
                output_backtrack right =
                    backtrack_M1s(r, k + 1, j - 1, evaluator, theta);
                right.shift(1); // as usual, multi-strand index shift
                out.merge(left);
                out.merge(right);
                return out;
            }
        }
    }

    // --- Case 5: external disconnect ---
    if (F(m, s, i + 1, r, j - 1, 0) == best) {
        output_backtrack inner =
            backtrack_F_multi(m, s, i + 1, r, j - 1, 0,
                              strands, C, M, M1_multi, F, evaluator, params, theta);
        out.merge(inner);
        return out;
    }

    std::cerr << "ERROR: backtrack_C_multi could not match any ClosedCase branch\n";
    return out;
}

output_backtrack backtrack_M_multi(
    int m, int s, int i, int r, int j,
    const std::unordered_map<int, std::string>& strands,
    Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
    RNAEnergyEvaluator& evaluator,
    vrna_param_t* params,
    int theta
) {
    output_backtrack out;

    float best = M(m, s, i, r, j);
    if (best >= inf_energy) return out;

    const std::string& seq_s = strands.at(s);
    const std::string& seq_r = strands.at(r);
    const auto& S1s = evaluator.get_S1_map();

    int len_s = int(seq_s.length()) - 1;
    int len_r = int(seq_r.length()) - 1;

    // --- Case 1 : i unpaired ---
    if (M(m, s, i + 1, r, j) + params->MLbase == best) {
        return backtrack_M_multi(m, s, i + 1, r, j,
                                 strands, C, M, M1_multi, F, evaluator, params, theta);
    }

    // --- Case 2 : i pairs within s ---
    for (int k = i + theta + 1; k <= len_s; ++k) {
        if (!can_pair(seq_s[i], seq_s[k])) continue;

        int type = vrna_get_ptype_md(
            S1s.at(s)[i], S1s.at(s)[k], &params->model_details);
        float cand = evaluator.get_C(s)[i][k] + params->MLintern[type]
                   + M(m, s, k + 1, r, j);
        if (cand == best) {
            output_backtrack left = backtrack_Cs(s, i, k, evaluator, theta);
            output_backtrack right =
                backtrack_M_multi(m, s, k + 1, r, j,
                                  strands, C, M, M1_multi, F, evaluator, params, theta);
            left.merge(right);
            return left;
        }
    }

    // --- Case 3 : i pairs with new strand t ---
    if (m >= 1) {
        for (int t = 1; t <= F.get_s_size(); ++t) {
            const std::string& seq_t = strands.at(t);
            int len_t = int(seq_t.length()) - 1;
            for (int m1 = 0; m1 < m; ++m1) {
                int m2 = m - m1 - 1;
                for (int k = 1; k <= len_t; ++k) {
                    if (!can_pair(seq_s[i], seq_t[k])) continue;

                    int type = vrna_get_ptype_md(
                        S1s.at(s)[i], S1s.at(t)[k], &params->model_details);

                    float cand = C(m1, s, i, t, k)
                               + params->MLintern[type]
                               + M(m2, t, k + 1, r, j);
                    if (cand == best) {
                        output_backtrack left =
                            backtrack_C_multi(m1, s, i, t, k,
                                              strands, C, M, M1_multi, F, evaluator, params, theta);
                        output_backtrack right =
                            backtrack_M_multi(m2, t, k + 1, r, j,
                                              strands, C, M, M1_multi, F, evaluator, params, theta);
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
            float cand = C(m, s, i, r, k) + params->MLintern[type];
            if (cand == best) {
                output_backtrack left =
                    backtrack_C_multi(m, s, i, r, k,
                                      strands, C, M, M1_multi, F, evaluator, params, theta);
                return left;
            }
        } else {
            float cand =
                params->MLclosing
              + M1Minimization(m, s, i + 1, r, k, strands, C, M1_multi, params, evaluator)
              + params->MLintern[type]
              + evaluator.get_M1(r)[k + 1][j];

            if (cand == best) {
                output_backtrack left =
                    backtrack_M1_multi(m, s, i + 1, r, k,
                                       strands, C, M, M1_multi, F, evaluator, params, theta);
                output_backtrack right =
                    backtrack_M1s(r, k + 1, j, evaluator, theta);
                right.shift(1);
                left.merge(right);
                return left;
            }
        }
    }

    std::cerr << "ERROR: backtrack_M_multi could not match any MultipleCase branch\n";
    return out;
}

output_backtrack backtrack_M1_multi(
    int m, int s, int i, int r, int j,
    const std::unordered_map<int, std::string>& strands,
    Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
    RNAEnergyEvaluator& evaluator,
    vrna_param_t* params,
    int theta
) {
    output_backtrack out;

    float best = M1_multi(m, s, i, r, j);
    if (best >= inf_energy) return out;

    const auto& S1s = evaluator.get_S1_map();

    int len_s = int(strands.at(s).length()) - 1;

    // --- Case 1: extend unpaired (i -> i+1) ---
    if (i < len_s) {
        float cand = M1_multi(m, s, i + 1, r, j) + params->MLbase;
        if (cand == best) {
            return backtrack_M1_multi(m, s, i + 1, r, j,
                                      strands, C, M, M1_multi, F, evaluator, params, theta);
        }
    }

    // --- Case 2: helix at (i,j) using C(m,s,i,r,j) ---
    int type = vrna_get_ptype_md(
        S1s.at(s)[i], S1s.at(r)[j], &params->model_details);

    float cand2 = C(m, s, i, r, j) + params->MLintern[type];
    if (cand2 == best) {
        output_backtrack inner =
            backtrack_C_multi(m, s, i, r, j,
                              strands, C, M, M1_multi, F, evaluator, params, theta);
        out.merge(inner);
        return out;
    }

    std::cerr << "ERROR: backtrack_M1_multi could not match any M1Minimization branch\n";
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

    std::vector<std::vector<int>> starting_points;
    try {
        starting_points = Find_all_start_backtrack(F);
    } catch (const std::runtime_error& e) {
        std::cerr << "Warning: " << e.what() << " Setting repartition to zero.\n";
        return secondary_structures; // empty
    }

    for (auto& sp : starting_points) {
        int m = sp[0];
        int s = sp[1];
        int i = sp[2];
        int r = sp[3];
        int j = sp[4];
        int c = sp[5];

        std::cout << "Starting point: F("
                  << m << "," << s << "," << i << ","
                  << r << "," << j << "," << c << ") = "
                  << F(m, s, i, r, j, c) << std::endl;

        output_backtrack secondary =
            backtrack_F_multi(m, s, i, r, j, c,
                              strands, C, M, M1_multi, F,
                              evaluator, params, theta);

        // As in old code: add r as last, s as first — your abstract sequence ordering
        secondary.add_sequence(r);
        secondary.add_sequence_front(s);
        secondary.print();

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
 * STEP 2,3 : Then for each combination, it generates the corresponding RNA strands and computes the Nussinov matrices.
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
         std::string bare1   = generate_triplet_repeat(triplet1, number_of_repeats);
         std::string strand1 = with_dollar(bare1); // add '$' for 1-based indexing
 
         internal_file   << triplet1 << ",";
         homogeneous_file << triplet1 << ",";
         heterogeneous_file << triplet1 << ",";
 
         for (const auto& triplet2 : triplets) {
             std::cout << "\n========== Triplet : " << triplet1
                       << " and " << triplet2 << " ==========" << std::endl;
 
             std::string bare2   = generate_triplet_repeat(triplet2, number_of_repeats);
             std::string strand2 = with_dollar(bare2); // add '$'
 
             // Create strands map (index -> sequence with $)
             std::unordered_map<int, std::string> strands;
             add_strand_if_unique(strands, strand1);
             add_strand_if_unique(strands, strand2);
             // NOTE: add_strand_if_unique is assumed not to add '$' itself
 
             // STEP 3: Build RNAEnergyEvaluator (computes single-strand F/C/M/M1)
             RNAEnergyEvaluator evaluator(strands);
             vrna_param_t* params = evaluator.get_params();
 
             // STEP 4: Compute multi-strand matrices
             int m_size = m_start - 1;                // as in old code
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
                         // same strand index in the abstract ordering
                         internal++;
                     } else {
                         // Using list_of_sequences as in old logic:
                         // strand_idx are 1-based indices in list_of_sequences.
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
 
             // Add a comma only if it's not the last column in this row
             if (triplet2 != triplets.back()) {
                 internal_file   << ",";
                 homogeneous_file << ",";
                 heterogeneous_file << ",";
             }
 
             // Update progress
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
//=============================================================================================//
int main() {
    std::cout << "\n==================== Strand Soup ====================" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    std::cout << "\n========== Setting the parameters ==========" << std::endl;
    int m_start = 4;                // number of strands in total
    int number_of_repeats = 10;     // triplet repeats per strand
    int sequence_length = number_of_repeats * 3; // real #bases
    std::cout << "  Number of strands : m = " << m_start << std::endl;
    std::cout << "  Length of the sequences : " << sequence_length << std::endl;
    std::cout << "  Theta : " << theta << std::endl;
    std::cout << "\n========== Generation of the strands ==========" << std::endl;
    std::unordered_map<int, std::string> strands;
    strands[1] = with_dollar(generate_triplet_repeat("GUU", number_of_repeats));
    strands[2] = with_dollar(generate_triplet_repeat("CAG", number_of_repeats));
    strands[3] = with_dollar(generate_triplet_repeat("ACG", number_of_repeats));
    std::cout << "Generated strands:" << std::endl;
    for (const auto& p : strands) {
        std::cout << "  Index " << p.first << " : " << p.second << std::endl;
    }

    std::cout << "\n========== Single-strand thermodynamic matrices ==========" << std::endl;
    RNAEnergyEvaluator evaluator(strands);
    vrna_param_t* params = evaluator.get_params();

    // If you want to print single-strand matrices for debugging:
    // for (const auto& p : strands) {
    //     int s = p.first;
    //     const auto& F_s = evaluator.get_F(s);
    //     const auto& C_s = evaluator.get_C(s);
    //     const auto& M_s = evaluator.get_M(s);
    //     const auto& M1_s = evaluator.get_M1(s);
    //     // print them as you like
    // }

    std::cout << "========== MainAuxiliaryMatrix (multi-strand) ==========" << std::endl;
    int m_size = m_start - 1;           // same logic as before
    int s_size = int(strands.size());   // number of distinct strand types
    int i_size = sequence_length;       // real #bases, indices 1..sequence_length
    int r_size = int(strands.size());
    int j_size = sequence_length;
    int c_size = 2;                     // connectivity bit

    Matrix5D C(m_size, s_size, i_size, r_size, j_size);
    Matrix5D M(m_size, s_size, i_size, r_size, j_size);
    Matrix5D M1_multi(m_size, s_size, i_size, r_size, j_size);
    Matrix6D F(m_size, s_size, i_size, r_size, j_size, c_size);

    MainAuxiliaryMatrix(strands, C, M, M1_multi, F, evaluator, params);

    std::cout << std::endl << "========== full_backtrack ==========" << std::endl;
    std::vector<output_backtrack> structures =
        full_backtrack(strands, C, M, M1_multi, F, evaluator, params, theta);

    // Just to show them:
    for (const auto& sec : structures) {
        sec.print();
    }

    std::cout << std::endl << "========== compute_all_triplet_repartition ==========" << std::endl;
    // Example: around 30 minutes for m_start = 2 and length = 10*3, etc.
    // compute_all_triplet_repartition(6, 4); // enable when you want the full sweep

std::cout << std::endl << "==================== End of the program ====================" << std::endl;
auto end_time = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed_time = end_time - start_time;
std::cout << "Computation completed in " << elapsed_time.count() << " seconds." << std::endl << std::endl;

    return 0;
}