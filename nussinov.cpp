/**
 * @file nussinov.cpp
 * @brief Implements the Nussinov algorithm for RNA secondary structure prediction.
 *
 * This file contains the implementation of the dynamic programming approach 
 * for RNA structure prediction, including matrix filling and backtracking.
 * 
 * The FillMatrix function is also used in the strand soup algorithm. (The backtrack for nussinov is reimplemented inside strand_soup.cpp with an other output format)
 * 
 * 
 * 
 * COMPILATION of this file for testing:
 * Go to the nussinov.ccp and Makefile directory and run the following command in the terminal:
 *     make nussinov
 * Execution:
 * After compilation, run the program with:
 *     ./nussinov.exe
 *
 * Dependencies:
 * - Requires "utilities.hpp" and "global_variables.hpp" for helper functions and global variables. The linking is done by the Makefile
 *
 * 
 * @date 2025-03
 */


#include "nussinov.hpp"
#include "utilities.hpp"
#include "global_variables.hpp"

#include <iostream>
#include <vector>
#include <string>





/**
 * @brief Fills the minimum free energy (MFE) matrix using the Nussinov algorithm.
 * 
 * 
 * @param sequence The RNA sequence (1-indexed; first character is `$`).
 * @param m The matrix to fill, representing minimum free energy values.
 * @param theta The minimum loop length constraint (default: global `theta`).
 * @param pair_energy The energy contribution of forming a base pair.
 */
void FillMatrix(const std::string& sequence, Matrix2D& m, int theta, float pair_energy){
    int len_seq = sequence.size()-1; // RNA sequence is 1-indexed

    // Fill the matrix using dynamic programming
    for (int i = len_seq; i >= 1; i--) {
        for (int j = i + theta + 1; j <= len_seq; j++) {
            // Case A: Position i remains unpaired
            m(i, j) = m(i + 1, j);

            // Case B: Positions i and j form a base pair
            if (can_pair(sequence[i], sequence[j])) {
                m(i, j) = std::min(m(i, j), m(i + 1, j - 1) + pair_energy);
            }

            // Case C: Position i pairs with another base k (i < k < j)
            for (int k = i + theta + 1; k < j; k++) {
                if (can_pair(sequence[i], sequence[k])) {
                    m(i, j) = std::min(m(i, j), m(i + 1, k - 1) + m(k + 1, j) + pair_energy);
                }
            }
        }
    }
}





/**
 * @brief Backtracks the MFE matrix to determine the RNA secondary structure.
 * 
 * This function extracts the optimal RNA base pairings by tracing back through 
 * the matrix filled by `FillMatrix()`. Returns one optimal structure among all optimals.
 * 
 * @param i Start index of the subsequence.
 * @param j End index of the subsequence.
 * @param m The matrix filled with free energy values.
 * @param sequence The RNA sequence.
 * @param S A vector storing base pairs (modified in-place). The pairs are stored as (i, j) where i and j are the indices of the bases in the sequence.
 * @param theta The minimum loop length constraint.
 */
void Backtrack(int i, int j, const Matrix2D& m, const std::string& sequence, std::vector<std::pair<int, int>>& S, int theta){
    const float epsilon = 1e-7; // Numerical tolerance. Useless here but might be usefull if the energy model is more complex.

    // Base case: Subsequence too short to form a pair
    if (j-i < theta) {
        return;
    }
    else {
        // Case A: Position i remains unpaired
        if (std::abs(m(i, j) - m(i + 1, j)) < epsilon) {
            Backtrack(i + 1, j, m, sequence, S, theta);
            return;
        }

        // Case B: Position i pairs with j
        else if (std::abs(m(i, j) - (m(i + 1, j - 1) + pair_energy)) < epsilon && can_pair(sequence[i], sequence[j])){
            S.push_back(std::make_pair(i,j));
            Backtrack(i + 1, j - 1, m, sequence, S, theta);
            return;
        }

        // Case C: Position i pairs with some k (i < k < j)
        else {
            for (int k = i + theta + 1; k < j; k++) {
                if (std::abs(m(i, j) - (m(i + 1, k - 1) + m(k + 1, j) + pair_energy)) < epsilon && can_pair(sequence[i], sequence[k])) {
                    S.push_back(std::make_pair(i,k));
                    Backtrack(i+1, k - 1, m, sequence, S, theta);
                    Backtrack(k + 1, j, m, sequence, S, theta);
                    break;
                }
            }
        }

    }
}





#ifdef NUSSINOV_TEST
/**
 * @brief Test function for the Nussinov algorithm.
 * 
 * This function runs a test case to validate the matrix filling and backtracking process.
 * 
 */
int main() {
    std::cout << std::endl << "========== Nussinov test ==========" << std::endl;
    std::string a_sequence = "$GCAACUGGCACAAAGGCCUCCUGG"; // Example RNA sequence //TODO : can be changed
    int len_seq = a_sequence.size()-1;

    std::cout << "Input Sequence: " << a_sequence << std::endl;
    std::cout << "Size of the sequence (without $): " << len_seq << std::endl;
    std::cout << "Theta: " << theta << std::endl;
    std::cout << "Pair energy: " << pair_energy << std::endl;

    std::cout << "========== Filling the energy matrix ==========" << std::endl;
    Matrix2D m(len_seq, len_seq);
    FillMatrix(a_sequence, m, theta, pair_energy);
    print_matrix(m, a_sequence);

    std::cout << "========== Backtracking the energy matrix... ==========" << std::endl;
    std::vector<std::pair<int, int>> S;
    Backtrack(1, len_seq, m, a_sequence, S, theta);

    std::cout << "========== Nussinnov results ==========" << std::endl;
    std::cout << "Secondary structure in base pair list representation: ";
    for (const auto& pair : S) {
        std::cout << "(" << pair.first << ", " << pair.second << ") ";
    }
    std::cout << std::endl;
    std::string dot_bracket_structure = bplist2db(S, len_seq);
    std::cout << "Secondary structure in dot-bracket representation: " << dot_bracket_structure << std::endl; 

    std::cout << "========== End Nussinov test ==========" << std::endl << std::endl;

    return 0;
}
#endif //NUSSINOV_TEST

