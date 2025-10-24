/**
 * @file nussinov.hpp
 * @brief Declares functions for RNA secondary structure prediction using the Nussinov algorithm.
 * The Nussinov algorithm and code is described in the README file with pseudo code.
 *
 * This file defines the core functions :
 * - for filling the dynamic programming matrix
 * - and performing backtracking to predict RNA secondary structures on single strands.
 *
 * @date 2025-03
 */

#ifndef NUSSINOV_HPP
#define NUSSINOV_HPP

#include <vector>
#include <string>
#include <utility>

#include "utilities.hpp"
#include "global_variables.hpp"



/**
 * @brief Fills the dynamic programming matrix of the Nussinov algorithm
 *
 * This function implements the Nussinov algorithm to compute the optimal
 * RNA secondary structure based on base pairing rules.
 *
 * @param sequence The RNA sequence.
 * @param m The 2D matrix used for dynamic programming.
 * @param theta Minimum loop length (default: global `theta`).
 * @param pair_energy Energy contribution of a base pair (default: global `pair_energy`).
 */
void FillMatrix(const std::string& sequence, Matrix2D& m, int theta= ::theta, float pair_energy = ::pair_energy);



/**
 * @brief Performs backtracking to determine the optimal RNA secondary structure.
 *
 * This function traces back through the computed matrix to extract the optimal
 * RNA base pairings according to the Nussinov algorithm.
 *
 * @param i Starting index.
 * @param j Ending index.
 * @param m The 2D matrix filled by `FillMatrix()`.
 * @param sequence The RNA sequence.
 * @param S A vector of base pair indices forming the predicted structure.
 * @param theta Minimum loop length (default: global `theta`).
 */
void Backtrack(int i, int j, const Matrix2D& m, const std::string& sequence, std::vector<std::pair<int, int>>& S, int theta = ::theta);



#endif // NUSSINOV_HPP