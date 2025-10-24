/**
 * @file global_variables.hpp
 * @brief Declares global constants for RNA structure calculations.
 *
 * This file contains global constants for RNA structure prediction, including
 * energy values and valid base pair combinations.
 *
 * @date 2025-03
 */

#ifndef GLOBAL_VARIABLES_HPP
#define GLOBAL_VARIABLES_HPP

#include <vector>
#include <utility>

/**
 * @brief Minimum distance between paired bases. (0 = no constraint)
 */
extern const int theta;


/**
 * @brief Energy of a base pair.
 */
extern const float pair_energy;

extern const float const_a;
extern const float const_b;
extern const float const_c;


/**
 * @brief Energy value representing infinity.
 */
extern const float inf_energy;


/**
 * @brief Valid base pair combinations.
 *
 * This vector contains pairs of characters representing valid base pairs in RNA.
 * Watson-Crick + Wobble pairs are included.
 */
extern const std::vector< std::pair<char, char> > possible_pairs;


#endif // GLOBAL_VARIABLES_HPP