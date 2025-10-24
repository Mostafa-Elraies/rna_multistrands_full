/**
 * @file global_variables.cpp
 * @brief Implements global constants for RNA structure calculations.
 *
 * This file defines constants for RNA secondary structure prediction.
 * These constants are used for energy calculations and base pairing.
 * 
 */

#include "global_variables.hpp"

const int theta = 3; // Usually set to 3

const float const_a{};
const float const_b{};
const float const_c{};


const float pair_energy = -1; // Can be changed to a dictionary pair/energy. + define a function to extract the energy of a specific pair.

const float inf_energy = std::numeric_limits<float>::infinity();

const std::vector< std::pair<char, char> > possible_pairs = {
    {'A','U'}, {'U','A'}, {'G','C'}, {'C','G'}, {'G','U'}, {'U','G'}

};


