/**
 * @file global_variables.cpp
 * @brief Implements global constants for RNA structure calculations.
 *
 * This file defines constants for RNA secondary structure prediction.
 * 
 */

#include "global_variables.hpp"

const int theta = 3; // Usually set to 3

const float inf_energy = std::numeric_limits<float>::infinity();

const std::vector< std::pair<char, char> > possible_pairs = {
    {'A','U'}, {'U','A'}, {'G','C'}, {'C','G'}, {'G','U'}, {'U','G'}

};


