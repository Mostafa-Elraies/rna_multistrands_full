/**
 * @file utilities.hpp
 * @brief Utility functions and classes for RNA secondary structure prediction.
 * 
 * This file defines helper functions and a `Matrix2D` class used in the Nussinov algorithm 
 * for RNA folding. The matrix class allows for 1-based indexing, unlike the standard 
 * 0-based indexing in C++. This choice is motivated by the 1-based syntax used in the paper of the strand soup algorithm.
 * 
 * Dependencies:
 * - Requires `global_variables.hpp` for shared constants.
 * 
 * @date 2025-03
 */

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "global_variables.hpp"

/**
 * @class Matrix2D
 * @brief A 2D matrix class with 1-based indexing.
 * 
 * @details This class is designed to store a 2D matrix of floats while allowing 
 *          access using a 1-based index. It provides methods for accessing and 
 *          modifying elements, as well as retrieving matrix dimensions.
 */
class Matrix2D{
    public:
    /**
     * @brief Default constructor (creates an empty matrix).
     */
    Matrix2D();

    /**
     * @brief Constructs a matrix with the given number of rows and columns.
     * @param rows Number of rows.
     * @param cols Number of columns.
     */
    Matrix2D(int rows, int cols);

    /**
     * @brief Constructs a matrix from an existing 2D vector.
     * @param matrix A 2D vector of floats.
     */
    Matrix2D(const std::vector<std::vector<float>>& matrix);

    /**
     * @brief Accesses an element of the matrix (modifiable).
     * @param i Row index (1-based).
     * @param j Column index (1-based).
     * @return Reference to the element at (i, j).
     */
    float& operator()(int i, int j);

    /**
     * @brief Accesses an element of the matrix (read-only).
     * @param i Row index (1-based).
     * @param j Column index (1-based).
     * @return Const reference to the element at (i, j).
     */
    const float& operator()(int i, int j) const;

    /**
     * @brief Returns the number of rows in the matrix.
     * @return Number of rows.
     */
    int get_rows() const;

    /**
     * @brief Returns the number of columns in the matrix.
     * @return Number of columns.
     */
    int get_cols() const;

    private:
    int rows, cols; ///< Number of rows in the matrix
    std::vector<std::vector<float>> data; ///< Matrix data storage.
};


/**
 * @brief Checks whether two RNA bases can form a valid pair.
 * @param base1 First RNA base.
 * @param base2 Second RNA base.
 * @return True if the bases form a valid pair, false otherwise.
 */
bool can_pair(char base1, char base2);


/**
 * @brief Prints the base pairings of an RNA sequence.
 * @param sequence The RNA sequence.
 * @param S A vector of paired indices.
 */
void print_pairing(const std::string& sequence, const std::vector<std::pair<int,int> >& S);


/**
 * @brief Converts a base pair list into a dot-bracket notation.
 * @param S A vector of paired indices.
 * @param size Length of the RNA sequence.
 * @return Dot-bracket representation of the secondary structure.
 */
std::string bplist2db(const std::vector<std::pair<int,int> >& S,int size);


/**
 * @brief Converts a dot-bracket notation into a base pair list.
 * @param structure Dot-bracket notation of the RNA secondary structure.
 * @return A pair consisting of a vector of base pairs and the structure length.
 */
std::pair<std::vector<std::pair<int,int> >,int> db2bplist(const std::string& structure);


/**
 * @brief Prints a formatted matrix with an associated sequence.
 * @param m The matrix to print.
 * @param sequence The RNA sequence.
 * @param cellWidth The width of each cell in the output (default: 3). This controls the padding for better alignment.
 */
void print_matrix(const Matrix2D& m, const std::string& sequence, int cellWidth = 3);


/**
 * @brief Adds a strand to the strands map if it is not already present.
 * 
 * @param strands A map of strand indices to their sequences.
 * @param strand The strand sequence to be added.
 */
void add_strand_if_unique(std::unordered_map<int, std::string>& strands, const std::string& strand);


/**
 * @brief Generates a random RNA sequence of a given length.
 * @param length The desired sequence length.
 * @return A randomly generated RNA sequence.
 */
std::string generate_random_sequence(int length);


/**
 * @brief Generates a Triplet Repeat.
 * @param sequence The sequence to be repeated.
 * @param number_of_repeats The number of times to repeat the sequence.
 * 
 * @return The triplet repeat
 */
std::string generate_triplet_repeat(std::string sequence, int number_of_repeats);



#endif //UTILITIES_HPP