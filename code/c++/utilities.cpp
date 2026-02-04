/**
 * @file utilities.cpp
 * @brief This file contains the implementation of utility functions and classes used for RNA structure manipulation.
 * 
 * The file includes:
 * - Matrix2D class with methods to manipulate 2D matrices.
 * - Functions for RNA base pairing and structure conversion (DotBracket (db) to BasePair list (bplist), and vice versa).
 * - A function to generate random RNA sequences.
 * - Helper functions for printing matrices and displaying results in a formatted way.
 * 
 * These utilities are used in the strand_soup
 * 
 * 
 * COMPILATION of this file for testing:
 * Go to the utilities.ccp directory where the Makefile should also be and run the following command in the terminal:
 *     make utilities
 * Execution:
 * After compilation, run the program with:
 *     ./utilities.exe
 *
 * Dependencies:
 * - Requires "global_variables.hpp" for global variables. The linking is done by the Makefile
 *
 * 
 * @date 2025-03
 */

#include "utilities.hpp"
#include <random>
#include <stdexcept>
#include <iostream>
#include <string>
#include <utility>
#include <iomanip>
#include <stack>


//==================== Matrix2D class implementation ====================

/**
 * @brief Default constructor for Matrix2D class.
 * 
 * Initializes the matrix with 0 rows and 0 columns.
 */
Matrix2D::Matrix2D() : rows(0), cols(0) {}

/**
 * @brief Parameterized constructor for Matrix2D class.
 * 
 * Initializes the matrix with the given number of rows and columns, setting all elements to 0.
 * 
 * @param rows (int): The number of rows in the matrix.
 * @param cols (int): The number of columns in the matrix.
 */
Matrix2D::Matrix2D(int rows, int cols) : rows(rows), cols(cols) {
    data = std::vector<std::vector<float>>(rows, std::vector<float>(cols, 0));}

/**
 * @brief Constructor for Matrix2D class that initializes it from an existing 2D vector.
 * 
 * Copies the data from the given matrix into the Matrix2D object.
 * 
 * @param matrix (std::vector<std::vector<float>>): The 2D vector to initialize the matrix from.
 */
Matrix2D::Matrix2D(const std::vector<std::vector<float>>& matrix) {
    rows = matrix.size();
    cols = matrix[0].size();
    data = matrix;
}

/**
 * @brief Accessor operator for the Matrix2D class.
 * 
 * Allows accessing matrix elements using (i, j) notation. Indices are 1-based.
 * 
 * @param i (int): The row index (1-based).
 * @param j (int): The column index (1-based).
 * 
 * @return float&: Reference to the element at position (i, j).
 */
float& Matrix2D::operator()(int i, int j) {
    return data[i - 1][j - 1];}

/**
 * @brief Const accessor operator for the Matrix2D class.
 * 
 * Allows accessing matrix elements using (i, j) notation, with read-only access. Indices are 1-based.
 * 
 * @param i (int): The row index (1-based).
 * @param j (int): The column index (1-based).
 * 
 * @return const float&: Constant reference to the element at position (i, j).
 */
const float& Matrix2D::operator()(int i, int j) const {
    return data[i - 1][j - 1];}

/**
 * @brief Get the number of rows in the matrix.
 * 
 * @return int: The number of rows.
 */
int Matrix2D::get_rows() const {
    return rows;
}

/**
 * @brief Get the number of columns in the matrix.
 * 
 * @return int: The number of columns.
 */
int Matrix2D::get_cols() const {
    return cols;
}


//==================== RNA structure and base pairing functions ====================


/**
 * @brief Check if two bases can pair.
 * 
 * This function checks if a given pair of RNA bases (base1, base2) can form a valid base pair.
 * 
 * @param base1 (char): The first RNA base.
 * @param base2 (char): The second RNA base.
 * 
 * @return bool: True if the bases can form a valid pair, false otherwise.
 */
bool can_pair(char base1, char base2){
    for (const auto& pair : possible_pairs){
        if (pair.first == base1 && pair.second == base2){
            return true;
        }
    }
    return false;
};




/**
 * @brief Print the base pairing for a given RNA sequence.
 * 
 * This function prints all the base pairs in the format (i, j), where i and j are the indices of the paired bases.
 * 
 * @param sequence (const std::string&): The RNA sequence.
 * @param S (const std::vector<std::pair<int,int>>&): A vector of pairs representing the base pairs (i, j).
 */
void print_pairing(const std::string& sequence, const std::vector<std::pair<int,int> >& S){
    for (const auto& pair : S){
        int i = pair.first;
        int j = pair.second;
        if (i>=0 && i<int(sequence.size()) && j>=0 && j<int(sequence.size())){
            std::cout << sequence[i] << " " << sequence[j] << std::endl;
        }
        else {
            std::cerr << "print_pairing - Invalid indices: (" << i << ", " << j << ")" << std::endl;

        }
    }
};




/**
 * @brief Convert a list of base pairs into a dot-bracket notation string.
 * 
 * This function converts a list of base pairs into a string representing the RNA secondary structure in dot-bracket notation.
 * 
 * @param S (const std::vector<std::pair<int, int>>&): A vector of base pairs.
 * @param size (int): The length of the RNA sequence.
 * 
 * @return std::string: A string representing the secondary structure in dot-bracket notation.
 */
std::string bplist2db(const std::vector<std::pair<int,int> >& S,int size){
    std::string structure(size,'.');

    for (const auto& pair : S){
        int i = pair.first;
        int j = pair.second;

        if (i>= 1 && i <= size && j >= 1 && j <= size){
            structure[i-1]='(';
            structure[j-1]=')';
        }
        else {
            std::cerr << "bplist2db - Invalid indices: (" << i << ", " << j << ")" << std::endl;
        }
    }
    return structure;
}



/**
 * @brief Convert a dot-bracket notation string into a list of base pairs.
 * 
 * This function parses a dot-bracket notation string and returns the base pairs as a list of pairs (i, j).
 * 
 * @param structure (const std::string&): The RNA secondary structure in dot-bracket notation.
 * 
 * @return std::pair<std::vector<std::pair<int, int>>, int>: A pair containing the list of base pairs and the length of the sequence.
 */
std::pair<std::vector<std::pair<int,int> >,int> db2bplist(const std::string& structure){
    std::vector<std::pair<int,int>> base_pairs;
    std::stack<int> open_brackets;
    int n = structure.size();
    for (int i=0; i<n; i++){
        if (structure[i] == '('){
            open_brackets.push(i+1);}
        else if (structure[i] == ')'){
            if (open_brackets.empty()){
                std::cerr << "db2bplist - Unmatched closing bracket at position " << i+1 << std::endl;
            }
            else {
                int j = open_brackets.top();
                open_brackets.pop();
                base_pairs.push_back(std::make_pair(j,i+1));
            }
        }
    }
    return {base_pairs, n};
};


/**
 * @brief Print a matrix of floats with formatted cell width.
 * 
 * This function prints a matrix with a given cell width for proper alignment. The matrix is displayed alongside an associated sequence.
 * 
 * @param m (const Matrix2D&): The matrix to be printed.
 * @param sequence (const std::string&): The RNA sequence associated with the matrix.
 * @param cellWidth (int): The width of each cell for formatting.
 */
void print_matrix(const Matrix2D& m, const std::string& sequence, int cellWidth){
    int n = m.get_rows();
    std::cout << std::string(cellWidth, ' '); // Leading spaces for alignment
    for (int i = 1; i <= n; ++i) {
        std::cout << std::setw(cellWidth) << sequence[i] << " "; 
    }
    std::cout << std::endl;
    for (int i = 1; i<=n; i++){
        std::cout << std::string(cellWidth, ' ');
        for (int j = 1; j<=n; j++){
            std::cout << std::setw(cellWidth) << m(i,j) << " ";
        }
        std::cout << std::endl;
    }
};



// Function to add a strand to the dictionary only if it is not already present
void add_strand_if_unique(std::unordered_map<int, std::string>& strands, const std::string& strand) {
    // Check if the strand already exists in the dictionary
    for (const auto& pair : strands) {
        if (pair.second == strand) {
            // Strand already exists, do not add it
            return;
        }
    }

    // If the strand is not found, add it with the next available index
    int next_index = strands.size() + 1;
    strands[next_index] = strand;
}




/**
 * @brief Generate a random RNA sequence of a given length.
 * 
 * This function generates a random RNA sequence consisting of the nucleotides 'A', 'C', 'G', and 'U'.
 * 
 * @param length (int): The length of the desired sequence.
 * 
 * @return std::string: The randomly generated RNA sequence.
 */
static thread_local std::mt19937 rng;

void seed_random(uint64_t seed) {
    rng.seed(static_cast<uint32_t>(seed));
}

std::string generate_random_sequence(int length) {
    if (length < 0) throw std::invalid_argument("length must be >= 0");
    static const char bases[4] = {'A','U','C','G'};
    std::uniform_int_distribution<int> dist(0, 3);

    std::string s;
    s.reserve((size_t)length);
    for (int i = 0; i < length; ++i) s.push_back(bases[dist(rng)]);
    return s;
}


/**
 * @brief Generate a triplet repeat sequence.
 * 
 * This function generates a sequence consisting of a specified number of repeats of a given triplet.
 * 
 * @param sequence (std::string): The triplet to be repeated.
 * @param number_of_repeats (int): The number of times the triplet should be repeated.
 * 
 * @return std::string: The generated triplet repeat sequence.
 */
std::string generate_triplet_repeat(std::string sequence, int number_of_repeats) {
    if (number_of_repeats < 0) {
        std::cerr << "Error: Length must be non-negative." << std::endl;
        return "";
    }
    std::string triplet = "$";
    for (int i = 0; i < number_of_repeats; ++i) {
        triplet += sequence;
    }
    return triplet;
}




#ifdef UTILITIES_TEST
int main() {
    std::cout << std::endl << "========== Test of the utilities functions ==========" << std::endl;
    std::string sequence = "$GCUAAAAGC";
    int size = sequence.size()-1;

    std::vector<std::pair<int, int>> S = {{1, 9}, {2, 8}, {3, 7}};

    std::cout << "  Input Sequence: " << sequence << std::endl;
    std::cout << "  Base-Pairs: ";
    for (const auto& pair : S) {
        std::cout << "(" << pair.first << ", " << pair.second << ") ";
    }
    std::cout << std::endl;
    std::cout << "  Size of the sequence (without $): " << size << std::endl;


    std::string structure = bplist2db(S, size);
    std::cout << "  Secondary structure in dot-bracket representation: " << structure << std::endl; 

    auto [base_pairs, length] = db2bplist(structure);
    std::cout << "  Secondary structure in base pair list: ";
    for (const auto& pair : base_pairs) {
        std::cout << "(" << pair.first << ", " << pair.second << ") ";
    }
    std::cout << std::endl;


    std::cout << "========== Test of the function print_matrix ==========" << std::endl;


    Matrix2D matrix({
        {0,0,1,1,1,1,1,2,3},
        {0,0,0,0,1,1,1,2,2},
        {0,0,0,0,1,1,1,1,1},
        {0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0},
    });

    print_matrix(matrix,sequence);

    std::cout << "========== Test of the function generate_random_sequence ==========" << std::endl;
    std::cout << "  Random sequence: " << generate_random_sequence(10) << std::endl;

    std::cout << "========== Test of the function generate_triplet_repeat ==========" << std::endl;
    std::string triplet = generate_triplet_repeat(3, "AAA");
    std::cout << "  Triplet repeat: " << triplet << std::endl;
    triplet = generate_triplet_repeat(3, "CGU");
    std::cout << "  Triplet repeat: " << triplet << std::endl;


    std::cout << "========== End of the utilities test ==========" << std::endl << std::endl;
    return 0;
}
#endif //UTILITIES_TEST