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
 * @param M The energy matrix.
 * @param nussinov_matrices A dictionary of Nussinov matrices.
 * @return float The minimum energy.
 */
float GeneralCaseMinimization (int m, int s, int i, int r, int j, int c, std::unordered_map<int, std::string> strands, Matrix6D& M,  std::unordered_map<std::string, Matrix2D> nussinov_matrices){
    float min_value = inf_energy;

    // Case 1 : i is left unpaired
    min_value = M(m,s,i+1,r,j,c);

    // Case 2 : i is paired to k of strand s
    for (int k=i+theta+1; k<=int(strands.at(s).length())-1; k++){
        if (can_pair(strands.at(s)[i],strands.at(s)[k])){
            float value = pair_energy + nussinov_matrices[strands.at(s)](i+1,k-1) + M(m,s,k+1,r,j,c);
            if (value < min_value){
                min_value = value;
                // std::cout << "CASE2 with k = " << k << std::endl;
            }
        }
    }

    // Case 3 : i is paired to k of new strand t
    if (m>=1){
        for (int t=1; t<=M.get_s_size(); t++){
            for (int m1=0; m1<m; m1++){
                int m2 = m - m1 - 1;
                for (int k=1; k<=int(strands.at(t).length())-1; k++){
                    if (can_pair(strands.at(s)[i],strands.at(t)[k])){
                        float value = pair_energy + M(m1,s,i+1,t,k-1,0) + M(m2,t,k+1,r,j,c);
                        if (value < min_value){
                            // std::cout << "CASE3 with k = " << k << std::endl;
                            min_value = value;
                        }
                    }
                }
            }
        }   
    }

    // Case 4 : i is paired to k of strand r
    for (int k=1; k<=j; k++){
        if (can_pair(strands.at(s)[i],strands.at(r)[k])){
            if (k == int(strands.at(r).length()-1)){
                float value = pair_energy + M(m,s,i+1,r,k-1,0);
                if (value < min_value){
                    // std::cout << "CASE4 with k = " << k << std::endl;
                    min_value = value;
                }}
            else {
                float value = pair_energy + M(m,s,i+1,r,k-1,0) + nussinov_matrices[strands.at(r)](k+1,j);
                if (value < min_value){min_value = value;
                    // std::cout << "CASE4 with k = " << k << std::endl;
                }}
        }
    }
    // if (min_value == M(m,s,i+1,r,j,c)){std::cout << "Case 1" << std::endl;}
    return min_value;
}





/**
 * @brief Computes the 6D matrix for the RNA strand soup problem
 * 
 * @param strands The dictionary of strands (index, sequence).
 * @param M The main matrix to be filled.
 * @param nussinov_matrices A dictionary of the nussinov matrices.
 */
void MainAuxiliaryMatrix(std::unordered_map<int, std::string> strands, Matrix6D& M, std::unordered_map<std::string, Matrix2D> nussinov_matrices){
    // std::cout << "MainAuxiliaryMatrix" << std::endl;
    // ================== Lexicographic Order (m,s,-i,r,j,c) ================== //
    for (int m=0; m<M.get_m_size(); m++){
        for (int s=1; s<= M.get_s_size(); s++){
            for (int i= M.get_i_size()+1; i>=0; i--){ //! we have added one column before and after for border effect
                for (int r=1; r<=M.get_r_size(); r++){
                    for (int j=0; j<=M.get_j_size()+1; j++){ //! we have added one column before and after for border effect
                        for (int c=0; c<=1; c++){
                            // std::cout << "M(" << m << "," << s << "," << i << "," << r << "," << j << "," << c << ")" << std::endl;
                            if (i > int(strands.at(s).length()-1)){ //! this is the case where the strand s is empty (border effect). Be careful the strands.length must be decreased by one because of the $ sign at the beginning of the strands!
                                // std::cout << "CASE s empty" << std::endl;
                                if(c==1){
                                    M(m,s,i,r,j,c) = inf_energy;}
                                else{
                                    if (m==0){
                                        if (j == 0) {
                                            M(m,s,i,r,j,c) = 0;} //additional condition in case both are empty
                                        else {
                                            // std::cout << "Nussinov  " << strands.at(r) << " " << 0 << " " << j;
                                            M(m,s,i,r,j,c) = nussinov_matrices[strands.at(r)](1,j);}
                                    }
                                    else{// we select the min M[m-1][t][0][r][j][1] over all t
                                        float min_value = inf_energy;
                                        for (int t=1; t<=M.get_s_size(); t++){
                                            if (M(m-1,t,1,r,j,1) < min_value){min_value = M(m-1,t,1,r,j,1);}}
                                        M(m,s,i,r,j,c) = min_value;
                                    }
                                }
                            }
                            else{
                                // std::cout << "CASE r empty" << std::endl;
                                if (j<1){//if r is empty
                                    if (c==1){
                                        M(m,s,i,r,j,c) = inf_energy;}
                                    else {
                                        if (m==0){
                                            if (i ==0){
                                                M(m,s,i,r,j,c) = 0;}
                                            else{
                                                M(m,s,i,r,j,c) = nussinov_matrices[strands.at(s)](i,strands.at(s).length()-1);}
                                            }
                                        else{// we select the min M[m-1][s][i][t][size(t)-1][1] over all t
                                            float min_value = inf_energy;
                                            for (int t=1; t<=M.get_s_size(); t++){
                                                if (M(m-1,s,i,t,strands[t].length()-1,1) < min_value){min_value = M(m-1,s,i,t,strands[t].length()-1,1);}} //! again be careful of the length of the strands real length = length -1 (-1 because of the $ sign).
                                            M(m,s,i,r,j,c) = min_value;
                                        }
                                    }
                                }
                                else {
                                    // std::cout << "CASE GENERAL" << std::endl;
                                    M(m,s,i,r,j,c) = GeneralCaseMinimization(m,s,i,r,j,c,strands,M,nussinov_matrices);
                                }
                                
                            }
                            //Displaying the values of the matrix with colors (red for borders, green for acceptable values (c==1), white for the rest)
                            // if (i == 0 || j== 0 || i == M.get_i_size() + 1 || j == M.get_j_size() + 1){
                            //     std::cout << RED << "M(" << m << "," << s << "," << i << "," << r << "," << j << "," << c << ") = " << M(m,s,i,r,j,c) << RESET << std::endl;
                            // }
                            // else{
                            //     if (c ==1){std::cout << GREEN << "M(" << m << "," << s << "," << i << "," << r << "," << j << "," << c << ") = " << M(m,s,i,r,j,c) << RESET << std::endl;}
                            //     else{std::cout << "M(" << m << "," << s << "," << i << "," << r << "," << j << "," << c << ") = " << M(m,s,i,r,j,c) << std::endl;}
                            // }
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
 * This function searches through the filled energy matrix (M) and identifies the minimum energy value. It returns all the coordinates (m,s,i,r,j,c) of the minimum energy, excluding the borders.
 * 
 * @param M The 6D energy matrix that stores the computed energy values for all possible configurations.
 * 
 * @return std::vector<std::vector<int>> A vector containing all the starting points (m,s,i,r,j,c) where the minimum energy is located.
 */
std::vector<std::vector<int>>  Find_all_start_backtrack(Matrix6D& M){
    std::vector<std::vector<int>> starting_points;
    float min_value = inf_energy;

    for (int s=1; s<=M.get_s_size(); s++){
        for (int r=1; r<= M.get_r_size(); r++){
            if (M(M.get_m_size()-1,s,1,r,M.get_j_size(),1) < min_value){
                min_value = M(M.get_m_size()-1,s,1,r,M.get_j_size(),1);
                starting_points.clear();
            }
            if (M(M.get_m_size()-1,s,1,r,M.get_j_size(),1) == min_value){
                starting_points.push_back({M.get_m_size()-1,s,1,r,M.get_j_size(),1});
            }
        }
    }
    if (min_value == inf_energy){
        throw std::runtime_error("No valid starting point found by Find_start_backtrack in the given energy matrix");
    }
    return starting_points;
}






// We need to define the square_backtrack header before the bubble_backtrack function because it is called by it.
output_backtrack square_backtrack(int m, int s, int i, int r, int j, int c, std::unordered_map<int, std::string> strands, Matrix6D& M, std::unordered_map<std::string, Matrix2D> nussinov_matrices);







/**
 * @brief Backtrack the minimum energy matrix of a single strand to find its secondary structure using the Nussinov algorithm (1D case).
 * 
 * This function recursively backtracks through the Nussinov matrix to find the optimal secondary structure of each single RNA strand.
 * It considers the three cases of pairing and returns the secondary structure as an output_backtrack so that it can be handled as the square_backtrack or the bubble_backtrack outputs.
 * The nussinov case corresponds to the spikes on the Fig7 of the paper AlMoB_submittedVersion.pdf (Schematic illustration of the dynamic programming scheme for the Strand Soup Interaction model)
 * 
 * @param s : the index of the strand
 * @param i : the first index of the subsequence 
 * @param j : the last index of the subsequence
 * @param strands : the dictionary of strands (int = index, string = sequence)
 * @param nussinov_matrices : the dictionary of Nussinov matrices (string = sequence, Matrix2D = energy matrix)
 * @param theta : the minimum distance between paired bases 
 * 
 * @return output_backtrack : the secondary structure
 */
output_backtrack nussinov_backtrack(int s, int i, int j,std::unordered_map<int, std::string> strands, std::unordered_map<std::string, Matrix2D> nussinov_matrices){
    // The cases are ordered according to the Fig7 of the paper AlMoB_submittedVersion.pdf (Schematic illustration of the dynamic programming scheme for the Strand Soup Interaction model)

    // std::cout << "nussinov_backtrack on " << "s="  << s << " i=" << i << " j=" << j << std::endl;

    if (j-i < theta){ // Base case: return an empty structure if subsequence length is below the threshold.
        return output_backtrack();
    }
    else {
        // Case A: Position i remains unpaired
        if (nussinov_matrices[strands.at(s)](i,j) == nussinov_matrices[strands.at(s)](i+1,j)){
            return nussinov_backtrack(s,i+1,j,strands,nussinov_matrices);
        }

        // Case B: Position i pairs with j
        else if (nussinov_matrices[strands.at(s)](i,j) == nussinov_matrices[strands.at(s)](i+1,j-1) + pair_energy && can_pair(strands.at(s)[i],strands.at(s)[j])){
            output_backtrack output = nussinov_backtrack(s,i+1,j-1,strands,nussinov_matrices);
            output.add_pair(1,i,1,j);
            return output;
        }

        // Case C: Position i pairs with some k (i < k < j)
        else {
            for (int k=i+theta+1; k<j; k++){
                if (nussinov_matrices[strands.at(s)](i,j) == nussinov_matrices[strands.at(s)](i+1,k-1) + nussinov_matrices[strands.at(s)](k+1,j) + pair_energy && can_pair(strands.at(s)[i],strands.at(s)[k])){
                    output_backtrack output1 = nussinov_backtrack(s,i+1,k-1,strands,nussinov_matrices);
                    output_backtrack output2 = nussinov_backtrack(s,k+1,j,strands,nussinov_matrices);
                    output1.merge(output2);
                    output1.add_pair(1,i,1,k);
                    return output1;
                }
            }
        }
    }
    std::cout << "ERROR : Unexpected behaviour in nussinov_backtrack" << std::endl; //The cases A,B and C should cover all the possibilities. The function should never reach this point.
    return output_backtrack();
}






/**
 * @brief Backtracks the bubble case.
 * 
 * This function recursively backtracks through the 6D matrix to find the optimal secondary structure for the bubble case. It considers the 4 cases and makes recursive calls via square_backtrack and nussinov_backtrack functions.
 * The bubble case corresponds to the non empty case of the strand soup problem described on the left on the Fig7 of the paper AlMoB_submittedVersion.pdf (Schematic illustration of the dynamic programming scheme for the Strand Soup Interaction model )
 * 
 * @param m : the number of remaining strands in the soup
 * @param s : the index of the starting strand (1-based)
 * @param i : the index of the starting nucleotide of the first strand (1-based)
 * @param r : the index of the ending strand (1-based)
 * @param j : the index of the ending nucleotide of the last second strand (1-based)
 * @param c : the connectivity bit
 * @param strands : the dictionary of strands (int = index, string = sequence)
 * @param M : the energy_matrix
 * @param nussinov_matrices : the dictionary of Nussinov matrices (string = sequence, Matrix2D = energy matrix)
 * 
 * @return output_backtrack : the secondary structure
 */
output_backtrack bubble_backtrack(int m, int s, int i, int r, int j, int c, std::unordered_map<int, std::string> strands, Matrix6D& M, std::unordered_map<std::string, Matrix2D> nussinov_matrices){
    // We go from the top to the bottom of the matrix according to the lexico order (m,s,i,r,j,c)
    // std::cout << "bubble_backtrack on " << "m=" << m << " s="  << s << " i=" << i << " r=" << r << " j=" << j << " c=" << c << std::endl;

    //CASE 1 : i is left unpaired
    if (M(m,s,i,r,j,c) == M(m,s,i+1,r,j,c)){
        // std::cout << "CASE1" << std::endl;
        return square_backtrack(m,s,i+1,r,j,c,strands,M,nussinov_matrices);
    }

    //CASE 2: i is paired to base k of strand s
    for (int k=i+theta+1; k <=int(strands.at(s).length())-1;k++){
        if (can_pair(strands.at(s)[i],strands.at(s)[k])){
            if (M(m,s,i,r,j,c) == pair_energy + nussinov_matrices[strands.at(s)](i+1,k-1) + M(m,s,k+1,r,j,c)){;
                // std::cout << "CASE2 " << "k = " << k << std::endl;
                output_backtrack output1 = nussinov_backtrack(s,i+1,k-1,strands,nussinov_matrices);
                output1.add_pair(1,i,1,k); // we add the pair between i and k of strand s 
                output_backtrack output2 = square_backtrack(m,s,k+1,r,j,c,strands,M,nussinov_matrices);
                output1.merge(output2);
                return output1;
            }
        }
    }

    //CASE 3: i is paired to base k of new strand t
    if (m>=1){
        // std::cout << "CASE3" << std::endl;
        for (int t=1; t<=M.get_s_size(); t++){
            for (int m1=0; m1<m; m1++){
                int m2 = m-m1-1;
                for (int k=1; k<=int(strands.at(t).length())-1;k++){
                    if (can_pair(strands.at(s)[i],strands.at(t)[k])){
                        if (M(m,s,i,r,j,c) == pair_energy + M(m1,s,i+1,t,k-1,0) + M(m2,t,k+1,r,j,c)){
                            output_backtrack output1 = square_backtrack(m1,s,i+1,t,k-1,0,strands,M,nussinov_matrices);
                            output1.add_pair(1,i,1+m1+1,k);
                            output1.add_sequence(t);
                            output_backtrack output2 = square_backtrack(m2,t,k+1,r,j,c,strands,M,nussinov_matrices); // there is a problem here
                            output2.shift(1+m1);
                            output1.merge(output2);
                            return output1;
                        }
                    }
                }
            }
        }
    }

    //CASE 4: i is paired to base k of strand r
    for (int k=1; k<=j;k++){
        if (can_pair(strands.at(s)[i],strands.at(r)[k])){
            if (k == int(strands.at(r).length()-1)){
                if (M(m,s,i,r,j,c) == pair_energy + M(m,s,i+1,r,k-1,0)){
                    // std::cout << "CASE4.1 i= "<< i << " is paired to k= "<< k << std::endl;
                    output_backtrack output = square_backtrack(m,s,i+1,r,k-1,0,strands,M,nussinov_matrices);
                    output.add_pair(1,i,1+m+1,k);
                    return output;
                }
            }
            else{
                if (M(m,s,i,r,j,c) == pair_energy + M(m,s,i+1,r,k-1,0) + nussinov_matrices[strands.at(r)](k+1,j)){
                    // std::cout << "CASE4.2 i= "<< i << " is paired to k= "<< k << std::endl;
                    output_backtrack output1 = square_backtrack(m,s,i+1,r,k-1,0,strands,M,nussinov_matrices);
                    output1.add_pair(1,i,1+m+1,k);
                    output_backtrack output2 = nussinov_backtrack(r,k+1,j,strands,nussinov_matrices);
                    output2.shift(1);
                    output1.merge(output2);
                    return output1;
                }
            }
        }
    }
    std::cout << "ERROR : Unexpected behaviour in bubble_backtrack" << std::endl;
    return output_backtrack();
}






/**
 * @brief Backtracks the square case of the RNA strand soup problem.
 * 
 * This function recursively backtracks through the 6D matrix to find the optimal secondary structure for the square case. It considers the cases where the first strand is empty or the second strand is empty and makes recursive calls via bubble_backtrack, square_backtrack or nussinov_backtrack functions.
 * The square case corresponds to the empty case of the strand soup problem described on the right on the Fig7 of the paper AlMoB_submittedVersion.pdf (Schematic illustration of the dynamic programming scheme for the Strand Soup Interaction model )
 * 
 * @param m : the number of remaining strands in the soup
 * @param s : the index of the starting strand (1-based)
 * @param i : the index of the starting nucleotide of the first strand (1-based)
 * @param r : the index of the ending strand (1-based)
 * @param j : the index of the ending nucleotide of the last second strand (1-based)
 * @param c : the connectivity bit
 * @param strands : the dictionary of strands (int = index, string = sequence)
 * @param M : the energy_matrix
 * @param nussinov_matrices : the dictionary of Nussinov matrices (string = sequence, Matrix2D = energy matrix)
 * 
 * @return output_backtrack : the secondary structure
 */
output_backtrack square_backtrack(int m, int s, int i, int r, int j, int c, std::unordered_map<int, std::string> strands, Matrix6D& M, std::unordered_map<std::string, Matrix2D> nussinov_matrices){
    // The cases are ordered according to the Fig7 of the paper AlMoB_submittedVersion.pdf (Schematic illustration of the dynamic programming scheme for the Strand Soup Interaction model )
    if (i>int(strands.at(s).length()-1)){ //! this is the case where the strand s is empty (border effect). Be careful the strands.length must be decreased by one because of the $ sign at the beginning of the strands!
        if (c==1){
            return output_backtrack(); 
        }
        else {
            if (m==0){
                return nussinov_backtrack(r,1,j,strands,nussinov_matrices);
            }
            else {
                for (int t=1; t<=M.get_s_size();t++){
                    if (M(m,s,i,r,j,c) == M(m-1,t,1,r,j,1)){
                        output_backtrack output = square_backtrack(m-1,t,1,r,j,1,strands,M,nussinov_matrices);
                        output.add_sequence(t);
                        output.shift(1);
                        return output;
                    }
                }
            }
        }
    }
    else {
        if (j<1){ //if r is empty (1-based)
            if (c==1){
                return output_backtrack();
            }
            else {
                if (m==0){
                    return nussinov_backtrack(s,i,strands.at(s).length()-1,strands,nussinov_matrices);
                    }
                else {
                    for (int t=1; t<=M.get_s_size();t++){
                        if (M(m,s,i,r,j,c) == M(m-1,s,i,t,strands.at(t).length()-1,1)){
                            output_backtrack output = square_backtrack(m-1,s,i,t,strands.at(t).length()-1,1,strands,M,nussinov_matrices);
                            output.add_sequence(t);
                            return output;
                        }
                    }
                }
            }
        }
        else{
            return bubble_backtrack(m,s,i,r,j,c,strands,M,nussinov_matrices);
        }
    }
    std::cout << "ERROR : Unexpected behaviour in square_backtrack" << std::endl;
    return output_backtrack();
}





/**
 * @brief Backtracks the 6D matrix to find the optimal secondary structure.
 * 
 * This function determines the starting point in the energy matrix and then calls the `square_backtrack` function to reconstruct the optimal secondary structure for the strand_soup problem.
 * 
 * @param strands A dictionary mapping strand indices (int) to their RNA sequences (std::string). [Input]
 * @param M The 6D energy matrix containing computed energy values. [Input]
 * @param nussinov_matrices A dictionary mapping RNA sequences (std::string) to their precomputed 2D Nussinov energy matrices (Matrix2D). [Input]
 * 
 * @throws std::runtime_error If the starting point cannot be determined.
 * 
 * @return void
 */
std::vector<output_backtrack> full_backtrack(const std::unordered_map<int, std::string> strands, Matrix6D& M, const std::unordered_map<std::string, Matrix2D> nussinov_matrices){

    std::vector<std::vector<int>> starting_points;
    try {
        starting_points = Find_all_start_backtrack(M);
    } catch (const std::runtime_error& e) {
        // Handle the case where no valid starting point is found
        std::cerr << "Warning: " << e.what() << " Setting repartition to zero." << std::endl;

        // Return an empty secondary structure
        return std::vector<output_backtrack>();
    }
    std::vector<output_backtrack> secondary_structures;
    for (auto& starting_point : starting_points){
        std::cout << "Starting point: M(" 
        << starting_point[0] << "," << starting_point[1] << "," 
        << starting_point[2] << "," << starting_point[3] << "," 
        << starting_point[4] << "," << starting_point[5] << ") = " 
        << M(starting_point[0], starting_point[1], starting_point[2], 
            starting_point[3], starting_point[4], starting_point[5]) 
        << std::endl;

        output_backtrack secondary_structure = square_backtrack(
            starting_point[0], starting_point[1], starting_point[2], 
            starting_point[3], starting_point[4], starting_point[5], 
            strands, M, nussinov_matrices);

        secondary_structure.add_sequence(starting_point[3]); // we need to add the starting strand
        secondary_structure.add_sequence_front(starting_point[1]); // and the ending strand
        secondary_structure.print();
        secondary_structures.push_back(secondary_structure);
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
void compute_all_triplet_repartition(int number_of_repeats,int m_start){
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

    std::string internal_filename = "../../results/internal_n=" + std::to_string(number_of_repeats*3) + "_m=" + std::to_string(m_start) + ".csv";
    std::string homogeneous_filename = "../../results/homogeneous_n=" + std::to_string(number_of_repeats*3) + "_m=" + std::to_string(m_start) + ".csv";
    std::string heterogeneous_filename = "../../results/heterogeneous_n=" + std::to_string(number_of_repeats*3) + "_m=" + std::to_string(m_start) + ".csv";


    std::ofstream internal_file(internal_filename);
    std::ofstream homogeneous_file(homogeneous_filename);
    std::ofstream heterogeneous_file(heterogeneous_filename);

    if (!internal_file.is_open() || !homogeneous_file.is_open() || !heterogeneous_file.is_open()) {
        std::cerr << "Error: Unable to open one or more files for writing." << std::endl;
        return;
    }

    internal_file << ",";
    homogeneous_file << ",";
    heterogeneous_file << ",";
    // Write the triplets as column headers without a trailing comma
    for (size_t i = 0; i < triplets.size(); ++i) {
        internal_file << triplets[i];
        homogeneous_file << triplets[i];
        heterogeneous_file << triplets[i];
        // Add a comma only if it's not the last triplet
        if (i != triplets.size() - 1) {
            internal_file << ",";
            homogeneous_file << ",";
            heterogeneous_file << ",";
        }
    }
    internal_file << "\n";
    homogeneous_file << "\n";
    heterogeneous_file << "\n";



    int total_combinations = triplets.size() * triplets.size();
    int current_combination = 0;
    auto start_time = std::chrono::high_resolution_clock::now();

    // STEP 2: Iterate over each triplet
    for (const auto& triplet1 : triplets) {
        std::string strand1 = generate_triplet_repeat(triplet1, number_of_repeats);
        internal_file << triplet1 << ",";
        homogeneous_file << triplet1 << ",";
        heterogeneous_file << triplet1 << ",";

        for (const auto& triplet2 : triplets){
            std::cout << "\n========== Triplet : " <<  triplet1 << " and " << triplet2 << "==========" << std::endl;
            std::string strand2 = generate_triplet_repeat(triplet2, number_of_repeats);
            
            // Create strands map
            std::unordered_map<int, std::string> strands;
            add_strand_if_unique(strands, strand1);
            add_strand_if_unique(strands, strand2);

            // STEP3 : Compute Nussinov matrices
            std::unordered_map<std::string, Matrix2D> nussinov_matrices;
            for (const auto& pair : strands) {
                const auto& seq = pair.second;
                Matrix2D energy_matrix(seq.length()-1, seq.length()-1);
                FillMatrix(seq, energy_matrix);
                nussinov_matrices[seq] = energy_matrix;
                // print_matrix(energy_matrix,seq);
            }

            // STEP 4: Compute the 6D matrix
            int m_size = m_start -1;
            int s_size = strands.size();
            int i_size = strand1.length()-1; // Length of the strand without the '$'
            int r_size = strands.size();
            int j_size = strand2.length()-1;
            int c_size = 2; // Connectivity bit: 0 or 1
            Matrix6D M(m_size, s_size, i_size, r_size, j_size, c_size);
            MainAuxiliaryMatrix(strands, M, nussinov_matrices);

            // STEP 5 : Backtrack the secondary structure
            std::vector<output_backtrack> secondary_structures = full_backtrack(strands, M, nussinov_matrices);

            // STEP 6 : Compute the repartition of pairs
            int internal = 0;
            int homogeneous = 0;
            int heterogeneous = 0;
            for (const auto& secondary_structure : secondary_structures){
                for (const auto& pair : secondary_structure.list_of_pairs){
                    int strand1 = pair[0];
                    int strand2 = pair[2];
                    if (strand1 == strand2){
                        internal++;
                    }
                    else {
                        if (secondary_structure.list_of_sequences[strand1-1] == secondary_structure.list_of_sequences[strand2-1]){
                            homogeneous++;
                        }
                        else {
                            heterogeneous++;
                        }
                    }
                }
            }
            int total = internal + homogeneous + heterogeneous;
            if (total > 0) {
                //STEP 7: Write the results to the respective CSV files
                internal_file << std::fixed << std::setprecision(3) << static_cast<float>(internal) / total;
                homogeneous_file << std::fixed << std::setprecision(3) << static_cast<float>(homogeneous) / total;
                heterogeneous_file << std::fixed << std::setprecision(3) << static_cast<float>(heterogeneous) / total;
                std::cout << "Internal : " << static_cast<float>(internal) / total << std::endl;
                std::cout << "Homogeneous : " << static_cast<float>(homogeneous) / total << std::endl;
                std::cout << "Heterogeneous : " << static_cast<float>(heterogeneous) / total << std::endl;
            } else {
                internal_file << "0";
                homogeneous_file << "0";
                heterogeneous_file << "0";
            }
            // Add a comma only if it's not the last column
            if (triplet2 != triplets.back()) {
                internal_file << ",";
                homogeneous_file << ",";
                heterogeneous_file << ",";
            }

            // Update progress
            current_combination++;
            float progress = (float(current_combination) / total_combinations) * 100;

            // Display progress bar
            std::cout << "\rProgress: [";
            int bar_width = 50;
            int pos = bar_width * progress / 100;
            for (int i = 0; i < bar_width; ++i) {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << std::fixed << std::setprecision(2) << progress << "%";
            std::cout.flush();
        }
        // End the current row in each file
        internal_file << "\n";
        homogeneous_file << "\n";
        heterogeneous_file << "\n";
    }
    // End timing
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end_time - start_time;

    // Display elapsed time
    std::cout << "\nComputation completed in " << elapsed_time.count() << " seconds." << std::endl;
    std::cout << "Results saved to internal.csv, homogeneous.csv, and heterogeneous.csv" << std::endl;
}







//=============================================================================================//








int main() {
    std::cout << "\n==================== Strand Soup ====================" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now(); // To check the time of the computation

    std::cout << "\n========== Setting the parameters ==========" << std::endl;
    int m_start = 2; // Number of sequences in total  //TODO : This parameter can be changed
    int number_of_repeats = 2; // Number of triplet repeats within each strand //TODO : This parameter can be changed
    int sequence_length = number_of_repeats * 3; // length of the sequences
    std::cout << "  Number of strands : m = " << m_start << std::endl;
    std::cout << "  Length of the sequences : " << sequence_length << std::endl;
    std::cout << "  Theta : " << theta << std::endl;
    std::cout << "  Pair energy : " << pair_energy << std::endl;


    std::cout << "\n========== Generation of the strands ==========" << std::endl;
    std::unordered_map<int, std::string> strands;
    strands[1] = generate_triplet_repeat("CCU", number_of_repeats); //TODO : The triplet can be changed
    strands[2] = generate_triplet_repeat("GGA", number_of_repeats); //TODO : The triplet can be changed
    std::cout << "Generated strands:" << std::endl;
    for (const auto& pair : strands){
        std::cout << "  Index " << pair.first << " : " << pair.second << std::endl;
    }


    std::cout << "\n========== Nussinov matrices ==========" << std::endl;
    std::unordered_map<std::string, Matrix2D> nussinov_matrices;
    for (const auto& pair : strands) {
        const auto& seq = pair.second;
        std::cout << "  Energy matrix for sequence : " << seq << std::endl;
        Matrix2D energy_matrix(sequence_length, sequence_length);
        FillMatrix(seq, energy_matrix);
        nussinov_matrices[seq] = energy_matrix;
        print_matrix(energy_matrix,seq);
        std::cout << std::endl;
    }

    std::cout << "========== MainAuxiliaryMatrix ==========" << std::endl;
    int m_size = m_start - 1 ; // Total number of dimensions for m (if we have 2 strands we only have the case m = 0 which is m_size = 1)
    int s_size = strands.size(); // Total number of different strands
    int i_size = sequence_length; //Max length of a strand
    int r_size = strands.size(); // We have the same set of strands as for s
    int j_size = sequence_length; // We have the same length for the strands
    int c_size = 2; // Connectivité : 0 ou 1
    Matrix6D M(m_size, s_size, i_size, r_size, j_size, c_size);

    MainAuxiliaryMatrix(strands, M, nussinov_matrices);


    std::cout << std::endl << "========== full_backtrack ==========" << std::endl;
    full_backtrack(strands, M, nussinov_matrices);



    std::cout << std::endl << "========== compute_all_triplet_repartition ==========" << std::endl;
    // around 30 minutes for m_start = 2 and strand_length = 10*3
    // 40s for m_start = 2 and strand_length = 2*3
    compute_all_triplet_repartition(2,2); //TODO : The Two parameters (number_of_repeats, m_start) can be changed




    std::cout << std::endl << "==================== End of the program ====================" << std::endl;
    // Display elapsed time
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end_time - start_time;
    std::cout << "Computation completed in " << elapsed_time.count() << " seconds." << std::endl << std::endl;

    return 0;
}