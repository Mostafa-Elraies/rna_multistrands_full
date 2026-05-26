#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include <unordered_map>
#include <vector>
 
 
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
 
 
 
  template<typename T>
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
        data = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>>(
            m_size, std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>(
            s_size, std::vector<std::vector<std::vector<std::vector<T>>>>(
            i_size + 2, std::vector<std::vector<std::vector<T>>>( //! we add one column before and after for border effect
            r_size, std::vector<std::vector<T>>(
            j_size + 2, std::vector<T>( //! we add one column before and after for border effect
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

T& operator()(int m, int s, int i, int r, int j, int c) {
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
    std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>> data;
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
 template<typename T>
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
        data = std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>(
            m_size, std::vector<std::vector<std::vector<std::vector<T>>>>(
            s_size, std::vector<std::vector<std::vector<T>>>(
            i_size + 2, std::vector<std::vector<T>>( //! we add one column before and after for border effect
            r_size, std::vector<T>(
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
    T& operator()(int m, int s, int i, int r, int j) {
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
    std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>> data;
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



