// strand_soup_core.hpp
#pragma once

#include <vector>
#include <unordered_map>
#include <string>

// forward declare evaluator + Vienna params
class RNAEnergyEvaluator;
struct vrna_param_s;
typedef struct vrna_param_s vrna_param_t;

extern const float inf_energy;
extern const int theta;

// ================= Matrix6D =================
class Matrix6D{
public:
    Matrix6D(int m_size, int s_size, int i_size, int r_size, int j_size, int c_size);
    float& operator()(int m, int s, int i, int r, int j, int c);

    int get_m_size() const;
    int get_s_size() const;
    int get_i_size() const;
    int get_r_size() const;
    int get_j_size() const;
    int get_c_size() const;

private:
    int m_size, s_size, i_size, r_size, j_size, c_size;
    std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>>> data;
};

// ================= Matrix5D =================
class Matrix5D{
public:
    Matrix5D(int m_size, int s_size, int i_size, int r_size, int j_size);
    float& operator()(int m, int s, int i, int r, int j);

    int get_m_size() const;
    int get_s_size() const;
    int get_i_size() const;
    int get_r_size() const;
    int get_j_size() const;

private:
    int m_size, s_size, i_size, r_size, j_size;
    std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>> data;
};

// ================= output_backtrack =================
class output_backtrack{
public:
    std::vector<int> list_of_sequences;
    std::vector<std::vector<int>> list_of_pairs;

    void add_sequence(int sequence);
    void add_sequence_front(int sequence);
    void add_pair(int x, int i, int y, int j);
    void shift(int shift);
    void merge(const output_backtrack& sub_output);
    void print() const;
};

// ================= API =================
void MainAuxiliaryMatrix(std::unordered_map<int, std::string> strands,
                         Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
                         RNAEnergyEvaluator& evaluator,
                         vrna_param_t* params);

std::vector<output_backtrack> full_backtrack(
    const std::unordered_map<int, std::string>& strands,
    Matrix5D& C, Matrix5D& M, Matrix5D& M1_multi, Matrix6D& F,
    RNAEnergyEvaluator& evaluator,
    vrna_param_t* params,
    int theta);
