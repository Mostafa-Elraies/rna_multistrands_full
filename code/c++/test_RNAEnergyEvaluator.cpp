#include <iostream>
#include <unordered_map>
#include "RNAEnergyEvaluator.hpp"

int main() {
    std::unordered_map<int, std::string> strands = {
        {1, "GCAUCUAGCUAUGC"},  // test sequence 1
        {2, "GGGAAAUCC"},       // test sequence 2
    };

    RNAEnergyEvaluator evaluator(strands);

    // Test matrices
    std::cout << "C[2][10] (strand 1): " << evaluator.get_C(1)[2][10] << std::endl;
    std::cout << "M[2][10] (strand 1): " << evaluator.get_M(1)[2][10] << std::endl;

    // Test loop energies
    std::cout << "Hairpin Energy (2-10): " << evaluator.hairpin_energy(2, 10, 1) << std::endl;

    std::cout << "Interior loop Energy (2,10,5,7): " 
              << evaluator.interior_loop_energy(2,10,5,7,1) << std::endl;

    std::cout << "Exterior loop Energy (1,9): " 
              << evaluator.exterior_loop_energy(1,9,1) << std::endl;

    std::cout << "Multiloop stem Energy (2,10): " 
              << evaluator.multiloop_stem_energy(2,10,1,false) << std::endl;

    return 0;
}
