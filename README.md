# RNA_Triplets

---

We aim at providing a C++ implementation of a dynamic programming scheme for a thermodynamical Turner energy model of a Strand soup Interaction algorithm. This model is found in the following paper : 

**RNA Triplet Repeats: Improved Algorithms for Structure Prediction and Interactions**  
Kimon Boehmer, Sarah J. Berkemer, Sebastian Will, Yann Ponty



---

## Project Overview

This project implements a loop based algorithm as in Turner thermodynamic model and extends it to the Strand Soup Interaction model for RNA secondary structure prediction.

___


# Strand soup algorithm incorporating turner energy model.

## The problem

 The goal is to find the structure that minimizez the thermodynamical Energy as in Tunerner energy model.
 The goal is to find the secondary structure that minimizes the energy. 


## Key Components

1. **Matrix6D Class**  
A 6-dimensional dynamic programming matrix used to store thermodynamic energy values for all multi-strand configurations of the Strand Soup model.

Each DP state is indexed as:
F(m, s, i, r, j, c)

where:
-  m is the number of additional strands involved (excluding endpoints),
-  s, i denote the left strand and position,
-  r, j denote the right strand and position,
-  c encodes whether the strands are connected or not.

2. **MainAuxiliaryMatrix**  
    The main dynamic programming routine that fills the 6D matrix using a thermodynamically consistent recurrence system.

### Key features:

-   Supports inter-strand base pairing
-   Explicit handling of multiloops, interior loops, hairpins, and external regions
-   Uses ViennaRNA low-level energy functions instead of heuristic constants
-   Ensures consistency with ViennaRNA’s loop decomposition and energy model
3. **RNAEnergyEvaluator**  
   central abstraction that encapsulates all thermodynamic energy calculations.

Responsibilities:

- Builds a ViennaRNA fold_compound for each strand
- Stores ViennaRNA sequence encodings (S1) for every strand
- Computes loop energies using ViennaRNA low-level APIs:
    hairpin loops
    interior loops (intra- and inter-strand)
    multiloop penalties (MLclosing, MLintern, MLbase)
    exterior stem energies
    
- Computes single-strand DP matrices (C, M, M1, F) consistently with ViennaRNA

    This class guarantees that cross-strand interactions use the same pairing types and penalties as single-strand folding.

4. **Backtracking System**  
    A unified backtracking framework that reconstructs optimal structures from the DP matrix.

### Key components:

- backtrack_F_multi_square: handles square configurations

- backtrack_F_multi_bubble: handles bubble configurations

- backtrack_C_multi, backtrack_M_multi, backtrack_M1_multi: reconstruct multiloop and helix contributions

5. **Utilities**
Helper functions and data structures used throughout the project, including:

- RNA sequence generation (random and triplet repeats)

- Base-pair validation

- Matrix helpers with 1-based indexing

- Structure visualization (dot-bracket conversion, pairing output)

- Random seeding and reproducibility helpers
## Compilation and Execution

1. **Compile the Strand Soup Algorithm**  
   Navigate to the `code/c++` directory and run:
   ```
    make strand_soup
    ./strand_soup.exe

   ```


2. **Clean Build Files**  
    To clean up build files, run:
    ```
    make clean
    ```

---

### Description of the file hierarchy

#### C++ Files
- **`strand_soup.cpp`**: Implements the multi-strand Strand Soup interaction model using a 6D thermodynamic dynamic programming matrix and unified backtracking.
- **`RNAEnergyEvaluator.cpp`**: Encapsulates all ViennaRNA-based energy evaluation logic for both single-strand and multi-strand contexts.
- **`utilities.cpp`**: Helper functions for RNA generation, pairing logic, matrix utilities, and structure formatting.
- **`global_variables.hpp`**: Defines global constants (e.g., `theta`) and shared variables used across the C++ implementation.
- **`Makefile`**: Build configuration for compiling the Strand Soup implementation and related utilities.



---

Please refer to the report in the resources for more details on the mathematical representation and the DP algorithm details.
---


### Implementation and choices
Most of the implementation (at least for the part associated to the algorithm) is done with a **1-based index**. Most of this re-indexing is hidden within access operator of the matrix classes. We consider it to be easier to debug as the math is 1-based.

The strands are 1-based strings. They are represented with a `$`sign in front of the sequence. Example : `$AUCAUCAUC`. Therefore, the real sequence length is sequence.length() - 1  ! 


The set $R$ is implemented using a `std::unordered_map<int, std::string>`. Each unique sequence is associated to an int.


# Concerning the `Matrix6D` : 
- It is 1-based. To access the elements, use the operator ().
- The dimensions of the matrix follow the lexicographic order and are therefore : m,s,i,r,j,c
- the size of the first dimension is m-1. To avoid confusion we use m_start for the initial number of strands in the soup and m_size for the size of the first dimension of the matrix. m_size = m_start-1. 
- the dimension of s and r corresponds to the number of different sequences in the set $R$
- i and j have been extended to handle border effects of the strand soup algorithm. Indeed, the square case is when we might call the function out of the sequence (s empty and/or r empty). Therefore we have chosen to add 2 columns one before at 0 and one at the end. The index are also 1-based here. We have to be careful when calling `Matrix6D.get_i_size()`as it returns the true length of range for i without the added borders.
- the dimension of c is 2.


### The backtracking is done with these functions : 
- `backtrack_F_multi_square`
- `bubble_backtrackbubble_backtrack`
- `backtrack_C_multi` 
- `backtrack_M_multi`
- `backtrack_M1_multi`
- `backtrack_Fs`
- `backtrack_Cs`
- `backtrack_M1s`
- `backtrack_Ms`

### The output of these backtracks is handled with a class : 
`output_backtrack` 

### This class contains two elements : 
- An ordered list of the type of strands involved in the secondary structure which is a sub-sequence of $\{t_1, \dots, t_m\}$.
- a list of all the pairings of this sub-structure : (x,i,y,j) with x the index of the first strand within the ordered list, i its base, y the second strand, j its base. 
  
There is a method shift that allows to shift the indices of the strands before merging them to the upper-problem. The idea is that when we merge the sub-problem, we add the oredered list somewhere in the ordered list of the upper case. Therefore the values x and y of the pairs involved in the subproblem must be shifted with respect to where the sub-rpbolem is within the upper problem. 

---

Thanks for reading 

---



















