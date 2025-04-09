# RNA_Triplets

---

We aim at providing a C++ implementation of a dynamic programming scheme for the Strand soup Interaction model. This model is found in the following paper : 

**RNA Triplet Repeats: Improved Algorithms for Structure Prediction and Interactions**  
Kimon Boehmer1, Sarah J. Berkemer1,2, Sebastian Will1, Yann Ponty1



---

## Project Overview

This project implements the Nussinov algorithm and extends it to the Strand Soup Interaction model for RNA secondary structure prediction.

___


## Nussinov algorithm

This algorithm aims at predicting the secondary structure of a single RNA strand. It does so by finding the secondary structure that maximizes the number of base pairs in a sequence.

### The problem

Given an RNA sequence S of length n. The goal is to find the structure that maximizes the number of base pairs.
It is represented via a Matrix M[i,j] :
- i and j represent indices of bases in the RNA sequence
- M[i,j] is the maximum number of base pairs in the sequence S[i,j]

### Solving

We use a dynamic programming approach. It is composed of two main functions :

1. **Filling the Minimum Energy Matrix**  
   Computes the minimum energy matrix using dynamic programming.

2. **Backtracking for the Minimum Energy Structure**  
   Extracts the optimal RNA secondary structure by tracing back through the computed matrix.


### Function 1 : Filling the minimum energy matrix

#### **Input**: 
- `ω` – RNA of size n

#### **Output**:
- `m` – Minimum energy matrix m

``` 
    Fonction FillMatrix (ω):
        m ← EmptyMatrix(n × n)
        ## Initialize with 0 all the values of the diagonal up to θ.
        for i ← 1 to n do
            for j ← i to min(i + θ, n) do
                mᵢ,ⱼ ← 0
        for i ← n to 1 do
            for j ← i + θ + 1 to n do
                Case A: Position i left without partner
                    mᵢ,ⱼ ← mᵢ₊₁,ⱼ
                Case B: Positions i and j form a base pair
                    mᵢ,ⱼ ← min(mᵢ,ⱼ, mᵢ₊₁,ⱼ₋₁ + E^ωᵢⱼ)
                Case C: Position i paired to k < j
                    for k ← i + θ + 1 to j − 1 do
                        mᵢ,ⱼ ← min(mᵢ,ⱼ, mᵢ₊₁,ₖ₋₁ + mₖ₊₁,ⱼ + E^ωᵢ,ₖ)
        return m
```

### Function 2: Backtracking for the Minimum Energy Structure

#### **Input**:
- `[i, j]` – Region under consideration  
- `m` – Dynamic programming matrix, previously computed by FillMatrix
- `ω` – RNA sequence of length \( n \)  

#### **Output**:
- $S^*$ – Structure minimizing free energy  

```
Function Backtrack(i, j, m, w):
    if j - i <= theta then 
        return *....* #The empty structure has min energy
    else
        **Case A: Position i left without partner**
        if mᵢ,ⱼ = mᵢ₊₁,ⱼ then
            S*ᵢ ← Backtrack(i+1, j, m, w)
            return • S*ᵢ
        Case B : Positions i and j form a base pair
        if mi,j = mi+1,j−1 + Eω then
            Sᵢ,ⱼ⋆ ← Backtrack(i+1, j−1, m, w)
            return S*ᵢⱼ
        Case C: Position i paired to k < j
        for k ← i + θ + 1 to j − 1 do
            if mᵢⱼ = mᵢ₊₁,ₖ₋₁ + mₖ₊₁,ⱼ + Eωᵢₖ then
                S*₁ ← Backtrack(i+1, k−1, m, w)
                S*₂ ← Backtrack(k+1, j, m, w)
                return (S1⋆)S2⋆
```

---

#### To run the c++ implementation for Nussinov
```
cd code/c++
make nussinov
./nussinov.exe
```

---

# Strand soup algorithm

### The problem

Instead of considering only one strand, here we consider a soup of strands. The goal is to find the secondary structure that minimizes the energy. Here we do so with a simplistic energy model which maximizes the base pair count.


Here is the schematic illustration of the dynamic programming scheme for the Strand Soup Interaction model (base pair-based energy model). The illustration is extracted from the reference paper figure 7.

![schematic illustration](ressources/Strand_soup_diagram.png)


### Key Components

1. **Matrix6D Class**  
   A 6-dimensional matrix used to store energy values for all possible configurations.

2. **MainAuxiliaryMatrix**  
    Fills the 6D matrix with energy values using a dynamic programming approach.

3. **Backtracking Functions**  
   - `square_backtrack`: Handles the square case of the Strand Soup problem.
   - `bubble_backtrack`: Handles the bubble case of the Strand Soup problem.
   - `nussinov_backtrack`: Handles single-strand backtracking.

4. **Utilities**  
   Helper functions for matrix manipulation, RNA sequence generation, and structure visualization.


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
- **`strand_soup.cpp`**: Implements the Strand Soup Interaction model using a 6D dynamic programming matrix.
- **`nussinov.cpp`**: Implements the Nussinov algorithm for single RNA strand secondary structure prediction.
- **`utilities.cpp`**: Contains helper functions for matrix manipulation, RNA sequence generation, and other utilities.
- **`global_variables.hpp`**: Defines global constants (e.g., `theta`, `pair_energy`) and shared variables used across the C++ implementation.
- **`Makefile`**: Provides build instructions for compiling the C++ code. Includes targets for `strand_soup`, `nussinov`, and `clean`.

#### Python Files
The python implementation was just an introduction to better understand the algorithm. It isn't usefull for the strand soup algorithm implemented in c++. The only exception being `Draw_matrices.ipynb`which is used to draw matrices from csv files.
- **`Draw_matrices.ipynb`** : Reads the csv tables produced by the strand_soup algorithm. Produces the heatmaps of the repartition of pair types with respect to the triplet repeats in the soup.
- **`main.ipynb`**: Python implementation of the Nussinov algorithm. Contains the functions FillMatrix, Backtrack and some visualisation via viennaRNA.
- **`utilities.py`**: Helper functions for matrix operations and RNA sequence handling in Python.
- **`global_variables.py`**: Defines global constants for the Python implementation.

---

## Details on the strand soup algorithm




---

## Mathematical Modeling of the Strand Soup Problem

The Strand Soup problem can be formalized as follows:

### **Input:**
- A set of sequences $R = \{r_1, r_2, \dots, r_p\}$, where each $r_i$ is an RNA sequence.
- An integer $m \in \mathbb{N}$, encoded in unary, representing the number of selected sequences.

### **Objective:**
Find a set of sequences $\{t_1, t_2, \dots, t_m\} \subseteq R$ and a secondary structure $S \in \Omega(\{t_1, \dots, t_m\})$ that minimizes the total energy $E$.

### **Total Energy Function:**
$$
\min_{t_1 \in R, \dots, t_m \in R} \, \min_{S \in \Omega(\{t_1, \dots, t_m\})} \, E(\{t_1, \dots, t_m\}, S)
$$

where:
- $\Omega(\{t_1, \dots, t_m\})$ is the set of all possible secondary structures for the sequences $\{t_1, \dots, t_m\}$.
- $E(\{t_1, \dots, t_m\}, S)$ is the energy associated with the structure $S$.

### **Constraints:**
1. The selected sequences $\{t_1, \dots, t_m\}$ must belong to the set $R$.
2. The secondary structure $S$ must follow RNA base-pairing rules (A-U, G-C, we also consider G-U).
3. Base pairs in $S$ must not cross.

--- 

### Implementation and choices
Most of the implementation (at least for the part associated to the algorithm) is done with a **1-based index**. Most of this re-indexing is hidden within access operator of the matrix classes. We consider it to be easier to debug as the math is 1-based.

The strands are 1-based strings. They are represented with a `$`sign in front of the sequence. Example : `$AUCAUCAUC`. Therefore, the real sequence length is sequence.length() - 1  ! 


The set $R$ is implemented using a `std::unordered_map<int, std::string>`. Each unique sequence is associated to an int.

For each of these sequences we compute the nussinov matrix before the strand_soup algorithm using `FillMatrix(sequence, energy_matrix)` with energy_matrix an empty `Matrix2D`. This class is described in `utilities.cpp`. This is usefull as we don't have to recompute the Nussinov algorithm each time inside the strand_soup.

Concerning the `Matrix6D` : 
- It is 1-based. To access the elements, use the operator ().
- The dimensions of the matrix follow the lexicographic order and are therefore : m,s,i,r,j,c
- the size of the first dimension is m-1. To avoid confusion we use m_start for the initial number of strands in the soup and m_size for the size of the first dimension of the matrix. m_size = m_start-1. 
- the dimension of s and r corresponds to the number of different sequences in the set $R$
- i and j have been extended to handle border effects of the strand soup algorithm. Indeed, the square case is when we might call the function out of the sequence (s empty and/or r empty). Therefore we have chosen to add 2 columns one before at 0 and one at the end. The index are also 1-based here. We have to be careful when calling `Matrix6D.get_i_size()`as it returns the true length of range for i without the added borders.
- the dimension of c is 2.

The filling of the matrix is done with `MainAuxiliaryMatrix`. At the begining we had 2 matrices : $M$ and $\bar{M}$ as in the paper. However, as suggested by Sebastian Will, we have decided to only keep the AuxiliaryMatrix $\bar{M}$. The other matrix is replaced by a function `GeneralCaseMinimization`which compares the 4 different cases and takes the min across them. (cf fig7 of the paper).

The backtracking is done with 3 functions : 
- `square_backtrack`
- `bubble_backtrack`
- `nussinov_backtrack` 

Which corresponds to the 3 different cases depicted in the figure 7. 

The output of these backtracks is handled with a class : 
`output_backtrack` 

This class contains two elements : 
- An ordered list of the type of strands involved in the secondary structure which is a sub-sequence of $\{t_1, \dots, t_m\}$.
- a list of all the pairings of this sub-structure : (x,i,y,j) with x the index of the first strand within the ordered list, i its base, y the second strand, j its base. 
  
There is a method shift that allows to shift the indices of the strands before merging them to the upper-problem. The idea is that when we merge the sub-problem, we add the oredered list somewhere in the ordered list of the upper case. Therefore the values x and y of the pairs involved in the subproblem must be shifted with respect to where the sub-rpbolem is within the upper problem. 


---

## Some results

### Example on the strand soup algorithm :

#### Parameter

- **Number of strands** : `m_start = 4`
- **Length of sequences** : `sequence_length = 10*3 = 30`
- **Theta** : `theta = 1`
- **Pair energy** : `pair_energy = -1`


#### Generated strands:
- Index 1 : $GUUGUUGUUGUUGUUGUUGUUGUUGUUGUU
- Index 2 : $CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG
- Index 3 : $ACGACGACGACGACGACGACGACGACGACG

$(GUU)^{10}$, $(CAG)^{10}$, $(ACG)^{10}$

#### Output :

Starting point: M(2,1,1,1,30,1) = -59

Ordered list of sequences: 1 2 2 1 

Pairs:  (1,30,2,2) (1,29,2,3) (1,28,2,4) (1,27,2,5) (1,26,2,6) (1,25,2,7) (1,24,2,8) (1,23,2,9) (1,22,2,10) (1,21,2,11) (1,20,2,12) (1,19,2,13) (1,18,2,14) (1,17,2,15) (1,16,2,16) (1,15,2,17) (1,14,2,18) (1,13,2,19) (1,12,2,20) (1,11,2,21) (1,10,2,22) (1,9,2,23) (1,8,2,24) (1,7,2,25) (1,6,2,26) (1,5,2,27) (1,4,2,28) (1,3,2,29) (1,2,2,30) (1,1,3,1) (3,30,4,2) (3,29,4,3) (3,28,4,4) (3,27,4,5) (3,26,4,6) (3,25,4,7) (3,24,4,8) (3,23,4,9) (3,22,4,10) (3,21,4,11) (3,20,4,12) (3,19,4,13) (3,18,4,14) (3,17,4,15) (3,16,4,16) (3,15,4,17) (3,14,4,18) (3,13,4,19) (3,12,4,20) (3,11,4,21) (3,10,4,22) (3,9,4,23) (3,8,4,24) (3,7,4,25) (3,6,4,26) (3,5,4,27) (3,4,4,28) (3,3,4,29) (3,2,4,30) 

---
Starting point: M(2,1,1,2,30,1) = -59

Ordered list of sequences: 1 1 2 2 

Pairs:  (1,30,2,1) (2,30,3,2) (2,29,3,3) (2,28,3,4) (2,27,3,5) (2,26,3,6) (2,25,3,7) (2,24,3,8) (2,23,3,9) (2,22,3,10) (2,21,3,11) (2,20,3,12) (2,19,3,13) (2,18,3,14) (2,17,3,15) (2,16,3,16) (2,15,3,17) (2,14,3,18) (2,13,3,19) (2,12,3,20) (2,11,3,21) (2,10,3,22) (2,9,3,23) (2,8,3,24) (2,7,3,25) (2,6,3,26) (2,5,3,27) (2,4,3,28) (2,3,3,29) (3,30,4,1) (2,2,4,2) (1,29,4,3) (1,28,4,4) (1,27,4,5) (1,26,4,6) (1,25,4,7) (1,24,4,8) (1,23,4,9) (1,22,4,10) (1,21,4,11) (1,20,4,12) (1,19,4,13) (1,18,4,14) (1,17,4,15) (1,16,4,16) (1,15,4,17) (1,14,4,18) (1,13,4,19) (1,12,4,20) (1,11,4,21) (1,10,4,22) (1,9,4,23) (1,8,4,24) (1,7,4,25) (1,6,4,26) (1,5,4,27) (1,4,4,28) (1,3,4,29) (1,2,4,30) 

---
Starting point: M(2,1,1,3,30,1) = -59

Ordered list of sequences: 1 3 1 3 

Pairs:  (1,2,1,4) (1,1,1,5) (2,30,3,1) (3,3,3,5) (3,2,3,6) (2,29,3,7) (2,28,3,8) (2,27,3,9) (2,26,3,10) (2,25,3,11) (2,24,3,12) (2,23,3,13) (2,22,3,14) (2,21,3,15) (2,20,3,16) (2,19,3,17) (2,18,3,18) (2,17,3,19) (2,16,3,20) (2,15,3,21) (2,14,3,22) (2,13,3,23) (2,12,3,24) (2,11,3,25) (2,10,3,26) (2,9,3,27) (2,8,3,28) (2,7,3,29) (2,6,3,30) (2,5,4,1) (2,4,4,2) (2,3,4,3) (2,2,4,4) (2,1,4,5) (1,30,4,6) (1,29,4,7) (1,28,4,8) (1,27,4,9) (1,26,4,10) (1,25,4,11) (1,24,4,12) (1,23,4,13) (1,22,4,14) (1,21,4,15) (1,20,4,16) (1,19,4,17) (1,18,4,18) (1,17,4,19) (1,16,4,20) (1,15,4,21) (1,14,4,22) (1,13,4,23) (1,12,4,24) (1,11,4,25) (1,10,4,26) (1,9,4,27) (1,8,4,28) (1,7,4,29) (1,6,4,30)

---
Starting point: M(2,2,1,1,30,1) = -59

Ordered list of sequences: 2 2 1 1 

Pairs:  (1,30,2,1) (2,30,3,2) (2,29,3,3) (2,28,3,4) (2,27,3,5) (2,26,3,6) (2,25,3,7) (2,24,3,8) (2,23,3,9) (2,22,3,10) (2,21,3,11) (2,20,3,12) (2,19,3,13) (2,18,3,14) (2,17,3,15) (2,16,3,16) (2,15,3,17) (2,14,3,18) (2,13,3,19) (2,12,3,20) (2,11,3,21) (2,10,3,22) (2,9,3,23) (2,8,3,24) (2,7,3,25) (2,6,3,26) (2,5,3,27) (2,4,3,28) (2,3,3,29) (3,30,4,1) (2,2,4,2) (1,29,4,3) (1,28,4,4) (1,27,4,5) (1,26,4,6) (1,25,4,7) (1,24,4,8) (1,23,4,9) (1,22,4,10) (1,21,4,11) (1,20,4,12) (1,19,4,13) (1,18,4,14) (1,17,4,15) (1,16,4,16) (1,15,4,17) (1,14,4,18) (1,13,4,19) (1,12,4,20) (1,11,4,21) (1,10,4,22) (1,9,4,23) (1,8,4,24) (1,7,4,25) (1,6,4,26) (1,5,4,27) (1,4,4,28) (1,3,4,29) (1,2,4,30) 

---
Starting point: M(2,2,1,2,30,1) = -59

Ordered list of sequences: 2 1 1 2 

Pairs:  (1,30,2,2) (1,29,2,3) (1,28,2,4) (1,27,2,5) (1,26,2,6) (1,25,2,7) (1,24,2,8) (1,23,2,9) (1,22,2,10) (1,21,2,11) (1,20,2,12) (1,19,2,13) (1,18,2,14) (1,17,2,15) (1,16,2,16) (1,15,2,17) (1,14,2,18) (1,13,2,19) (1,12,2,20) (1,11,2,21) (1,10,2,22) (1,9,2,23) (1,8,2,24) (1,7,2,25) (1,6,2,26) (1,5,2,27) (1,4,2,28) (1,3,2,29) (1,2,2,30) (1,1,3,1) (3,30,4,2) (3,29,4,3) (3,28,4,4) (3,27,4,5) (3,26,4,6) (3,25,4,7) (3,24,4,8) (3,23,4,9) (3,22,4,10) (3,21,4,11) (3,20,4,12) (3,19,4,13) (3,18,4,14) (3,17,4,15) (3,16,4,16) (3,15,4,17) (3,14,4,18) (3,13,4,19) (3,12,4,20) (3,11,4,21) (3,10,4,22) (3,9,4,23) (3,8,4,24) (3,7,4,25) (3,6,4,26) (3,5,4,27) (3,4,4,28) (3,3,4,29) (3,2,4,30) 

---
Starting point: M(2,2,1,3,30,1) = -59

Ordered list of sequences: 2 1 1 3 

Pairs:  (1,30,2,2) (1,29,2,3) (1,28,2,4) (1,27,2,5) (1,26,2,6) (1,25,2,7) (1,24,2,8) (1,23,2,9) (1,22,2,10) (1,21,2,11) (1,20,2,12) (1,19,2,13) (1,18,2,14) (1,17,2,15) (1,16,2,16) (1,15,2,17) (1,14,2,18) (1,13,2,19) (1,12,2,20) (1,11,2,21) (1,10,2,22) (1,9,2,23) (1,8,2,24) (1,7,2,25) (1,6,2,26) (1,5,2,27) (1,4,2,28) (1,3,2,29) (1,2,2,30) (3,30,4,1) (3,28,4,2) (3,27,4,3) (3,26,4,4) (3,25,4,5) (3,24,4,6) (3,23,4,7) (3,22,4,8) (3,21,4,9) (3,20,4,10) (3,19,4,11) (3,18,4,12) (3,17,4,13) (3,16,4,14) (3,15,4,15) (3,14,4,16) (3,13,4,17) (3,12,4,18) (3,11,4,19) (3,10,4,20) (3,9,4,21) (3,8,4,22) (3,7,4,23) (3,6,4,24) (3,5,4,25) (3,4,4,26) (3,3,4,27) (3,2,4,28) (3,1,4,29) (1,1,4,30) 

---
Starting point: M(2,3,1,1,30,1) = -59

Ordered list of sequences: 3 2 1 1 

Pairs:  (1,30,2,1) (2,30,3,2) (2,29,3,3) (2,28,3,4) (2,27,3,5) (2,26,3,6) (2,25,3,7) (2,24,3,8) (2,23,3,9) (2,22,3,10) (2,21,3,11) (2,20,3,12) (2,19,3,13) (2,18,3,14) (2,17,3,15) (2,16,3,16) (2,15,3,17) (2,14,3,18) (2,13,3,19) (2,12,3,20) (2,11,3,21) (2,10,3,22) (2,9,3,23) (2,8,3,24) (2,7,3,25) (2,6,3,26) (2,5,3,27) (2,4,3,28) (2,3,3,29) (2,2,3,30) (1,29,4,1) (1,28,4,2) (1,27,4,3) (1,26,4,4) (1,25,4,5) (1,24,4,6) (1,23,4,7) (1,22,4,8) (1,21,4,9) (1,20,4,10) (1,19,4,11) (1,18,4,12) (1,17,4,13) (1,16,4,14) (1,15,4,15) (1,14,4,16) (1,13,4,17) (1,12,4,18) (1,11,4,19) (1,10,4,20) (1,9,4,21) (1,8,4,22) (1,7,4,23) (1,6,4,24) (1,5,4,25) (1,4,4,26) (1,3,4,27) (1,2,4,28) (1,1,4,29) 

---
Starting point: M(2,3,1,2,30,1) = -59

Ordered list of sequences: 3 1 1 2 

Pairs:  (1,29,2,1) (1,28,2,2) (1,27,2,3) (1,26,2,4) (1,25,2,5) (1,24,2,6) (1,23,2,7) (1,22,2,8) (1,21,2,9) (1,20,2,10) (1,19,2,11) (1,18,2,12) (1,17,2,13) (1,16,2,14) (1,15,2,15) (1,14,2,16) (1,13,2,17) (1,12,2,18) (1,11,2,19) (1,10,2,20) (1,9,2,21) (1,8,2,22) (1,7,2,23) (1,6,2,24) (1,5,2,25) (1,4,2,26) (1,3,2,27) (1,2,2,28) (1,1,2,29) (2,30,3,1) (3,30,4,2) (3,29,4,3) (3,28,4,4) (3,27,4,5) (3,26,4,6) (3,25,4,7) (3,24,4,8) (3,23,4,9) (3,22,4,10) (3,21,4,11) (3,20,4,12) (3,19,4,13) (3,18,4,14) (3,17,4,15) (3,16,4,16) (3,15,4,17) (3,14,4,18) (3,13,4,19) (3,12,4,20) (3,11,4,21) (3,10,4,22) (3,9,4,23) (3,8,4,24) (3,7,4,25) (3,6,4,26) (3,5,4,27) (3,4,4,28) (3,3,4,29) (3,2,4,30) 

---
Starting point: M(2,3,1,3,30,1) = -59

Ordered list of sequences: 3 1 1 3 

Pairs:  (1,29,2,1) (1,28,2,2) (1,27,2,3) (1,26,2,4) (1,25,2,5) (1,24,2,6) (1,23,2,7) (1,22,2,8) (1,21,2,9) (1,20,2,10) (1,19,2,11) (1,18,2,12) (1,17,2,13) (1,16,2,14) (1,15,2,15) (1,14,2,16) (1,13,2,17) (1,12,2,18) (1,11,2,19) (1,10,2,20) (1,9,2,21) (1,8,2,22) (1,7,2,23) (1,6,2,24) (1,5,2,25) (1,4,2,26) (1,3,2,27) (1,2,2,28) (1,1,2,29) (3,30,4,1) (3,28,4,2) (3,27,4,3) (3,26,4,4) (3,25,4,5) (3,24,4,6) (3,23,4,7) (3,22,4,8) (3,21,4,9) (3,20,4,10) (3,19,4,11) (3,18,4,12) (3,17,4,13) (3,16,4,14) (3,15,4,15) (3,14,4,16) (3,13,4,17) (3,12,4,18) (3,11,4,19) (3,10,4,20) (3,9,4,21) (3,8,4,22) (3,7,4,23) (3,6,4,24) (3,5,4,25) (3,4,4,26) (3,3,4,27) (3,2,4,28) (3,1,4,29) (2,30,4,30) 


---

Heatmaps for the repartition of base pairs :


n = 9*3 and m = 3
Computation completed in 4844.99 seconds.























