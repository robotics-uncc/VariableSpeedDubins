// MathToolsArma_Combinatornics.h

#ifndef MATH_TOOLS_ARMA_COMBINATRONICS_H
#define MATH_TOOLS_ARMA_COMBINATRONICS_H

#include<armadillo>
#include<MathToolsArma_Combinatronics.h>
#include<MathToolsArma_Polytopes.h>

namespace MathTools {

// return all n! permutations of vector x of dimension n, as rows of a matrix
arma::mat permute(arma::vec x);

// Heap's recursive permutation algorithm
void heapPermute(int n, arma::vec &x, int &j, arma::mat &xperm);

// returns matrix with rows of all permutations of integer pairs that sum to n 
// i.e. {(n,0), (n-1,1) , ... , (0,n)}
arma::mat constantSumPairs(int n);

// returns a matrix with rows that contain all possible sequences of m integers 
// (including zero) that sum to n 
arma::mat constantSumSubsets(int n, int m);

// recursively populates the seqMat matrix with the desired sequences 
// (used in conjuction with the function above)
void constantSumSubsets(int n, int &m, int &curDepth, int &curRow, 
                        arma::mat &seqMat);

} // namespace

#endif
