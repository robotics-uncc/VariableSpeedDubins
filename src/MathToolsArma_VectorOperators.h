// MathToolsArma_VectorOperators.h
// Provides some basic vector operations.

#ifndef MATH_TOOLS_ARMA_VECTOR_OPERATORS_H
#define MATH_TOOLS_ARMA_VECTOR_OPERATORS_H

// standard C++ headers
#include<vector> // std::vector

// external
#include<armadillo>
//#include<MathTools_VectorOperators.h>

namespace MathTools {

// -----------------------------------------------
// COMPUTING VECTOR PROPERTIES
// -----------------------------------------------

// Returns true of each element of x satisfies x_L <= x <= x_U
bool checkVectorBounds(arma::vec x, arma::vec x_L, arma::vec x_U);

// Finds minimum element(s) in a vector given a tolerance
std::vector<int> minIndicesWithTolerance(arma::vec testVector, 
                                         double tolerance);

// -----------------------------------------------
// GENERATING VECTORS
// -----------------------------------------------

// Generate a sequence of evenly spaced points along an interval defining an 
// angle, with a give direction where (direction = +1) is counter-clockwise, 
// and (direction = -1) is clockwise
arma::colvec polarspace_arma(double angleInitial, double angleFinal, 
							 int npts, int direction);

// -----------------------------------------------
// MODIFYING VECTORS / VECTOR OPERATIONS
// -----------------------------------------------

// Given a vector x, and a scalar t, returns an augmented vector z = [x ; t]^T
arma::colvec augmentVector(double t, arma::colvec x);

// Return cartesian product of m input basis vectors of equal size (n).
// Rows of the (m x n) basisVectorsMatrix are the input vectors
// i.e. 
// basisVectorsMatrix = [   v1 ;   ---> an n-dimensional row vector
//						              v2 ;						   
//						              ..
//						              vm ];
//
// the cartesian product is a (n^m, n) matrix  
// cartesianProductMatrix = v1 x v2 x ... x vm
//
// where each row of cartesianProductMatrix corresponds to an n-dimensional 
// vector that is an element of the cartesian product set.
arma::mat cartesianProduct(arma::mat basisVectorsMatrix);

// replace vector element x(i) with x(j) and vice versa
void swap(arma::vec &x, int i, int j);

// -----------------------------------------------
// MODIFYING MATRICES
// -----------------------------------------------

// addVectorToMatrixColoumnWises
arma::mat addVectorToMatrixColoumnWise(arma::rowvec a, arma::mat B);

// Join coloumn vectors/matrices into a matrix
arma::mat join_horiz(arma::vec a, arma::vec b, arma::vec c);
arma::mat join_horiz(arma::mat a, arma::mat b, arma::mat c);
arma::mat join_horiz(arma::vec a, arma::vec b, arma::vec c, arma::vec d);
arma::mat join_horiz(arma::mat a, arma::mat b, arma::mat c, arma::mat d);

// Join row vectors/matrices into a matrix
arma::mat join_vert(arma::vec a, arma::vec b, arma::vec c);
arma::mat join_vert(arma::mat a, arma::mat b, arma::mat c);
arma::mat join_vert(arma::vec a, arma::vec b, arma::vec c, arma::vec d);
arma::mat join_vert(arma::mat a, arma::mat b, arma::mat c, arma::mat d);


// -----------------------------------------------
// MODIFYING MATRICES
// -----------------------------------------------

// Convert vectors
arma::vec stdToArmaVector(std::vector<double> inputVector);

} // namespace

#endif
