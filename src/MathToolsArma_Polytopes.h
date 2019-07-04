// MathTools_Polytopes.h

#ifndef MATH_TOOLS_POLYTOPES_H
#define MATH_TOOLS_POLYTOPES_H

#include<armadillo>

namespace MathTools {

// -----------------------------------------------
// VERTEX ENUMERATION
// -----------------------------------------------

// Returns a matrix with rows corresponding to the vertices of a polytope 
// Polytope: xl <= x  <= xu,  
//           gl <= Cx <= gu
arma::mat vertexEnumeration(arma::vec xl, arma::vec xu, arma::vec gl, 
                            arma::vec gu, arma::mat C);

// Polytope: xl <= x  <= xu,  
//                 Cx <= gu
arma::mat vertexEnumeration(arma::vec xl, arma::vec xu, arma::vec gu, 
                            arma::mat C);

// Polytope: xl <= x  <= xu
arma::mat vertexEnumeration(arma::vec xl, arma::vec xu);

// An interface to cddlib written by Komei Fukuda for vertex enumeration.
// The algorithm implements the Double Description Method of Motzkin et al. 
// https://www.inf.ethz.ch/personal/fukudak/cdd_home/
// The polytope must be specified in H-format:    b - Ax <= 0
// where arma::mat [bmAx] has the coeffs of the concatanated matrix [b,-A]
arma::mat runCddlib(int numConstraints, int numVars, arma::mat &bmAx);

// -----------------------------------------------
// REJECTION SAMPLING
// -----------------------------------------------

// Returns samples from a polytope generated using the rejection sampling method
// Polytope: xl <= x  <= xu,  
//           gl <= Cx <= gu
// Inputs: 
// - numSamplesMax, an upper bound on the number of samples to be returned
// - xl, xu, are vectors of size N that define the lower and upper bounds of x
// - gl, gu, are vectors of size M define the lower and upper bounds of Cx
// - C is a matrix of size (M x N) defining the "halfspaces"/constraints      
// Outputs:
// - a matrix with rows corresponding to samples inside the polytope
arma::mat uniformRejectionSample(int numSamplesMax, arma::vec xl, 
								                 arma::vec xu, arma::vec gl, arma::vec gu,  
                                 arma::mat C);

// Returns the number of samples per dimension that will generate a sample size
// that satsifies numSamplesMaxw with numDims number of dimensions
int uniformSamplesPerDim(int numSamplesMax, int numDims);

// -----------------------------------------------
// BARYCENTRIC SAMPLING
// -----------------------------------------------

// Returns samples from a polytope generated using a weighted barycentric
// coordinate method
// Polytope: xl <= x  <= xu,  
//           gl <= Cx <= gu
arma::mat barycentricSample(int numSamplesMax, arma::vec xl, 
                            arma::vec xu, arma::vec gl, arma::vec gu, 
                            arma::mat A);
// Polytope: xl <= x  <= xu,  
//                Cx  <= gu
arma::mat barycentricSample(int numSamplesMax, arma::vec xl,  arma::vec xu, 
                            arma::vec gu, arma::mat A);

// Converts a set of vertices defining a polytope in n-dimensional space
// into barycentric samples with an upper bound of numSamplesMax on the 
// number of samples
arma::mat verticesToBarycentricSamples(int &numSamplesMax, arma::mat &vertices,
                                       int &n);

// Returns the number of samples per vertex that will generate a sample size
// that satsifies numSamplesMax with numVertices number of vertices
int barycentricSamplesPerVertex(int numSamplesMax, int numVertices);


} // namespace

#endif
