// VSDUtils.h
// Provides utilities for the VariableSpeedDubinsSolver
// Last Modefied: 18-Jan-2016, Artur Wolek

#ifndef VSDUTILS_H
#define VSDUTILS_H

// C++ headers
#include<string>

// external headers
#include<armadillo>

//custom
#include<MathTools.h>
#include<MathToolsArma.h>


namespace VarSpeedDubins {

// ----------------------------------------------------------------------------
// VCSDUtils_pathGeneration.cpp
// ----------------------------------------------------------------------------

// returns the displacement vector for a "S" segment
arma::vec S(double psi, double L);

// returns the displacement vector for a "BCB" segment
arma::vec BCB(double aUnsigned, double bUnsigned, double gUnsigned, double psi0, 
	            double sgnK, double R, double r);

// returns a matrix defining the path for a "S" segment
arma::mat Spath(double psi0, double L, double nomSpacing);

// returns a matrix defining the path for a "BCB" segment
arma::mat BCBpath(double aUnsigned, double bUnsigned, double gUnsigned, 
	                double psi0, double sgnK, double R, double r, 
                  double nomSpacing);

// returns endpoint for a given pathClass and orientation
arma::vec pathEndpoint(arma::vec paramsShort, std::string pathClass, 
	                     std::string pathType, std::string orientation, 
                       double R, double r);
	
// returns history for a given pathClass and orientation
void pathHistory(arma::mat &pathHistory, arma::vec paramsShort, 
	              std::string pathClass, std::string pathType, 
                std::string orientation, double R, double r, double nomSpacing);

// extremalSwitches
arma::mat extremalSwitches(arma::vec paramsShort, std::string pathClass, 
                           std::string pathType, std::string orientation, 
	                         double R, double r);

// ----------------------------------------------------------------------------
// VCSDUtils_startPoint.cpp
// ----------------------------------------------------------------------------
arma::vec guessStartPoint(std::string algorithm, std::string pathClass, 
                          std::string pathType, std::string orientation, 
                          arma::vec xl, arma::vec xu, arma::vec bl, 
                          arma::vec bu, arma::mat A, double R, double r, int n,
                          int m, arma::vec stateFinal, int numSamplesPerVar);


double normErrorToEndpoint(arma::vec xTest, arma::vec xTrue);

arma::mat barycentricSample(int numSamplesPerVar, arma::vec xl, 
                                       arma::vec xu, arma::vec gl, arma::vec gu, 
                                       arma::mat A);
// ----------------------------------------------------------------------------
// VCSDUtils_displacementsToPath.cpp
// ----------------------------------------------------------------------------

// Suppose a state variable (such as x-position) for a particular path is 
// defined by a sequence of three unique curves/line segments. Assume that the 
// displacement along each segment, relative to the start point of that 
// segment, is given by the vectors seg1, seg2 and seg3 where seg1 is the time 
// history of the state variable along the first segment in the sequence, and 
// seg2 and seg3 correspond to the second and third segments. Then we may join
// the three segments to form a continuous state history for the state variable
// of interest by concatonating the vectors and shifting them appropriately.
// The function below returns the continuous state history as described above. 
void displacementsToPath(arma::vec &segJoined, arma::vec &seg1, 
                         arma::vec &seg2);
void displacementsToPath(arma::vec &segJoined, arma::vec &seg1, arma::vec &seg2, 
	                       arma::vec &seg3);

// similarily for a matrix of many states, we can compute the continuous state
// histories using the following functions
void displacementsToPath(arma::mat &segJoined, arma::mat &seg1, 
                         arma::mat &seg2);
void displacementsToPath(arma::mat &segJoined, arma::mat &seg1, arma::mat &seg2, 
	                       arma::mat &seg3);
void displacementsToPath(arma::mat &segJoined, arma::mat &seg1, arma::mat &seg2, 
	                       arma::mat &seg3, arma::mat &seg4);

// ----------------------------------------------------------------------------
// VCSDUtils_
// ----------------------------------------------------------------------------


// returns pLong given pShort and pathClass, pathType
// VCSDUtils_expandParamsShort.cpp
arma::vec expandParamsShort(arma::vec paramsShort, std::string pathClass, 
	                          std::string pathType);

// cost function

int compareCosts(double costRef, double costNew, double tol);

// VCSDUtils_costFunction.cpp
double costFunction(arma::vec paramsLong, std::string pathClass, 
	                  std::string pathType, std::string orientation, double R);

// cost function gradient
// VCSDUtils_costFunction.cpp
arma::vec costFunctionGradient(double R, std::string pathClass, 
                               std::string pathType, std::string orientation);

// linear inequality constraint
// VCSDUtils_linearInequalityConstraints.cpp
void linearConstraints(std::string pathClass, std::string pathType, 
		                   arma::vec &x_l, arma::vec &x_u, arma::vec &g_l, 
                       arma::vec &g_u, arma::mat &A, double R, double r, 
                       double asubopt_pi, double Lmax);

// non-linear in-equality constraint 
arma::vec nonlinearInequalityConstraints(arma::vec paramsShort, 
	                                       std::string pathClass, 
                                         std::string pathType, 
                                         std::string orientation, 
	                                       double r, double R);

double asubopt_func(double &r, double &R, double &beta);
double dcineq_dbeta_func(double &r, double &R, double &beta);

// non-linear equality constraint gradient
// VCSDUtils_nonlinearEqualityConstraintJacobian.cpp
arma::mat nonlinearEqualityConstraintJacobian(arma::vec paramsShort, 
	                                            std::string pathClass, 
                                              std::string orientation, 
                                              std::string pathType, 
	                                            double r, double R);

// ----------------------------------------------------------------------------
// VCSDUtils_pathProperties.cpp
// ----------------------------------------------------------------------------

// returns dimensions of the nlp problem for a given pathClass
void NLPdimensions(std::string pathClass, std::string pathType,
		               int &numParams, int &numConstraints, int &numElementsJacob);

// sets the curvature params k1, k2 given the pathClass and orientation
void determineCurvatureParams(std::string pathClass, std::string orientation, 
	                            double &k1, double &k2);

// returns a string of the orientation given the initial and final 
// curvature (k1, k2), and the pathClass
std::string determineOrientation(double k1, double k2, std::string pathClass);

// list of all possible path types/ orientations/ classes
std::vector< std::vector<std::string> > candidateList();

// ----------------------------------------------------------------------------
// VCSDUtils_suboptimalityConditions.cpp
// ----------------------------------------------------------------------------

// checkSuboptimality, returns true if suboptimal
bool checkSuboptimality(arma::vec paramsShort, std::string pathClass, 
                        std::string pathType, double R, double r);


// checkBTBSuboptimality
// check if a three or four turn sequences contains a B-BCB-B extremal with 
// equal length B arcs. returns true if suboptimal.
bool checkBTBSuboptimality(arma::vec paramsShort, std::string pathClass, 
                           std::string pathType, double R, double r);


// checkAlphaSuboptimality
// check if an extremal sequence violates the alphaSubopt condition.
// returns true if suboptimal.
bool checkAlphaSuboptimality(arma::vec paramsShort, std::string pathClass, 
	                           std::string pathType, double R, double r);
// alphaSubopt
// return value for alphaSubopt
// (i.e., for a BCB sequence, if the B arcs subtend an angle alphaSubopt
// then the sequence is suboptimal)
double alphaSubopt(double beta, double R, double r);

arma::vec  randomParam(std::string pathClass, std::string pathType, double R, 
                       double r, double Lmax);


} // namespace 

//arma::vec solveBCs(std::string pathClass, std::string pathType, 
//                   std::string orientation, double R, double r,
//                   double xf, double yf, double hf, 
//                   arma::vec paramGuess);

#endif
