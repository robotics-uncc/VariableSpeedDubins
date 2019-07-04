// testSolveSinglePath.cpp
// This code generates a VSDpath and attempts passes the endpoint to the solver.
// Last Modified: 02-Feb-2015, Artur Wolek

// c++ standard
#include<string>

// external
#include<armadillo>

#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"

// custom
#include<VSDUtils.h>
#include<VSDPath.h>
#include<VSDProblem.h>
#include<VSDSolver.h>
#include<VSDNLP.h>
#include<MathTools.h>
#include<MathToolsArma.h>


int main(){


	// ------------------------------------------------------------------
	// Generate a Variable Speed Dubins Path with Known Parameters
	// ------------------------------------------------------------------
	// define a path	
	VarSpeedDubins::Path vsdpath;
	// define a parameter vector
	arma::vec paramsShort = {0.2, 0.5, 0.4, 0.1};
	// set the vehicle turn "C" and "B" turn radii
	vsdpath.set_turnRadii(0.1, 1.0);	
	// nominal spacing of the resulting path
	vsdpath.set_nominalSpacing(0.001);
	// set file to write path data and octave plotting commands 
	std::string fileName = "vsdpath.m";
	// temporary variable holding the x and y variable names
	std::string xVarName, yVarName;
	// set path parameters	
	std::vector< std::vector<std::string> > candList =  VarSpeedDubins::candidateList();
	
	
	int numPaths = 64;
	arma::vec solutionStatus(numPaths);
	arma::vec optimalityStatus(numPaths);

	std::vector<std::string> cand;


	double costRef, costOpt;

	int i = 0;
//	for (int i = 0; i < numPaths; i++){
		std::cout << "candidate : " << i << std::endl;
		cand = candList[i];
		std::string pathClass = cand[0];
		std::string pathType = cand[1];
		std::string pathOrientation = cand[2];
		vsdpath.set_pathClass(pathClass);
		vsdpath.set_pathType(pathType);
		vsdpath.set_pathOrientation(pathOrientation);		
		vsdpath.set_paramsShort(paramsShort);		
		// compute the path
		vsdpath.computePathHistory();
		// dislay path properties
		// create x,y variables 
		xVarName = "x";
		yVarName = "y";
		// write the path data and plotting commands
		vsdpath.writePathPlotCommands(fileName,1,xVarName,yVarName,"b");
		// compute endpoint
		vsdpath.computeEndpoint();
		xVarName = "xe";
		yVarName = "ye";
		// write the endpoint data and plotting commands
		vsdpath.writeEndpointPlotCommands(fileName,1,xVarName,yVarName,"r");
		//vsdpath.print();
		costRef = vsdpath.get_cost();

		// ------------------------------------------------------------------
		// Solve Variable Speed Dubins Problem
		// ------------------------------------------------------------------

		// create a problem statement
		VarSpeedDubins::ProblemStatement vsdprob;
		vsdprob.set_xFinal(vsdpath.get_xFinal());
		vsdprob.set_yFinal(vsdpath.get_yFinal());
		vsdprob.set_hFinal(vsdpath.get_hFinal());
		vsdprob.set_R(vsdpath.get_R());
		vsdprob.set_r(vsdpath.get_r());

		VarSpeedDubins::Solver vsdsolver;
    vsdsolver.set_problemStatement(vsdprob);
    arma::vec solveStatusAll = vsdsolver.get_solveStatusAll();
    solveStatusAll.print("solveStatusAll");

	  return 0;

}



