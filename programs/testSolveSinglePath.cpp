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


//  int numSamplesPerVar = 5;
//  std::string samplingAlgorithm = "barycentric";

  int numSamplesMax = 1000;
  //std::string samplingAlgorithm = "rejection-sample";
  std::string samplingAlgorithm = "barycentric";

	// ------------------------------------------------------------------
	// Generate a Variable Speed Dubins Path with Known Parameters
	// ------------------------------------------------------------------
	// define a path	
	VarSpeedDubins::Path vsdpath;
  int candidateID = 63;	// max 75
	// define a parameter vector
  arma::vec paramsShort = {0.2,   0.5,   0.7,   0.3};

	// set the vehicle turn "C" and "B" turn radii
	vsdpath.set_turnRadii(0.3, 1.0);	
	// nominal spacing of the resulting path
	vsdpath.set_nominalSpacing(0.001);
	// set file to write path data and octave plotting commands 
	std::string fileName = "vsdpath.m";
	// temporary variable holding the x and y variable names
	std::string xVarName, yVarName;
  std::cout << " get cand list " << std::endl;
	// set path parameters	
	std::vector< std::vector<std::string> > candList =  VarSpeedDubins::candidateList();
	


	std::vector<std::string> cand = candList[candidateID];

	std::string pathClass = cand[0];
	std::string pathType = cand[1];
	std::string pathOrientation = cand[2];
	vsdpath.set_pathClass(pathClass);
	vsdpath.set_pathType(pathType);
	vsdpath.set_pathOrientation(pathOrientation);		
	vsdpath.set_paramsShort(paramsShort);		

  std::cout << " compute the path " << std::endl;
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
	vsdpath.writeSwitchingPtPlotCommands(fileName,1,xVarName,yVarName,"ro");
  std::cout << " REFERENCE PATH" << std::endl;
  vsdpath.print();

	// ------------------------------------------------------------------
	// Solve Variable Speed Dubins Problem
	// ------------------------------------------------------------------
  std::cout << " solve the problem " << std::endl;
	// create a problem statement
	VarSpeedDubins::ProblemStatement vsdprob;
	vsdprob.set_xFinal(vsdpath.get_xFinal());
	vsdprob.set_yFinal(vsdpath.get_yFinal());
	vsdprob.set_hFinal(vsdpath.get_hFinal());
//  double xFinal = 0.0; // vsdpath.get_xFinal(); // 
//  double yFinal = 2.0; // vsdpath.get_yFinal(); // 
//  double hFinal = 2.0*M_PI/3.0; // vsdpath.get_hFinal(); // 
//  vsdprob.set_stateFinal(xFinal, yFinal, hFinal);
	vsdprob.set_R(vsdpath.get_R());
	vsdprob.set_r(vsdpath.get_r());

	// create a candidate 
	VarSpeedDubins::Path vsdcand;
	vsdcand.set_turnRadii(0.3, 1.0);
	vsdcand.set_pathClass(pathClass);
	vsdcand.set_pathType(pathType);
	vsdcand.set_pathOrientation(pathOrientation);	

  // create a VSDNLP from the problem statement + candidate
  SmartPtr<VSDNLP> vsdnlp = new VSDNLP();
	vsdnlp->set_problemStatement(vsdprob);
	vsdnlp->set_candidate(vsdcand);
	vsdnlp->set_Lmax(10.0);
  vsdnlp->set_samplingAlgorithm(samplingAlgorithm);
  //arma::vec paramGuess = {0, 0, 0, 0};   
	//vsdnlp->set_startingPoint(paramsShort);

	std::cout << "VarSpeedDubins::Solver vsdsolver;" <<std::endl;
	// pass nlp to solver
	VarSpeedDubins::Solver vsdsolver;
	std::cout << "vsdsolver.set_vsdnlp(vsdnlp); " <<std::endl;
	vsdsolver.set_vsdnlp(vsdnlp);
  vsdsolver.set_numSamplesMax(numSamplesMax);
  vsdsolver.set_samplingAlgorithm(samplingAlgorithm);
	std::cout << "vsdsolver.solve(); " <<std::endl;
	vsdsolver.solveGivenNLP();	
  int solveStatus = vsdsolver.get_solveStatus();
	// ------------------------------------------------------------------
	// Plot the solution
	// ------------------------------------------------------------------
  std::cout << "get_optSoln" << std::endl;
	arma::vec optSoln = vsdnlp->get_optSoln();
	vsdcand.set_paramsShort(optSoln);
	vsdcand.computePathHistory();
	vsdcand.computeEndpoint();	

	vsdcand.writePathPlotCommands(fileName,1,"x","y","g--");
	vsdcand.writeSwitchingPtPlotCommands(fileName,1,xVarName,yVarName,"ko");
  vsdcand.print();
	double costRef = vsdpath.get_cost();	
	double costOpt = vsdcand.get_cost();
	std::cout << "costRef : " << costRef << std::endl;
	std::cout << "costOpt : " << costOpt << std::endl;
  double tol = 0.0001;
	int optimalityStatus = VarSpeedDubins::compareCosts(costRef, costOpt, tol);
  std::cout << "solveSatus : " << solveStatus << std::endl;
	std::cout << "optimalityStatus : " << optimalityStatus << std::endl;
	// run the script in octave
//
	MathTools::runOctaveScript(fileName, "--persist");
	return 0;


}



