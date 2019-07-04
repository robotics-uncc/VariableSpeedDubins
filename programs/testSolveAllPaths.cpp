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
	
	
	int numPaths = 76;
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

		//vsdprob.print();

		// create a candidate 
//		VarSpeedDubins::Path vsdcand;
//		vsdcand.set_turnRadii(0.1, 1.0);
//		vsdcand.set_pathClass(pathClass);
//		vsdcand.set_pathType(pathType);
//		vsdcand.set_pathOrientation(pathOrientation);	

		// create a VSDNLP from the problem statement + candidate
    //	  SmartPtr<VSDNLP> vsdnlp = new VSDNLP();
    //		vsdnlp->set_problemStatement(vsdprob);
    //		vsdnlp->set_candidate(vsdcand);
    //		vsdnlp->set_Lmax(50.0);
    //		vsdnlp->set_startingPoint(paramsShort);
		//std::cout << "vsdnlp->print(); " <<std::endl;
		//vsdnlp->print();
		//std::cout << "VarSpeedDubins::Solver vsdsolver;" <<std::endl;
		// pass nlp to solver

//		//std::cout << "vsdsolver.set_vsdnlp(vsdnlp); " <<std::endl;
//		vsdsolver.set_vsdnlp(vsdnlp);
//		//std::cout << "vsdsolver.solve(); " <<std::endl;
//		vsdsolver.solveGivenNLP();	
//		solutionStatus(i) = vsdsolver.get_solveStatus();
//		// ------------------------------------------------------------------
//		// Plot the solution
//		// ------------------------------------------------------------------

//		arma::vec optSoln = vsdnlp->get_optSoln();
//		//optSoln.print("optSoln");
//		vsdcand.set_paramsShort(optSoln);
//		vsdcand.computePathHistory();
//		vsdcand.computeEndpoint();	
//		//vsdcand.print();
//		costOpt = vsdcand.get_cost();
//		vsdcand.writePathPlotCommands(fileName,1,"x","y","g");

//    double costTol = 0.0001;
//		if (costOpt < costRef){
//			optimalityStatus(i) = 1;
//		}
//		else if (MathTools::checkBounds(costOpt, costRef-costTol, costRef+costTol)){
//			optimalityStatus(i) = 2;
//		}
//		else {
//			optimalityStatus(i) = 0;
//		}


//	}
//  



  //	for (int i = 0; i < numPaths; i++){
  //		cand = candList[i];
  //		std::string pathClass = cand[0];
  //		std::string pathType = cand[1];
  //		std::string pathOrientation = cand[2];
  //		std::cout << i << ")  " << pathClass << "," << pathType << "," << pathOrientation << ","
  //		<< MathTools::numberToString(solutionStatus(i)) << "," 
  //		<< MathTools::numberToString(optimalityStatus(i)) << std::endl;
  //	}

	//solutionStatus.print("solutionStatus:");
	// run the script in octave
//	MathTools::runOctaveScript(fileName);
  std::cout << "script complete." << std::endl;
	return 0;


}



