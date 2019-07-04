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
  std::string fileName = "vsdpath.m";
  // create a problem statement
  VarSpeedDubins::ProblemStatement vsdprob;

  // intial guess selection algorithm options
  std::string algorithm = "barycentric";
  int numSamplesMax = 1000;

  double umax = 15*M_PI/180;
  double Vmax = 20;
  double Vmin = 10;
  double R = Vmax/umax;
  double r = Vmin/umax;

  double xFinal = 150; // vsdpath.get_xFinal(); // 
  double yFinal = -150; // vsdpath.get_yFinal(); // 
  double hFinal = M_PI/3.0; // vsdpath.get_hFinal(); // 
  vsdprob.set_stateFinal(xFinal, yFinal, hFinal);



  vsdprob.set_R(R);
  vsdprob.set_r(r);
  vsdprob.print();
  VarSpeedDubins::Solver vsdsolver;
    vsdsolver.set_numSamplesMax(numSamplesMax);
    vsdsolver.set_samplingAlgorithm(algorithm);
  vsdsolver.set_problemStatement(vsdprob);
  // solve
  vsdsolver.solveAll();
  vsdsolver.printAll();
  
  // outputs related to optimal soln


  arma::colvec solveStatus = vsdsolver.get_solveStatusAll();
  arma::colvec candCosts =  vsdsolver.get_candCosts();
  double multipleSolnFlag = vsdsolver.get_multipleSolnFlag();
  solveStatus.print("solveStatus");
  candCosts.print("candCostS");
  std::cout << "multipleSolnFlag" << multipleSolnFlag << std::endl;

  std::vector<int> optPathTypes = vsdsolver.get_optPathTypes();
  arma::mat optPathParam = vsdsolver.get_optPathParams();
  arma::vec optPathCosts = vsdsolver.get_optPathCosts();
  std::cout << "optPathTypes" << std::endl;
  MathTools::debugPrintVector(optPathTypes);
  optPathParam.print("optPathParam");
  optPathCosts.print("optPathCosts");

  for (int i = 0; i < optPathTypes.size(); i++){
	  // ------------------------------------------------------------------
	  // Generate a Variable Speed Dubins Path with Known Parameters
	  // ------------------------------------------------------------------
	  // define a path	
	  VarSpeedDubins::Path vsdpath;
    int candidateID = optPathTypes.at(i);
	  // define a parameter vector
    arma::vec paramsShort = optPathParam.col(i);
	  // set the vehicle turn "C" and "B" turn radii
	  vsdpath.set_turnRadii(r, R);	
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
	  xVarName = "x_optimal" + MathTools::numberToString(i);
	  yVarName = "y_optimal" + MathTools::numberToString(i);
	  // write the path data and plotting commands
	  vsdpath.writePathPlotCommands(fileName,3,xVarName,yVarName,"b");
    vsdpath.print();
  }


  vsdsolver.plotAll(fileName);
  MathTools::runOctaveScript(fileName);
  return 0;
}

