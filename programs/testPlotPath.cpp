// testPlotPath.cpp
// This code generates and plots a user defined VSDPath
// Last Modified: 02-Feb-2015, Artur Wolek

// c++ standard
#include<string>

// external
#include<armadillo>

// custom
#include<VSDUtils.h>
#include<VSDPath.h>
#include<VSDProblem.h>
#include<MathTools.h>
#include<MathToolsArma.h>

int main(){

	// define a path	
	VarSpeedDubins::Path vsdpath;
	// define a parameter vector
	//arma::vec paramsShort = {0.2, 0.5, 0.4, 0.1};
  arma::vec paramsShort = {0.7630, 1.3114, 0.7994, 0.010, 0.010};
  paramsShort.print("paramsShort");
	// set the vehicle turn "C" and "B" turn radii
	vsdpath.set_turnRadii(0.1, 1.0);
	// nominal spacing of the resulting path
	vsdpath.set_nominalSpacing(0.001);
	// set file to write path data and octave plotting commands 
	std::string fileName = "vsdpath.m";
	// temporary variable holding the x and y variable names
	std::string xVarName, yVarName;
	// set path parameters	
	vsdpath.set_pathClass("TST");
	vsdpath.set_pathType("CB-S-BC");
	vsdpath.set_pathOrientation("LSL");		
	vsdpath.set_paramsShort(paramsShort);		
	// compute the path
	vsdpath.computePathHistory();
	vsdpath.computeEndpoint();
	double cost = vsdpath.get_cost();
	// dislay path properties
	vsdpath.print();
	// create x,y variables 
	xVarName = "x";
	yVarName = "y";
	// write the path data and plotting commands
	vsdpath.writePathPlotCommands(fileName,1,xVarName,yVarName,"b");
	// compute endpoint
	xVarName = "xe";
	yVarName = "ye";
	// write the endpoint data and plotting commands
	vsdpath.writeEndpointPlotCommands(fileName,1,xVarName,yVarName,"r");
	// run the script in octave
	//MathTools::runOctaveScript(fileName);


  // test a generic BCB-S-BCB sequence
  arma::mat seg1 = VarSpeedDubins::BCBpath(0, 0.7630, 1.3114, 0, 1.0, 1.0, 0.1, 0.001); 
  double angle2 = 0.7630 + 1.3114;
  arma::mat seg2 = VarSpeedDubins::Spath(angle2, 0.7994, 0.001);
  arma::mat seg3 = VarSpeedDubins::BCBpath(0.010, 0.010, 0, 0.7630, 1.0, 1.0, 0.1, 0.001); 


  std::string testFn = "test.m";
  MathTools::writeMatrixData(testFn, seg1, "seg1");
  MathTools::writeMatrixData(testFn, seg2, "seg2");
  MathTools::writeMatrixData(testFn, seg3, "seg3");
  MathTools::runOctaveScript(testFn);
	return 0;
}

