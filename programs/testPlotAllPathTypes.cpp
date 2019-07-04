// testPlotAllPathTypes.cpp
// This code generates and plots all 72 VSD path types, given a 4-element 
// parameter vector.
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
	arma::vec paramsShort = {0.2, 0.3, 0.5, 0.3};
	// set the vehicle turn "C" and "B" turn radii
	vsdpath.set_turnRadii(0.1, 1.0);
	// nominal spacing of the resulting path
	vsdpath.set_nominalSpacing(0.001);
	// set file to write path data and octave plotting commands 
	std::string fileName = "vsdpath.m";
	// temporary variable holding the x and y variable names
	std::string xVarName, yVarName;
	// list of all candidates
	std::vector< std::vector<std::string> > candidateList = 
		VarSpeedDubins::candidateList();
	// the current candidate variable 
	std::vector<std::string> candidate;
	// begin cycling through all candidates
	std::cout << "loop started..." << std::endl;
	// Note: only 64 paths from the candidate list are considered for now
	for (int i = 0; i < 64; i++){
		// set current candidate
		candidate = candidateList[i]; 
		// set path parameters	
		vsdpath.set_pathClass(candidate[0]);
		vsdpath.set_pathType(candidate[1]);
		vsdpath.set_pathOrientation(candidate[2]);		
		vsdpath.set_paramsShort(paramsShort);		
		// compute the path
		vsdpath.computePathHistory();
		// display path properties
		vsdpath.print();
		// create x,y variables with index i appended so that they are unique
		xVarName = "x" + MathTools::numberToString(i);
		yVarName = "y" + MathTools::numberToString(i);
		// write the path data and plotting commands
		vsdpath.writePathPlotCommands(fileName,1,xVarName,yVarName,"b");
		// compute endpoint
		vsdpath.computeEndpoint();
		xVarName = "xe" + MathTools::numberToString(i);
		yVarName = "ye" + MathTools::numberToString(i);
		// write the endpoint data and plotting commands
		vsdpath.writeEndpointPlotCommands(fileName,1,xVarName,yVarName,"r");
	}
	// the above function plots all paths on one figure, we want to 
	// view each path in a seperate window.
	
	// open the file
	std::ofstream outfile(fileName, std::ofstream::app);
	// write plot commands
	outfile << "figure(2); \n";
	for (int i=0; i < 64; i++){
		outfile << "axis(\"equal\"); plot(x" << i << ",y" << i << 
		",'linewidth',2); hold on; \n";
		outfile << "plot(xe" << i << ",ye" << i << 
		",'ro','linewidth',2); pause(); disp(\"" << i << "\"); \n";
	}
	// close file
    outfile.close();
	// run the script in octave
	MathTools::runOctaveScript(fileName);
	return 0;
}

