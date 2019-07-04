// MathToolsArma_OctavePlotting.h
// Provides a set of plotting tools for use with GNU Octave.

#ifndef MATH_TOOLS_ARMA_PLOTTING_TOOLS_H
#define MATH_TOOLS_ARMA_PLOTTING_TOOLS_H

// standard headers
#include<vector>
#include<string>

#include<MathToolsArma_OctavePlotting.h>

namespace MathTools {

// Generates octave commands to write vector data to .m file
void writeVectorData(std::string fileName, arma::colvec data, 
	                   std::string varName);
void writeVectorData(std::string fileName, arma::rowvec data, 
	                   std::string varName);

// Generate octave commands to write matrix data to .m file
void writeMatrixData(std::string fileName, arma::mat &data, 
                     std::string varName);

// writeVetorArmVecToCell
void writeVectorArmaVecToCell(std::string fileName, 
                              std::vector< arma::vec > data,
                              std::string varName);

void writeVectorArmaMatToCell(std::string fileName, 
                              std::vector< arma::mat > data,
                              std::string varName);

// Generates octave commands to plot a curve (with std::vector, 
// arma::vec passed in as arguments)
void plot1DCurve(std::string fileName, int figureNumber, arma::vec xData, 
	               arma::vec yData, std::string xVarName, std::string yVarName, 
	               std::string style);

// Generates octave commands to plot a series of data points from a matrix
void plotMatrixRows(std::string fileName, int figureNumber, arma::mat data, 
	                  std::string dataName, std::string colorStyle);

// Generate octave commands to plot a 3D point
void plotPoints3D(std::string fileName, int figureNumber, arma::vec x, 
				          arma::vec y, arma::vec z, std::string xVarName, 
				          std::string yVarName, std::string zVarName, 
				          std::string style);

// Generates octave commands to plot a 2D arrow
void plotArrow2D(std::string fileName, int figureNumber, arma::vec start, 
				         arma::vec end, std::string arrowName, std::string colorStyle);

// Generates octave commands to plot an ellipse
void plotEllipse2D(std::string fileName, int figureNumber, 
				           arma::vec ellipseParams, std::string ellipseName, 
				           std::string colorStyle);

}

#endif
