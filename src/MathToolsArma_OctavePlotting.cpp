// MathTools_OctavePlotting.cpp

// standard headers
#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<vector>

// third-party
#include<armadillo>

// custom
#include<MathToolsArma_OctavePlotting.h>
#include<MathToolsArma_Geometry.h>
#include<MathToolsArma_VectorOperators.h>

#include<MathTools.h>


// writeMatrixData
void MathTools::writeMatrixData(std::string fileName, arma::mat &data, 
                                std::string varName){
	std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << varName << "= [ ";
	// assume that data is
	std::vector< int > currentRow;
  for (int i = 0; i < data.n_rows; i++){ // for each row
		for (int j = 0; j < data.n_cols; j++){
			outfile << " , " << data(i,j);
		}
        outfile << " ; "; // end of row
  }	
	outfile << "]; \n";
  outfile.close();	   	// close file
}

// writeVectorData for arm::colvec data
void MathTools::writeVectorData(std::string fileName, arma::colvec data, 
		                            std::string varName){
	std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << varName << "= [ " << data(0); // write coordinates
    for (int i = 1; i < data.n_elem; i++)
    {
    	outfile << " , " << data(i);
    }	
	outfile << "]; \n";
  outfile.close();	   	// close file
	return;
}

// writeVectorData for arm::rowvec data
void MathTools::writeVectorData(std::string fileName, arma::rowvec data, 
		                            std::string varName){
	arma::colvec colData = data.t();
	MathTools::writeVectorData(fileName, colData, varName);
	return;
}

// writeVectorIntToCell
void MathTools::writeVectorArmaVecToCell(std::string fileName, 
                                         std::vector< arma::vec > data,  
                                         std::string varName){
	std::ofstream outfile(fileName, std::ofstream::app); // open for writing 
  outfile << varName << "= cell \n";
  for (int i = 0; i < data.size(); i++){
    outfile << varName << "{" << i+1 << "} = [";
    for (int j = 0; j < (data.at(i)).n_elem; j++){
      outfile << "," << (data.at(i))(j);
    }
	  outfile << "]; \n";
  }	
  outfile.close();	   	// close file  
}

// writeVectorArmaMatToCell
void MathTools::writeVectorArmaMatToCell(std::string fileName, 
                                         std::vector< arma::mat > data,  
                                         std::string varName){
	std::ofstream outfile(fileName, std::ofstream::app); // open for writing 
  outfile << varName << "= cell \n";
  for (int i = 0; i < data.size(); i++){
    outfile << varName << "{" << i+1 << "} = [";
    for (int j = 0; j < (data.at(i)).n_rows; j++){
      for (int k = 0; k < (data.at(i)).n_cols; k++){
        outfile << "," << (data.at(i))(j,k);
      }
	    outfile << "; \n";
    }
	  outfile << "]; \n";
  }	
  outfile.close();	   	// close file  
}

// plot1DCurve for arm::vec data
void MathTools::plot1DCurve(std::string fileName, int figureNumber, 
		                        arma::vec xData, arma::vec yData, 
                            std::string xVarName, std::string yVarName, 
                            std::string colorStyle){	
	MathTools::writeVectorData(fileName, xData, xVarName);
	MathTools::writeVectorData(fileName, yData, yVarName);
	MathTools::plot1DCurve(fileName,figureNumber,xVarName,yVarName,colorStyle);
}

// plot1DCurve for arm::vec data
void MathTools::plotMatrixRows(std::string fileName, int figureNumber, 
		                           arma::mat data, std::string dataName, 
                               std::string colorStyle){
	// assume cols of matrix are indpendent samples, and rows of matrix 
	// correspond to particular variables
	// the number of rows must be less than <=3
	int numRows = data.n_rows;
	int numCols = data.n_cols;
	if (numRows < numCols){ // assume the rows are data pts and the cols are 
						   // indpe variables
		if (numRows == 2){
			arma::rowvec xDataRow = data.row(0);
			arma::rowvec yDataRow = data.row(1);
			arma::colvec xData = xDataRow.t();
			arma::colvec yData = yDataRow.t();
			std::string xDataString = dataName + "_x";
			std::string yDataString = dataName + "_y";
			MathTools::plot1DCurve(fileName, figureNumber, xData, yData, 
				xDataString, yDataString, colorStyle);
		}
		else if (numRows == 3){
			arma::vec xData = data.row(0).t();
			arma::vec yData = data.row(1).t();
			arma::vec zData = data.row(2).t();
			std::string xDataString = dataName + "_x";
			std::string yDataString = dataName + "_y";
			std::string zDataString = dataName + "_z";
			MathTools::plotPoints3D(fileName, figureNumber, xData, yData, zData, 
								   xDataString, yDataString, zDataString, 
								   colorStyle);
		}
		else {
			std::cout << "MathTools: plotMatrixRows(): Error: "
					  << "Invalid dimension of input matrix: (" << numRows 
					  << " x " << numCols << ")" << std::endl;
		}
	}
	else {
		if (numCols == 2){
			arma::vec xData = data.col(0);
			arma::vec yData = data.col(1);
			std::string xDataString = dataName + "_x";
			std::string yDataString = dataName + "_y";
			MathTools::plot1DCurve(fileName, figureNumber, xData, yData, 
				xDataString, yDataString, colorStyle);
		}
		else if (numCols == 3){
			arma::vec xData = data.col(0);
			arma::vec yData = data.col(1);
			arma::vec zData = data.col(2);
			std::string xDataString = dataName + "_x";
			std::string yDataString = dataName + "_y";
			std::string zDataString = dataName + "_z";
			MathTools::plotPoints3D(fileName, figureNumber, xData, yData, zData, 
								   xDataString, yDataString, zDataString, 
								   colorStyle);
		}
		else {
			std::cout << "MathTools: plotMatrixRows(): Error: "
					  << "Invalid dimension of input matrix: (" << numRows 
					  << " x " << numCols << ")" << std::endl;
		}
	}
}

// plotPoint3D
void MathTools::plotPoints3D(std::string fileName, int figureNumber, 
							               arma::vec x, arma::vec y, arma::vec z, 
							               std::string xVarName, std::string yVarName, 
							               std::string zVarName, std::string style){
	MathTools::writeVectorData(fileName, x, xVarName);
	MathTools::writeVectorData(fileName, y, yVarName);
	MathTools::writeVectorData(fileName, z, zVarName);
	std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << "figure " << figureNumber << " \n ";
	outfile << "plot3("<< xVarName << " , " << yVarName << " , " 
			<< zVarName <<",'" << style << "', 'linewidth',3); hold on; \n";
  // close file
  outfile.close();
}

// plotArrow2D
void MathTools::plotArrow2D(std::string fileName, int figureNumber, 
		arma::vec startPt, arma::vec endPt, std::string arrowName, 
		std::string colorStyle){
	// trunk of the arrow
	arma::vec xTrunk = {startPt(0), endPt(0)};
	arma::vec yTrunk = {startPt(1), endPt(1)};
	MathTools::plot1DCurve(fileName, figureNumber, xTrunk, yTrunk, arrowName + 
		"_x", arrowName + "_y", colorStyle);
}

// plotEllipse2D
void MathTools::plotEllipse2D(std::string fileName, int figureNumber, 
		arma::vec ellipseParams, std::string ellipseName, 
		std::string style){
	arma::mat pts = MathTools::generateEllipse(ellipseParams, 100);
	std::cout << "pts : " << pts.size() << std::endl;
	plot1DCurve(fileName, figureNumber, pts.col(0), pts.col(1), ellipseName + 
		"_x", ellipseName + "_y", style);
}

