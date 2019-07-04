#ifndef MATH_TOOLS_VECTOR_OPERATORS_H
#define MATH_TOOLS_VECTOR_OPERATORS_H

// standard C++ headers
#include<vector> // std::vector
#include<cmath> // floor, cos, M_PI
#include<algorithm> //std::min_element
#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<iomanip>      // std::setprecision
#include<sstream>

namespace MathTools {

// return the sign of a number
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}


// Finds minimum element(s) in a vector given a tolerance
std::vector<int> minIndicesWithTolerance(std::vector<double> testVector, 
										                     double tolerance);

// determine (smallest) angular angular distance between two angles (rad)
double polarDistance(double a, double b);

// Computes Euclidean distance between two vectors by taking the p-norm
double distance(std::vector<double> referenceVector, 
				std::vector<double> testVector);

// return the modulus of a number, e.g., mod(3pi/2, pi) = pi/2, forces positive
double        mod(double number, double base); 

// Modify a std::vector by appending a second vector to its end
void append(std::vector<double> &baseVector, std::vector<double> &appendVector);

// Returns the index of the minimum element in a vector
int minElement(std::vector<double> testVector);

// Returns the value of the minimum element in a vector
double minValue(std::vector<double> testVector);

// Generate a sequence of evenly spaced points within prescribed bounds
// if desired skip nSkip several points from the begining
std::vector<double> linspace(double xmin, double xmax, double spacing, 
							               int nSkip);
std::vector<double> linspace(double xmin, double xmax, double spacing);

// check if a double value is an integer
bool isInteger(double value);

// Generates octave commands to write vector data to .m file
void writeVectorData(std::string fileName, std::vector<double> data, 
	                   std::string varName);	

// Generates octave commands to plot a curve (with reference to existing vectrs)
void plot1DCurve(std::string fileName, int figureNumber, std::string xVarName, 
	               std::string yVarName, std::string style);

// Set axis properties e.g. "equal" , "tight"
void axisProperty(std::string fileName, std::string property);

// Generate octave commands to write matrix data to .m file
void writeMatrixData(std::string fileName, 
		                 std::vector< std::vector<double> > &data, 
                     std::string varName);
void writeMatrixData(std::string fileName, 
		                 std::vector< std::vector<int> > &data, 
                     std::string varName);

// Runs a (.m) script in octave
void runOctaveScript(std::string fileName);
void runOctaveScript(std::string fileName, std::string options);

// Add command
void addCommand(std::string fileName, std::string command);

// check if an interval [a,b] intersects [c,d] at any point
bool checkIntervalIntersect(double a, double b, double c, double d);

// Generates octave commands to plot a 2D point
void plotPoint(std::string fileName, int figureNumber, double x, double y, 
			         std::string xVarName, std::string yVarName, std::string style);

// returns an integer number of points that divides arcLength by the 
// nominal spacing. Handles case when requested spacing is smaller than arcleng
int safeCurveSpacing(double arcLength, double nominalSpacing);

// returns the n-th row and k-th entry in Pascal's triangle
// Note: there are n+1 entires per row.
// i.e.   n = 1             1     ( k = 1 entries )
//        n = 2           1 2 1
//        n = 3          1 3 3 1 
//        n = 4         1 4 6 4 1  ( k <= 5 entries)
int pascalsTriangle(int n, int k);

// return binomial coefficient "n choose k"
// returns bionomial "n choose k"
// The binomial(n,k) is equal to the n-th row and k-th entry of Pascal's 
// triangle. This approach does not suffer from the overflow that occurs with
// factorial values 13! or greater. 
int binomial(int n, int k);

// check if x is within a given interval [a,b]
bool checkBounds(double x, double a, double b);

// Generate octave commands to write a single variable to .m file
void writeVariable(std::string fileName, double x, std::string varName);

// returns the factorial: x! 
  // Warning: 13! >  2,147,483,647 ( the largest admissible int value )
  // thus overflow will occur if n >= 13
int factorial(int x);

// check if an integer is even
bool isEven(int x);

// PI is given as M_PI
// angles and conversions
const double RAD2DEG = 180.0/M_PI;
const double DEG2RAD = M_PI/180.0;

// Print a vector to standard output
void debugPrintVector(std::vector<double> varVector, bool debugFlag);
void debugPrintVector(std::vector<double> varVector);
void debugPrintVector(std::vector<int> varVector, bool debugFlag);
void debugPrintVector(std::vector<int> varVector);

// numberToString
// convert a number to a string
template <typename T>
std::string numberToString ( T Number ){
	std::stringstream ss;
	ss << Number;
	return ss.str();
}
} // namespace

#endif
