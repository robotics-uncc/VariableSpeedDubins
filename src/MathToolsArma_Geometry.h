// MathToolsArma_Geometry.h

#ifndef MATH_TOOLS_ARMA_GEOMETRY_H
#define MATH_TOOLS_ARMA_GEOMETRY_H

// C++ standard
#include<vector>

// external
#include<armadillo>
#include<MathToolsArma_Geometry.h>

namespace MathTools{


// calculate the arc-length of a curve given as y = f(x)
void arclength(arma::vec &x, arma::vec &y, arma::vec &arclength);

// Generates (x,y) coordinates defining a two dimensional ellipse
	// ellipseParams is a vector with the following fields:
	// ellipseParams(0) = xc, x coord of center of ellipse
	// ellipseParams(1) = yc, y coord of center of ellipse
	// ellipseParams(2) = semiMajor, axis of ellipse
	// ellipseParams(3) = semiMinor, axis of ellipse
	// ellipseParams(4) = angle, relative to horizontal (+ CCW) of the semi-Major
arma::mat generateEllipse(arma::colvec ellipseParams, int numPts);

// Generates (x,y) coordinates defining a circualr arc
	// arcParams is a vector with the following fields:
	// arcParams(0) = xc, x coord of center of arc
	// arcParams(1) = yc, y coord of center of arc
	// arcParams(2) = R, radius of arc
	// arcParams(3) = angleInit, initial angle relative to horizontal
	// arcParams(4) = angleFinal, final angle relative to horizontal
	// numPts = points defining the arc
	// direction = 1, CCW is -1 CW
arma::mat generateCircularArc(arma::colvec arcParams, int numPts, int direction);

// Generates (x,y) coordinates defining a circle
	// circleParams is a vector with the following fields:
	// circleParams(0) = xc, x coord of center of circle
	// circleParams(1) = yc, y coord of center of circle
	// circleParams(2) = R, radius of circle
arma::mat generateCircle(arma::colvec circleParams, int numPts);

// rotationMatrix
// returns the matrix R that rotates points in the x-y plane counter-clockwise 
// by an angle theta. i.e., x' = Rx where x is the original coordinate and 
// x' is the transformed point
arma::mat rotationMatrix(double theta);

// rotateCoordinates2D
// rotates each row of coords by a rotation matrix defined by R(theta) and 
// returns a list of rotated coordinates
arma::mat rotateCoordinates2D(arma::mat coords, double theta);

}

#endif
