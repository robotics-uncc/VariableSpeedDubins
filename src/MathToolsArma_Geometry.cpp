// MathTools_Geometry.cpp

// custom
#include<MathToolsArma_Geometry.h>
#include<MathTools.h>
#include<MathToolsArma_VectorOperators.h>
//#include<MathToolsArma_NumericalIntegration.h>

// arclength
void MathTools::arclength(arma::vec &x, arma::vec &y, arma::vec &arclength){
  // the formula for arc-length of a curve y(x)
  // s = \Integral_(a,b)  Sqrt( 1 + (dy/dx)^2 ) dx
  // numerically take derivative
  arma::vec dy_dx = diff(y)/diff(x); // result is 1 elem shorter than x,y
  arma::vec dx = diff(x);
  // define the integrand
  arma::vec integrand = arma::sqrt(1.0 + (dy_dx % dy_dx));
  // simple integration
  arclength.set_size(x.size());
  arclength(0) = 0.0; // initial condition
  for (int i = 1; i < x.size(); i++){
    arclength(i) = integrand(i-1)*dx(i-1) + arclength(i-1);
  }  
  // add to initial condition
} 

// generateEllipse (arma::vec)
arma::mat MathTools::generateEllipse(arma::vec ellipseParams, int numPts){
	// The general parametric equation of an ellipse, centered around (xc, yc),
	// inclined at angle phi relative to the x-axis
	// and with semi major axis a, and semi-minor axis b is:
	//
	//	x(t) = xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi)
	//	y(t) = yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi)		
	// 	
	//	where: 		0 <= t <= 2*pi
	//
	double xc = ellipseParams(0);
	double yc = ellipseParams(1);
	double semiMajor = ellipseParams(2);
	double semiMinor = ellipseParams(3);
	double angle = ellipseParams(4);
	arma:: vec t = arma::linspace(0.0, 2.0*M_PI, numPts);
	arma:: vec xData = xc + semiMajor*cos(t)*cos(angle) - 
		semiMinor*sin(t)*sin(angle);
	arma:: vec yData = yc + semiMajor*cos(t)*sin(angle) + 
		semiMinor*sin(t)*cos(angle);
	arma::mat ellipsePoints(numPts, 2);
	ellipsePoints.col(0) = xData;
	ellipsePoints.col(1) = yData;
	return ellipsePoints;
}


// generateCircularArc
arma::mat MathTools::generateCircularArc(arma::vec arcParams, int numPts, 
		int direction){
	// The general paramtric equation for a full circle centered around (xc, yc) 
	// with radius R is given by:
	//
	// x(t) = R*cos(t) + xc		where: 		0 <= t <= 2*pi
	// y(t) = R*sin(t) + yc
	// 	
	// This circle will begin with t = 0 at it's "east-most" point: 
	// (x,y) = (R + xc, yc) 
	// To generate an arc we need to modify the range of t appropriately. 
	double xc = arcParams(0);
	double yc = arcParams(1);
	double R = arcParams(2);
	double tInit = arcParams(3);
	double tFinal = arcParams(4);
	// linear spacing in angle is equivalent to linear spacing in arc-length 
	// for a circle
	arma::vec t = MathTools::polarspace_arma(tInit, tFinal, numPts, direction);
	arma::mat arcPts;
	arcPts.col(0) = R*cos(t) + xc;
	arcPts.col(1) = R*sin(t) + yc;
	return arcPts;
}



// Generates (x,y) coordinates defining a two dimensional circle
arma::mat MathTools::generateCircle(arma::colvec circleParams, int numPts){
	// add parameters necessary for arc definition
	circleParams(3) = 0;
	circleParams(4) = 2.0*M_PI;
	return MathTools::generateCircularArc(circleParams, numPts, 1.0); 
	// default is CCW circle
}


// rotationMatrix
arma::mat MathTools::rotationMatrix(double theta){
  arma::mat R(2,2);
  R = { {cos(theta), -sin(theta)},
        {sin(theta), cos(theta)} };
  return R;
}

// rotateCoordinates2D
arma::mat MathTools::rotateCoordinates2D(arma::mat coords, double theta){
  // create rotation matrix
  arma::mat R = MathTools::rotationMatrix(theta);
  // assume that coords has npts rows and 2 cols
  int npts = coords.n_rows;
  // temporary vectors
  arma::colvec x(2), xtransformed(2);  
  arma::mat rotatedCoords(npts,2);
  // go through each row
  for (int i = 0; i < npts; i++){
    x = coords.row(i).t(); // input point column vector
    xtransformed = R*x; // transformed point column vector
    rotatedCoords.row(i) = xtransformed.t(); // save to output
  }  
  return rotatedCoords;
}


