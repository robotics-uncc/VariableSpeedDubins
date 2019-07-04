// MathTools_VectorOperators.cpp

// standard C++ headers
#include<iostream> // std::cout
#include<cmath> // std::sqrt
#include<algorithm> //std::min_element

// custom
#include<MathTools.h>
#include<MathToolsArma_VectorOperators.h>
//#include<MathTools_BasicOperators.h>
//#include<MathTools_Constants.h>
//#include<MathTools_DisplayDebug.h>

// -----------------------------------------------
// COMPUTING VECTOR PROPERTIES
// -----------------------------------------------

//checkVectorBounds
bool MathTools::checkVectorBounds(arma::vec x, arma::vec x_L, arma::vec x_U){
  // check if x greater than x_L, element-wise
	arma::umat vec1 = ((x - x_L) >= 0);
  // check if x less than x_U, element-wise
	arma::umat vec2 = ((x_U - x) >= 0);
	int prod = 1; // default true
	for (int i = 0; i < vec1.n_elem; i++){
    // if there is even one zero in vec1 or vec2 this product will be zero
		prod = vec1(i)*vec2(i)*prod; 
	}
	if (prod == 1){ // check result
		return true;
	}
	else {
		return false;
	}
}

// minIndicesWithTolerance
std::vector<int> MathTools::minIndicesWithTolerance(
		                                            arma::vec testVector, 
                                                double tolerance){
	// initialize vector of min indices
	std::vector<int> minIndices;
	// store index of smallest value in testVector
  arma::uword minInd;
	// store smallest value for reference
	double minValueRef = testVector.min(minInd);
	minIndices.push_back((int) minInd); 
	// set this element to largest possible value since it has been recorded
	testVector(minInd) = std::numeric_limits<double>::max(); 
	int i = 1;
	// check if there are additional minimum values within the specified tol
	while ( std::abs(testVector.min(minInd) - minValueRef) <= tolerance ){
		minIndices.push_back((int) minInd); // store additional index
		// set to largest possible value	
		testVector(minInd) = std::numeric_limits<double>::max(); 
		i++;
	}
  //MathTools::debugPrintVector(minIndices,true);
	return minIndices;
}

// -----------------------------------------------
// GENERATING VECTORS
// -----------------------------------------------

// polarspace
arma::colvec MathTools::polarspace_arma(double angleInitial, double angleFinal, 
		int npts, int direction){
	// Generate a sequnec of evenly spaced points along an interval defining an 
	// angle, with a give direction
	// assumes less than one revolution (i.e. the interval requested is < 2PI)
	double angleInitialMod = MathTools::mod(angleInitial, 2.0*M_PI);
	double angleFinalMod;
	if (direction == 1) { //counter-clockwise
		angleFinalMod = MathTools::mod(angleFinal, 2.0*M_PI);
		if (angleFinalMod > angleInitialMod){ // interval excludes 2*PI
		}
		else { // the interval includes 2*PI, so we make the final angle > 2*PI
			angleFinalMod = angleFinalMod + 2.0*M_PI; 
		}
	}
	else if (direction == -1){ //clockwise
		// by making the final angle negative, we force the clockwise direction
		angleFinalMod = MathTools::mod(angleFinal, 2.0*M_PI) - 2.0*M_PI;
	}
	else {
		std::cout << "MathTools::polarspace, "
			 	      << "Error: direction specified is ambiguous" << std::endl;
	}
	return arma::linspace(angleInitialMod, angleFinalMod, npts);
}

// -----------------------------------------------
// MODIFYING VECTORS / VECTOR OPERATIONS
// -----------------------------------------------

// augmentVector
arma::colvec MathTools::augmentVector(double t, arma::colvec x){
	arma::colvec z(x.n_elem + 1);
	z(arma::span(0,x.n_elem - 1)) = x; // set all but last element of z to x
	z(x.n_elem) = t; // set last element of z to t
	return z;

}

// caretesianProduct
arma::mat MathTools::cartesianProduct(arma::mat basisVectorsMatrix){

	int m = basisVectorsMatrix.n_cols; // number of elements in dimension
	int n = basisVectorsMatrix.n_rows; // number of dimensions
	
	int cpNumPts = pow(m,n);
	arma::mat cpMat(cpNumPts,n);
	arma::mat indMat = arma::ones<arma::mat>(cpNumPts,n);
	indMat = indMat*-1;
	arma::vec indVec = arma::zeros<arma::vec>(n); // gives current index 
  int i = 0;
	int j = 1;
	int curElem = 0;
	bool incremented = false;
	while (i < cpNumPts){
		indVec(0) = curElem;				
		indMat.row(i) = indVec.t(); // indices of new pt of cp matrix
		// add new entry to cpMat according to indices of indVec
		for (int k = 0; k < n; k++){ 
			cpMat(i,k) = basisVectorsMatrix(k,indVec(k));
		}
		curElem++;
		i++;
		if (curElem >= m && (i < cpNumPts) ){
			curElem = 0;
			while(incremented == false){
				if (indVec(j) < m - 1){
					indVec(j) = indVec(j) + 1; 
					incremented = true;
					j = 1;
				}
				else { 
					if (j < n - 1) {// move on to next dimension
						indVec(j) = 0; 
						j = j + 1;
					}
					else {
						j = 1;
					}
				}
			}
			incremented = false;	
		}
	}
	return cpMat;
}

// swap
void MathTools::swap(arma::vec &x, int i, int j){
  double ival = x(i); 
  x(i) = x(j);
  x(j) = ival;
  return;
}

// -----------------------------------------------
// MODIFYING MATRICES
// -----------------------------------------------

// addVectorToMatrixColoumnWise
arma::mat MathTools::addVectorToMatrixColoumnWise(arma::rowvec a, arma::mat B){
	arma::mat Bshifted(size(B));
	for (int i = 0; i < B.n_cols - 1; i++){
		Bshifted.col(i) = B.col(i) + a(i);
	}
	return Bshifted;
}

// join_horiz
arma::mat MathTools::join_horiz(arma::colvec a, arma::colvec b, arma::colvec c){
	arma::mat a_b = arma::join_horiz(a,b);
	arma::mat a_b_c = arma::join_horiz(a_b,c);
	return a_b_c;
}

// join_horiz
arma::mat MathTools::join_horiz(arma::mat a, arma::mat b, arma::mat c){
	arma::mat a_b = arma::join_horiz(a,b);
	arma::mat a_b_c = arma::join_horiz(a_b,c);
	return a_b_c;
}

// join_horiz
arma::mat MathTools::join_horiz(arma::vec a, arma::vec b, arma::vec c, 
								arma::vec d){
	arma::mat a_b = arma::join_horiz(a,b);
	arma::mat a_b_c = arma::join_horiz(a_b,c);
	arma::mat a_b_c_d = arma::join_horiz(a_b_c,d);
	return a_b_c_d;
}

// join_horiz
arma::mat MathTools::join_horiz(arma::mat a, arma::mat b, arma::mat c, 
								arma::mat d){
	arma::mat a_b = arma::join_horiz(a,b);
	arma::mat a_b_c = arma::join_horiz(a_b,c);
	arma::mat a_b_c_d = arma::join_horiz(a_b_c,d);
	return a_b_c_d;
}

// join_vert
arma::mat MathTools::join_vert(arma::vec a, arma::vec b, arma::vec c){
	arma::mat a_b = arma::join_vert(a,b);
	arma::mat a_b_c = arma::join_vert(a_b,c);
	return a_b_c;
}

// join_vert
arma::mat MathTools::join_vert(arma::mat a, arma::mat b, arma::mat c){
	arma::mat a_b = arma::join_vert(a,b);
	arma::mat a_b_c = arma::join_vert(a_b,c);
	return a_b_c;
}

// join_vert
arma::mat MathTools::join_vert(arma::vec a, arma::vec b, arma::vec c, 
							   arma::vec d){
	arma::mat a_b = arma::join_vert(a,b);
	arma::mat a_b_c = arma::join_vert(a_b,c);
	arma::mat a_b_c_d = arma::join_vert(a_b_c,d);
	return a_b_c_d;
}

// join_vert
arma::mat MathTools::join_vert(arma::mat a, arma::mat b, arma::mat c, 
							   arma::mat d){
	arma::mat a_b = arma::join_vert(a,b);
	arma::mat a_b_c = arma::join_vert(a_b,c);
	arma::mat a_b_c_d = arma::join_vert(a_b_c,d);
	return a_b_c_d;
}

// stdToArmaVector
arma::vec MathTools::stdToArmaVector(std::vector<double> inputVector){
  //   
  arma::vec outputVector(inputVector.size());
  for (int i = 0; i < inputVector.size(); i++){
    outputVector(i) = inputVector.at(i);
  }
  return outputVector;
}
