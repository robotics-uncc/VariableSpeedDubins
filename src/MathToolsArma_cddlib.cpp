
#include<MathToolsArma_Polytopes.h> // functions prototyped here
#include<MathToolsArma_VectorOperators.h>
#include<armadillo>

// cddlib
#include "setoper.h"
#include "cdd.h"

// -----------------------------------------------
// VERTEX ENUMERATION
// -----------------------------------------------


// vertexEnumeration
arma::mat MathTools::vertexEnumeration(arma::vec xl, arma::vec xu, arma::vec gl, 
                                       arma::vec gu, arma::mat C){
  // get size of problem
	int numVars = xl.n_elem;
  int numConstraints = 2*numVars + 2*gu.n_elem;
  // input is of the form:
  // xl <= x  <= xu
  // gl <= Cx <= gu
  // but we need to conver this to one big system of the form
  // Dx <= b
  // the first n constraints are lower bounds
  // the next n constraints are upper bounds
  arma::mat I(numVars,numVars,arma::fill::eye); // identity matrix
  arma::vec b = MathTools::join_vert(xl, xu, gl, gu);
  arma::mat D = MathTools::join_vert(-I, I, -C, C); 
  arma::mat bmAx = arma::join_horiz(b, -D); // b-Ax >= 0
  return MathTools::runCddlib(numConstraints, numVars, bmAx);
}

// vertexEnumeration
arma::mat MathTools::vertexEnumeration(arma::vec xl, arma::vec xu, arma::vec gu, 
                                       arma::mat C){
	int numVars = xl.n_elem;
  int numConstraints = 2*numVars + gu.n_elem;
  arma::mat I(numVars,numVars,arma::fill::eye); // identity matrix
  arma::vec b = MathTools::join_vert(xl, xu, gu);
  arma::mat D = MathTools::join_vert(-I, I, C); 
  arma::mat bmAx = arma::join_horiz(b, -D); // b-Ax >= 0
  return MathTools::runCddlib(numConstraints, numVars, bmAx);
}

// vertexEnumeration
arma::mat MathTools::vertexEnumeration(arma::vec xl, arma::vec xu){
	int numVars = xl.n_elem;
  int numConstraints = 2*numVars;
  arma::mat I(numVars,numVars,arma::fill::eye); // identity matrix
  arma::vec b = arma::join_vert(xl, xu);
  arma::mat D = arma::join_vert(-I, I); 
  arma::mat bmAx = arma::join_horiz(b, -D); // b-Ax >= 0
  return MathTools::runCddlib(numConstraints, numVars, bmAx);
}

// runCddlib
arma::mat MathTools::runCddlib(int numConstraints, int numVars, 
                               arma::mat &bmAx){
  // the following is modified from a cddlib example to accept arma:mat
  dd_PolyhedraPtr poly;
  dd_MatrixPtr A, B, G;
  dd_rowrange m; 
  dd_colrange d;
  dd_ErrorType err;
  m=numConstraints; 
  d=numVars+1;
  dd_set_global_constants();  /* First, this must be called to use cddlib. */
  A=dd_CreateMatrix(m,d);
  // define the A matrix
  for (int i = 0; i < numConstraints; i++){
    for (int j = 0; j < numVars+1; j++){
      dd_set_d(A->matrix[i][j],bmAx(i,j));
    }
  }
  A->representation=dd_Inequality;
  poly=dd_DDMatrix2Poly(A, &err);  /* compute the second (generator) rep */
  G=dd_CopyGenerators(poly);
  // Optional: Print matrices
  //dd_WriteMatrix(stdout,A);  printf("\n");
  //dd_WriteMatrix(stdout,G);

  // extract vertex values
  int rowsG = G->rowsize;
  int colsG = G->colsize;
  arma::mat vertices(rowsG, colsG-1);
  for (int i = 0; i < rowsG; i++){
    for (int j = 1; j < colsG; j++){
      vertices(i,j-1) = (G->matrix[i][j])[0];
    }
  }
  dd_FreeMatrix(A);
  dd_FreeMatrix(G);
  return vertices;
}

