#include<MathTools.h>
#include<MathToolsArma_Polytopes.h>
#include<MathToolsArma_Combinatronics.h>
#include<MathToolsArma_VectorOperators.h>

// cddlib
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
//


// -----------------------------------------------
// REJECTION SAMPLING
// -----------------------------------------------

// uniformRejectionSample
// x is a N dimensional vector and each dimension of x is partitioned into an 
// equal number S of values. 
// e.g. if the first dimension x1 varies from xmin, xmax then the partitioning
// is [xmin, xmin + (xmax-xmin)/S, ... , xmax];
// This is repeated for every dimension. The cartesian product of the 
// resulting set is computed (it containts S^N number of samples). Every 
// sample in this resulting cartesian product set is tested to see if it 
// satisfies the constraints. Only the samples that satisfy the constraints are 
// returned. The numSamplesMax input specifies an upper bound on the number of
// samples that will be generated, i.e. the largest admissible S is chosen 
// such that S^N <= numSamplesMax
arma::mat MathTools::uniformRejectionSample(int numSamplesMax, arma::vec xl, 
								                            arma::vec xu, arma::vec gl, 
                                            arma::vec gu, arma::mat A){
	int n = xl.n_elem; // get the dimension of x
  int numSamplesPerDim = uniformSamplesPerDim(numSamplesMax, n);
	arma::mat basisVectorMatrix(n,numSamplesPerDim); 
	arma::vec basisVectorCurrent;
	for (int i = 0; i < n; i++){
		basisVectorCurrent = arma::linspace<arma::vec>(xl(i), xu(i), 
                                                   numSamplesPerDim);
		basisVectorMatrix.row(i) = basisVectorCurrent.t();
	}
	// generate cartesian product 
	arma::mat cp = MathTools::cartesianProduct(basisVectorMatrix);
  int numSamplesCandidates = pow(numSamplesPerDim,n);
	arma::mat validSamples = arma::zeros<arma::mat>(numSamplesCandidates, n);
	arma::vec testVector;
	arma::vec testResult;
	// check if constraints are satisfied, if so add to 
  int k = 0;
	for (int j = 0; j < cp.n_rows; j++){
		testVector = cp.row(j).t();
		testResult = A*(testVector);
		if (MathTools::checkVectorBounds(testResult, gl, gu)){
			validSamples.row(k) = testVector.t();
			k++;
		}
	}
	validSamples.shed_rows(k-1,numSamplesCandidates-1);
	return validSamples;
}

// uniformSamplesPerDim
int MathTools::uniformSamplesPerDim(int numSamplesMax, int numDims){	
  int numSamplesPerDim = 1;
  int numSamples = 1;
  // increment numSamplesPerDim until the max is reached
  // recall that the number of samples in numSamplesPerDim^numDims
  while (numSamples <= numSamplesMax){    
    numSamplesPerDim++;
    numSamples = pow(numSamplesPerDim,numDims);
  }
  return numSamplesPerDim-1;
}

// -----------------------------------------------
// BARYCENTRIC SAMPLING
// -----------------------------------------------


// barycentricSample
arma::mat MathTools::barycentricSample(int numSamplesMax, arma::vec xl, 
                                       arma::vec xu, arma::vec gl, arma::vec gu, 
                                       arma::mat A){
  arma::mat vertices = MathTools::vertexEnumeration(xl, xu, gl, gu, A);
  int n = xl.n_elem;
  arma::mat samples = MathTools::verticesToBarycentricSamples(numSamplesMax, 
                                                              vertices, n);
  return samples;
}

// barycentricSample
arma::mat MathTools::barycentricSample(int numSamplesMax, arma::vec xl, 
                                       arma::vec xu, arma::vec gu, arma::mat A){
  arma::mat vertices = MathTools::vertexEnumeration(xl, xu, gu, A);
  int n = xl.n_elem;
  arma::mat samples = MathTools::verticesToBarycentricSamples(numSamplesMax, 
                                                              vertices, n);
  return samples;
}

// verticesToBarycentricSamples
arma::mat MathTools::verticesToBarycentricSamples(int &numSamplesMax, 
                                                  arma::mat &vertices,
                                                  int &n){
  int numVertices = vertices.n_rows;
  int numSamplesPerVertex = MathTools::barycentricSamplesPerVertex(
                                                                 numSamplesMax, 
                                                                 numVertices);
  arma::mat weightsInt = MathTools::constantSumSubsets(numSamplesPerVertex, 
                                                       numVertices);
  arma::mat weights = weightsInt / numSamplesPerVertex;
  arma::mat samples(weights.n_rows, n);
  for (int i = 0; i < weights.n_rows; i++){ 
    arma::vec sample = arma::zeros<arma::vec>(n);
    for (int j = 0; j < numVertices; j++){
      sample += weights(i,j) * vertices.row(j).t();
    }
    samples.row(i) = sample.t();
  }
  return samples;
}

// returns the number of samples per vertex required such that the resulting 
// number of barycentric samples does not exceed numSamplesMax
int MathTools::barycentricSamplesPerVertex(int numSamplesMax, int numVertices){
  int numSamplesPerVertex = 1;
  int numSamples = 1;
  while (numSamples <= numSamplesMax){    
    numSamplesPerVertex++;
    numSamples = MathTools::binomial(numSamplesPerVertex + numVertices - 1, 
                                     numVertices - 1 ); 
  }
  return numSamplesPerVertex-1;
}

