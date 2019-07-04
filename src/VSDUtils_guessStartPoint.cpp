#include<VSDUtils.h>

// guessStartPoint
arma::vec VarSpeedDubins::guessStartPoint(std::string algorithm, 
                          std::string pathClass, std::string pathType, 
                          std::string orientation, arma::vec xl, arma::vec xu, 
                          arma::vec bl, arma::vec bu, arma::mat A, double R, 
                          double r, int n, int m, arma::vec stateFinal, 
                          int numSamplesMax){
  // generate a list of valid samples
  arma::mat samples;
  if (algorithm.compare("rejection-sample")==0){
    //std::cout << "rejection-sample" << std::endl;
    samples = MathTools::uniformRejectionSample(numSamplesMax, xl, xu, bl, 
                                                bu, A);
  }
  else if (algorithm.compare("barycentric")==0){
    //std::cout << "barycentric" << std::endl;
    samples = MathTools::barycentricSample(numSamplesMax, xl, xu, bu, A);  
  }
  //std::cout << " sample complete" << std::endl;
  arma::mat endpoints(samples.n_rows,3);
  arma::vec sampleNormError(samples.n_rows);
  // cycle through all samples
  for (int i = 0; i < samples.n_rows; i++){
    // temporary colvec for holding params of the current sample
    arma::vec params = samples.row(i).t();
    // compute the endpoint corresponding to the current sample
    arma::vec endpoint = VarSpeedDubins::pathEndpoint(params, pathClass, 
                                                      pathType, orientation, 
                                                      R, r); 
    endpoints.row(i) = endpoint.t();
    // compute normalized error of current endpoint
    sampleNormError(i) = normErrorToEndpoint(endpoint, stateFinal);
  }

  // check for lowest cost endpoint
  arma::uword index;
  double min_val = sampleNormError.min(index);
  arma::vec bestGuessParam = samples.row(index).t();
  return bestGuessParam;
}

// normErrorToEndpoint
double VarSpeedDubins::normErrorToEndpoint(arma::vec xTest, 
                                           arma::vec xTrue){
  double Rtrue = sqrt(xTrue(0)*xTrue(0) + xTrue(1)*xTrue(1));
  double xerror = (xTest(0) - xTrue(0))*(xTest(0) - xTrue(0));
  double yerror = (xTest(1) - xTrue(1))*(xTest(1) - xTrue(1));
  double herror = abs(xTest(2) - xTrue(2));
  return sqrt(xerror + yerror)/Rtrue + herror/(2.0*M_PI);
}
