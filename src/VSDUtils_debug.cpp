#include<VSDUtils.h>


// randomParam
arma::vec VarSpeedDubins::randomParam(std::string pathClass, 
                                      std::string pathType, double R, double r, 
                                      double Lmax){
	double asubopt_pi = VarSpeedDubins::alphaSubopt(M_PI, R, r);
  // get constraints
  arma::vec x_l, x_u;
  arma::vec g_l, g_u;
  arma::mat A;
  VarSpeedDubins::linearConstraints(pathClass, pathType, x_l, x_u, g_l, g_u, 
                                    A, R, r, asubopt_pi, Lmax);
  bool paramAdmissible = false;
  arma::vec rv(4);
  // generate a random vector, check if admissible, continue until 
  // an admissible random vector is found
  rv = rv.randu() % x_u; // multiple element-wise to scale
  while (paramAdmissible == false){
    rv = rv.randu() % x_u;
    // check if sample is within the polytope
    // return true if sample is valid
    bool cond1 = MathTools::checkVectorBounds(A*rv, g_l, g_u); 
    // check if suboptimality conditions are violated
    // return true if sample is suboptimal (invalid)
    bool cond2 = VarSpeedDubins::checkSuboptimality(rv, pathClass, pathType, 
                                                    R, r);
    paramAdmissible = (cond1==true && cond2==false);
  }
  return rv;
}
