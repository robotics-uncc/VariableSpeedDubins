// VCSDUtils_suboptimalityConditions.cpp
#include<VSDUtils.h>

// checkSuboptimality
bool VarSpeedDubins::checkSuboptimality(arma::vec paramsShort, 
		                                         std::string pathClass, 
                                             std::string pathType, 
                                             double R, double r){
  // check both conditions
  //bool conditionBTB = VarSpeedDubins::checkBTBSuboptimality(paramsShort, 
  //                                                          pathClass, pathType, 
  //                                                          R, r);
  bool conditionAlpha = VarSpeedDubins::checkAlphaSuboptimality(paramsShort, 
                                                            pathClass, pathType, 
                                                            R, r);
//  std::cout << "conditionBTB" << conditionBTB << std::endl;
//  std::cout << "conditionAlpha" << conditionAlpha << std::endl;
  // if either condition is true, then the candidate is suboptimal and
  // a true boolean is returned
  //return (conditionBTB || conditionAlpha);
  return (conditionAlpha);
}

bool VarSpeedDubins::checkBTBSuboptimality(arma::vec paramsShort, 
		                                         std::string pathClass, 
                                             std::string pathType, 
                                             double R, double r){
  // only applies to three or four turns
  if ( (pathClass.compare("TT") == 0) || (pathClass.compare("TST") == 0) ) {
    return false;
  }
 	// expand paramsShort into long form
	arma::vec paramsLong = VarSpeedDubins::expandParamsShort(paramsShort, 
                                                           pathClass, pathType);
 	double a1 = paramsLong(0);
	double b1 = paramsLong(1);
	double g1 = paramsLong(2);
	//double L  = paramsLong(3);
	double a2 = paramsLong(4);
	double b2 = paramsLong(5);
	double g2 = paramsLong(6);
	double a3 = paramsLong(7);
	double b3 = paramsLong(8);
	double g3 = paramsLong(9);
	double a4 = paramsLong(10);
	double b4 = paramsLong(11);
	double g4 = paramsLong(12);

  // specify alpha (a) and beta (b) as appropriate to pathClass
  double a, b;  
  if ( (pathClass.compare("TTTT") == 0) ) {
    if (pathType.compare("B-TT-B") == 0){
      //std::cout << "g1 : " << g1 << std::endl;
      //std::cout << "a4 : " << a4 << std::endl;
      a = std::max(g1,a4);
      b = b2;
      //std::cout << "a : " << a << std::endl;
      //std::cout << "b : " << b << std::endl;
    }
  }
  else if ( (pathClass.compare("TTT") == 0) ) {
    if (pathType.compare("B-T-BCB") == 0){
      a = g1;
      b = b2;
    }
    else if (pathType.compare("B-T-BC") == 0){
      a = g1;
      b = b2;
    }
    else if (pathType.compare("B-T-B") == 0){
      a = std::min(g1,a3);
      b = b2;
    }
    else if (pathType.compare("BCB-T-B") == 0){
      a = a3;
      b = b2;
    }
    else if (pathType.compare("CB-T-B") == 0){
      a = a3;
      b = b2;
    }
  }

  // required definitions
	double pi = M_PI;
	double delx = 2.0*((R-r)*sin(b/2.0)-R*sin(a+b/2)); // actuall "d"
	double xc = 2.0*(R-r)*sin(b/2.0);

  // th1 and th2 must be defined using the smallest angle (set above)
	double th1 = a + b/2.0 - pi/2.0;
	double th2 = b/2.0 - pi/2.0;
  //std::cout << "th1 : " << th1 << std::endl;
  //std::cout << "th2 : " << th2 << std::endl;

	// Case A
  bool casea = false;
  bool caseb = false;
	if ( (0.0 <= abs(xc)) && (abs(xc) <= 2.0*(R+r)) ){
		double th_a_min = acos(delx/(2.0*(R+r)));
    double th_a_max = acos(delx/(4.0*(R)));
    //std::cout << "th_a_min : " << th_a_min << std::endl;
    //std::cout << "th_a_max : " << th_a_max << std::endl;
    // if the interval [th2, th1] intersects [th_a_min, th_a_max] then
    // the Case A arc exists "above" the circles
    bool casea1 = MathTools::checkIntervalIntersect(th2,th1,th_a_min,th_a_max);
    // if the interval [th2, th1] intersects [-th_a_max, -th_a_min] then
    // the Case A arc exists "above" the circles
    bool casea2 = MathTools::checkIntervalIntersect(th2,th1,-th_a_max,-th_a_min);
    casea = (casea1 || casea2); // if either is true then arc exists
	}
  
	// Case B
	if ( (0.0 <= abs(xc)) && (abs(xc) <= 2.0*(R-r)) ){
    // the Case B arc exists on the interval [-th_b_max, th_b_max]
    double th_b_max = b/2.0 - pi/2.0;
//    std::cout << "th_b_max : " << th_b_max << std::endl;
    // if the interval [th2, th1] intersects [-th_b_max, th_b_max] then
    // the Case B arc exists
    caseb = MathTools::checkIntervalIntersect(th2,th1,-th_b_max,th_b_max);
  }

  //std::cout << " case a " << casea << std::endl;
  //std::cout << " case b " << caseb << std::endl;

  return (casea || caseb);
}


//// checkBTBSuboptimality
//bool VarSpeedDubins::checkBTBSuboptimality(arma::vec paramsShort, 
//		                                       std::string pathClass, 
//                                           std::string pathType){
// 	// expand paramsShort into long form
//	arma::vec paramsLong = VarSpeedDubins::expandParamsShort(paramsShort, 
//                                                           pathClass, pathType);
// 	double a1 = paramsLong(0);
//	double b1 = paramsLong(1);
//	double g1 = paramsLong(2);
//	//double L  = paramsLong(3);
//	double a2 = paramsLong(4);
//	double b2 = paramsLong(5);
//	double g2 = paramsLong(6);
//	double a3 = paramsLong(7);
//	double b3 = paramsLong(8);
//	double g3 = paramsLong(9);
//	double a4 = paramsLong(10);
//	double b4 = paramsLong(11);
//	double g4 = paramsLong(12);
//  // check conditions
//  if ( (pathClass.compare("TTT") == 0) ) {
//    if (pathType.compare("B-BCB-BCB-B") == 0){
//      if (std::min(g1,a4) >= a2){
//        return true;
//      }
//    }
//  }
//  else if ( (pathClass.compare("TTT") == 0) ) {
//    if (pathType.compare("B-BCB-BCB") == 0){
//      if (g1 >= a2){
//        return true;
//      }
//    }
//    else if (pathType.compare("B-BCB-BC") == 0){
//      if (g1 >= a2){
//        return true;
//      }
//    }
//    else if (pathType.compare("B-BCB-B") == 0){
//      if (std::min(g1,a3) >= a2){
//        return true;
//      }
//    }
//    else if (pathType.compare("BCB-BCB-B") == 0){
//      if (a3 >= g2){
//        return true;
//      }
//    }
//    else if (pathType.compare("CB-BCB-B") == 0){
//      if (a3 >= g2){
//        return true;
//      }
//    }
//  }
//  return false;
//}

// checkAlphaSuboptimality
bool VarSpeedDubins::checkAlphaSuboptimality(arma::vec paramsShort, 
		                                         std::string pathClass, 
                                             std::string pathType, 
                                             double R, double r){
	// expand paramsShort into long form
	arma::vec paramsLong = VarSpeedDubins::expandParamsShort(paramsShort, 
                                                           pathClass, pathType);
  // expand long form into individual parameters of a BCB-S-BCB curve
	double a1 = paramsLong(0);
	double b1 = paramsLong(1);
	double g1 = paramsLong(2);
  // paramsLong(3) always corresponds to L even if it is unused 
	double a2 = paramsLong(4);
	double b2 = paramsLong(5);
	double g2 = paramsLong(6);
  // intiialize variables
	double alpha, beta, asubopt, asubopt1, asubopt2;
  // compare each path class
	if ( (pathClass.compare("TT") == 0) ) {
		// check first turn
		alpha = std::min(a1,g1);
		beta = b1;
		asubopt1 = VarSpeedDubins::alphaSubopt(beta,R,r);
		// check second turn
		alpha = std::min(a2,g2);
		beta = b2;
		asubopt2 = VarSpeedDubins::alphaSubopt(beta,R,r);
		return ( (alpha >= asubopt1) || (alpha >= asubopt2) );
	}	
	else if ( (pathClass.compare("TTT") == 0) || (pathClass.compare("TTTT") == 0) ){
    // the boundary turns are constrained to be shorter than the internal turns
    // thus it is sufficient to just check an internal turn segment for 
    // suboptimality
		beta = b2;
		alpha = a2;
		asubopt = VarSpeedDubins::alphaSubopt(beta,R,r);
		return (alpha >= asubopt);
	}
	else if (pathClass.compare("TST") == 0){
		// check first turn
		alpha = std::min(a1,g1);
		beta = b1;
		asubopt1 = VarSpeedDubins::alphaSubopt(beta,R,r);
		// check second turn
		alpha = std::min(a2,g2);
		beta = b2;
		asubopt2 = VarSpeedDubins::alphaSubopt(beta,R,r);
		return ( (alpha >= asubopt1) || (alpha >= asubopt2) );
	}	
}

// alphaSubopt
double VarSpeedDubins::alphaSubopt(double beta, double R, double r){
	return M_PI - beta/2.0 -asin(sin(beta/2.0)*(R-r)/(R+r));
}
