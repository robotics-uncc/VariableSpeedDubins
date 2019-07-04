// VCSDUtils_costFunction.cpp
// Last Modefied: 18-Jan-2016, Artur Wolek

#include<VSDUtils.h>

// compareCosts
int VarSpeedDubins::compareCosts(double costRef, double costNew, double tol){
  double deltaCost = costRef - costNew;
  //std::cout << "deltaCost : " <<  deltaCost << std::endl;
  //std::cout << "abs(deltaCost) : " << std::abs(deltaCost) << std::endl;
  //std::cout << "tol " << tol << std::endl;
  if ( std::abs(deltaCost) <= tol){
    return 2; // costNew is equal to costRef
  }
  else if ( deltaCost > tol ){
    return 1; // costNew is lower than costRef 
  }
  else {
    return 0; // costNew is higher than costRef
  }
}

double VarSpeedDubins::costFunction(arma::vec paramsShort, 
		                                std::string pathClass, std::string pathType, 
                                    std::string orientation, double R){
	arma::vec paramsLong = VarSpeedDubins::expandParamsShort(paramsShort, 
		                                                       pathClass, pathType);
	double a1 = paramsLong(0);
	double b1 = paramsLong(1);
	double g1 = paramsLong(2);
	double L = paramsLong(3);
	double a2 = paramsLong(4);
	double b2 = paramsLong(5);
	double g2 = paramsLong(6);
	double a3 = paramsLong(7);
	double b3 = paramsLong(8);
	double g3 = paramsLong(9);
	double a4 = paramsLong(10);
	double b4 = paramsLong(11);
	double g4 = paramsLong(12);


	double cost;
	if (pathClass.compare("TST") == 0){
		cost = R*(a1+b1+g1+a2+b2+g2)+L;
	}
	else if (pathClass.compare("TT") == 0){
    cost = R*(a1+b1+g1+a2+b2+g2);
	}
	else if (pathClass.compare("TTT") == 0){
		cost = R*(a1+b1+g1+a2+b2+g2+a3+b3+g3+a4+b4+g4);
	}
	else if (pathClass.compare("TTTT") == 0){
		cost = R*(a1+b1+g1+a2+b2+g2+a3+b3+g3+a4+b4+g4);
	}
	return cost;
}


arma::vec VarSpeedDubins::costFunctionGradient(double R, std::string pathClass, 
		std::string pathType, std::string orientation){
	arma::vec gradient;
  #if DEBUG_LEVEL > 3
    std::cout << "VarSpeedDubins::costFunctionGradient" << std::endl;
  #endif 
	if ( pathClass.compare("TST") == 0 ){
			if (pathType.compare("BCB-S-BCB") == 0){
   	     		// p = [a1 g1 L g2];
   	     		gradient = { R, 2.0*R, 1.0, R};
			}
			else if (pathType.compare("BCB-S-BC") == 0){
   	     		// p = [a1 g1 L b2];
   	     		gradient = { R, 2.0*R, 1.0, R};	
			}		
			else if (pathType.compare("BCB-S-B") == 0){
   	     		// p = [a1 g1 L a2];
   	     		gradient = { R, R, 1.0, R};	
			}		
			else if (pathType.compare("CB-S-BCB") == 0){
   	     		// p = [b1 g1 L g2];
   	     		gradient = { R, 2.0*R, 1.0, R};	
			}					
			else if (pathType.compare("CB-S-BC") == 0){
   	     		// p = [b1 g1 L b2];
   	     		gradient = { R, 2.0*R, 1.0, R};	
			}				
			else if (pathType.compare("CB-S-B") == 0){
   	     		// p = [b1 g1 L a2];
   	     		gradient = { R, R, 1.0, R};	
			}		
			else if (pathType.compare("B-S-BCB") == 0){
   	     		// p = [g1 L a2 g2];
   	     		gradient = { R, 1.0, R, R};	
			}		
			else if (pathType.compare("B-S-BC") == 0){
   	     		// p = [g1 L a2 b2];
   	     		gradient = { R, 1.0, R, R};	
			}								
	}
	else if ( (pathClass.compare("TT") == 0) ){
			if (pathType.compare("BCB-BCB") == 0){
   	     		gradient = { R, 2.0*R, 2.0*R, R};
			}
			else if (pathType.compare("BCB-BC") == 0){
   	     		gradient = { R, R, 2.0*R, R};
			}			
			else if (pathType.compare("BCB-B") == 0){
   	     		gradient = { R, R, R, R};
			}			
			else if (pathType.compare("CB-BCB") == 0){
   	     		gradient = { R, 2.0*R, R, R};
			}	
			else if (pathType.compare("CB-BC") == 0){
   	     		gradient = { R, 2.0*R, R};
			}		
			else if (pathType.compare("CB-B") == 0){
   	     		gradient = { R, R, R};
			}	
			else if (pathType.compare("B-BCB") == 0){
   	     		gradient = { R, R, R, R};
			}	
			else if (pathType.compare("B-BC") == 0){
   	     		gradient = { R, R, R};
			}											
	}
	else if ( (pathClass.compare("TTT") == 0)){
			if (pathType.compare("BCB-T-BCB") == 0){
   	     		gradient = { R, 3.0*R, 4.0*R, R};
			}	
			else if (pathType.compare("BCB-T-BC") == 0){
   	     		gradient = { R, 2.0*R, 4.0*R, R};
			}	
			else if (pathType.compare("BCB-T-B") == 0){
   	     		gradient = { R, 2.0*R, 3.0*R, R};
			}	
			else if (pathType.compare("CB-T-BCB") == 0){
   	     		gradient = { R, 4.0*R, 2.0*R, R};
			}	
			else if (pathType.compare("CB-T-BC") == 0){
   	     		gradient = { R, 4.0*R, R, R};
			}		
			else if (pathType.compare("CB-T-B") == 0){
   	     		gradient = { R, 3.0*R, R, R};
			}		
			else if (pathType.compare("B-T-BCB") == 0){
   	     		gradient = { R, 3.0*R, 2.0*R, R};
			}		
			else if (pathType.compare("B-T-BC") == 0){
   	     		gradient = { R, 3.0*R, R, R};
			}	
			else if (pathType.compare("B-T-B") == 0){
   	     		gradient = { R, 2.0*R, R, R};
			}																
	}
	else if ( (pathClass.compare("TTTT") == 0)){
			if (pathType.compare("CB-TT-BC") == 0){
   	     		gradient = { R, 6.0*R, 2.0*R, R};
			}				
			else if (pathType.compare("CB-TT-B") == 0){
   	     		gradient = { R, 5.0*R, 2.0*R, R};
			}	
			else if (pathType.compare("B-TT-BC") == 0){
   	     		gradient = { R, 5.0*R, 2.0*R, R};
			}	
			else if (pathType.compare("B-TT-B") == 0){
   	     		gradient = { R, 4.0*R, 2.0*R, R};
			}														
	}	
	return gradient;
}


