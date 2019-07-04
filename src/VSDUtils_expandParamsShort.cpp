// VCSDUtils_expandParamsShort.cpp
// Last Modefied: 18-Jan-2016, Artur Wolek

#include<VSDUtils.h>

arma::vec VarSpeedDubins::expandParamsShort(arma::vec p, std::string pathClass, std::string pathType){
	// the input is a (short) list of at most four parameters
	// we want to define the following exhaustive list of parameters
	double a1 = 0;
	double b1 = 0;
	double g1 = 0;
	double L = 0;
	double a2 = 0;
	double b2 = 0;
	double g2 = 0;
	double a3 = 0;
	double b3 = 0;
	double g3 = 0;
	double a4 = 0;
	double b4 = 0;
	double g4 = 0;
	// for convenience 
	double pi = M_PI;
	// check pathClass
	if ( pathClass.compare("TST") == 0 ){
			// check pathType
			if (pathType.compare("BCB-S-BCB") == 0){
				// change non-zero parameters
				a1 = p(0);
				b1 = pi;
				g1 = p(1);
				L = p(2);
				a2 = g1;
				b2 = pi;
				g2 = p(3);
			}
			else if (pathType.compare("BCB-S-BC") == 0){
				a1 = p(0);
				b1 = pi;
				g1 = p(1);
				L = p(2);
				a2 = g1;
				b2 = p(3);
				g2 = 0;
			}		
			else if (pathType.compare("BCB-S-B") == 0){
				a1 = p(0);
				b1 = pi;
				g1 = p(1);
				L = p(2);
				a2 = p(3);
				b2 = 0;
				g2 = 0;
			}		
			else if (pathType.compare("CB-S-BCB") == 0){
				a1 = 0;
				b1 = p(0);
				g1 = p(1);
				L = p(2);
				a2 = g1;
				b2 = pi;
				g2 = p(3);
			}					
			else if (pathType.compare("CB-S-BC") == 0){
				a1 = 0;
				b1 = p(0);
				g1 = p(1);
				L = p(2);
				a2 = g1;
				b2 = p(3);
				g2 = 0;
			}				
			else if (pathType.compare("CB-S-B") == 0){
				a1 = 0;
				b1 = p(0);
				g1 = p(1);
				L = p(2);
				a2 = p(3);
				b2 = 0;
				g2 = 0;
			}		
			else if (pathType.compare("B-S-BCB") == 0){
				a1 = 0;
				b1 = 0;
				g1 = p(0);
				L = p(1);
				a2 = p(2);
				b2 = pi;
				g2 = p(3);
			}		
			else if (pathType.compare("B-S-BC") == 0){
				a1 = 0;
				b1 = 0;
				g1 = p(0);
				L = p(1);
				a2 = p(2);
				b2 = p(3);
				g2 = 0;
			}								
	}
	else if ( (pathClass.compare("TT") == 0) ){
			if (pathType.compare("BCB-BCB") == 0){
				a1 = p(0);
				b1 = p(1);
				g1 = p(2);
				a2 = g1;
				b2 = b1;
				g2 = p(3);
			}
			else if (pathType.compare("BCB-BC") == 0){
				a1 = p(0);
				b1 = p(1);
				g1 = p(2);
				a2 = g1;
				b2 = p(3);
				g2 = 0;
			}			
			else if (pathType.compare("BCB-B") == 0){
				a1 = p(0);
				b1 = p(1);
				g1 = p(2);
				a2 = p(3);
				b2 = 0;
				g2 = 0;
			}				
			else if (pathType.compare("CB-BCB") == 0){
				a1 = 0;
				b1 = p(0);
				g1 = p(1);
				a2 = g1;
				b2 = p(2);
				g2 = p(3);
			}	
			else if (pathType.compare("CB-BC") == 0){
				a1 = 0;
				b1 = p(0);
				g1 = p(1);
				a2 = g1;
				b2 = p(2);
				g2 = 0;
			}		
			else if (pathType.compare("CB-B") == 0){
				a1 = 0;
				b1 = p(0);
				g1 = p(1);
				a2 = p(2);
				b2 = 0;
				g2 = 0;
			}	
			else if (pathType.compare("B-BCB") == 0){
				a1 = 0;
				b1 = 0;
				g1 = p(0);
				a2 = p(1);
				b2 = p(2);
				g2 = p(3);
			}	
			else if (pathType.compare("B-BC") == 0){
				a1 = 0;
				b1 = 0;
				g1 = p(0);
				a2 = p(1);
				b2 = p(2);
				g2 = 0;  
			}												
	}
	else if ( (pathClass.compare("TTT") == 0)){
			if (pathType.compare("BCB-T-BCB") == 0){
				a1 = p(0);
				b1 = p(1);
				g1 = p(2);
				a2 = g1;
				b2 = b1;
				g2 = a2;
				a3 = a2;
				b3 = b1;
				g3 = p(3);
			}	
			else if (pathType.compare("BCB-T-BC") == 0){
				a1 = p(0);
				b1 = p(1);
				g1 = p(2);
				a2 = g1;
				b2 = b1;
				g2 = a2;
				a3 = a2;
				b3 = p(3);
				g3 = 0;
			}	
			else if (pathType.compare("BCB-T-B") == 0){
				a1 = p(0);
				b1 = p(1);
				g1 = p(2);
				a2 = g1;
				b2 = b1;
				g2 = a2;
				a3 = p(3);
				b3 = 0;
				g3 = 0;
			}	
			else if (pathType.compare("CB-T-BCB") == 0){
				a1 = 0;
				b1 = p(0);
				g1 = p(1);
				a2 = g1;
				b2 = p(2);
				g2 = a2;
				a3 = a2;
				b3 = b2;
				g3 = p(3);
			}	
			else if (pathType.compare("CB-T-BC") == 0){
				a1 = 0;
				b1 = p(0);
				g1 = p(1);
				a2 = g1;
				b2 = p(2);
				g2 = a2;
				a3 = a2;
				b3 = p(3);
				g3 = 0;
			}		
			else if (pathType.compare("CB-T-B") == 0){
				a1 = 0;
				b1 = p(0);
				g1 = p(1);
				a2 = g1;
				b2 = p(2);
				g2 = a2;
				a3 = p(3);
				b3 = 0;
				g3 = 0;
			}		
			else if (pathType.compare("B-T-BCB") == 0){
				a1 = 0;
				b1 = 0;
				g1 = p(0);
				a2 = p(1);
				b2 = p(2);
				g2 = a2;
				a3 = a2;
				b3 = b2;
				g3 = p(3);
			}		
			else if (pathType.compare("B-T-BC") == 0){
			    a1 = 0;
				b1 = 0;
				g1 = p(0);
				a2 = p(1);
				b2 = p(2);
				g2 = a2;
				a3 = a2;
				b3 = p(3);
				g3 = 0;
			}	
			else if (pathType.compare("B-T-B") == 0){
				a1 = 0;
				b1 = 0;
				g1 = p(0);
				a2 = p(1);
				b2 = p(2);
				g2 = a2;
				a3 = p(3);
				b3 = 0;
				g3 = 0;
			}																
	}
	else if ( (pathClass.compare("TTTT") == 0)){
			if (pathType.compare("BCB-TT-BCB") == 0){
				a1 = p(0);
				b1 = p(1);
				g1 = p(2);
				a2 = g1;
				b2 = b1;
				g2 = g1;
				a4 = g1;
				b4 = b1;
				g4 = p(3);
			}	
			else if (pathType.compare("BCB-TT-BC") == 0){
				a1 = p(0);
				b1 = p(1);
				g1 = p(2);
				a2 = g1;
				b2 = b1;
				g2 = a2;
				a4 = g1;
				b4 = p(3);
				g4 = 0;
			}		
			else if (pathType.compare("BCB-TT-B") == 0){
				a1 = p(0);
				b1 = p(1);
				g1 = p(2);
				a2 = g1;
				b2 = b1;
				g2 = a2;
				a4 = p(3);
				b4 = 0;
				g4 = 0;
			}		
			else if (pathType.compare("CB-TT-BCB") == 0){
				a1 = 0;
				b1 = p(0);
				g1 = p(1);
				a2 = g1;
				b2 = p(2);
				g2 = a2;
				a4 = g1;
				b4 = b2;
				g4 = p(3);
			}		
			else if (pathType.compare("CB-TT-BC") == 0){
				a1 = 0;
				b1 = p(0);
				g1 = p(1);
				a2 = g1;
				b2 = p(2);
				g2 = a2;
				a4 = g1;
				b4 = p(3);
				g4 = 0;
			}		
			else if (pathType.compare("CB-TT-B") == 0){
				a1 = 0;
				b1 = p(0);
				g1 = p(1);
				a2 = g1;
				b2 = p(2);
				g2 = a2;
				a4 = p(3);
				b4 = 0;
				g4 = 0;
			}		
			else if (pathType.compare("B-TT-BCB") == 0){
				a1 = 0;
				b1 = 0;
				g1 = p(0);
				a2 = p(1);
				b2 = p(2);
				g2 = a2;
				a4 = a2;
				b4 = b2;
				g4 = p(3);
			}				
			else if (pathType.compare("B-TT-BC") == 0){
				a1 = 0;
				b1 = 0;
				g1 = p(0);
				a2 = p(1);
				b2 = p(2);
				g2 = a2;
				a4 = a2;
				b4 = p(3);
				g4 = 0;
			}		
			else if (pathType.compare("B-TT-B") == 0){
				a1 = 0;
				b1 = 0;
				g1 = p(0);
				a2 = p(1);
				b2 = p(2);
				g2 = a2;
				a4 = p(3);
				b4 = 0;
				g4 = 0;
			}		
			a3 = a2;
    	b3 = b2;
    	g3 = a2;															
	}	
	// group parameters into pLong vector
	arma::vec pLong = {a1,b1,g1,L,a2,b2,g2,a3,b3,g3,a4,b4,g4};
	return pLong;
}



