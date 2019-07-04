// VCSDUtils_linearConstraints.cpp
// Last Modefied: 13-Feb-2016, Artur Wolek

#include<VSDUtils.h>

void VarSpeedDubins::linearConstraints(std::string pathClass, 
		                                   std::string pathType, arma::vec &x_l, 
                                       arma::vec &x_u, arma::vec &g_l, 
		                                   arma::vec &g_u, arma::mat &A, double R, 
                                       double r, double asubopt_M_PI, 
                                       double Lmax){
	arma::colvec b;
	if ( pathClass.compare("TST") == 0 ){
			if (pathType.compare("BCB-S-BCB") == 0){
   	     		// p = [a1 g1 L g2];
				// p(i) > 0 for i = 1, to 4 (an implicit constraint)
				// 0 <= a1 <= asubopt_M_PI
				// 0 <= g1 <= 2*M_PI
				// 0 <= L  <= Lmax
				// 0 <= g2 <= subopt_M_PI
        A = { {1, -1, 0, 0},            	
       			  {0, -1, 0, 1} };
				x_l = {0, 0, 0, 0};
				x_u = {asubopt_M_PI, M_PI/2.0, Lmax, asubopt_M_PI};
				g_l = {-M_PI, -M_PI};
        		g_u = {0, 0};
				
			}
			else if (pathType.compare("BCB-S-BC") == 0){
   	     		// p = [a1 g1 L b2];
				x_l = {0, 0, 0, 0};
				x_u = {asubopt_M_PI, M_PI/2.0, Lmax, M_PI};
        A = { {1, -1, 0, 0} };
				g_l = {-M_PI};
        g_u = {0};
			}		
			else if (pathType.compare("BCB-S-B") == 0){
   	     		// p = [a1 g1 L a2];
				x_l = {0, 0, 0, 0};
				x_u = {asubopt_M_PI, M_PI/2.0, Lmax, M_PI};
          A = { {1, -1, 0, 0},
                {0, -1, 0, 1} };
				g_l = {-M_PI, -M_PI};
        		g_u = {0, 0};
			}		
			else if (pathType.compare("CB-S-BCB") == 0){
        // p = [b1 g1 L g2];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, M_PI/2.0, Lmax, asubopt_M_PI};
          A = { {0, -1, 0, 1}};
				g_l = {-M_PI};
        		g_u = {0};
			}					
			else if (pathType.compare("CB-S-BC") == 0){
   	     		// p = [b1 g1 L b2];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, M_PI, Lmax, M_PI};
          A = { {1, 2, 0, 1},
                {0, 2, 0, 1}};
				g_l = {0, 0};
        		g_u = {2*M_PI, 2*M_PI};
			}				
			else if (pathType.compare("CB-S-B") == 0){
   	     		// p = [b1 g1 L a2];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, M_PI, Lmax, M_PI};
          A = { {1, 2, 0, 1}, 
					      {0, -1, 0, 1}};
				g_l = {0, -M_PI};
        		g_u = {2.0*M_PI, 0};
			}		
			else if (pathType.compare("B-S-BCB") == 0){
   	     		// p = [g1 L a2 g2];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, Lmax, M_PI, asubopt_M_PI};
          A = { {1, 0, -1, 0},
                {0, 0, -1, 1}}; 
				g_l = {-M_PI, -M_PI};
        		g_u = {0, 0};
			}		
			else if (pathType.compare("B-S-BC") == 0){
   	     		// p = [g1 L a2 b2];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, Lmax, M_PI, M_PI};
          A = { {1, 0, -1, 0},
                {0, 0, 2, 1}};
				g_l = {-M_PI, 0};
        		g_u = {0, 2.0*M_PI};
			}								
	}
	else if ( (pathClass.compare("TT") == 0) ){
			if (pathType.compare("BCB-BCB") == 0){
				// p = [a1 b1 g1 g2];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, 2.0*M_PI, M_PI, M_PI};
          A = { {1,  0, -1, 0},
                {0,  0,-1, 1},
                {0,  1, 2, 0} };
				g_l = {-M_PI, -M_PI, 0};
        		g_u = {0, 0, 2.0*M_PI};
			}
			else if (pathType.compare("BCB-BC") == 0){
				// p = [a1 b1 g1 b2];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, 2.0*M_PI, M_PI, 2.0*M_PI};
        		A = { {1,  0, -1, 0},
            		  {0, -1, 0, 1},
            		  {0,  1, 2, 0} };
				g_l = {-M_PI, -2.0*M_PI, 0};
        		g_u = {0, 0, 2.0*M_PI};
			}			
			else if (pathType.compare("BCB-B") == 0){
				// p = [a1 b1 g1 a2];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, 2.0*M_PI, M_PI, M_PI};
        		A = { {1,  0, -1, 0},
            		  {0,  0, -1, 1},
            		  {0,  1, 2, 0} };
				g_l = {-M_PI, -M_PI, 0};
        		g_u = {0, 0, 2.0*M_PI};
			}			
			else if (pathType.compare("CB-BCB") == 0){
				// p = [b1 g1 b2 g2];
				x_l = {0, 0, 0, 0};
				x_u = {2.0*M_PI, M_PI, 2.0*M_PI, M_PI};
        A = { {1,  0, -1, 0},
              {0,  -1, 0, 1},
              {0,  2,  1, 0} };
				g_l = {-2.0*M_PI, -M_PI, 0};
        g_u = {0, 0, 2.0*M_PI};
			}	
			else if (pathType.compare("B-BCB") == 0){
				// p = [g1 a2 b2 g2];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, Lmax, M_PI, M_PI};
        		A = { {1,  -1, 0, 0},
            		  {0,  -1, 0, 1},
            		  {0,  2, 1, 0} };
				g_l = {-M_PI, -M_PI, 0};
        		g_u = {0, 0, 2.0*M_PI};
			}			
	}	
	else if ( (pathClass.compare("TTTT") == 0)){
		if (pathType.compare("CB-TT-BC") == 0){
				// p = [b1 g1 b2 b4];
				x_l = {0, 0, 0, 0};
				x_u = {2.0*M_PI, M_PI, 2.0*M_PI, 2.0*M_PI};
        		A = { {1, 0, -1, 0},
            		  {0, 0, -1, 1},
            		  {0, 2, 1,  0} };
				g_l = {-2.0*M_PI, -2.0*M_PI, 0};
        g_u = {0, 0, 2.0*M_PI};
		}		
		else if (pathType.compare("CB-TT-B") == 0){
				// p = [b1 g1 b2 a4];
				x_l = {0, 0, 0, 0};
				x_u = {2.0*M_PI, M_PI, 2.0*M_PI, M_PI};
        		A = { {1, 0, -1, 0},
            		  {0, -1, 0, 1},
            		  {0, 2, 1,  0} };
				g_l = {-2.0*M_PI, -M_PI, 0};
        g_u = {0, 0, 2.0*M_PI};
		}	
		else if (pathType.compare("B-TT-BC") == 0){
				// p = [g1 a2 b2 b4];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, M_PI, 2.0*M_PI, 2.0*M_PI};
        		A = { {1, -1, 0, 0},
            		  {0,  0, -1, 1},
            		  {0,  2, 1, 0} };
				g_l = {-M_PI, -2.0*M_PI, 0};
        g_u = {0, 0, 2.0*M_PI};
		}	
		else if (pathType.compare("B-TT-B") == 0){
				// p = [g1 a2 b2 a4];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, M_PI, 2.0*M_PI, M_PI};
        		A = { {1, -1, 0, 0},
            		  {0, -1, 0, 1},
            		  {0,  2, 1, 0} };
				g_l = {-M_PI, -M_PI, 0};
        g_u = {0, 0, 2.0*M_PI};
		}																	
	}	
	else if ( (pathClass.compare("TTT") == 0) ){
			if ( pathType.compare("CB-T-BCB") == 0 ){
				// p = [b1 g1 b2 g3];
				x_l = {0, 0, 0, 0};
				x_u = {2.0*M_PI, M_PI, 2.0*M_PI, M_PI};
        		A = { {1,  0, -1, 0},
            		  {0, -1, 0, 1},
            		  {0,  2, 1, 0} };
				g_l = {-2.0*M_PI, -M_PI, 0};
        g_u = {0, 0, 2.0*M_PI};
			}		
			else if ( pathType.compare("CB-T-BC") == 0 ){
				// p = [b1 g1 b2 b3];
				x_l = {0, 0, 0, 0};
				x_u = {2.0*M_PI, M_PI, 2.0*M_PI, 2.0*M_PI};
        		A = { {1,  0, -1, 0},
            		  {0,  0, -1, 1},
            		  {0,  2, 1, 0} };
				g_l = {-2.0*M_PI, -2.0*M_PI, 0};
        g_u = {0, 0, 2.0*M_PI};
			}		
			else if (pathType.compare("CB-T-B") == 0){
				// p = [b1 g1 b2 a3];
				x_l = {0, 0, 0, 0};
				x_u = {2.0*M_PI, M_PI, 2.0*M_PI, M_PI};
        		A = { {1,  0, -1, 0},
            		  {0,  -1, 0, 1},
            		  {0,  2,  1, 0} };
				g_l = {-2.0*M_PI, -M_PI, 0};
        		g_u = {0, 0, 2.0*M_PI};
			}	
			else if ( pathType.compare("B-T-BCB") == 0 ){
				// p = [g1 a2 b2 g3];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, M_PI, 2.0*M_PI, M_PI};
        		A = { {1, -1, 0, 0},
            		  {0, -1, 0, 1},
            		  {0,  2, 1, 0} };
				g_l = {-M_PI, -M_PI, 0};
        		g_u = {0, 0, 2.0*M_PI};
			}		
			else if (pathType.compare("B-T-BC") == 0){
				// p = [g1 a2 b2 b3];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, M_PI, 2.0*M_PI, 2.0*M_PI};
        		A = { {1, -1, 0, 0},
            		  {0,  0, -1, 1},
            		  {0,  2, 1, 0} };
				g_l = {-M_PI, -2.0*M_PI, 0};
        		g_u = {0, 0, 2.0*M_PI};
			}	
			else if (pathType.compare("B-T-B") == 0){
				// p = [g1 a2 b2 a3];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, M_PI, 2.0*M_PI, M_PI};
        		A = { {1, -1, 0, 0},
            		  {0, -1, 0, 1},
            		  {0,  2, 1, 0} };
				g_l = {-M_PI, -M_PI, 0};
        		g_u = {0, 0, 2.0*M_PI};
			}			
			else if (pathType.compare("BCB-T-BC") == 0){
				// p = [a2 b2 g2 b4];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, 2.0*M_PI, M_PI, 2.0*M_PI};
        		A = { {1, -1, 0, 0},
            		  {0, -1, 0, 1},
            		  {0,  2, 1, 0} };
				g_l = {-M_PI, -2.0*M_PI, 0};
        g_u = {0, 0, 2.0*M_PI};
			}			
			else if (pathType.compare("BCB-T-B") == 0){
				// p = [a2 b2 g2 a4];
				x_l = {0, 0, 0, 0};
				x_u = {M_PI, 2.0*M_PI, M_PI, M_PI};
        		A = { {1,  0, -1, 0},
            		  {0,  0, -1, 1},
            		  {0,  2, 1, 0} };
				g_l = {-M_PI, -M_PI, 0};
        		g_u = {0, 0, 2.0*M_PI};
			}																				
	}	
}

