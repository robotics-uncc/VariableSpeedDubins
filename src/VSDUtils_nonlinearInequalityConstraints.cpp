// VSDUtils_nonlinearInequalityConstraints.cpp
// Last Modefied: 18-Jan-2016, Artur Wolek

#include<VSDUtils.h>


double VarSpeedDubins::asubopt_func(double &r, double &R, double &beta){
	return (M_PI - beta/2.0 -asin(sin(beta/2.0)*(R-r)/(R+r)));
}

double VarSpeedDubins::dcineq_dbeta_func(double &r, double &R, double &beta){
	double denom = (1-1*(-r+R)*(-r+R)/((r+R)*(r+R))*sin(0.5*beta)*sin(0.5*beta));
	return (0.5+0.5*(-r+R)/(r+R)*cos(0.5*beta)/(denom*denom));
}

arma::vec VarSpeedDubins::nonlinearInequalityConstraints(arma::vec paramsShort, std::string pathClass, std::string pathType, std::string orientation, double r, double R){
	// returns cineq,  cineq_grad

	// a function that sets k1, and k2
	double k1, k2;
	VarSpeedDubins::determineCurvatureParams(pathClass, orientation, k1, k2);
	// a function that sets these params
	arma::vec paramsLong = VarSpeedDubins::expandParamsShort(paramsShort, pathClass, pathType);
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

	double alpha, beta, asubopt, dcineq_dbeta, denom;
	double pi = M_PI;

	arma::vec cineq;
	arma::mat cineq_grad;

	if ( pathClass.compare("TTTT") == 0 ){
		beta = b2;
		alpha = a2;	
		cineq(0) = alpha - VarSpeedDubins::asubopt_func(r,R,beta);
		// dcineq_dalpha = 1;
		// dcineq_dp1
		cineq_grad(0,0) = 0;
		cineq_grad(0,1) = 1;
		cineq_grad(0,2) = VarSpeedDubins::dcineq_dbeta_func(r, R, beta);
		cineq_grad(0,3) = 0;
	}
	else if ( pathClass.compare("TTT") == 0 ){			
		if ( (pathType.compare("B-T-BCB") == 0) || (pathType.compare("B-T-BC") == 0) || (pathType.compare("B-T-B") == 0) ){
			beta = b2;
			alpha = a2;	
			cineq(0) = alpha - VarSpeedDubins::asubopt_func(r,R,beta);
			cineq_grad(0,0) = 0;
			cineq_grad(0,1) = 1;
			cineq_grad(0,2) = VarSpeedDubins::dcineq_dbeta_func(r, R, beta);
			cineq_grad(0,3) = 0;
		}
		else if ( (pathType.compare("BCB-T-B") == 0)){
			beta = b1;
			alpha = g1;	
			cineq(0) = alpha - VarSpeedDubins::asubopt_func(r,R,beta);
			cineq_grad(0,0) = 0;
			cineq_grad(0,1) = VarSpeedDubins::dcineq_dbeta_func(r, R, beta);
			cineq_grad(0,2) = 1;
			cineq_grad(0,3) = 0;
		}
		else if ( (pathType.compare("CB-T-B") == 0)){
			beta = b2;
			alpha = g1;	
			cineq(0) = alpha - VarSpeedDubins::asubopt_func(r,R,beta);
			cineq_grad(0,0) = 0;
			cineq_grad(0,1) = 1;
			cineq_grad(0,2) = VarSpeedDubins::dcineq_dbeta_func(r, R, beta);
			cineq_grad(0,3) = 0;
		}
	}
	else if ( (pathClass.compare("TT") == 0) || (pathClass.compare("UMAX") == 0) ){		
		if (pathType.compare("B-BCB") == 0){
			beta = b2;
			alpha = g2;	
			cineq(0) = alpha - VarSpeedDubins::asubopt_func(r,R,beta);
			cineq_grad(0,0) = 0;
			cineq_grad(0,1) = 0;
			cineq_grad(0,2) = VarSpeedDubins::dcineq_dbeta_func(r, R, beta);
			cineq_grad(0,3) = 1;
		}
		else if ( (pathType.compare("B-BC") == 0) || (pathType.compare("CB-BC") == 0) || (pathType.compare("CB-B") == 0) ){
			cineq(0) = 0;
			cineq_grad(0,0) = 0;
			cineq_grad(0,1) = 0;
			cineq_grad(0,2) = 0;
		}
		else if (pathType.compare("BCB-BCB") == 0){
			beta = b1;
			alpha = a1;	
			cineq(0) = alpha - VarSpeedDubins::asubopt_func(r,R,beta);
			cineq_grad(0,0) = 0;
			cineq_grad(0,1) = 1;
			cineq_grad(0,2) = VarSpeedDubins::dcineq_dbeta_func(r, R, beta);
			cineq_grad(0,3) = 0;
			//////////////////////////////////////
			beta = b1;
			alpha = g2;	
			cineq(1) = alpha - VarSpeedDubins::asubopt_func(r,R,beta);
			cineq_grad(1,0) = 0;
			cineq_grad(1,1) = 1;
			cineq_grad(1,2) = VarSpeedDubins::dcineq_dbeta_func(r, R, beta);
			cineq_grad(1,3) = 0;
		}
		else if (pathType.compare("CB-BCB") == 0){
			beta = b2;
			alpha = g2;	
			cineq(0) = alpha - VarSpeedDubins::asubopt_func(r,R,beta);
			cineq_grad(0,0) = 0;
			cineq_grad(0,1) = 0;
			cineq_grad(0,2) = VarSpeedDubins::dcineq_dbeta_func(r, R, beta);
			cineq_grad(0,3) = 1;
		}
	}
	else if (pathClass.compare("TST") == 0){
		cineq(0) = 0;
		cineq_grad(0) = 0;
		cineq_grad(1) = 0;
		cineq_grad(2) = 0;
		cineq_grad(4) = 0;
	}
}

