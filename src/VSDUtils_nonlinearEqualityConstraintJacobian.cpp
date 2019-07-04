// VCSDUtils_nonlinearEqualityConstraintJacobian.cpp
// Last Modefied: 18-Jan-2016, Artur Wolek

#include<VSDUtils.h>

arma::mat VarSpeedDubins::nonlinearEqualityConstraintJacobian(arma::vec 
		paramsShort, std::string pathClass, std::string pathType, 
		std::string orientation, double r, double R){
	
	// a function that sets k1, and k2
	double k1, k2;
	determineCurvatureParams(pathClass, orientation, k1, k2);
	// a function that sets these params
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
	
	double pi = M_PI;

	arma::mat grad_out;
	if ( pathClass.compare("TST") == 0 ){
		grad_out.set_size(4,3);
		if (pathType.compare("BCB-S-BCB") == 0){
			////////////////////////////////////////////////
			// dx_da1
			grad_out(0,0) = ((R-r)*(k1*cos(g1*k2+k1*(pi+a1+g1))-k1*cos(a1*k1+(pi+g1)*(k1+k2)))+R*k1*cos(k1*(pi+a1+g1)+k2*(pi+g1+g2))-R*k1*cos(k1*(pi+a1+g1)))/k2-((R-r)*(k1*cos(k1*(pi+a1))-k1*cos(a1*k1))-R*k1*cos(k1*(pi+a1+g1)))/k1-L*k1*sin(k1*(pi+a1+g1));
			// dx_dg1
			grad_out(1,0) = R*cos(k1*(pi+a1+g1))-((cos(a1*k1+(pi+g1)*(k1+k2))*(k1+k2)-cos(g1*k2+k1*(pi+a1+g1))*(k1+k2))*(R-r)-R*cos(k1*(pi+a1+g1)+k2*(pi+g1+g2))*(k1+k2)+R*k1*cos(k1*(pi+a1+g1)))/k2-L*k1*sin(k1*(pi+a1+g1));
			// dx_dL
			grad_out(2,0) = cos(k1*(a1+g1+pi));
			// dx_dg2
			grad_out(3,0) = R*cos(k1*(a1+g1+pi)+k2*(g1+g2+pi));
			////////////////////////////////////////////////
			// dy_da1
			grad_out(0,1) = -((R-r)*(k1*sin(k1*(pi+a1))-k1*sin(a1*k1))-R*k1*sin(k1*(pi+a1+g1)))/k1+((R-r)*(k1*sin(g1*k2+k1*(pi+a1+g1))-k1*sin(a1*k1+(pi+g1)*(k1+k2)))+R*k1*sin(k1*(pi+a1+g1)+k2*(pi+g1+g2))-R*k1*sin(k1*(pi+a1+g1)))/k2+L*k1*cos(k1*(pi+a1+g1));
			// dy_dg1
			grad_out(1,1) = R*sin(k1*(pi+a1+g1))-((sin(a1*k1+(pi+g1)*(k1+k2))*(k1+k2)-sin(g1*k2+k1*(pi+a1+g1))*(k1+k2))*(R-r)-R*sin(k1*(pi+a1+g1)+k2*(pi+g1+g2))*(k1+k2)+R*k1*sin(k1*(pi+a1+g1)))/k2+L*k1*cos(k1*(pi+a1+g1));
			// dy_dL
			grad_out(2,1) = sin(k1*(a1+g1+pi));
			// dy_dg2
			grad_out(3,1) = R*sin(k1*(a1+g1+pi)+k2*(g1+g2+pi));
			////////////////////////////////////////////////
			// dpsi_da1
			grad_out(0,2) = 1.0/k1;
			// dpsi_dg1
			grad_out(1,2) = 1.0/k1+1.0/k2;
			// dpsi_dL
			grad_out(2,2) = 0.0;
			// dpsi_dg2
			grad_out(3,2) = 1.0/k2;
        }
		else if (pathType.compare("BCB-S-BC") == 0){
			////////////////////////////////////////////////
			// dx_da1
			grad_out(0,0) = (k1*r*cos(k1*(pi+a1+g1)+k2*(b2+g1))+k1*cos(g1*k2+k1*(pi+a1+g1))*(R-r)-R*k1*cos(k1*(pi+a1+g1)))/k2-((R-r)*(k1*cos(k1*(pi+a1))-k1*cos(a1*k1))-R*k1*cos(k1*(pi+a1+g1)))/k1-L*k1*sin(k1*(pi+a1+g1));
			// dx_dg1
			grad_out(1,0) = R*cos(k1*(pi+a1+g1))+(cos(g1*k2+k1*(pi+a1+g1))*(k1+k2)*(R-r)-R*k1*cos(k1*(pi+a1+g1))+r*cos(k1*(pi+a1+g1)+k2*(b2+g1))*(k1+k2))/k2-L*k1*sin(k1*(pi+a1+g1));
			// dx_dL
			grad_out(2,0) = cos(k1*(a1+g1+pi));
			// dx_db2
			 grad_out(3,0) = r*cos((b2+g1)*k2+k1*(a1+g1+pi));
			////////////////////////////////////////////////
			// dy_da1
			grad_out(0,1) = -((R-r)*(k1*sin(k1*(pi+a1))-k1*sin(a1*k1))-R*k1*sin(k1*(pi+a1+g1)))/k1+(k1*r*sin(k1*(pi+a1+g1)+k2*(b2+g1))+k1*sin(g1*k2+k1*(pi+a1+g1))*(R-r)-R*k1*sin(k1*(pi+a1+g1)))/k2+L*k1*cos(k1*(pi+a1+g1));
			// dy_dg1
			grad_out(1,1) = R*sin(k1*(pi+a1+g1))+(r*sin(k1*(pi+a1+g1)+k2*(b2+g1))*(k1+k2)+sin(g1*k2+k1*(pi+a1+g1))*(k1+k2)*(R-r)-R*k1*sin(k1*(pi+a1+g1)))/k2+L*k1*cos(k1*(pi+a1+g1));
			// dy_dL
			grad_out(2,1) = sin(k1*(a1+g1+pi));
			// dy_db2
			grad_out(3,1) = r*sin((b2+g1)*k2+k1*(a1+g1+pi));
			////////////////////////////////////////////////
			// dpsi_da1
			grad_out(0,2) = 1.0/k1;
			// dpsi_dg1
			grad_out(1,2) = 1.0/k1+1.0/k2;
			// dpsi_dL
			grad_out(2,2) = 0;
			// dpsi_db2
			grad_out(3,2) = 1.0/k2;
		}
		else if (pathType.compare("BCB-S-B") == 0){		
			////////////////////////////////////////////////
			// dx_da1
			grad_out(0,0) = -((R-r)*(k1*cos(k1*(pi+a1))-k1*cos(a1*k1))-R*k1*cos(k1*(pi+a1+g1)))/k1-L*k1*sin(k1*(pi+a1+g1))-(R*k1*sin(a2*k2*(1.0/2.0)+k1*(pi+a1+g1))*sin(a2*k2*(1.0/2.0))*2.0)/k2;
			// dx_dg1
			grad_out(1,0) = R*cos(k1*(pi+a1+g1))-L*k1*sin(k1*(pi+a1+g1))-(R*k1*sin(a2*k2*(1.0/2.0)+k1*(pi+a1+g1))*sin(a2*k2*(1.0/2.0))*2.0)/k2;
			// dx_dL
			grad_out(2,0) = cos(k1*(a1+g1+pi));
			// dx_da2
			grad_out(3,0) = R*cos(a2*k2*(1.0/2.0)+k1*(pi+a1+g1))*cos(a2*k2*(1.0/2.0))-R*sin(a2*k2*(1.0/2.0)+k1*(pi+a1+g1))*sin(a2*k2*(1.0/2.0));
			////////////////////////////////////////////////
			// dy_da1
			grad_out(0,1) = -((R-r)*(k1*sin(k1*(pi+a1))-k1*sin(a1*k1))-R*k1*sin(k1*(pi+a1+g1)))/k1+L*k1*cos(k1*(pi+a1+g1))+(R*k1*cos(a2*k2*(1.0/2.0)+k1*(pi+a1+g1))*sin(a2*k2*(1.0/2.0))*2.0)/k2;
			// dy_dg1
			grad_out(1,1)  = R*sin(k1*(pi+a1+g1))+L*k1*cos(k1*(pi+a1+g1))+(R*k1*cos(a2*k2*(1.0/2.0)+k1*(pi+a1+g1))*sin(a2*k2*(1.0/2.0))*2.0)/k2;
			//dy_dL
			grad_out(2,1) = sin(k1*(a1+g1+pi));
			// dy_da2
			grad_out(3,1) = R*cos((1.0/2.0)*a2*k2+k1*(a1+g1+pi))*sin((1.0/2.0)*a2*k2)+R*cos((1.0/2.0)*a2*k2)*sin((1.0/2.0)*a2*k2+k1*(a1+g1+pi));
			////////////////////////////////////////////////
			// dpsi_da1
			grad_out(0,2) = 1.0/k1;
			// dpsi_dg1
			grad_out(1,2) = 1.0/k1;
			// dpsi_dL
			grad_out(2,2) = 0.0;
			// dpsi_da2
			grad_out(3,2) = 1.0/k2;
		}
		else if (pathType.compare("CB-S-BCB") == 0){
			////////////////////////////////////////////////
			// dx_db1
			grad_out(0,0) = -L*k1*sin(k1*(b1+g1))-(-R*(k1*k1)*cos(k2*(pi+g1+g2)+k1*(b1+g1))+k1*(k1*cos(k2*(pi+g1)+k1*(b1+g1))-k1*cos(b1*k1+g1*(k1+k2)))*(R-r)+R*k1*cos(k1*(b1+g1))*(k1-k2)+k1*k2*cos(b1*k1)*(R-r))/(k1*k2);
			// dx_dg1
			grad_out(1,0) = -(k1*(cos(k2*(pi+g1)+k1*(b1+g1))*(k1+k2)-cos(b1*k1+g1*(k1+k2))*(k1+k2))*(R-r)+R*k1*cos(k1*(b1+g1))*(k1-k2)-R*k1*cos(k2*(pi+g1+g2)+k1*(b1+g1))*(k1+k2))/(k1*k2)-L*k1*sin(k1*(b1+g1));
			// dx_dL
			grad_out(2,0) = cos((b1+g1)*k1);
			// dx_dg2
			grad_out(3,0) = R*cos((b1+g1)*k1+k2*(g1+g2+pi));
			////////////////////////////////////////////////
			// dy_db1
			grad_out(0,1) = -((R-r)*(k1*sin(k2*(pi+g1)+k1*(b1+g1))-k1*sin(b1*k1+g1*(k1+k2)))-R*k1*sin(k2*(pi+g1+g2)+k1*(b1+g1))+R*k1*sin(k1*(b1+g1)))/k2+(R*k1*sin(k1*(b1+g1))-k1*sin(b1*k1)*(R-r))/k1+L*k1*cos(k1*(b1+g1));
			// dy_dg1
			grad_out(1,1) = -((sin(k2*(pi+g1)+k1*(b1+g1))*(k1+k2)-sin(b1*k1+g1*(k1+k2))*(k1+k2))*(R-r)+R*k1*sin(k1*(b1+g1))-R*sin(k2*(pi+g1+g2)+k1*(b1+g1))*(k1+k2))/k2+R*sin(k1*(b1+g1))+L*k1*cos(k1*(b1+g1));
			// dy_dL
			grad_out(2,1) = sin((b1+g1)*k1);
			// dy_dg2
			grad_out(3,1) = R*sin((b1+g1)*k1+k2*(g1+g2+pi));
			////////////////////////////////////////////////
			// dpsi_db1
			grad_out(0,2) = 1.0/k1;
			// dpsi_dg1
			grad_out(1,2) = 1.0/k1+1.0/k2;
			// dpsi_dL
			grad_out(2,2) = 0.0;
			// dpsi_dg2
			grad_out(3,2) = 1.0/k2;
		}	
		else if (pathType.compare("CB-S-BC") == 0){
			////////////////////////////////////////////////
			// dx_db1
			grad_out(0,0) = -(-(k1*k1)*r*cos(k1*(b1+g1)+k2*(b2+g1))-R*(k1*k1)*cos(g1*k2+k1*(b1+g1))+(k1*k1)*r*cos(g1*k2+k1*(b1+g1))+R*(k1*k1)*cos(k1*(b1+g1))+L*(k1*k1)*k2*sin(k1*(b1+g1))+k1*k2*cos(b1*k1)*(R-r)-R*k1*k2*cos(k1*(b1+g1)))/(k1*k2);
			// dx_dg1
			grad_out(1,0) = -(R*(k1*k1)*cos(k1*(b1+g1))-k1*r*cos(k1*(b1+g1)+k2*(b2+g1))*(k1+k2)-R*k1*cos(g1*k2+k1*(b1+g1))*(k1+k2)+L*(k1*k1)*k2*sin(k1*(b1+g1))+k1*r*cos(g1*k2+k1*(b1+g1))*(k1+k2)-R*k1*k2*cos(k1*(b1+g1)))/(k1*k2);
			// dx_dL
			grad_out(2,0) = cos((b1+g1)*k1);
			// dx_db2
			grad_out(3,0) = r*cos((b1+g1)*k1+(b2+g1)*k2);
			////////////////////////////////////////////////
			// dy_db1
			grad_out(0,1) = -(-k1*(k1*sin(g1*k2+k1*(b1+g1))*(R-r)+k1*r*sin(k1*(b1+g1)+k2*(b2+g1))+L*k1*k2*cos(k1*(b1+g1)))+R*k1*sin(k1*(b1+g1))*(k1-k2)+k1*k2*sin(b1*k1)*(R-r))/(k1*k2);
			// dy_dg1
			grad_out(1,1) = (k1*(sin(g1*k2+k1*(b1+g1))*(k1+k2)*(R-r)+r*sin(k1*(b1+g1)+k2*(b2+g1))*(k1+k2)+L*k1*k2*cos(k1*(b1+g1)))-R*k1*sin(k1*(b1+g1))*(k1-k2))/(k1*k2);
			// dy_dL
			grad_out(2,1) = sin((b1+g1)*k1);
			// dy_db2
			grad_out(3,1) = r*sin((b1+g1)*k1+(b2+g1)*k2);
			////////////////////////////////////////////////
			// dpsi_db1
			grad_out(0,2) = 1.0/k1;
			// dpsi_dg1
			grad_out(1,2) = 1.0/k1+1.0/k2;
			// dpsi_dL
			grad_out(2,2) = 0.0;
			// dpsi_db2
			grad_out(3,2) = 1.0/k2;
		}	
		else if (pathType.compare("CB-S-B") == 0){
			////////////////////////////////////////////////
			// dx_db1
			grad_out(0,0) = -(-R*(k1*k1)*cos(a2*k2+k1*(b1+g1))+R*k1*cos(k1*(b1+g1))*(k1-k2)+L*(k1*k1)*k2*sin(k1*(b1+g1))+k1*k2*cos(b1*k1)*(R-r))/(k1*k2);
			// dx_dg1
			grad_out(1,0) = -(-R*(k1*k1)*cos(a2*k2+k1*(b1+g1))+R*k1*cos(k1*(b1+g1))*(k1-k2)+L*(k1*k1)*k2*sin(k1*(b1+g1)))/(k1*k2);
			// dx_dL
			grad_out(2,0) = cos((b1+g1)*k1);
			// dx_da2
			grad_out(3,0) = R*cos((b1+g1)*k1+a2*k2);
			////////////////////////////////////////////////
			// dy_db1
			grad_out(0,1) = (R*(k1*k1)*sin(a2*k2+k1*(b1+g1))-R*k1*sin(k1*(b1+g1))*(k1-k2)+L*(k1*k1)*k2*cos(k1*(b1+g1))-k1*k2*sin(b1*k1)*(R-r))/(k1*k2);
			// dy_dg1
			grad_out(1,1) = (R*(k1*k1)*sin(a2*k2+k1*(b1+g1))-R*k1*sin(k1*(b1+g1))*(k1-k2)+L*(k1*k1)*k2*cos(k1*(b1+g1)))/(k1*k2);
			// dy_dL
			grad_out(2,1) = sin((b1+g1)*k1);
			// dy_da2
			grad_out(3,1) = R*sin((b1+g1)*k1+a2*k2);
			////////////////////////////////////////////////
			// dpsi_db1
			grad_out(0,2) = 1.0/k1;
			// dpsi_dg1
			grad_out(1,2) = 1.0/k1;
			// dpsi_dL
			grad_out(2,2) = 0.0;
			// dpsi_da2
			grad_out(3,2) = 1.0/k2;
		}		
		else if (pathType.compare("B-S-BCB") == 0){
			////////////////////////////////////////////////
			// dx_dg1
			grad_out(0,0) = R*cos(g1*k1)+(R*(k1*cos(g1*k1+k2*(pi+a2+g2))-k1*cos(g1*k1))-(R-r)*(k1*cos(g1*k1+k2*(pi+a2))-k1*cos(a2*k2+g1*k1)))/k2-L*k1*sin(g1*k1);
			// dx_dL
			grad_out(1,0) = cos(g1*k1);
			// dx_da2
			grad_out(2,0) = -((R-r)*(k2*cos(g1*k1+k2*(pi+a2))-k2*cos(a2*k2+g1*k1))-R*k2*cos(g1*k1+k2*(pi+a2+g2)))/k2;
			// dx_dg2
			grad_out(3,0) = R*cos(g1*k1+k2*(a2+g2+pi));
			////////////////////////////////////////////////
			// dy_dg1
			grad_out(0,1) = R*sin(g1*k1)-((k1*sin(g1*k1+k2*(pi+a2))-k1*sin(a2*k2+g1*k1))*(R-r)-R*k1*sin(g1*k1+k2*(pi+a2+g2))+R*k1*sin(g1*k1))/k2+L*k1*cos(g1*k1);
			// dy_dL
			grad_out(1,1) = sin(g1*k1);
			// dy_da2
			grad_out(2,1) = -((k2*sin(g1*k1+k2*(pi+a2))-k2*sin(a2*k2+g1*k1))*(R-r)-R*k2*sin(g1*k1+k2*(pi+a2+g2)))/k2;
			// dy_dg2
			grad_out(3,1) = R*sin(g1*k1+k2*(a2+g2+pi));
			////////////////////////////////////////////////
			// dpsi_dg1
			grad_out(0,2) = 1.0/k1;
			// dpsi_dL
			grad_out(1,2) = 0.0;
			// dpsi_da2
			grad_out(2,2) = 1.0/k2;
			// dpsi_dg2
			grad_out(3,2) = 1.0/k2;
		}	
		else if (pathType.compare("B-S-BC") == 0){
			//////////////////////////////////////////////
			// dx_dg1
			grad_out(0,0) = ((k1*k1)*r*cos(g1*k1+k2*(a2+b2))+(k1*k1)*cos(a2*k2+g1*k1)*(R-r)-R*k1*cos(g1*k1)*(k1-k2)-L*(k1*k1)*k2*sin(g1*k1))/(k1*k2);
			// dx_dL
			grad_out(1,0) = cos(g1*k1);
			// dx_da2
			grad_out(2,0) = (k1*k2*r*cos(g1*k1+k2*(a2+b2))+k1*k2*cos(a2*k2+g1*k1)*(R-r))/(k1*k2);
			// dx_db2
			grad_out(3,0) = r*cos(g1*k1+(a2+b2)*k2);
			////////////////////////////////////////////////
			// dy_dg1
			grad_out(0,1) = ((k1*k1)*r*sin(g1*k1+k2*(a2+b2))+(k1*k1)*sin(a2*k2+g1*k1)*(R-r)-R*k1*sin(g1*k1)*(k1-k2)+L*(k1*k1)*k2*cos(g1*k1))/(k1*k2);
			// dy_dL
			grad_out(1,1) = sin(g1*k1);
			// dy_da2
			grad_out(2,1) = (k1*k2*r*sin(g1*k1+k2*(a2+b2))+k1*k2*sin(a2*k2+g1*k1)*(R-r))/(k1*k2);
			// dy_db2
			grad_out(3,1) = r*sin(g1*k1+(a2+b2)*k2);
			////////////////////////////////////////////////
			// dpsi_dg1
			grad_out(0,2) = 1.0/k1;
			// dpsi_dL
			grad_out(1,2) = 0.0;
			// dpsi_da2
			grad_out(2,2) = 1.0/k2;
			//dpsi_dg2
			grad_out(3,2) = 1.0/k2;
		}				
	}
////////////////////////////////////////////////////////////
	else if ( (pathClass.compare("TT") == 0) ){	
		if (pathType.compare("BCB-BCB") == 0){
			grad_out.set_size(4,3);
			// dx_da1
			grad_out(0,0) = ((-r+R)*(k1*cos(a1*k1)-k1*cos((a1+b1)*k1))+k1*R*cos((a1+b1+g1)*k1))/k1+
                      (R*(-(k1*cos((a1+b1+g1)*k1))+k1*cos((a1+b1+g1)*k1+(b1+g1+g2)*k2))+
                      (-r+R)*(k1*cos((a1+b1+g1)*k1+g1*k2)-k1*cos(a1*k1+(b1+g1)*(k1+k2))))/k2;
			grad_out(1,0) = (-(k1*(-r+R)*cos((a1+b1)*k1))+k1*R*cos((a1+b1+g1)*k1))/k1+
                      (R*(-(k1*cos((a1+b1+g1)*k1))+(k1+k2)*cos((a1+b1+g1)*k1+(b1+g1+g2)*k2))+
                      (-r+R)*(k1*cos((a1+b1+g1)*k1+g1*k2)-(k1+k2)*cos(a1*k1+(b1+g1)*(k1+k2))))/k2;
			grad_out(2,0) = R*cos((a1+b1+g1)*k1)+(R*(-(k1*cos((a1+b1+g1)*k1))+(k1+k2)*cos((a1+b1+g1)*k1+(b1+g1+g2)*k2))+
                      (-r+R)*((k1+k2)*cos((a1+b1+g1)*k1+g1*k2)-(k1+k2)*cos(a1*k1+(b1+g1)*(k1+k2))))/k2;
			grad_out(3,0) = R*cos((a1+b1+g1)*k1+(b1+g1+g2)*k2);
			// dy_da1
			grad_out(0,1) = ((r-R)*(-(k1*sin(a1*k1))+k1*sin((a1+b1)*k1))+k1*R*sin((a1+b1+g1)*k1))/k1+
                      (R*(-(k1*sin((a1+b1+g1)*k1))+k1*sin((a1+b1+g1)*k1+(b1+g1+g2)*k2))+
                      (r-R)*(-(k1*sin((a1+b1+g1)*k1+g1*k2))+k1*sin(a1*k1+(b1+g1)*(k1+k2))))/k2;
			grad_out(1,1) = (k1*(r-R)*sin((a1+b1)*k1)+k1*R*sin((a1+b1+g1)*k1))/k1+
                      (R*(-(k1*sin((a1+b1+g1)*k1))+(k1+k2)*sin((a1+b1+g1)*k1+(b1+g1+g2)*k2))+
                      (r-R)*(-(k1*sin((a1+b1+g1)*k1+g1*k2))+(k1+k2)*sin(a1*k1+(b1+g1)*(k1+k2))))/k2;
			grad_out(2,1) = R*sin((a1+b1+g1)*k1)+(R*(-(k1*sin((a1+b1+g1)*k1))+(k1+k2)*sin((a1+b1+g1)*k1+(b1+g1+g2)*k2))+
                      (r-R)*(-((k1+k2)*sin((a1+b1+g1)*k1+g1*k2))+(k1+k2)*sin(a1*k1+(b1+g1)*(k1+k2))))/k2;
			grad_out(3,1) = R*sin((a1+b1+g1)*k1+(b1+g1+g2)*k2);
			// dpsi_da1
			grad_out(0,2) = 1/k1;
			grad_out(1,2) = 1/k1+1/k2;
			grad_out(2,2) = 1/k1+1/k2;
			grad_out(3,2) = 1/k2;
		}
		else if (pathType.compare("BCB-BC") == 0){
			grad_out.set_size(4,3);
			//dx_da1
			grad_out(0,0) = ((-r+R)*(k1*cos(a1*k1)-k1*cos((a1+b1)*k1))+k1*R*cos((a1+b1+g1)*k1))/k1+
                      (-(k1*R*cos((a1+b1+g1)*k1))+k1*(-r+R)*cos((a1+b1+g1)*k1+g1*k2)+
                      k1*r*cos((a1+b1+g1)*k1+(b2+g1)*k2))/k2; 
			grad_out(1,0) = (-(k1*(-r+R)*cos((a1+b1)*k1))+k1*R*cos((a1+b1+g1)*k1))/k1+
                      (-(k1*R*cos((a1+b1+g1)*k1))+k1*(-r+R)*cos((a1+b1+g1)*k1+g1*k2)+
                      k1*r*cos((a1+b1+g1)*k1+(b2+g1)*k2))/k2;
			grad_out(2,0) = R*cos((a1+b1+g1)*k1)+(-(k1*R*cos((a1+b1+g1)*k1))+(k1+k2)*(-r+R)*cos((a1+b1+g1)*k1+g1*k2)+
                      (k1+k2)*r*cos((a1+b1+g1)*k1+(b2+g1)*k2))/k2;
			grad_out(3,0) = r*cos((a1+b1+g1)*k1+(b2+g1)*k2);
			// dy_da1
			grad_out(0,1) = ((r-R)*(-(k1*sin(a1*k1))+k1*sin((a1+b1)*k1))+k1*R*sin((a1+b1+g1)*k1))/k1+
                      (-(k1*R*sin((a1+b1+g1)*k1))-k1*(r-R)*sin((a1+b1+g1)*k1+g1*k2)+
                      k1*r*sin((a1+b1+g1)*k1+(b2+g1)*k2))/k2;
			grad_out(1,1) = (k1*(r-R)*sin((a1+b1)*k1)+k1*R*sin((a1+b1+g1)*k1))/k1+
                      (-(k1*R*sin((a1+b1+g1)*k1))-k1*(r-R)*sin((a1+b1+g1)*k1+g1*k2)+
                      k1*r*sin((a1+b1+g1)*k1+(b2+g1)*k2))/k2;
			grad_out(2,1) = R*sin((a1 + b1 + g1)*k1) + (-(k1*R*sin((a1 + b1 + g1)*k1)) - 
                      (k1 + k2)*(r - R)*sin((a1 + b1 + g1)*k1 + g1*k2) + 
                      (k1 + k2)*r*sin((a1 + b1 + g1)*k1 + (b2 + g1)*k2))/k2;
			grad_out(3,1) = r*sin((a1+b1+g1)*k1+(b2+g1)*k2);
			// dpsi_da1
			grad_out(0,2) = 1/k1;
			grad_out(1,2) = 1/k1;
			grad_out(2,2) = 1/k1+1/k2;
			grad_out(3,2) = 1/k2;
		}	
		else if (pathType.compare("BCB-B") == 0){
			grad_out.set_size(4,3);
			// dx_da1
			grad_out(0,0) = ((-r+R)*(k1*cos(a1*k1)-k1*cos((a1+b1)*k1))+k1*R*cos((a1+b1+g1)*k1))/k1-
                      (2*k1*R*sin((a2*k2)/2.)*sin((a1+b1+g1)*k1+(a2*k2)/2.))/k2;
			grad_out(1,0) = (-(k1*(-r+R)*cos((a1+b1)*k1))+k1*R*cos((a1+b1+g1)*k1))/k1-
                      (2*k1*R*sin((a2*k2)/2.)*sin((a1+b1+g1)*k1+(a2*k2)/2.))/k2;
			grad_out(2,0) = R*cos((a1+b1+g1)*k1)-(2*k1*R*sin((a2*k2)/2.)*sin((a1+b1+g1)*k1+(a2*k2)/2.))/k2;
			grad_out(3,0) = R*cos((a2*k2)/2.)*cos((a1+b1+g1)*k1+(a2*k2)/2.)-R*sin((a2*k2)/2.)*sin((a1+b1+g1)*k1+(a2*k2)/2.);
			// dy_da1
			grad_out(0,1) = ((r-R)*(-(k1*sin(a1*k1))+k1*sin((a1+b1)*k1))+k1*R*sin((a1+b1+g1)*k1))/k1+
                      (2*k1*R*cos((a1+b1+g1)*k1+(a2*k2)/2.)*sin((a2*k2)/2.))/k2;
			grad_out(1,1) = (k1*(r-R)*sin((a1+b1)*k1)+k1*R*sin((a1+b1+g1)*k1))/k1+
                      (2*k1*R*cos((a1+b1+g1)*k1+(a2*k2)/2.)*sin((a2*k2)/2.))/k2;
			grad_out(2,1) = R*sin((a1+b1+g1)*k1)+(2*k1*R*cos((a1+b1+g1)*k1+(a2*k2)/2.)*sin((a2*k2)/2.))/k2;
			grad_out(3,1) = R*cos((a1+b1+g1)*k1+(a2*k2)/2.)*sin((a2*k2)/2.)+R*cos((a2*k2)/2.)*sin((a1+b1+g1)*k1+(a2*k2)/2.);
			// dpsi_da1
			grad_out(0,2) = 1/k1;
			grad_out(1,2) = 1/k1;
			grad_out(2,2) = 1/k1;
			grad_out(3,2) = 1/k2;
		}		
		else if (pathType.compare("CB-BCB") == 0){
			grad_out.set_size(4,3);
			// dx_da1
			grad_out(0,0) = (k1*(r-R)*cos(b1*k1)+k1*R*cos((b1+g1)*k1))/k1+
                      (R*(-(k1*cos((b1+g1)*k1))+k1*cos((b1+g1)*k1+(b2+g1+g2)*k2))+
                      (r-R)*(k1*cos((b1+g1)*k1+(b2+g1)*k2)-k1*cos(b1*k1+g1*(k1+k2))))/k2;
			grad_out(1,0) = R*cos((b1+g1)*k1)+(R*(-(k1*cos((b1+g1)*k1))+(k1+k2)*cos((b1+g1)*k1+(b2+g1+g2)*k2))+
                      (r-R)*((k1+k2)*cos((b1+g1)*k1+(b2+g1)*k2)-(k1+k2)*cos(b1*k1+g1*(k1+k2))))/k2;
			grad_out(2,0) = (k2*(r-R)*cos((b1+g1)*k1+(b2+g1)*k2)+k2*R*cos((b1+g1)*k1+(b2+g1+g2)*k2))/k2;
			grad_out(3,0) = R*cos((b1+g1)*k1+(b2+g1+g2)*k2);
			// dy_da1
			grad_out(0,1) = (-(k1*(-r+R)*sin(b1*k1))+k1*R*sin((b1+g1)*k1))/k1+
                      (R*(-(k1*sin((b1+g1)*k1))+k1*sin((b1+g1)*k1+(b2+g1+g2)*k2))+
                      (r-R)*(k1*sin((b1+g1)*k1+(b2+g1)*k2)-k1*sin(b1*k1+g1*(k1+k2))))/k2;
			grad_out(1,1) = R*sin((b1+g1)*k1)+(R*(-(k1*sin((b1+g1)*k1))+(k1+k2)*sin((b1+g1)*k1+(b2+g1+g2)*k2))+
                      (r-R)*((k1+k2)*sin((b1+g1)*k1+(b2+g1)*k2)-(k1+k2)*sin(b1*k1+g1*(k1+k2))))/k2;
			grad_out(2,1) = (k2*(r-R)*sin((b1+g1)*k1+(b2+g1)*k2)+k2*R*sin((b1+g1)*k1+(b2+g1+g2)*k2))/k2;
			grad_out(3,1) = R*sin((b1+g1)*k1+(b2+g1+g2)*k2);
			// dpsi_da1 
			grad_out(0,2) = 1/k1;
			grad_out(1,2) = 1/k1+1/k2;
			grad_out(2,2) = 1/k2;
			grad_out(3,2) = 1/k2;
		}	
		else if (pathType.compare("B-BCB") == 0){
			grad_out.set_size(4,3);
			// dx_da1
			grad_out(0,0) = R*cos(g1*k1)+((-r+R)*(k1*cos(g1*k1+a2*k2)-k1*cos(g1*k1+(a2+b2)*k2))+
                      R*(-(k1*cos(g1*k1))+k1*cos(g1*k1+(a2+b2+g2)*k2)))/k2;
			grad_out(1,0) = ((-r+R)*(k2*cos(g1*k1+a2*k2)-k2*cos(g1*k1+(a2+b2)*k2))+k2*R*cos(g1*k1+(a2+b2+g2)*k2))/k2;
			grad_out(2,0) = (-(k2*(-r+R)*cos(g1*k1+(a2+b2)*k2))+k2*R*cos(g1*k1+(a2+b2+g2)*k2))/k2;
			grad_out(3,0) = R*cos(g1*k1+(a2+b2+g2)*k2);
			// dy_da1
			grad_out(0,1) = (-(k1*(k1-k2)*R*sin(g1*k1))+k1*(r-R)*(-(k1*sin(g1*k1+a2*k2))+k1*sin(g1*k1+(a2+b2)*k2))+
                      k1*k1*R*sin(g1*k1+(a2+b2+g2)*k2))/(k1*k2);
			grad_out(1,1) = (k1*(r-R)*(-(k2*sin(g1*k1+a2*k2))+k2*sin(g1*k1+(a2+b2)*k2))+k1*k2*R*sin(g1*k1+(a2+b2+g2)*k2))/(k1*k2);
			grad_out(2,1) = (k1*k2*(r-R)*sin(g1*k1+(a2+b2)*k2)+k1*k2*R*sin(g1*k1+(a2+b2+g2)*k2))/(k1*k2);
			grad_out(3,1) = R*sin(g1*k1+(a2+b2+g2)*k2);
			// dpsi_da1
			grad_out(0,2) = 1/k1;
			grad_out(1,2) = 1/k2;
			grad_out(2,2) = 1/k2;
			grad_out(3,2) = 1/k2;
		}		
	}
////////////////////////////////////////////////////////////
	else if ( pathClass.compare("TTT") == 0 ){	
		grad_out.set_size(4,3);
		if (pathType.compare("BCB-T-BCB") == 0){
			grad_out(0,0) = (R*(k1*cos(k1*(a1-g1))*-2.0+k1*cos(k1*(a1+b1+g1))*2.0+k1*cos(k1*(a1+b1+g3)))-k1*cos(k1*(a1+b1))*(R-r)*3.0+k1*cos(a1*k1)*(R-r)*3.0)/k1;
			grad_out(1,0) = (R*(k1*cos(k1*(a1+b1+g1))*2.0+k1*cos(k1*(a1+b1+g3)))-k1*cos(k1*(a1+b1))*(R-r)*3.0)/k1;
			grad_out(2,0) = (R*(k1*cos(k1*(a1-g1))*2.0+k1*cos(k1*(a1+b1+g1))*2.0))/k1;
			grad_out(3,0) = R*cos((a1+b1+g3)*k1);
			grad_out(0,1) = (R*(k1*sin(k1*(a1-g1))*-2.0+k1*sin(k1*(a1+b1+g1))*2.0+k1*sin(k1*(a1+b1+g3)))-k1*sin(k1*(a1+b1))*(R-r)*3.0+k1*sin(a1*k1)*(R-r)*3.0)/k1;
			grad_out(1,1) = (R*(k1*sin(k1*(a1+b1+g1))*2.0+k1*sin(k1*(a1+b1+g3)))-k1*sin(k1*(a1+b1))*(R-r)*3.0)/k1;
			grad_out(2,1) = (R*(k1*sin(k1*(a1-g1))*2.0+k1*sin(k1*(a1+b1+g1))*2.0))/k1;
			grad_out(3,1) = R*sin((a1+b1+g3)*k1);
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = 1.0/k1;
			grad_out(2,2) = 0.0;
			grad_out(3,2) = 1.0/k1;
        }
        else if (pathType.compare("BCB-T-BC") == 0){
			grad_out(0,0) = (k1*r*cos(k1*(a1+b3))-k1*cos(k1*(a1+b1))*(R-r)*2.0-R*k1*cos(k1*(a1-g1))*2.0+R*k1*cos(k1*(a1+b1+g1))*2.0+k1*cos(a1*k1)*(R-r)*3.0)/k1;
			grad_out(1,0) = -(k1*cos(k1*(a1+b1))*(R-r)*2.0-R*k1*cos(k1*(a1+b1+g1))*2.0)/k1;
			grad_out(2,0) = (R*k1*cos(k1*(a1-g1))*2.0+R*k1*cos(k1*(a1+b1+g1))*2.0)/k1;
			grad_out(3,0) = r*cos((a1+b3)*k1);
			grad_out(0,1) = (k1*r*sin(k1*(a1+b3))-k1*sin(k1*(a1+b1))*(R-r)*2.0-R*k1*sin(k1*(a1-g1))*2.0+R*k1*sin(k1*(a1+b1+g1))*2.0+k1*sin(a1*k1)*(R-r)*3.0)/k1;
			grad_out(1,1) = -(k1*sin(k1*(a1+b1))*(R-r)*2.0-R*k1*sin(k1*(a1+b1+g1))*2.0)/k1;
			grad_out(2,1) = (R*k1*sin(k1*(a1-g1))*2.0+R*k1*sin(k1*(a1+b1+g1))*2.0)/k1;
			grad_out(3,1) = r*sin((a1+b3)*k1);
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = 0.0;
			grad_out(2,2) = 0.0;
			grad_out(3,2) = 1.0/k1;     
		}   
        else if (pathType.compare("BCB-T-B") == 0){
			grad_out(0,0) = (R*(k1*cos(k1*(a1-g1))*-2.0+k1*cos(k1*(a1+b1+g1))*2.0+k1*cos(k1*(a1+a3-g1)))-k1*cos(k1*(a1+b1))*(R-r)*2.0+k1*cos(a1*k1)*(R-r)*2.0)/k1;
			grad_out(1,0) = -(k1*cos(k1*(a1+b1))*(R-r)*2.0-R*k1*cos(k1*(a1+b1+g1))*2.0)/k1;
			grad_out(2,0) = (R*(k1*cos(k1*(a1-g1))*2.0+k1*cos(k1*(a1+b1+g1))*2.0-k1*cos(k1*(a1+a3-g1))))/k1;
			grad_out(3,0) = R*cos((a1+a3+(-1)*g1)*k1);
			grad_out(0,1) = (R*(k1*sin(k1*(a1-g1))*-2.0+k1*sin(k1*(a1+b1+g1))*2.0+k1*sin(k1*(a1+a3-g1)))-k1*sin(k1*(a1+b1))*(R-r)*2.0+k1*sin(a1*k1)*(R-r)*2.0)/k1;
			grad_out(1,1) = -(k1*sin(k1*(a1+b1))*(R-r)*2.0-R*k1*sin(k1*(a1+b1+g1))*2.0)/k1;
			grad_out(2,1) = (R*(k1*sin(k1*(a1-g1))*2.0+k1*sin(k1*(a1+b1+g1))*2.0-k1*sin(k1*(a1+a3-g1))))/k1;
			grad_out(3,1) = R*sin((a1+a3+(-1)*g1)*k1);
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = 0.0;
			grad_out(2,2) = -1.0/k1;
			grad_out(3,2) = 1.0/k1;
        }	
        else if (pathType.compare("CB-T-BCB") == 0){
			grad_out(0,0)  = (R*(k1*cos(k1*(b1+g1))*2.0+k1*cos(k1*(b1+g3))-k1*cos(k1*(-b1+b2+g1))*2.0)-k1*cos(b1*k1)*(R-r)*3.0+k1*cos(k1*(b1-b2))*(R-r)*2.0)/k1;
			grad_out(1,0) = (R*(k1*cos(k1*(b1+g1))*2.0+k1*cos(k1*(-b1+b2+g1))*2.0))/k1;
			grad_out(2,0) = (R*k1*cos(k1*(-b1+b2+g1))*2.0-k1*cos(k1*(b1-b2))*(R-r)*2.0)/k1;
			grad_out(3,0) = R*cos((b1+g3)*k1);
			grad_out(0,1) = (k1*sin(k1*(b1-b2))*(R-r)*2.0+R*k1*sin(k1*(-b1+b2+g1))*2.0+R*k1*sin(k1*(b1+g1))*2.0+R*k1*sin(k1*(b1+g3))-k1*sin(b1*k1)*(R-r)*3.0)/k1;
			grad_out(1,1) = -(R*k1*sin(k1*(-b1+b2+g1))*2.0-R*k1*sin(k1*(b1+g1))*2.0)/k1;
			grad_out(2,1) = -(k1*sin(k1*(b1-b2))*(R-r)*2.0+R*k1*sin(k1*(-b1+b2+g1))*2.0)/k1;
			grad_out(3,1) = R*sin((b1+g3)*k1);
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = 0.0;
			grad_out(2,2) = 0.0;
			grad_out(3,2) = 1.0/k1;
        }       
        else if (pathType.compare("CB-T-BC") == 0){
			grad_out(0,0) = (R*k1*cos(k1*(-b1+b2+g1))*-2.0+k1*r*cos(k1*(b1-b2+b3))+R*k1*cos(k1*(b1+g1))*2.0-k1*cos(b1*k1)*(R-r)*2.0+k1*cos(k1*(b1-b2))*(R-r)*2.0)/k1;
			grad_out(1,0) = (R*k1*cos(k1*(-b1+b2+g1))*2.0+R*k1*cos(k1*(b1+g1))*2.0)/k1;
			grad_out(2,0) = -(R*k1*cos(k1*(-b1+b2+g1))*-2.0+k1*r*cos(k1*(b1-b2+b3))+k1*cos(k1*(b1-b2))*(R-r)*2.0)/k1;
			grad_out(3,0) = r*cos((b1+(-1)*b2+b3)*k1);
			grad_out(0,1) = (k1*sin(k1*(b1-b2))*(R-r)*2.0+R*k1*sin(k1*(-b1+b2+g1))*2.0+k1*r*sin(k1*(b1-b2+b3))+R*k1*sin(k1*(b1+g1))*2.0-k1*sin(b1*k1)*(R-r)*2.0)/k1;
			grad_out(1,1) = -(R*k1*sin(k1*(-b1+b2+g1))*2.0-R*k1*sin(k1*(b1+g1))*2.0)/k1;
			grad_out(2,1) = -(k1*sin(k1*(b1-b2))*(R-r)*2.0+R*k1*sin(k1*(-b1+b2+g1))*2.0+k1*r*sin(k1*(b1-b2+b3)))/k1;
			grad_out(3,1) = r*sin((b1+(-1)*b2+b3)*k1);
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = 0.0;
			grad_out(2,2) = -1.0/k1;
			grad_out(3,2) = 1.0/k1;
        }     
        else if (pathType.compare("CB-T-B") == 0){
			grad_out(0,0) = (R*(k1*cos(k1*(b1+g1))*2.0+k1*cos(k1*(a3+b1-b2-g1))-k1*cos(k1*(-b1+b2+g1))*2.0)-k1*cos(b1*k1)*(R-r)*2.0+k1*cos(k1*(b1-b2))*(R-r))/k1;
			grad_out(1,0) = (R*(k1*cos(k1*(b1+g1))*2.0-k1*cos(k1*(a3+b1-b2-g1))+k1*cos(k1*(-b1+b2+g1))*2.0))/k1;
			grad_out(2,0) = -(R*(k1*cos(k1*(a3+b1-b2-g1))-k1*cos(k1*(-b1+b2+g1))*2.0)+k1*cos(k1*(b1-b2))*(R-r))/k1;
			grad_out(3,0) = R*cos((a3+b1+(-1)*b2+(-1)*g1)*k1);
			grad_out(0,1) = (R*k1*sin(k1*(a3+b1-b2-g1))+k1*sin(k1*(b1-b2))*(R-r)+R*k1*sin(k1*(-b1+b2+g1))*2.0+R*k1*sin(k1*(b1+g1))*2.0-k1*sin(b1*k1)*(R-r)*2.0)/k1;
			grad_out(1,1) = -(R*k1*sin(k1*(a3+b1-b2-g1))+R*k1*sin(k1*(-b1+b2+g1))*2.0-R*k1*sin(k1*(b1+g1))*2.0)/k1;
			grad_out(2,1) = -(R*k1*sin(k1*(a3+b1-b2-g1))+k1*sin(k1*(b1-b2))*(R-r)+R*k1*sin(k1*(-b1+b2+g1))*2.0)/k1;
			grad_out(3,1) = R*sin((a3+b1+(-1)*b2+(-1)*g1)*k1);
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = -1.0/k1;
			grad_out(2,2) = -1.0/k1;
			grad_out(3,2) = 1.0/k1;
        }
        else if (pathType.compare("B-T-BCB") == 0){
			grad_out(0,0) = -(R*(k1*cos(k1*(a2*2.0+b2-g1))-k1*cos(g1*k1))*2.0+k1*cos(k1*(a2-g1))*(R-r)*2.0-R*k1*cos(k1*(-a2+g1+g3))-k1*cos(k1*(a2+b2-g1))*(R-r)*2.0)/k1;
			grad_out(1,0) = (k1*cos(k1*(a2-g1))*(R-r)*2.0-R*k1*cos(k1*(-a2+g1+g3))-k1*cos(k1*(a2+b2-g1))*(R-r)*2.0+R*k1*cos(k1*(a2*2.0+b2-g1))*4.0)/k1;
			grad_out(2,0) = -(k1*cos(k1*(a2+b2-g1))*(R-r)*2.0-R*k1*cos(k1*(a2*2.0+b2-g1))*2.0)/k1;
			grad_out(3,0) = R*cos((a2+(-1)*g1+(-1)*g3)*k1);
			grad_out(0,1)  = (R*(k1*sin(k1*(a2*2.0+b2-g1))*2.0+k1*sin(g1*k1)*2.0+k1*sin(k1*(-a2+g1+g3)))+k1*sin(k1*(a2-g1))*(R-r)*2.0-k1*sin(k1*(a2+b2-g1))*(R-r)*2.0)/k1;
			grad_out(1,1) = -(R*(k1*sin(k1*(a2*2.0+b2-g1))*4.0+k1*sin(k1*(-a2+g1+g3)))+k1*sin(k1*(a2-g1))*(R-r)*2.0-k1*sin(k1*(a2+b2-g1))*(R-r)*2.0)/k1;
			grad_out(2,1) = (k1*sin(k1*(a2+b2-g1))*(R-r)*2.0-R*k1*sin(k1*(a2*2.0+b2-g1))*2.0)/k1;
			grad_out(3,1) = (-1)*R*sin((a2+(-1)*g1+(-1)*g3)*k1);
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = -1.0/k1;
			grad_out(2,2) = 0.0;
			grad_out(3,2) = 1.0/k1;			
        }  
        else if (pathType.compare("B-T-BC") == 0){
			grad_out(0,0) = -(R*(k1*cos(k1*(a2*2.0+b2-g1))-k1*cos(g1*k1))*2.0+k1*cos(k1*(a2-g1))*(R-r)-k1*r*cos(k1*(a2+b2-b3-g1))-k1*cos(k1*(a2+b2-g1))*(R-r)*2.0)/k1;
			grad_out(1,0) = (k1*cos(k1*(a2-g1))*(R-r)-k1*r*cos(k1*(a2+b2-b3-g1))-k1*cos(k1*(a2+b2-g1))*(R-r)*2.0+R*k1*cos(k1*(a2*2.0+b2-g1))*4.0)/k1;
			grad_out(2,0) = -(k1*r*cos(k1*(a2+b2-b3-g1))+k1*cos(k1*(a2+b2-g1))*(R-r)*2.0-R*k1*cos(k1*(a2*2.0+b2-g1))*2.0)/k1;
			grad_out(3,0) = r*cos((a2+b2+(-1)*b3+(-1)*g1)*k1);
			grad_out(0,1) = (k1*sin(k1*(a2-g1))*(R-r)-k1*r*sin(k1*(a2+b2-b3-g1))-k1*sin(k1*(a2+b2-g1))*(R-r)*2.0+R*k1*sin(k1*(a2*2.0+b2-g1))*2.0+R*k1*sin(g1*k1)*2.0)/k1;
			grad_out(1,1) = -(k1*sin(k1*(a2-g1))*(R-r)-k1*r*sin(k1*(a2+b2-b3-g1))-k1*sin(k1*(a2+b2-g1))*(R-r)*2.0+R*k1*sin(k1*(a2*2.0+b2-g1))*4.0)/k1;
			grad_out(2,1) = (k1*r*sin(k1*(a2+b2-b3-g1))+k1*sin(k1*(a2+b2-g1))*(R-r)*2.0-R*k1*sin(k1*(a2*2.0+b2-g1))*2.0)/k1;
			grad_out(3,1) = (-1)*r*sin((a2+b2+(-1)*b3+(-1)*g1)*k1);
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = -1.0/k1;
			grad_out(2,2) = -1.0/k1;
			grad_out(3,2) = 1.0/k1;
        }                  	
        else if (pathType.compare("B-T-B") == 0){
			grad_out(0,0) = (R*(k1*cos(k1*(a2*2.0-a3+b2-g1))-k1*cos(k1*(a2*2.0+b2-g1))*2.0+k1*cos(g1*k1)*2.0)-k1*cos(k1*(a2-g1))*(R-r)+k1*cos(k1*(a2+b2-g1))*(R-r))/k1;
			grad_out(1,0) = -(R*(k1*cos(k1*(a2*2.0-a3+b2-g1))*2.0-k1*cos(k1*(a2*2.0+b2-g1))*4.0)-k1*cos(k1*(a2-g1))*(R-r)+k1*cos(k1*(a2+b2-g1))*(R-r))/k1;
			grad_out(2,0) = -(R*(k1*cos(k1*(a2*2.0-a3+b2-g1))-k1*cos(k1*(a2*2.0+b2-g1))*2.0)+k1*cos(k1*(a2+b2-g1))*(R-r))/k1;
			grad_out(3,0) = R*cos(((-2)*a2+a3+(-1)*b2+g1)*k1);
			grad_out(0,1) = (R*(-k1*sin(k1*(a2*2.0-a3+b2-g1))+k1*sin(k1*(a2*2.0+b2-g1))*2.0+k1*sin(g1*k1)*2.0)+k1*sin(k1*(a2-g1))*(R-r)-k1*sin(k1*(a2+b2-g1))*(R-r))/k1;
			grad_out(1,1) = (R*(k1*sin(k1*(a2*2.0-a3+b2-g1))*2.0-k1*sin(k1*(a2*2.0+b2-g1))*4.0)-k1*sin(k1*(a2-g1))*(R-r)+k1*sin(k1*(a2+b2-g1))*(R-r))/k1;
			grad_out(2,1) = (R*(k1*sin(k1*(a2*2.0-a3+b2-g1))-k1*sin(k1*(a2*2.0+b2-g1))*2.0)+k1*sin(k1*(a2+b2-g1))*(R-r))/k1;
			grad_out(3,1) = R*sin(((-2)*a2+a3+(-1)*b2+g1)*k1);
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = (-2)*1.0/k1;
			grad_out(2,2) = -1.0/k1;
			grad_out(3,2) = 1.0/k1;
        }        
	}
////////////////////////////////////////////////////////////
	else if ( pathClass.compare("TTTT") == 0 ){	
		grad_out.set_size(4,3);
		if (pathType.compare("BCB-TT-BCB") == 0){
			////
			grad_out(0,0) = -(R*(k1*cos(k1*(a1-g1))*2.0+k1*cos(k1*(a1-g4))-k1*cos(k1*(a1+b1+g1))*4.0)+k1*cos(k1*(a1+b1))*(R-r)*4.0-k1*cos(a1*k1)*(R-r)*4.0)/k1;
			grad_out(1,0) = -(k1*cos(k1*(a1+b1))*(R-r)*4.0-R*k1*cos(k1*(a1+b1+g1))*4.0)/k1;
			grad_out(2,0) = (R*(k1*cos(k1*(a1-g1))*2.0+k1*cos(k1*(a1+b1+g1))*4.0))/k1;
			grad_out(3,0) = R*cos((a1+(-1)*g4)*k1);
			////
			grad_out(0,1) = -(R*(k1*sin(k1*(a1-g1))*2.0+k1*sin(k1*(a1-g4))-k1*sin(k1*(a1+b1+g1))*4.0)+k1*sin(k1*(a1+b1))*(R-r)*4.0-k1*sin(a1*k1)*(R-r)*4.0)/k1;
			grad_out(1,1) = -(k1*sin(k1*(a1+b1))*(R-r)*4.0-R*k1*sin(k1*(a1+b1+g1))*4.0)/k1;
			grad_out(2,1) = (R*(k1*sin(k1*(a1-g1))*2.0+k1*sin(k1*(a1+b1+g1))*4.0))/k1;
			grad_out(3,1) = R*sin((a1+(-1)*g4)*k1);
			////
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = 0.0;
			grad_out(2,2) = 0.0;
			grad_out(3,2) = -1.0/k1;
			////
		}
		else if (pathType.compare("BCB-TT-BC") == 0){
			////
			grad_out(0,0) = -(k1*r*cos(k1*(a1+b1-b4))+k1*cos(k1*(a1+b1))*(R-r)*4.0+R*k1*cos(k1*(a1-g1))*2.0-R*k1*cos(k1*(a1+b1+g1))*4.0-k1*cos(a1*k1)*(R-r)*3.0)/k1;
			grad_out(1,0) = -(k1*r*cos(k1*(a1+b1-b4))+k1*cos(k1*(a1+b1))*(R-r)*4.0-R*k1*cos(k1*(a1+b1+g1))*4.0)/k1;
			grad_out(2,0) = (R*k1*cos(k1*(a1-g1))*2.0+R*k1*cos(k1*(a1+b1+g1))*4.0)/k1;
			grad_out(3,0) = r*cos((a1+b1+(-1)*b4)*k1);
			////
			grad_out(0,1) = -(k1*r*sin(k1*(a1+b1-b4))+k1*sin(k1*(a1+b1))*(R-r)*4.0+R*k1*sin(k1*(a1-g1))*2.0-R*k1*sin(k1*(a1+b1+g1))*4.0-k1*sin(a1*k1)*(R-r)*3.0)/k1;
			grad_out(1,1) = -(k1*r*sin(k1*(a1+b1-b4))+k1*sin(k1*(a1+b1))*(R-r)*4.0-R*k1*sin(k1*(a1+b1+g1))*4.0)/k1;
			grad_out(2,1) = (R*k1*sin(k1*(a1-g1))*2.0+R*k1*sin(k1*(a1+b1+g1))*4.0)/k1;
			grad_out(3,1) = r*sin((a1+b1+(-1)*b4)*k1);
			////
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = 1.0/k1;
			grad_out(2,2) = 0.0;
			grad_out(3,2) = -1.0/k1;
			////
		}	
		else if (pathType.compare("BCB-TT-B") == 0){
			////
			grad_out(0,0) = -(R*(k1*cos(a4*k1-k1*(a1+b1+g1))+k1*cos(k1*(a1-g1))*2.0-k1*cos(k1*(a1+b1+g1))*4.0)+k1*cos(k1*(a1+b1))*(R-r)*3.0-k1*cos(a1*k1)*(R-r)*3.0)/k1;
			grad_out(1,0) = -(R*(k1*cos(a4*k1-k1*(a1+b1+g1))-k1*cos(k1*(a1+b1+g1))*4.0)+k1*cos(k1*(a1+b1))*(R-r)*3.0)/k1;
			grad_out(2,0) = (R*(-k1*cos(a4*k1-k1*(a1+b1+g1))+k1*cos(k1*(a1-g1))*2.0+k1*cos(k1*(a1+b1+g1))*4.0))/k1;
			grad_out(3,0) = R*cos(a4*k1+(-1)*(a1+b1+g1)*k1);
			////
			grad_out(0,1) = (R*(k1*sin(a4*k1-k1*(a1+b1+g1))-k1*sin(k1*(a1-g1))*2.0+k1*sin(k1*(a1+b1+g1))*4.0)-k1*sin(k1*(a1+b1))*(R-r)*3.0+k1*sin(a1*k1)*(R-r)*3.0)/k1;
			grad_out(1,1) = (R*(k1*sin(a4*k1-k1*(a1+b1+g1))+k1*sin(k1*(a1+b1+g1))*4.0)-k1*sin(k1*(a1+b1))*(R-r)*3.0)/k1;
			grad_out(2,1) = (R*(k1*sin(a4*k1-k1*(a1+b1+g1))+k1*sin(k1*(a1-g1))*2.0+k1*sin(k1*(a1+b1+g1))*4.0))/k1;
			grad_out(3,1) = (-1)*R*sin(a4*k1+(-1)*(a1+b1+g1)*k1);
			////
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = 1.0/k1;
			grad_out(2,2) = 1.0/k1;
			grad_out(3,2) = -1.0/k1;
			////
		}	
		else if (pathType.compare("CB-TT-BCB") == 0){
			////
			grad_out(0,0) = -(R*(k1*cos(k1*(b1+g1))*-4.0+k1*cos(k1*(-b1+b2+g1))*2.0+k1*cos(k1*(-b1+b2+g4)))+k1*cos(b1*k1)*(R-r)*4.0-k1*cos(k1*(b1-b2))*(R-r)*3.0)/k1;
			grad_out(1,0) = (R*(k1*cos(k1*(b1+g1))*4.0+k1*cos(k1*(-b1+b2+g1))*2.0))/k1;
			grad_out(2,0) = (R*(k1*cos(k1*(-b1+b2+g1))*2.0+k1*cos(k1*(-b1+b2+g4)))-k1*cos(k1*(b1-b2))*(R-r)*3.0)/k1;
			grad_out(3,0) = R*cos((b1+(-1)*b2+(-1)*g4)*k1);
			////
			grad_out(0,1) = (R*(k1*sin(k1*(b1+g1))*4.0+k1*sin(k1*(-b1+b2+g1))*2.0+k1*sin(k1*(-b1+b2+g4)))+k1*sin(k1*(b1-b2))*(R-r)*3.0-k1*sin(b1*k1)*(R-r)*4.0)/k1;
			grad_out(1,1) = (R*(k1*sin(k1*(b1+g1))*4.0-k1*sin(k1*(-b1+b2+g1))*2.0))/k1;
			grad_out(2,1) = -(R*(k1*sin(k1*(-b1+b2+g1))*2.0+k1*sin(k1*(-b1+b2+g4)))+k1*sin(k1*(b1-b2))*(R-r)*3.0)/k1;
			grad_out(3,1) = R*sin((b1+(-1)*b2+(-1)*g4)*k1);
			////
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = 0.0;
			grad_out(2,2) = -1.0/k1;
			grad_out(3,2) = -1.0/k1;
			////
		}		
		else if (pathType.compare("CB-TT-BC") == 0){
			////
			grad_out(0,0) = -(R*k1*cos(k1*(-b1+b2+g1))*2.0-R*k1*cos(k1*(b1+g1))*4.0+k1*cos(b1*k1)*(R-r)*4.0+k1*r*cos(k1*(b1-b4))-k1*cos(k1*(b1-b2))*(R-r)*2.0)/k1;
			grad_out(1,0) = (R*k1*cos(k1*(-b1+b2+g1))*2.0+R*k1*cos(k1*(b1+g1))*4.0)/k1;
			grad_out(2,0) = (R*k1*cos(k1*(-b1+b2+g1))*2.0-k1*cos(k1*(b1-b2))*(R-r)*2.0)/k1;
			grad_out(3,0) = r*cos((b1+(-1)*b4)*k1);
			////
			grad_out(0,1) = (k1*sin(k1*(b1-b2))*(R-r)*2.0+R*k1*sin(k1*(-b1+b2+g1))*2.0+R*k1*sin(k1*(b1+g1))*4.0-k1*sin(b1*k1)*(R-r)*4.0-k1*r*sin(k1*(b1-b4)))/k1;
			grad_out(1,1) = -(R*k1*sin(k1*(-b1+b2+g1))*2.0-R*k1*sin(k1*(b1+g1))*4.0)/k1;
			grad_out(2,1) = -(k1*sin(k1*(b1-b2))*(R-r)*2.0+R*k1*sin(k1*(-b1+b2+g1))*2.0)/k1;
			grad_out(3,1) = r*sin((b1+(-1)*b4)*k1);
			////
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = 0.0;
			grad_out(2,2) = 0.0;
			grad_out(3,2) = -1.0/k1;
			//// 
		}
		else if (pathType.compare("CB-TT-B") == 0){
			////
			grad_out(0,0) = -(R*(k1*cos(k1*(b1+g1))*-4.0+k1*cos(a4*k1-k1*(b1+g1))+k1*cos(k1*(-b1+b2+g1))*2.0)+k1*cos(b1*k1)*(R-r)*3.0-k1*cos(k1*(b1-b2))*(R-r)*2.0)/k1;
			grad_out(1,0) = (R*(k1*cos(k1*(b1+g1))*4.0-k1*cos(a4*k1-k1*(b1+g1))+k1*cos(k1*(-b1+b2+g1))*2.0))/k1;
			grad_out(2,0) = (R*k1*cos(k1*(-b1+b2+g1))*2.0-k1*cos(k1*(b1-b2))*(R-r)*2.0)/k1;
			grad_out(3,0) = R*cos(a4*k1+(-1)*(b1+g1)*k1);
			////
			grad_out(0,1) = (R*(k1*sin(k1*(b1+g1))*4.0+k1*sin(a4*k1-k1*(b1+g1))+k1*sin(k1*(-b1+b2+g1))*2.0)+k1*sin(k1*(b1-b2))*(R-r)*2.0-k1*sin(b1*k1)*(R-r)*3.0)/k1;
			grad_out(1,1) = (R*(k1*sin(k1*(b1+g1))*4.0+k1*sin(a4*k1-k1*(b1+g1))-k1*sin(k1*(-b1+b2+g1))*2.0))/k1;
			grad_out(2,1) = -(k1*sin(k1*(b1-b2))*(R-r)*2.0+R*k1*sin(k1*(-b1+b2+g1))*2.0)/k1;
			grad_out(3,1) = (-1)*R*sin(a4*k1+(-1)*(b1+g1)*k1);
			////
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = 1.0/k1;
			grad_out(2,2) = 0.0;
			grad_out(3,2) = -1.0/k1;
			////
		}	
		else if (pathType.compare("B-TT-BCB") == 0){
			////
			grad_out(0,0) = -(R*(k1*cos(k1*(a2*2.0+b2-g1))*2.0+k1*cos(k1*(a2+b2-g1+g4))-k1*cos(g1*k1)*4.0)+k1*cos(k1*(a2-g1))*(R-r)*3.0-k1*cos(k1*(a2+b2-g1))*(R-r)*3.0)/k1;
			grad_out(1,0) = (R*(k1*cos(k1*(a2*2.0+b2-g1))*4.0+k1*cos(k1*(a2+b2-g1+g4)))+k1*cos(k1*(a2-g1))*(R-r)*3.0-k1*cos(k1*(a2+b2-g1))*(R-r)*3.0)/k1;
			grad_out(2,0) = (R*(k1*cos(k1*(a2*2.0+b2-g1))*2.0+k1*cos(k1*(a2+b2-g1+g4)))-k1*cos(k1*(a2+b2-g1))*(R-r)*3.0)/k1;
			grad_out(3,0) = R*cos((a2+b2+(-1)*g1+g4)*k1);
			////
			grad_out(0,1) = (R*(-k1*sin(g1*k1-k1*(a2+b2+g4))+k1*sin(k1*(a2*2.0+b2-g1))*2.0+k1*sin(g1*k1)*4.0)+k1*sin(k1*(a2-g1))*(R-r)*3.0-k1*sin(k1*(a2+b2-g1))*(R-r)*3.0)/k1;
			grad_out(1,1) = (R*(k1*sin(g1*k1-k1*(a2+b2+g4))-k1*sin(k1*(a2*2.0+b2-g1))*4.0)-k1*sin(k1*(a2-g1))*(R-r)*3.0+k1*sin(k1*(a2+b2-g1))*(R-r)*3.0)/k1;
			grad_out(2,1) = (R*(k1*sin(g1*k1-k1*(a2+b2+g4))-k1*sin(k1*(a2*2.0+b2-g1))*2.0)+k1*sin(k1*(a2+b2-g1))*(R-r)*3.0)/k1;
			grad_out(3,1) = R*sin(g1*k1+(-1)*(a2+b2+g4)*k1);
			////
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = -1.0/k1;
			grad_out(2,2) = -1.0/k1;
			grad_out(3,2) = -1.0/k1;
			////
		}	
		else if (pathType.compare("B-TT-BC") == 0){
			////
			grad_out(0,0) = -(k1*cos(k1*(a2-g1))*(R-r)*3.0+k1*r*cos(k1*(a2+b4-g1))-k1*cos(k1*(a2+b2-g1))*(R-r)*2.0+R*k1*cos(k1*(a2*2.0+b2-g1))*2.0-R*k1*cos(g1*k1)*4.0)/k1;
			grad_out(1,0) = (k1*cos(k1*(a2-g1))*(R-r)*3.0+k1*r*cos(k1*(a2+b4-g1))-k1*cos(k1*(a2+b2-g1))*(R-r)*2.0+R*k1*cos(k1*(a2*2.0+b2-g1))*4.0)/k1;
			grad_out(2,0) = -(k1*cos(k1*(a2+b2-g1))*(R-r)*2.0-R*k1*cos(k1*(a2*2.0+b2-g1))*2.0)/k1;
			grad_out(3,0) = r*cos((a2+b4+(-1)*g1)*k1);
			////
			grad_out(0,1) = (k1*sin(k1*(a2-g1))*(R-r)*3.0+k1*r*sin(k1*(a2+b4-g1))-k1*sin(k1*(a2+b2-g1))*(R-r)*2.0+R*k1*sin(k1*(a2*2.0+b2-g1))*2.0+R*k1*sin(g1*k1)*4.0)/k1;
			grad_out(1,1) = -(k1*sin(k1*(a2-g1))*(R-r)*3.0+k1*r*sin(k1*(a2+b4-g1))-k1*sin(k1*(a2+b2-g1))*(R-r)*2.0+R*k1*sin(k1*(a2*2.0+b2-g1))*4.0)/k1;
			grad_out(2,1) = (k1*sin(k1*(a2+b2-g1))*(R-r)*2.0-R*k1*sin(k1*(a2*2.0+b2-g1))*2.0)/k1;
			grad_out(3,1) = (-1)*r*sin((a2+b4+(-1)*g1)*k1);
			////
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = -1.0/k1;
			grad_out(2,2) = 0.0;
			grad_out(3,2) = -1.0/k1;
			////
		}		
		else if (pathType.compare("B-TT-B") == 0){
			////
			grad_out(0,0) = -(R*(k1*cos(k1*(a2*2.0+b2-g1))*2.0-k1*cos(g1*k1)*4.0+k1*cos(k1*(a4-g1))-k1*cos(k1*(a2+b2-g1))*2.0)+k1*cos(k1*(a2-g1))*(R-r)*2.0+k1*r*cos(k1*(a2+b2-g1))*2.0)/k1;
			grad_out(1,0) = (R*(k1*cos(k1*(a2*2.0+b2-g1))*4.0-k1*cos(k1*(a2+b2-g1))*2.0)+k1*cos(k1*(a2-g1))*(R-r)*2.0+k1*r*cos(k1*(a2+b2-g1))*2.0)/k1;
			grad_out(2,0) = (R*(k1*cos(k1*(a2*2.0+b2-g1))*2.0-k1*cos(k1*(a2+b2-g1))*2.0)+k1*r*cos(k1*(a2+b2-g1))*2.0)/k1;
			grad_out(3,0) = R*cos((a4+(-1)*g1)*k1);
			////
			grad_out(0,1) = (R*(k1*sin(k1*(a2*2.0+b2-g1))*2.0+k1*sin(g1*k1)*4.0+k1*sin(k1*(a4-g1))-k1*sin(k1*(a2+b2-g1))*2.0)+k1*sin(k1*(a2-g1))*(R-r)*2.0+k1*r*sin(k1*(a2+b2-g1))*2.0)/k1;
			grad_out(1,1) = -(R*(k1*sin(k1*(a2*2.0+b2-g1))*4.0-k1*sin(k1*(a2+b2-g1))*2.0)+k1*sin(k1*(a2-g1))*(R-r)*2.0+k1*r*sin(k1*(a2+b2-g1))*2.0)/k1;
			grad_out(2,1) = -(R*(k1*sin(k1*(a2*2.0+b2-g1))*2.0-k1*sin(k1*(a2+b2-g1))*2.0)+k1*r*sin(k1*(a2+b2-g1))*2.0)/k1;
			grad_out(3,1) = (-1)*R*sin((a4+(-1)*g1)*k1);
			////
			grad_out(0,2) = 1.0/k1;
			grad_out(1,2) = 0.0;
			grad_out(2,2) = 0.0;
			grad_out(3,2) = -1.0/k1;
			////
		}				
	}	
	//grad_out.print("grad_out");
	return grad_out;
}





