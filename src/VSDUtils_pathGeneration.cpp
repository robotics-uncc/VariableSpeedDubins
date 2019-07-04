// VSDUtils_pathGeneration.cpp
// Last Modefied: 18-Jan-2016, Artur Wolek

// custom
#include<VSDUtils.h>
#include<MathTools.h>
#include<MathToolsArma.h>

// S
arma::vec VarSpeedDubins::S(double psi, double L){
	double delx = L*cos(psi);
	double dely = L*sin(psi);
	double delpsi = 0;
	arma::vec delVector = {delx, dely, delpsi}; 
	return delVector;
}

// BCB
arma::vec VarSpeedDubins::BCB(double aUnsigned, double bUnsigned, 
                              double gUnsigned, double psi0, double sgnK, 
                              double R, double r){
	// set sense of turn
	double a = aUnsigned*sgnK;
	double b = bUnsigned*sgnK;
	double g = gUnsigned*sgnK;
	// compute displacements
	double delx = sgnK*(R*(sin(a+b+g+psi0)-sin(psi0))
                +(R-r)*(sin(a+psi0)-sin(psi0+a+b)));
	double dely = sgnK*(R*(cos(psi0)-cos(a+b+g+psi0))
                +(R-r)*(cos(psi0+a+b)-cos(a+psi0)));
	double delpsi = (a+b+g);
	arma::vec delVector = {delx, dely, delpsi}; 
	return delVector;
}

// Spath
arma::mat VarSpeedDubins::Spath(double psi0, double L, double nomSpacing){
	int numPts = MathTools::safeCurveSpacing(L, nomSpacing);
	// generate points
	arma::vec Lvec = arma::linspace<arma::vec>(0,L,numPts);
	arma::vec delx = Lvec * cos(psi0);
	arma::vec dely = Lvec * sin(psi0);
	arma::vec delpsi = arma::zeros<arma::vec>(arma::size(Lvec));
	return MathTools::join_horiz(delx, dely, delpsi);
}

// BCBpath
arma::mat VarSpeedDubins::BCBpath(double aUnsigned, 
	double bUnsigned, double gUnsigned, double psi0, double sgnK, double R, 
	double r, double nomSpacing){
	
	// set sense of turn (alpha, beta, gamma)
	double a = aUnsigned*sgnK;
	double b = bUnsigned*sgnK;
	double g = gUnsigned*sgnK;

	// given the nominal spacing, compute number of pts per segment
	int numPts_a = MathTools::safeCurveSpacing(aUnsigned, nomSpacing);
	int numPts_b = MathTools::safeCurveSpacing(bUnsigned, nomSpacing);
	int numPts_g = MathTools::safeCurveSpacing(gUnsigned, nomSpacing);

	// initialize
	arma::vec aVec, bVec, gVec;
	arma::vec delxa, delxb, delxg;
	arma::vec delya, delyb, delyg;
	arma::vec delpsia, delpsib, delpsig;
	
	// compute displacement for alpha segment
	if (numPts_a > 1){ // only for number of pts greater than 1
		aVec = arma::linspace<arma::vec>(0,a,numPts_a);
		delxa = sgnK*R*(sin(aVec+psi0)-sin(psi0));
		delya = sgnK*R*(cos(psi0)-cos(aVec+psi0));
		delpsia = aVec;	
	}
	
	// beta segment
	if (numPts_b > 1){
		double psiBeta = psi0 + a;
		bVec = arma::linspace<arma::vec>(0,b,numPts_b);
		delxb = sgnK*r*(sin(bVec+psiBeta)-sin(psiBeta));
		delyb = sgnK*r*(cos(psiBeta)-cos(bVec+psiBeta));
		delpsib = bVec;
	}
	// gamma segment
	if (numPts_g > 1){
		double psiGamma = psi0 + a + b;
		gVec = arma::linspace<arma::vec>(0,g,numPts_g);
		delxg = sgnK*R*(sin(gVec+psiGamma)-sin(psiGamma));
		delyg = sgnK*R*(cos(psiGamma)-cos(gVec+psiGamma));
		delpsig = gVec;
	}
	// compile path history
	arma::vec xhist, yhist, psihist;
	//std::cout << "computing displacements" << std::endl;
	displacementsToPath(xhist, delxa, delxb, delxg);
	displacementsToPath(yhist, delya, delyb, delyg);
	displacementsToPath(psihist, delpsia, delpsib, delpsig);
	return MathTools::join_horiz(xhist, yhist, psihist);
}



// pathEndpoint
arma::vec VarSpeedDubins::pathEndpoint(arma::vec paramsShort, 
	std::string pathClass, std::string pathType, std::string orientation, 
	double R, double r){
	// initialize curvature parameters
	double k1, k2;
	// set curvature parameters
	determineCurvatureParams(pathClass, orientation, k1, k2);
	// expand to long parameter form
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
	// endpoint
	arma::vec endpoint;
	// find endpoint
	if ( pathClass.compare("TST") == 0 ){
		arma::vec firstTurn = BCB(a1, b1, g1, 0.0, k1, R, r);
		arma::vec straightSegment = S(k1*(a1+b1+g1),L);
		arma::vec secondTurn = BCB(a2, b2, g2, k1*(a1+b1+g1), k2, R, r);
		endpoint =  firstTurn + straightSegment + secondTurn;
	}
	else if ((pathClass.compare("TT") == 0)){
		arma::vec firstTurn = BCB(a1, b1, g1, 0.0, k1, R, r);
		arma::vec secondTurn = BCB(a2, b2, g2, k1*(a1+b1+g1), k2, R, r);
		endpoint = firstTurn + secondTurn;
	}
	else if ( (pathClass.compare("TTT") == 0) ){
		arma::vec firstTurn = BCB(a1, b1, g1, 0.0, k1, R, r);
		double psiTemp = k1*(a1+b1+g1);
		arma::vec secondTurn = BCB(a2, b2, g2, psiTemp, -k1, R, r);
		psiTemp = psiTemp - k1*(a2+b2+g2);
		arma::vec thirdTurn = BCB(a3, b3, g3, psiTemp, k1, R, r);
		endpoint =  firstTurn + secondTurn + thirdTurn;
	}		
	else if ( (pathClass.compare("TTTT") == 0) ){
		arma::vec firstTurn = BCB(a1, b1, g1, 0.0, k1, R, r);
		double psiTemp = k1*(a1+b1+g1);
		arma::vec secondTurn = BCB(a2, b2, g2, psiTemp, -k1, R, r);
		psiTemp = psiTemp - k1*(a2+b2+g2);
		arma::vec thirdTurn = BCB(a3, b3, g3, psiTemp, k1, R, r);
		psiTemp = psiTemp + k1*(a3+b3+g3);
		arma::vec fourthTurn = BCB(a4, b4, g4, psiTemp, -k1, R, r);
		endpoint = firstTurn + secondTurn + thirdTurn + fourthTurn;
	}	
	// wrap angle
	endpoint(2) = MathTools::mod(endpoint(2), 2.0*M_PI);
	return endpoint;
}

// pathHistory
void VarSpeedDubins::pathHistory(arma::mat &pathHistory, 
	                               arma::vec paramsShort, std::string pathClass, 
                                 std::string pathType, std::string orientation, 
                                 double R, double r, double nomSpacing){
	// initialize curvature parameters
	double k1, k2;
	// set curvature parameters
	determineCurvatureParams(pathClass, orientation, k1, k2);
	// expand to long parameter form
	arma::vec paramsLong = VarSpeedDubins::expandParamsShort(paramsShort, 
                                                           pathClass, 
                                                           pathType);
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

	// endpoint
	double arcLength = arma::sum(paramsLong);
	// find endpoint
	arma::mat firstTurn, secondTurn, thirdTurn, fourthTurn, straightSegment;
	if ( pathClass.compare("TST") == 0 ){
		firstTurn = BCBpath(a1, b1, g1, 0.0, k1, R, r, nomSpacing);
		straightSegment = Spath(k1*(a1+b1+g1),L, nomSpacing);
		secondTurn = BCBpath(a2, b2, g2, k1*(a1+b1+g1), k2, R, r, nomSpacing);
		displacementsToPath(pathHistory, firstTurn, straightSegment, secondTurn); 
	}
	else if ((pathClass.compare("TT") == 0)){
		firstTurn = BCBpath(a1, b1, g1, 0.0, k1, R, r, nomSpacing);
		secondTurn = BCBpath(a2, b2, g2, k1*(a1+b1+g1), k2, R, r, nomSpacing);
		displacementsToPath(pathHistory, firstTurn, secondTurn); 	
	}	
	else if ( (pathClass.compare("TTT") == 0) ){			
		firstTurn = BCBpath(a1, b1, g1, 0.0, k1, R, r, nomSpacing);
		double psiTemp = k1*(a1+b1+g1);
		secondTurn = BCBpath(a2, b2, g2, psiTemp, -k1, R, r, nomSpacing);
		psiTemp = psiTemp - k1*(a2+b2+g2);
		thirdTurn = BCBpath(a3, b3, g3, psiTemp, k1, R, r, nomSpacing);
		displacementsToPath(pathHistory, firstTurn, secondTurn, thirdTurn);
	}		
	else if ( (pathClass.compare("TTTT") == 0) ){	
		firstTurn = BCBpath(a1, b1, g1, 0.0, k1, R, r, nomSpacing);
		double psiTemp = k1*(a1+b1+g1);
		secondTurn = BCBpath(a2, b2, g2, psiTemp, -k1, R, r, nomSpacing);
		psiTemp = psiTemp - k1*(a2+b2+g2);
		thirdTurn = BCBpath(a3, b3, g3, psiTemp, k1, R, r, nomSpacing);
		psiTemp = psiTemp + k1*(a3+b3+g3);
		fourthTurn = BCBpath(a4, b4, g4, psiTemp,-k1,R,r, nomSpacing);
		displacementsToPath(pathHistory, firstTurn, secondTurn, 
			                  thirdTurn, fourthTurn);
	}	
	return;
}


// extremalSwitches
arma::mat VarSpeedDubins::extremalSwitches(arma::vec paramsShort, 
	                                         std::string pathClass, 
                                           std::string pathType, 
                                           std::string orientation, 
	                                         double R, double r){
  //std::cout << "extremalSwitches start" << std::endl;
	// initialize curvature parameters
	double k1, k2;
	// set curvature parameters
	determineCurvatureParams(pathClass, orientation, k1, k2);
	// expand to long parameter form
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
	// extremalSwitches points
  arma::mat extremalSwitches;
	if ( pathClass.compare("TST") == 0 ){
    extremalSwitches.set_size(3,7);
    extremalSwitches.col(0) = BCB(a1, 0.0, 0.0, 0.0, k1, R, r);
    extremalSwitches.col(1) = BCB(a1,  b1, 0.0, 0.0, k1, R, r);
    extremalSwitches.col(2) = BCB(a1,  b1,  g1, 0.0, k1, R, r);
		arma::vec firstTurn =     BCB(a1,  b1,  g1, 0.0, k1, R, r);
		arma::vec straightSegment = S(k1*(a1+b1+g1),L);
    extremalSwitches.col(3) = firstTurn + straightSegment;
    extremalSwitches.col(4) = BCB(a2, 0.0, 0.0, k1*(a1+b1+g1), k2, R, r) 
                              + firstTurn + straightSegment;
    extremalSwitches.col(5) = BCB(a2,  b2, 0.0, k1*(a1+b1+g1), k2, R, r) 
                              + firstTurn + straightSegment;
    extremalSwitches.col(6) = BCB(a2,  b2,  g2, k1*(a1+b1+g1), k2, R, r)  
                              + firstTurn + straightSegment;
	}
	else if ((pathClass.compare("TT") == 0)){
    extremalSwitches.set_size(3,6);
    extremalSwitches.col(0) = BCB(a1, 0.0, 0.0, 0.0, k1, R, r);
    extremalSwitches.col(1) = BCB(a1,  b1, 0.0, 0.0, k1, R, r);
    extremalSwitches.col(2) = BCB(a1,  b1,  g1, 0.0, k1, R, r);
		arma::vec firstTurn =     BCB(a1,  b1,  g1, 0.0, k1, R, r);
    extremalSwitches.col(3) = BCB(a2, 0.0, 0.0, k1*(a1+b1+g1), k2, R, r) 
                              + firstTurn;
    extremalSwitches.col(4) = BCB(a2,  b2, 0.0, k1*(a1+b1+g1), k2, R, r) 
                              + firstTurn;
    extremalSwitches.col(5) = BCB(a2,  b2,  g2, k1*(a1+b1+g1), k2, R, r) 
                              + firstTurn;
	}
	else if ( (pathClass.compare("TTT") == 0) ){
    extremalSwitches.set_size(3,9);
    extremalSwitches.col(0) = BCB(a1, 0.0, 0.0, 0.0, k1, R, r);
    extremalSwitches.col(1) = BCB(a1,  b1, 0.0, 0.0, k1, R, r);
    extremalSwitches.col(2) = BCB(a1,  b1,  g1, 0.0, k1, R, r);
		arma::vec firstTurn =     BCB(a1,  b1,  g1, 0.0, k1, R, r);
    extremalSwitches.col(3) = BCB(a2, 0.0, 0.0, k1*(a1+b1+g1), -k1, R, r) 
                              + firstTurn;
    extremalSwitches.col(4) = BCB(a2,  b2, 0.0, k1*(a1+b1+g1), -k1, R, r) 
                              + firstTurn;
    extremalSwitches.col(5) = BCB(a2,  b2,  g2, k1*(a1+b1+g1), -k1, R, r) 
                              + firstTurn;
		arma::vec secondTurn =    BCB(a2,  b2,  g2, k1*(a1+b1+g1), -k1, R, r);
		double psiTemp = k1*(a1+b1+g1) - k1*(a2+b2+g2);
    extremalSwitches.col(6) = BCB(a3, 0.0, 0.0, psiTemp, k1, R, r)
                              + firstTurn + secondTurn;
    extremalSwitches.col(7) = BCB(a3,  b3, 0.0, psiTemp, k1, R, r)
                              + firstTurn + secondTurn;
    extremalSwitches.col(8) = BCB(a3,  b3,  g3, psiTemp, k1, R, r)
                              + firstTurn + secondTurn;
	}		
	else if ( (pathClass.compare("TTTT") == 0) ){
    extremalSwitches.set_size(3,12);
    // first turn
    extremalSwitches.col(0) = BCB(a1, 0.0, 0.0, 0.0, k1, R, r);
    extremalSwitches.col(1) = BCB(a1,  b1, 0.0, 0.0, k1, R, r);
    extremalSwitches.col(2) = BCB(a1,  b1,  g1, 0.0, k1, R, r);
		arma::vec firstTurn =     BCB(a1,  b1,  g1, 0.0, k1, R, r);
    // second turn
    extremalSwitches.col(3) = BCB(a2, 0.0, 0.0, k1*(a1+b1+g1), -k1, R, r) 
                              + firstTurn;
    // 
    extremalSwitches.col(4) = BCB(a2,  b2, 0.0, k1*(a1+b1+g1), -k1, R, r) 
                              + firstTurn;
    //
    extremalSwitches.col(5) = BCB(a2,  b2,  g2, k1*(a1+b1+g1), -k1, R, r) 
                              + firstTurn;
    //
		arma::vec secondTurn =    BCB(a2,  b2,  g2, k1*(a1+b1+g1), -k1, R, r);
    // third turn
		double psiTemp = k1*(a1+b1+g1) - k1*(a2+b2+g2);
    extremalSwitches.col(6) = BCB(a3, 0.0, 0.0, psiTemp, k1, R, r)
                              + firstTurn + secondTurn;
    extremalSwitches.col(7) = BCB(a3,  b3, 0.0, psiTemp, k1, R, r)
                              + firstTurn + secondTurn;
    extremalSwitches.col(8) = BCB(a3,  b3,  g3, psiTemp, k1, R, r)
                              + firstTurn + secondTurn;
		arma::vec thirdTurn =     BCB(a3,  b3,  g3, psiTemp, k1, R, r);
    // fourth turn
		psiTemp = k1*(a1+b1+g1) - k1*(a2+b2+g2) + k1*(a3+b3+g3);
    extremalSwitches.col(9) =  BCB(a4, 0.0, 0.0, psiTemp, -k1, R, r)
                              + firstTurn + secondTurn + thirdTurn;
    extremalSwitches.col(10) = BCB(a4,  b4, 0.0, psiTemp, -k1, R, r)
                              + firstTurn + secondTurn + thirdTurn;
    extremalSwitches.col(11) = BCB(a4,  b4,  g4, psiTemp, -k1, R, r)
                              + firstTurn + secondTurn + thirdTurn;
	}	
  //std::cout << "extremalSwitches end" << std::endl;
	return extremalSwitches;
}

