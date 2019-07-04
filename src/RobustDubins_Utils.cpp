#include<RobustDubins_Utils.h>

// custom
#include<RobustDubins_Problem.h>
#include<RobustDubins_Solver.h>
#include<RobustDubins_Path.h>
#include<MathTools.h>


void RobustDubins::DubinsChain( double turnRadius, double nomSpacing, vd xPts, vd yPts, vd hRadPts, 
                                vd & xPath, vd & yPath, vd & hRadPath ){
  int numTargs = xPts.size();
  // iterate for each target, adding waypoints to the output vectors
  for (int i=1; i < numTargs; i++){

	  // create the problem
	  RobustDubins::Problem problemStatement;
    problemStatement.set_minTurningRadius(turnRadius);

	  // create the solver
	  RobustDubins::Solver rds;
    std::vector<double> xCurrentPath, yCurrentPath, hCurrentPath; // 

    // update problem statement
    problemStatement.set_stateInitial(xPts.at(i-1), yPts.at(i-1), 
                                      hRadPts.at(i-1));
    problemStatement.set_stateFinal(xPts.at(i),yPts.at(i), 
                                    hRadPts.at(i));
    problemStatement.print();

    // update solver
    rds.set_problemStatement(problemStatement);
	  rds.solve();

    // print results
    rds.print();

    // get waypoints/states of the current path
    rds.get_optimalWaypointsSetSpacing(xCurrentPath, yCurrentPath, hCurrentPath, nomSpacing);
   
    if ( i > 1 ){ // erase first pointers 
      xCurrentPath.erase( xCurrentPath.begin() );
      yCurrentPath.erase( yCurrentPath.begin() );
      hCurrentPath.erase( hCurrentPath.begin() );
    }

    // append to overall path
    MathTools::append(xPath,xCurrentPath);
    MathTools::append(yPath,yCurrentPath);
    MathTools::append(hRadPath,hCurrentPath);

    RobustDubins::Path optimalPath = rds.get_optimalPath();
    //optimalPath.print(); // prints info about this object
  }

}



std::string RobustDubins::DubinsPathClass(const std::string & pathType){
  std::string pathClass;
	if ( (pathType.compare("LRL") == 0) || (pathType.compare("RLR") == 0) ){
		pathClass = "BBB";
	}
	else if ( (pathType.compare("LSL") == 0) || (pathType.compare("LSR") == 0) 
         || (pathType.compare("RSL") == 0) || (pathType.compare("RSR") == 0)  ){
		pathClass = "BSB";
	}
  #ifndef NDEBUG  
  else {
    throw std::runtime_error("RobustDubins::DubinsPathClass: incorrect path type.");
  }
  #endif
	return pathClass;
}

bool RobustDubins::DubinsPathCurvatureSigns(const std::string & pathType, 
                                            double & ki, double & km, 
                                            double & kf){
	if (pathType.compare("LRL") == 0){
		ki = 1.0;
		km = -1.0;
		kf = 1.0;
	}
	else if (pathType.compare("RLR") == 0){
		ki = -1.0;
		km = 1.0;
		kf = -1.0;
	}
	else if (pathType.compare("LSL") == 0){
		ki = 1.0;
		km = 1.0;
		kf = 1.0;
	}
	else if (pathType.compare("LSR") == 0){
		ki = 1.0;
		km = 1.0;
		kf = -1.0;
	}
	else if (pathType.compare("RSL") == 0){
		ki = -1.0;
		km = 1.0;
		kf = 1.0;
	}
	else if (pathType.compare("RSR") == 0){
		ki = -1.0;
		km = 1.0;
		kf = -1.0;
	}
  #ifndef NDEBUG
  else {
    throw std::runtime_error("RobustDubins::DubinsPathCurvatureSigns: incorrect path type.");
    return false;
  }
  #endif
  return true;
}

bool RobustDubins::computeDubinsPoint(const std::string & pathType, 
		                                     const double & aUnsigned, 
                                         const double & bUnsigned, 
                                         const double & cUnsigned,
                                         const double & R,
                                         const double & xInit,
                                         const double & yInit,
                                         const double & hInit,
                                         const double & arclength,
                                         vd & endpoint){
  // error-checking 
  double tol = 0.001*R;
  #ifndef NDEBUG
  assert(arclength >= 0);
  assert(arclength <= aUnsigned + cUnsigned + bUnsigned + tol);
  assert(R > 0);
  #endif
	std::string pathClass = RobustDubins::DubinsPathClass(pathType);

  // here we determine if the arclength give corresponds to the the first,
  // second, or third segment of the Dubins path. 
  if ( arclength <= aUnsigned ){ // point along first arc either BSB, BBB
    return RobustDubins::computeDubinsEndpoint(pathType, 
                                               arclength, 0.0, 0.0, 
                                               xInit, yInit, hInit, 
                                               R, endpoint);
  }
  if ( pathClass.compare("BBB")==0 ){
    if ( arclength <= aUnsigned + bUnsigned ){
      // point along second arc if BBB
      return RobustDubins::computeDubinsEndpoint(pathType, aUnsigned, 
                                                 arclength - aUnsigned, 
                                                 0.0, xInit, yInit, hInit, 
                                                 R, endpoint);
    }
    else{
      // point along third arc if BBB
      return RobustDubins::computeDubinsEndpoint(pathType, aUnsigned, 
                                       bUnsigned, 
                                       arclength - aUnsigned - bUnsigned, 
                                       xInit, yInit, hInit, R, endpoint);
    }
  }
  else { 
    if ( arclength <= aUnsigned + bUnsigned ){
      // point along second (straight) segment if BSB 
      return RobustDubins::computeDubinsEndpoint(pathType, aUnsigned, 
                                                 arclength - aUnsigned, 
                                                 0.0, xInit, yInit, hInit, 
                                                 R, endpoint);
    }
    else {
      // point along second arc  if BSB 
      return RobustDubins::computeDubinsEndpoint(pathType, aUnsigned, 
                                       bUnsigned, 
                                       arclength - bUnsigned - aUnsigned, 
                                       xInit, yInit, hInit, R, endpoint);
    }
  }
}

bool RobustDubins::computeDubinsEndpoint(const std::string & pathType,
		                                     const double & aUnsigned, 
                                         const double & bUnsigned, 
                                         const double & cUnsigned,
                                         const double & xInit,
                                         const double & yInit,
                                         const double & hInit,
                                         const double & R,
                                         vd & endpoint){
  // error checking 
  #ifndef NDEBUG
  assert(aUnsigned >= 0);
  assert(bUnsigned >= 0);
  assert(cUnsigned >= 0);
  assert(R > 0);
  assert(endpoint.size()==3);
  #endif 
  std::string pathClass = RobustDubins::DubinsPathClass(pathType);
  // compute the endpoint 
	if ( pathClass.compare("BBB") == 0 ){
		RobustDubins::computeDubinsBBBendpoint( pathType, aUnsigned, bUnsigned,
			                                      cUnsigned, xInit, yInit, hInit, 
                                            R, endpoint );
	}
	else if ( pathClass.compare("BSB") == 0 ){
		RobustDubins::computeDubinsBSBendpoint( pathType, aUnsigned, bUnsigned,
			                                      cUnsigned, xInit, yInit, hInit,  
                                            R, endpoint );
	}
  return true;
}

bool RobustDubins::computeDubinsBBBendpoint(const std::string & pathType, 
		                                        const double & aUnsigned, 
                                            const double & bUnsigned, 
                                            const double & cUnsigned,
                                            const double & xInit,
                                            const double & yInit,
                                            const double & hInit,
                                            const double & R,
                                            vd & endpoint){
	double ki, km, kf;
  RobustDubins::DubinsPathCurvatureSigns(pathType, ki, km, kf);
  // here the a,b,c parameters are signed angles not arc-lengths
  double aSigned = ki*aUnsigned/R;
  double bSigned = km*bUnsigned/R;
  double cSigned = kf*cUnsigned/R;
  // analytical expressions 
	endpoint[0] = xInit + R*( 2*ki*sin( aSigned + hInit ) - ki*sin(hInit)
                            - 2*ki*sin( aSigned + bSigned + hInit)  
                            + ki*sin( aSigned + bSigned + cSigned + hInit) );
	endpoint[1] = yInit + R*( ki*cos(hInit) - 2*ki*cos( aSigned + hInit) 
                               + 2*ki*cos( aSigned + bSigned + hInit) 
                               - ki*cos( aSigned + bSigned + cSigned + hInit ));
	endpoint[2] = MathTools::mod( hInit + aSigned + bSigned + cSigned, 2.0*M_PI );
  return true;
}

bool RobustDubins::computeDubinsBSBendpoint(const std::string & pathType, 
		                                        const double & aUnsigned, 
                                            const double & bUnsigned, 
                                            const double & cUnsigned,
                                            const double & xInit,
                                            const double & yInit,
                                            const double & hInit,
                                            const double & R,
                                            vd & endpoint){
	double ki, km, kf;
  RobustDubins::DubinsPathCurvatureSigns(pathType, ki, km, kf);
  // here the a,c parameters are signed angles not arc-lengths
  double aSigned = ki*aUnsigned/R;
  double cSigned = kf*cUnsigned/R;
  // analytical expressions 
	endpoint[0] = xInit + R*ki*sin( aSigned + hInit ) - R*ki*sin(hInit)
                      + bUnsigned*cos( aSigned + hInit ) 
                      + R*kf*( sin( aSigned + cSigned + hInit) 
                               - sin( aSigned + hInit) );
	endpoint[1] = yInit + R*ki*(-cos( aSigned + hInit ) + cos(hInit) ) 
                      + bUnsigned*sin(aSigned + hInit) 
                      + R*kf*(  - cos( aSigned + cSigned + hInit) 
                                + cos( aSigned + hInit) );
	endpoint[2] = MathTools::mod( hInit + aSigned + cSigned , 2.0*M_PI);
	return true;
}
