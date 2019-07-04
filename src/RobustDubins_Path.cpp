#include<limits> 
#include<RobustDubins_Path.h>
#include<RobustDubins_Solver.h>
#include<MathTools.h>

// constructor
RobustDubins::Path::Path(){
  m_xInitial = 0.0;
  m_yInitial = 0.0;
  m_hInitial = 0.0;
  m_xFinal = 0.0;
  m_yFinal = 0.0;
  m_hFinal = 0.0;
  m_pathType = "N/A";
  m_pathClass = "N/A";
  m_aParamUnsigned = 0.0;
  m_bParamUnsigned = 0.0;
  m_cParamUnsigned = 0.0;
  m_cost = std::numeric_limits<double>::infinity();
  m_ki = 0.0;
  m_km = 0.0;
  m_kf = 0.0;
  m_num_pts = 15;
  m_minTurnRadius = 1.0;
  m_userInputNumWpts = false;
  m_userInputSpacing = false;
  m_pathHistoryComputedFlag = false;
  m_endpointComputedFlag = false;
  m_solutionStatus = "N/A";
}

// set functions
void RobustDubins::Path::set_pathType(const std::string & pathType){
  #ifndef NDEBUG
  if ( !pathType.compare("LSL")==0 && !pathType.compare("LSR")==0 && 
       !pathType.compare("RSL")==0 && !pathType.compare("RSR")==0 &&
       !pathType.compare("LRL")==0 && !pathType.compare("RLR")==0    ){
    throw std::runtime_error("RobustDubins::Path::set_pathType: Invalid entry");
  }
  #endif
  m_pathType = pathType;
  if ( pathType.compare("LSL")==0 || pathType.compare("LSR")==0 || 
       pathType.compare("RSL")==0 || pathType.compare("RSR")==0 ){
    m_pathClass = "BSB";
  }
  else {
    m_pathClass = "BBB";
  }
}

void RobustDubins::Path::set_aParamUnsigned(const double & aParamUnsigned){
  #ifndef NDEBUG
  if ( aParamUnsigned < 0 ){
    throw std::runtime_error("RobustDubins::Path::set_aParamUnsigned. Supplied parameter is negative.");
  }
  #endif
  m_aParamUnsigned = aParamUnsigned;
}

void RobustDubins::Path::set_bParamUnsigned(const double & bParamUnsigned){
  #ifndef NDEBUG
  if ( bParamUnsigned < 0 ){
    throw std::runtime_error("RobustDubins::Path::set_bParamUnsigned. Supplied parameter is negative.");
  }
  #endif
  m_bParamUnsigned = bParamUnsigned;
}

void RobustDubins::Path::set_cParamUnsigned(const double & cParamUnsigned){
  #ifndef NDEBUG
  if ( cParamUnsigned < 0 ){
    throw std::runtime_error("RobustDubins::Path::set_cParamUnsigned. Supplied parameter is negative.");
  }
  #endif 
  m_cParamUnsigned = cParamUnsigned;
}

void RobustDubins::Path::set_abcParamVectorUnsigned(
        const vd & abcParamVectorUnsigned){
	set_aParamUnsigned(abcParamVectorUnsigned[0]);
	set_bParamUnsigned(abcParamVectorUnsigned[1]);
	set_cParamUnsigned(abcParamVectorUnsigned[2]);
  set_cost();
}

void RobustDubins::Path::set_abcParamVectorUnsigned( 
                                                 const double & aParamUnsigned, 
		                                             const double & bParamUnsigned, 
                                                 const double & cParamUnsigned){
	set_aParamUnsigned(aParamUnsigned);
	set_bParamUnsigned(bParamUnsigned);
	set_cParamUnsigned(cParamUnsigned);
  set_cost();
}

void RobustDubins::Path::set_initialState( const double & xInit, 
                                           const double & yInit, 
                                           const double & hInit){
  m_xInitial = xInit;
  m_yInitial = yInit;
  #ifndef NDEBUG
  if ( hInit < 0 || hInit > 2.0 *M_PI ){
    throw std::runtime_error("RobustDubins::Path::set_initialState: 0 <= hInit <= 2pi is required (radians) ");
  }
  #endif
  m_hInitial = hInit;
}

void RobustDubins::Path::set_solutionStatus( const std::string & solutionStatus){
  #ifndef NDEBUG
  if ( solutionStatus.compare("infeasible")==1 && 
       solutionStatus.compare("feasible")==1 && 
       solutionStatus.compare("optimal")==1  ){
    printf("RobustDubins::Path::set_solutionStatus received: %s\n",solutionStatus.c_str());
    throw std::runtime_error("RobustDubins::Path::set_solutionStatus: solutionStatus invalid.");
  }
  #endif
  m_solutionStatus = solutionStatus;
};

void RobustDubins::Path::set_minTurnRadius( const double & turnRadius){
  #ifndef NDEBUG
  if ( turnRadius <= 0 ){
    throw std::runtime_error("RobustDubins::Path::set_minTurnRadius: turnRadius > 0 is required");
  }
  #endif
  m_minTurnRadius = turnRadius;
};

void RobustDubins::Path::set_num_pts( const int & num_pts ){
  #ifndef NDEBUG
  if ( num_pts <= 0 ){
    throw std::runtime_error("RobustDubins::Path::set_num_pts: num_pts > 0 is required");
  }
  #endif
  m_num_pts = num_pts;
  m_userInputNumWpts = true;
  m_userInputSpacing = false; 
};

void RobustDubins::Path::set_spacing( const double & spacing ){
  #ifndef NDEBUG
  if ( spacing <= 0 ){
    throw std::runtime_error("RobustDubins::Path::set_spacing: spacing > 0 is required");
  }
  #endif
  m_nom_spacing = spacing;
  m_userInputNumWpts = false;
  m_userInputSpacing = true; 
};


void RobustDubins::Path::set_cost(){
	if (m_solutionStatus.compare("feasible") == 0 || 
			m_solutionStatus.compare("optimal") == 0 ){
		m_cost = m_aParamUnsigned + m_bParamUnsigned + m_cParamUnsigned; 
	}
	else {
		m_cost = std::numeric_limits<double>::infinity();	
	}
  return;
}

// main functions 
void RobustDubins::Path::computeEndpoint(){
  if ( !m_endpointComputedFlag ){
    vd endPt(3); 
    RobustDubins::computeDubinsEndpoint( m_pathType, m_aParamUnsigned,
                                         m_bParamUnsigned, m_cParamUnsigned, 
                                         m_xInitial, m_yInitial, m_hInitial, 
                                         m_minTurnRadius, endPt );
    m_xFinal = endPt.at(0);
    m_yFinal = endPt.at(1);
    m_hFinal = endPt.at(2);
  }
}

void RobustDubins::Path::computePathHistory(){
  // to avoid re-computing 
  if ( m_pathHistoryComputedFlag == false ){
    double arc_length_tot = m_aParamUnsigned + m_bParamUnsigned 
                            + m_cParamUnsigned;
    if ( !m_userInputNumWpts && !m_userInputSpacing ){
      double nom_spacing = m_minTurnRadius/10;
      m_num_pts = std::floor(arc_length_tot/nom_spacing);
    }
    else if ( m_userInputSpacing ){ 
      m_num_pts = std::floor(arc_length_tot/m_nom_spacing); 
    }
    vd arcPts = MathTools::linspace(0, arc_length_tot, m_num_pts);
    vd currentPoint(3);
    for (int i = 0; i < m_num_pts; i++ ){
			  RobustDubins::computeDubinsPoint( m_pathType, 
		                                      m_aParamUnsigned, 
                                          m_bParamUnsigned, 
                                          m_cParamUnsigned,
                                          m_minTurnRadius,
                                          m_xInitial,
                                          m_yInitial,
                                          m_hInitial,
                                          arcPts[i],
                                          currentPoint );
        pushBackStateHistory(currentPoint);
    }
    computeEndpoint();
    m_pathHistoryComputedFlag = true;
  }
}

void RobustDubins::Path::pushBackStateHistory(vd 
		currentPoint){
	m_xHistory.push_back(currentPoint[0]);
	m_yHistory.push_back(currentPoint[1]);
	m_hHistory.push_back(currentPoint[2]);
}

void RobustDubins::Path::print(){
	printf("==============================================\n");
	printf("pathType: %s\n", m_pathType.c_str());
  printf("solutionStatus: %s\n", m_solutionStatus.c_str() );
  printf("turnRadius : %3.3f\n", m_minTurnRadius);
	printf("aParamUnsigned: %3.3f\n", m_aParamUnsigned );
	printf("bParamUnsigned: %3.3f\n", m_bParamUnsigned );
	printf("cParamUnsigned: %3.3f\n", m_cParamUnsigned );
	printf("xInitial: %3.3f\n", m_xInitial );
	printf("yInitial: %3.3f\n", m_yInitial );
	printf("hInitial (deg): %3.3f\n", m_hInitial*180/M_PI );
	printf("xFinal: %3.3f\n", m_xFinal );
	printf("yFinal: %3.3f\n", m_yFinal );
	printf("hFinal (deg): %3.3f\n", m_hFinal*180/M_PI );
	printf("cost: %3.3f\n", m_cost);
	printf("==============================================\n");
}

void RobustDubins::Path::writePathOctavePlotCommands( 
                                                   const std::string & fileName, 
                                                   const int & figureNumber, 
                                                   const std::string & xVarName, 
                                                   const std::string & yVarName, 
                                                   const std::string & style){
  computePathHistory();
  // write path plotting commands
	MathTools::writeVectorData(fileName, m_xHistory, xVarName);
	MathTools::writeVectorData(fileName, m_yHistory, yVarName);
	MathTools::plot1DCurve(fileName, figureNumber, xVarName, yVarName, style);
  MathTools::axisProperty(fileName, "equal");
}

// get functions 
double RobustDubins::Path::get_xFinal(){
  if ( !m_endpointComputedFlag ){
    computeEndpoint();
  }
  return m_xFinal;
};

double RobustDubins::Path::get_yFinal(){
  if ( !m_endpointComputedFlag ){
    computeEndpoint();
  }
  return m_yFinal;
};

double RobustDubins::Path::get_hFinal(){
  if ( !m_endpointComputedFlag ){
    computeEndpoint();
  }
  return m_hFinal;
};

vd RobustDubins::Path::get_xHistory(){
  if ( !m_pathHistoryComputedFlag ){
    computePathHistory();
  }   
  return m_xHistory;
};

vd RobustDubins::Path::get_yHistory(){
  if ( !m_pathHistoryComputedFlag ){
    computePathHistory();
  }   
  return m_yHistory;
};

vd RobustDubins::Path::get_hHistory(){
  if ( !m_pathHistoryComputedFlag ){
    computePathHistory();
  }   
  return m_hHistory;
};

//get_orbitParams
void RobustDubins::Path::get_orbitParams(vd & orbitCenterX, 
                                         vd & orbitCenterY,
                                         vd & entryAngles,
                                         vd & orbitDir){
  // x,y,h coordinates along start/end of segments 
  vd segX, segY, segH;
  get_segPts(segX,segY,segH);
  double ki, km, kf;
  DubinsPathCurvatureSigns(m_pathType, ki, km, kf);
  orbitCenterX.push_back(segX[0] + m_minTurnRadius*cos(segH[0] + ki*M_PI/2.0));
  orbitCenterY.push_back(segY[0] + m_minTurnRadius*sin(segH[0] + ki*M_PI/2.0));
  entryAngles.push_back(segH[0] - ki*M_PI/2.0);
  orbitDir.push_back(ki);
  if ( m_pathClass.compare("BBB")==0 ){
    orbitCenterX.push_back(segX[1]+m_minTurnRadius*cos(segH[1] + km*M_PI/2.0));
    orbitCenterY.push_back(segY[1]+m_minTurnRadius*sin(segH[1] + km*M_PI/2.0));
    entryAngles.push_back(segH[1] - km*M_PI/2.0);
    orbitDir.push_back(km);
  }
  orbitCenterX.push_back(segX[2] + m_minTurnRadius*cos(segH[2] + kf*M_PI/2.0));
  orbitCenterY.push_back(segY[2] + m_minTurnRadius*sin(segH[2] + kf*M_PI/2.0));
  entryAngles.push_back(segH[2] - kf*M_PI/2.0);
  orbitDir.push_back(kf);
}
//get_segPts
void RobustDubins::Path::get_segPts( vd & segX, vd & segY, vd & segH ){
  // start point
  segX.push_back(m_xInitial);
  segY.push_back(m_yInitial);
  segH.push_back(m_hInitial);
  // end of first
  double ki, km, kf;
  DubinsPathCurvatureSigns(m_pathType, ki, km, kf);
  segX.push_back( m_xInitial + delxTurn(m_hInitial, ki*m_aParamUnsigned) );
  segY.push_back( m_yInitial + delyTurn(m_hInitial, ki*m_aParamUnsigned) );
  segH.push_back( m_hInitial + ki*m_aParamUnsigned );
  // end of second
  if ( m_pathClass.compare("BBB")==0 ){
    segX.push_back( segX[1] + delxTurn( segH[1] , km*m_bParamUnsigned) );
    segY.push_back( segY[1] + delyTurn( segH[1] , km*m_bParamUnsigned) );
    segH.push_back( segH[1] + km*m_bParamUnsigned );
  }
  else {
    segX.push_back( segX[1] + delxStraight( segH[1] , km*m_bParamUnsigned) );
    segY.push_back( segY[1] + delyStraight( segH[1] , km*m_bParamUnsigned) );
    segH.push_back( segH[1] );
  }
  // end of third
  segX.push_back( m_xFinal );
  segY.push_back( m_yFinal );
  segH.push_back( m_hFinal );
}


// delxTurn
double RobustDubins::Path::delxTurn( const double & psi0, 
                                     const double & signedAngle){
  return (double)MathTools::sign(signedAngle)*m_minTurnRadius*(sin(psi0+signedAngle)-sin(psi0));
}

// delyTurn
double RobustDubins::Path::delyTurn( const double & psi0, 
                                     const double & signedAngle ){
  return (double)MathTools::sign(signedAngle)*m_minTurnRadius*(-cos(psi0+signedAngle)+cos(psi0));
}

// delhTurn
double RobustDubins::Path::delhTurn( const double & psi0, 
                                     const double & signedAngle ){
  return psi0 + signedAngle;
}

// delxStraight
double RobustDubins::Path::delxStraight( const double & psi0, 
                                         const double & length ){
  return length*cos(psi0);
}

// delyStraight
double RobustDubins::Path::delyStraight( const double & psi0, 
                                         const double & length ){
  return length*sin(psi0);
}




