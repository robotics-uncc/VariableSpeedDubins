#include<RobustDubins_Problem.h>
#include<assert.h>
#include<cmath>

// constructor
RobustDubins::Problem::Problem(){
  initialize();
}

void RobustDubins::Problem::initialize(){
  m_xInitial = 0.0;
  m_yInitial = 0.0;
  m_hInitial = 0.0;
  m_xFinal = 0.0;
  m_yFinal = 0.0;
  m_hFinal = 0.0;
  m_minTurningRadius = 1.0;
  m_startPtInputFlag = false;
  m_minTurningRadiusInputFlag = false;
}

// set functions
void RobustDubins::Problem::set_stateFinal( const double & xFinal, 
                                            const double & yFinal, 
                                            const double & hFinal ){
	set_xFinal(xFinal);
	set_yFinal(yFinal);
	set_hFinal(hFinal);
}

void RobustDubins::Problem::set_stateFinal( const vd & stateFinal ){
	set_xFinal(stateFinal[0]);
	set_yFinal(stateFinal[1]);
	set_hFinal(stateFinal[2]);
}

void RobustDubins::Problem::set_xFinal( const double & xFinal ){
  m_xFinal = xFinal;
  m_endPtDefined = true;
}

void RobustDubins::Problem::set_yFinal( const double & yFinal ){
  m_yFinal = yFinal;
  m_endPtDefined = true;
}

void RobustDubins::Problem::set_hFinal( const double & hFinal ){
  #ifndef NDEBUG
  if ( (hFinal > 2.0*M_PI) || (hFinal < 0.0 ) ){
    throw std::runtime_error("RobustDubins::Problem::set_hFinal: Invalid heading.");
  }
  #endif
  m_hFinal = hFinal;
  m_endPtDefined = true;
}

void RobustDubins::Problem::set_stateInitial(const double & xInitial, 
                                             const double & yInitial, 
                                             const double & hInitial ){
	set_xInitial(xInitial);
	set_yInitial(yInitial);
	set_hInitial(hInitial); 
}

void RobustDubins::Problem::set_stateInitial( const vd & stateInitial ){
	set_xInitial(stateInitial[0]);
	set_yInitial(stateInitial[1]);
	set_hInitial(stateInitial[2]);
}

void RobustDubins::Problem::set_xInitial( const double & xInitial ){
  m_xInitial = xInitial;
  m_startPtInputFlag = true;
}; 

void RobustDubins::Problem::set_yInitial( const double & yInitial ){
  m_yInitial = yInitial;
  m_startPtInputFlag = true;
};

void RobustDubins::Problem::set_hInitial( const double & hInitial ){
  #ifndef NDEBUG
  if ( (hInitial > 2.0*M_PI) || (hInitial < 0.0 ) ){
    throw std::runtime_error("RobustDubins::Problem::set_hInitial: Invalid heading.");
  }
  #endif
  m_hInitial = hInitial;
  m_startPtInputFlag = true;
};

void RobustDubins::Problem::set_minTurningRadius( const double & minTurningRadius){
  #ifndef NDEBUG
  if ( minTurningRadius < 0 ){
    throw std::runtime_error("RobustDubins::Problem::set_minTurningRadius: Invalid entry.");
  }
  #endif 
	m_minTurningRadius = minTurningRadius;
  m_minTurningRadiusInputFlag = true;
};

// main functions 
bool RobustDubins::Problem::isDefined(){
  return m_endPtDefined;
}

void RobustDubins::Problem::print(){
	printf("=========================================\n");
	printf("Dubins Problem Statement \n");
	printf("=========================================\n");
  printf("(xInit, yInit, hInitDeg) : \t(%12.12f, %12.12f, %12.12f) \n",
            m_xInitial, m_yInitial, m_hInitial*180/M_PI);
  printf("(xFinal, yFinal, hFinalDeg) : \t(%12.12f, %12.12f, %12.12f) \n",m_xFinal, 
            m_yFinal, m_hFinal*180/M_PI);
  printf("Min. Turning Radius : \t\t%3.3f \n", m_minTurningRadius);
}

// get functions
vd RobustDubins::Problem::get_stateInitial(){
  vd stateInitial = {m_xInitial, m_yInitial, m_hInitial};
  return stateInitial;
}

double RobustDubins::Problem::get_xFinal(){
  #ifndef NDEBUG
  if ( !m_endPtDefined ){
    throw std::runtime_error("RobustDubins::Problem::get_xFinal: Final state not defined.");
  }
  #endif
  return m_xFinal;
};

double RobustDubins::Problem::get_yFinal(){
  #ifndef NDEBUG
  if ( !m_endPtDefined ){
    throw std::runtime_error("RobustDubins::Problem::get_yFinal: Final state not defined.");
  }
  #endif 
  return m_yFinal;
};

double RobustDubins::Problem::get_hFinal(){
  #ifndef NDEBUG
  if ( !m_endPtDefined ){
    throw std::runtime_error("RobustDubins::Problem::get_hFinal: Final state not defined.");
  }
  #endif 
  return m_hFinal;
};

vd RobustDubins::Problem::get_stateFinal(){
  #ifndef NDEBUG
  if ( !m_endPtDefined ){
    throw std::runtime_error("RobustDubins::Problem::get_stateFinal: Final state not defined.");
  }
  #endif 
  vd stateFinal = {m_xFinal, m_yFinal, m_hFinal};
  return stateFinal;
}
