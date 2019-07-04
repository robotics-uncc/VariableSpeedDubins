// Problem.h
// Defines a planar path planning problem.
// Last Modefied: 31-Jan-2015, Artur Wolek

#ifndef VSDPROBLEM_H
#define VSDPROBLEM_H

// C++ headers
#include<vector>
#include<string>

// external headers
#include<armadillo>

// custom
#include<VSDPath.h>
#include<VSDUtils.h>

namespace VarSpeedDubins {

class ProblemStatement {

private:

	// inputs
	double m_xInitial;
	double m_yInitial;
	double m_hInitial;

	double m_xFinal; // Final x position
	double m_yFinal; // Final y position
	double m_hFinal; // Final theta (heading) angle

	// inputs 
//	double m_minSpeed;
//	double m_maxSpeed;
//	double m_maxTurnRate;
	double m_R;
	double m_r;

public:

	// constructor and destructor
	ProblemStatement();	
	virtual ~ProblemStatement();

	
	// set functions
	void set_xInitial(double xInitial){m_xInitial = xInitial;}; 
	void set_yInitial(double yInitial){m_yInitial = yInitial;};
	void set_hInitial(double hInitial){m_hInitial = hInitial;};
	void set_stateInitial(arma::colvec stateInitial);
	void set_stateInitial(double xInitial, double yInitial, double hInitial);

	void set_xFinal(double xFinal){m_xFinal = xFinal;}; 
	void set_yFinal(double yFinal){m_yFinal = yFinal;};
	void set_hFinal(double hFinal){m_hFinal = hFinal;};
	void set_stateFinal(arma::colvec stateFinal);
	void set_stateFinal(double xFinal, double yFinal, double hFinal);

//	void set_minSpeed(double minSpeed){m_minSpeed = minSpeed;};
//	void set_maxSpeed(double maxSpeed){m_maxSpeed = maxSpeed;};	
//	void set_maxTurnRate(double maxTurnRate){m_maxTurnRate = maxTurnRate;};
	void set_R(double R){m_R = R;};
	void set_r(double r){m_r = r;};
  void set_turnRadii(double R, double r){m_R = R; m_r = r;};

	// normalize problem functions
	void print();

	// get functions
	double get_xInitial(){return m_xInitial;}; /// User-input final x position
	double get_yInitial(){return m_yInitial;};
	double get_hInitial(){return m_hInitial;};
	arma::colvec get_stateInitial();

	double get_xFinal(){return m_xFinal;}; /// User-input final x position
	double get_yFinal(){return m_yFinal;};
	double get_hFinal(){return m_hFinal;};
	arma::colvec get_stateFinal();

	double get_R(){return m_R;};
	double get_r(){return m_r;};
//	double get_minSpeed(){return m_minSpeed;};
//	double get_maxSpeed(){return m_maxSpeed;};
//	double get_maxTurnRate(){return m_maxTurnRate;};
}; // class

} //namespace 

#endif


