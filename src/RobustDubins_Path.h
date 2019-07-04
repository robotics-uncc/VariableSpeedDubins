#ifndef ROBUST_DUBINS_PATH_H
#define ROBUST_DUBINS_PATH_H

#include<vector>
#include<string>

#include<RobustDubins_Problem.h>

namespace RobustDubins {

// The RobustDubins::Path object can be used on its own (e.g., to draw a Dubins
// path with known parameters, or it can be the output of 
// the RobustDubins::Solver

class Path {

private: 

	// Consider the following model of a planar mobile robot with position (x,y) 
	// and heading, h.
	// 
	// 	dx/dt = v cos h
	// 	dy/dt = v sin h
	// 	dh/dt = u		where  |u| <= v/R is the turn rate control input
	//					v is the (fixed) speed
	//					R is the minimum turning radius
	//
	// A Dubins path is the shortest path that joins two states (x0, y0, h0) and 
	// (xFinal, yFinal, hFinal) while satisfying the above motion model. 

  // the initial state
  double m_xInitial;
  double m_yInitial;
  double m_hInitial;

	// the final state
	double m_xFinal;
	double m_yFinal;
	double m_hFinal;

  double m_minTurnRadius;

	// The Dubins path was first described in the following paper:
	// 
	// L.E. Dubins "On Curves of Minimal Length with a Constraint on Average 
	// Curvature, and with Prescribed Initial and Final Tangents" American 
	// Journal of Mathematics 79(3), pp.497-516, 1957.
	//	 
	// Dubins showed that any minimum-length path joining the two states 
	// (discussed above) is a member of the following set: 
	// {LRL, RLR, LSL, LSR, RSL, RSR}, where "L" represents a left turning arc, 
	// "R" represents a right turning arc, and "S" represents a straight segment

	// a string representing an element of {LRL, RLR, LSL, LSR, RSL, RSR}
	std::string m_pathType;	

	// The paths of type : LRL and RLR consist of a turn-turn-turn sequence
	// We refer to this as a "BBB" class path where B stands for "bang" control 
	// input since the turn rate is saturated.
	// The paths of the type : LSL, LSR, RSL, RSR consist of turn-straight-turn 
	// sequence. We refer to this as a "BSB" class path where S stands for 
	// "singular" or "straight"

  // a string representing the "class" of the path
	std::string m_pathClass;	 
				
	// To solve the path planning problem we construct all (feasible) 
	// candidates from this set and compare their respective costs. 
	// The lowest cost candidate(s) gives the optimal path.  
	//
	// Let the "a" parameter be the arc-length of the first arc in a Dubins path
	// Let the "b" parameter be the arc-length of the second arc 
	// (or straight segment) in a Dubins path.
	// Let the "c" parameter be the arc-length of the third arc in a Dubins path
  // Note: if the turn radius is unity, then the arc-length of a,b is equal
  // to the angle of the arc subtended (in radians) but in general this is not 
  // the case 

	double m_aParamUnsigned; 
  double m_bParamUnsigned;  
  double m_cParamUnsigned; 
	double m_cost; // cost of the path is the total arc-length

	// to simplify the computation we denote the sign of each arc with the 
	// following symbols
	double m_ki; // initial arc curvature sign, (+) is a left/CCW turn 
	double m_km; // middle arc curvature sign, (+) for straight line
	double m_kf; // final arc curvature sign

  // the solution to the problem is a path type and triple (a,b,c) however
  // some users may prefer the planner to return a series of waypoints
	int m_num_pts; // total number of points to construct the waypoint path
  double m_nom_spacing;

	bool m_userInputNumWpts; // true if user uses set_num_pts() 
  bool m_userInputSpacing; // true if user uses 
	void computePathSpacing(); // computes path spacing if user gives no input 

	// adds the currentState to the state histories
	void pushBackStateHistory(vd currentState); 

	// a candidate Dubins path of a specific pathType can have three possible 
	// statuses:
	// "infeasible" - there is no solution for that pathType and boundary cond.
	// "feasible" - a solution exists, but it is not optimal
	// "optimal" - a solution exists, and it is optimal
	std::string m_solutionStatus;

	// the time and state histories along the path are output if requested
	vd m_xHistory;
	vd m_yHistory;
	vd m_hHistory; 
	
  bool m_pathHistoryComputedFlag; // true if computePathHistory() called 
  bool m_endpointComputedFlag; // true if computeEndpoint() called 

  // utility functions computing x,y,h displacement from a turn or straight seg.
  double delxTurn(const double & psi0, const double & signedAngle);
  double delyTurn(const double & psi0, const double & signedAngle);
  double delhTurn(const double & psi0, const double & signedAngle);
  double delxStraight(const double & psi0, const double & length);
  double delyStraight(const double & psi0, const double & length);

public: 
  // constructor and destruct
	Path();
	virtual ~Path(){};

	// set functions
  // ---------------------------------------------------------------------------
  // when defining a path manually (not output from solver), the following 
  // functions must be set 
	void set_pathType(const std::string & pathType);

	void set_aParamUnsigned( const double & aParamUnsigned );
	void set_bParamUnsigned( const double & bParamUnsigned );
	void set_cParamUnsigned( const double & cParamUnsigned );
	void set_abcParamVectorUnsigned( const vd & abcParamVectorUnsigned );
	void set_abcParamVectorUnsigned( const double & aParamUnsigned, 
                                       const double & bParamUnsigned, 
									                     const double & cParamUnsigned );

  // optional: set one of these (but not both) if a waypoint path is desired 
  void set_num_pts( const int & num_pts );
  void set_spacing( const double & spacing );

  // optional 
	void set_solutionStatus( const std::string & solutionStatus );
  void set_initialState( const double & xInit, 
                         const double & yInit, 
                         const double & hInit );
  void set_minTurnRadius( const double & turnRadius );

	// main functions
  // ---------------------------------------------------------------------------
  void computeEndpoint(); // compute the endpoint given a,b,c parameters 
	void computePathHistory(); // compute states along path given a,b,c parameters
	void print(); // display path properties to screen 
  
  // write commands for plotting 
	void writePathOctavePlotCommands( const std::string & fileName, 
                                    const int & figureNumber, 
										                const std::string & xVarName, 
                                    const std::string & yVarName, 
										                const std::string & style);
  void set_cost(); // should becalled to update cost once a,b,c params defined 



	// get functions
  // ---------------------------------------------------------------------------
	std::string get_pathType(){return m_pathType;};
  std::string get_pathClass(){return m_pathClass;}; 

	std::string get_solutionStatus(){return m_solutionStatus;};
	double get_aParamUnsigned(){return m_aParamUnsigned;};
	double get_bParamUnsigned(){return m_bParamUnsigned;}; 
	double get_cParamUnsigned(){return m_cParamUnsigned;};

	double get_xInitial(){return m_xInitial;};
	double get_yInitial(){return m_yInitial;};
	double get_hInitial(){return m_hInitial;};
	double get_xFinal();
	double get_yFinal();
	double get_hFinal();
	double get_cost(){return m_cost;};


	vd get_timeHistory();
	vd get_xHistory();
	vd get_yHistory();
	vd get_hHistory();

  // additional properties of each segment
  void get_segPts( vd & segX, 
                   vd & segY, 
                   vd & segH );
  void get_orbitParams(vd & orbitCenterX, 
                       vd & orbitCenterY,
                       vd & entryAngles,
                       vd & orbitDir );

}; // class

} // namespace
#endif
