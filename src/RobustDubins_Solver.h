#ifndef ROBUST_DUBINS_SOLVER_H
#define ROBUST_DUBINS_SOLVER_H

#include<vector>
#include<string>

#include<RobustDubins_Problem.h>
#include<RobustDubins_Path.h>
#include<RobustDubins_Utils.h>

namespace RobustDubins {

class Solver{

private:
	// problem statement
  RobustDubins::Problem m_problemStatementInput; // may not be scaled
	RobustDubins::Problem m_problemStatementNorm; // normalized problem 

	// solver tolerances
	double m_distanceErrorTolerance;
	double m_thetaErrorTolerance;
	double m_costTolerance;

	// private variables
	vd m_costVector; // cost of all candidate paths 
	int m_numMinCostPaths; // number of min cost paths (there can be more than 1)
	vi m_minCostPaths; // ID values of the min cost paths 
	int m_optimalSolutionID; // assume that: 0-LSL, 1-LSR, 2-RSL, 3-RSR, 4-LRL, 
                           //              5-RLR, 6-LSR/RSR, 7-LRL/RLR
	std::string m_optimalSolutionType; // the string corresponding to the above ID


	double m_cth; // cos (theta) for conveneince 
	double m_sth;
  double m_optCost; // cost of the optimal path

  // flags
  bool m_psNormalizedFlag; // true if the prob. statement req. normalization

	// private functions
  void normalizeProblem(); // defines m_problemStatementNorm if the input 
                           // problem m_problemStatementInput is non-standard
	void solveBBB(); // solve for the BBB type paths 
	void solveBSB();	// solve for the BSB type paths 

  // during the solution procedure we check if a candidate path reaches the goal
	bool checkCandidateEndpoint(const std::string & pathType, 
                              const double & aCandidate, 
			                        const double & bCandidate, 
                              const double & cCandidate );
	bool compareEndpoints( const vd & testState, 
			                   const vd & trueState );

  // check if a BBB candidate satisfies necessary conditions for optimality 
	bool checkBBBconditions( const double & aCandidate, 
                           const double & bCandidate, 
			                     const double & cCandidate);

  void transformSolutions(); // the solve procedures operate on the normalized 
                             // problem. Once that is complete, we convert
                             // path parameters to the true problem. 

	void compareCandidates(); // finally, the feasible candidates are compared 
                            // and the lowest cost solution is declared optimal 

	void determineSolutionType(); // defines m_optimalSolutionID 


	// solver output
	RobustDubins::Path m_LSL, m_LSR, m_RSL, m_RSR, m_LRL, m_RLR;

	// vector of pointers to above listed paths
	std::vector< RobustDubins::Path* > m_solnPtrs = { &m_LSL, &m_LSR, &m_RSL, 
                                                    &m_RSR, &m_LRL, &m_RLR };
	
public:

	// constructor and destructor
	Solver();
	virtual ~Solver(){};

	// required: set functions
	void set_problemStatement( RobustDubins::Problem & problemStatement );

  // optional: set tolerances
	void set_distanceErrorTolerance(const double & distanceErrorTolerance);
	void set_thetaErrorTolerance(const double & thetaErrorTolerance);
	void set_costTolerance(const double & costTolerance);

	// main functions
	void solve();
  void print();
	void writeOctaveCommandsToPlotSolution( const std::string & fileName, 
                                          const int & figNum);

	// get functions
	RobustDubins::Path get_RLR(){return m_RLR;};
	RobustDubins::Path get_LRL(){return m_LRL;};
	RobustDubins::Path get_LSL(){return m_LSL;};
	RobustDubins::Path get_LSR(){return m_LSR;};
	RobustDubins::Path get_RSL(){return m_RSL;};
	RobustDubins::Path get_RSR(){return m_RSR;};
	RobustDubins::Path get_optimalPath();
  double get_optCost(){return m_optCost;};

	int get_optimalSolutionID(){return m_optimalSolutionID;};
  void get_optimalWaypoints(vd &x, vd &y, vd &h);
  void get_optimalWaypointsSetSpacing( vd &x, vd &y, vd &h, double spacing);

}; // class

} // namespace

#endif
