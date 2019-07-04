#ifndef VSD_SOLVER
#define VSD_SOLVER

// c++ standard
#include<string>

// external
#include<armadillo>

#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"

// custom
#include<VSDUtils.h>
#include<VSDPath.h>
#include<VSDProblem.h>
#include<VSDNLP.h>

// custom other projects
#include<RobustDubins.h>
#include<MathTools.h>
#include<MathToolsArma.h>

// standard
#include <cassert>

// custom

namespace VarSpeedDubins {

class Solver {

private: 

  // Dubins-type paths
  VarSpeedDubins::ProblemStatement m_vsdprob;
  RobustDubins::Problem m_rdprob_bang;
  RobustDubins::Problem m_rdprob_corn;
  RobustDubins::Solver m_rdsolver_bang;
	RobustDubins::Solver m_rdsolver_corn;
  int m_numPathsDubins;
  std::vector< RobustDubins::Path > m_candPathsDubins;

  // solver options
  int m_numSamplesMax;
  double m_Lmax;
  double m_equalCostTolerance;

  // initial guess options
  std::string m_samplingAlgorithm;
  arma::vec m_paramGuess;
  bool m_paramGuessSuppliedFlag;

  // ipopt nlp
	Ipopt::SmartPtr<VSDNLP> m_vsdnlp;
  void set_NLP_options(SmartPtr<IpoptApplication> &app);

  // solver internal variables
  int m_numPaths;
  std::vector< std::vector<std::string> > m_candList;
  std::vector< VarSpeedDubins::Path > m_candPaths;

  // solver functions
  void sortFeasiblePaths();

  // solver output
	int m_solveStatus; // -1 unsolved, 0 failed, 1 solved
  arma::vec m_candCosts;
  arma::vec m_solveStatusAll;
  // acommodate multipel or  single "optimal" solution
  std::vector<int> m_optPathTypes;
  arma::mat m_optPathParams;
  arma::vec m_optPathCosts;
  double m_multipleSolnFlag;



public: 
  // constructor and destructor
	Solver();
	virtual ~Solver(){};
  
	// main functions
	int solveGivenNLP(); // solve m_vsdnlp
  void solveAll();	
  void printAll();
  void plotAll(std::string fileName);

	// set functions
  void set_problemStatement(VarSpeedDubins::ProblemStatement vsdprob);
  void set_Lmax(double Lmax){m_Lmax = Lmax;};
  void set_numSamplesMax(int numSamplesMax){
    m_numSamplesMax = numSamplesMax;};
  void set_samplingAlgorithm(std::string samplingAlgorithm){
    m_samplingAlgorithm = samplingAlgorithm;};

  // paramGuess is only supplied if the algorithm should not auto-generate
  // it's own initial guess from one of the available sampling methods
  void set_paramGuess(arma::vec paramGuess){
    m_paramGuess = paramGuess;
    m_paramGuessSuppliedFlag = true;};

  // a NLP should be provided if solveGivenNLP() will be used. 
  // the default is for the solver to solve for all candidate paths
	void set_vsdnlp(SmartPtr<VSDNLP> vsdnlp){m_vsdnlp = vsdnlp;};

	// get functions
	int get_solveStatus(){return m_solveStatus;};
  arma::vec get_candCosts(){return m_candCosts;};

  // outputs related to optimal soln
  std::vector<int> get_optPathTypes(){return m_optPathTypes;};
  arma::mat get_optPathParams(){return m_optPathParams;};
  arma::vec get_optPathCosts(){return m_optPathCosts;};

  double get_multipleSolnFlag(){return m_multipleSolnFlag;};

  arma::vec get_solveStatusAll(){return m_solveStatusAll;};
  std::vector< VarSpeedDubins::Path > get_candPaths(){return m_candPaths;};

}; // class

} // namespace

#endif
