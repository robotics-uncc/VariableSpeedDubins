// Copyright (C) 2005, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: hs071_nlp.hpp 1861 2010-12-21 21:34:47Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-09

#ifndef VSD_NLP
#define VSD_NLP

// third-party
#include "IpTNLP.hpp"

// standard
#include <iostream>

// custom
#include<VSDUtils.h>
#include<VSDProblem.h>
#include<VSDPath.h>
#include<armadillo>
#include<MathTools.h>
#include<MathToolsArma.h>


// Source: http://www.coin-or.org/Ipopt/documentation/
//
// IPOPT (Interior Point Optimizer, pronounced ``Eye-Pea-Opt'') is an open 
// source software package for large-scale nonlinear optimization. It can be 
// used to solve general nonlinear programming problems of the form 
//
// 			min f(x)	objective function
//
//		s.t. g^L <= g(x) <= g^U  	non-linear constraints
// 		s.t. x^L <= x 	 <= x^U
//
// Requires f and g to be twice continuously differentiable
//
// In order to solve a problem, IPOPT needs more information than just the 
// problem definition (for example, the derivative information).
// 
// The following information is required by IPOPT: 
//
// 	1. Problem dimensions 
// 		- number of variables
//		- number of constraints 
//	2. Problem bounds
//    	- variable bounds
//    	- constraint bounds
//	3. Initial starting point
//    	- Initial values for the primal x variables
//    	- Initial values for the multipliers (required for warm start option)
//	4. Problem Structure
//    	- number of nonzeros in the Jacobian of the constraints
//    	- number of nonzeros in the Hessian of the Lagrangian function
//    	- sparsity structure of the Jacobian of the constraints
//    	- sparsity structure of the Hessian of the Lagrangian function
//	5. Evaluation of Problem Functions
//		- Information evaluated using a given point ( $ x, \lambda, \sigma_f$ 
//		  coming from IPOPT)
//    	- Objective function, $ f(x)$
//    	- Gradient of the objective $ \nabla f(x)$
//    	- Constraint function values, $ g(x)$
//    	- Jacobian of the constraints, $ \nabla g(x)^T$
//    	- Hessian of the Lagrangian function, 
//		 \sigma_f \nabla^2 f(x) + \sum_{i=1}^m\lambda_i\nabla^2 g_i(x)$
//    	-(this is not required if a quasi-Newton options is chosen to 
//		approximate the second derivatives)


using namespace Ipopt;

// A Variable Speed Dubins (VSD) Nonlinear Problem
class VSDNLP : public TNLP
{
private:
	VarSpeedDubins::ProblemStatement m_vsdprob;
	VarSpeedDubins::Path m_vsdcand;

	int m_numParams, m_numConstraints, m_numElemJacobian;

	arma::vec m_xl, m_xu, m_gl, m_gu;
	arma::vec m_startingPoint;
	bool m_startingPointSuppliedFlag;
	arma::mat m_A;
  int m_numSamplesMax;
	double m_Lmax, m_asubopt_pi;
	arma::vec m_optSoln;
  std::string m_samplingAlgorithm;

	void startingPoint(Number* x);
	void objectiveFunction(Number& obj_value, const Number* x);
	void objectiveFunctionGradient(const Number* x, Number* grad_f);
	void evaluateConstraints(const Number* x, Number* g);
	void evaluateConstraintsJacobian(const Number* x, Number* values, Index n, 
			Index m);
public: 
	void set_problemStatement(VarSpeedDubins::ProblemStatement vsdprob){
		      m_vsdprob = vsdprob;};
	void set_candidate(VarSpeedDubins::Path vsdcand){m_vsdcand = vsdcand;};
	void set_startingPoint(arma::vec startingPoint){
		      m_startingPoint = startingPoint;
		      m_startingPointSuppliedFlag = "true";};
  void set_numSamplesMax(int numSamplesMax){
          m_numSamplesMax = numSamplesMax;};
	void print();
  void initializeCandidate();
	void set_Lmax(double Lmax){m_Lmax = Lmax;};
	arma::vec get_optSoln(){return m_optSoln;};
  void set_samplingAlgorithm(std::string samplingAlgorithm){
    m_samplingAlgorithm = samplingAlgorithm;};

// Do not change function prototypes below this line
//---------------------------------------------------------------------------
public:
  /** default constructor */
  VSDNLP();

  /** default destructor */
  virtual ~VSDNLP();

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

//  /** Method to return:
//   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
//   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
//   */
//  virtual bool eval_h(Index n, const Number* x, bool new_x,
//                      Number obj_factor, Index m, const Number* lambda,
//                      bool new_lambda, Index nele_hess, Index* iRow,
//                      Index* jCol, Number* values);

//  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
				 const IpoptData* ip_data,
				 IpoptCalculatedQuantities* ip_cq);
  //@}

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *  
   */
  //@{
  //  HS071_NLP();
  VSDNLP(const VSDNLP&);
  VSDNLP& operator=(const VSDNLP&);
  //@}
}; //end class


#endif
