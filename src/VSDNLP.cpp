// VSDNLP.cpp
// Last Modified: Artur Wolek, 31-Jan-2016

// Derived from the following example: hs071_nlp.cpp
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "VSDNLP.h"

#include <cassert>
#include <iostream>

using namespace Ipopt;

// constructor
VSDNLP::VSDNLP(){
  m_startingPointSuppliedFlag = false;
  m_numSamplesMax = 5;
  m_samplingAlgorithm = "barycentric";
}

//destructor
VSDNLP::~VSDNLP()
{}

void VSDNLP::initializeCandidate(){
  //std::cout << "NLPdimensions" << std::endl;

	VarSpeedDubins::NLPdimensions(m_vsdcand.get_pathClass(), 
		                            m_vsdcand.get_pathType(), m_numParams, 
                                m_numConstraints, m_numElemJacobian);
  //std::cout << "m_asubopt_pi" << std::endl;
	m_asubopt_pi = VarSpeedDubins::alphaSubopt(M_PI, m_vsdcand.get_R(), 
		                                         m_vsdcand.get_r());
  //std::cout << "linearConstraints" << std::endl;
	VarSpeedDubins::linearConstraints(m_vsdcand.get_pathClass(), 
		                                m_vsdcand.get_pathType(), m_xl, m_xu, m_gl, 
                                    m_gu, m_A, m_vsdcand.get_R(), 
                                    m_vsdcand.get_r(), m_asubopt_pi, m_Lmax);
  // if starting point not supplied, must compute a guess
  if (m_startingPointSuppliedFlag == false){
    //std::cout << "computing a guess..." << std::endl;
    m_startingPoint = VarSpeedDubins::guessStartPoint(m_samplingAlgorithm,
                                          m_vsdcand.get_pathClass(),  
                                          m_vsdcand.get_pathType(), 
                                          m_vsdcand.get_orientation(), 
                                          m_xl, m_xu, m_gl, m_gu, m_A, 
                                          m_vsdcand.get_R(), m_vsdcand.get_r(), 
                                          m_numParams, m_numConstraints,
                                          m_vsdprob.get_stateFinal(), 
                                          m_numSamplesMax);
    //std::cout << "computing a guess complete..." << std::endl;  
  }
}

void VSDNLP::print(){
	std::cout << "-----------------------------------------------" << std::endl;
	std::cout << "                  VSD NLP                      " << std::endl;
	std::cout << "-----------------------------------------------" << std::endl;
	m_vsdprob.print();
	m_vsdcand.printCandidate();
}

// returns the size of the problem
bool VSDNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  #if DEBUG_LEVEL > 1
    std::cout << "VSDNLP::get_nlp_info start" << std::endl;
  #endif

	m_xl.set_size(m_numParams);	
	m_xu.set_size(m_numParams);	
	m_gl.set_size(m_numConstraints);	
	m_gu.set_size(m_numConstraints);
	m_A.set_size(m_numConstraints-3,m_numParams);



  	// The dimension of the parameter vector (typically n = 4)
  	n = m_numParams;

  	// includes both equality and inequality constraints (combined into g(x))
	// BCB-S-BCB : 4 linear constraints, 3 B.C.'s
  	m = m_numConstraints;

  	// in this example the jacobian of the constraints is dense
  	nnz_jac_g = m_numElemJacobian;

	// IPOPT has an option to approximate the Hessian of the Lagrangian by a 
	// limited-memory quasi-Newton method (L-BFGS). You can use this feature by 
	// setting the hessian_approximation option to the value limited-memory. 
	// In this case, it is not necessary to implement the Hessian computation 
	// method eval_h in TNLP.
	// http://www.coin-or.org/Ipopt/documentation/node52.html

  	// Since we are using an approximation, we comment the hessian nnz:
  	// nnz_h_lag = 10;


  	// use the C style indexing (0-based)
  	index_style = TNLP::C_STYLE;

  #if DEBUG_LEVEL > 1
    std::cout << "VSDNLP::get_nlp_info complete" << std::endl;
    std::cout << "-> m_numParams : " << m_numParams << std::endl;
    std::cout << "-> m_numConstraints : " << m_numConstraints << std::endl;
    std::cout << "-> m_numElemJacobian : " << m_numElemJacobian << std::endl;
  #endif

  	return true;
}

// returns the variable bounds
bool VSDNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  #if DEBUG_LEVEL > 1
  std::cout << "VSDNLP::get_bounds_info start" << std::endl;
  #endif 
	// get constraint bounds
	VarSpeedDubins::linearConstraints(m_vsdcand.get_pathClass(), 
		                                m_vsdcand.get_pathType(), m_xl, m_xu, m_gl, 
                                    m_gu, m_A, m_vsdcand.get_R(), 
                                    m_vsdcand.get_r(), m_asubopt_pi, m_Lmax);
  //std::cout << " linearConstraints passed" << std::endl;
	// set parameter bounds
	for (int i = 0; i < n; i++){
		x_u[i] = m_xu(i);
		x_l[i] = m_xl(i);
	}
  //std::cout << " set parameter bounds passed" << std::endl;
  //std::cout << " m : " << m << std::endl;
	// set linear constraint bounds
	for (int i = 0; i < m - 3; i++){
		g_u[i] = m_gu(i);
		g_l[i] = m_gl(i);
	}
  //std::cout << " set linear constraint bounds passed" << std::endl;
	// add nonlinear constraint bounds (boundary conditions)
	g_l[m-3] = m_vsdprob.get_xFinal();
	g_l[m-2] = m_vsdprob.get_yFinal();
	g_l[m-1] = m_vsdprob.get_hFinal();
	g_u[m-3] = g_l[m-3];
	g_u[m-2] = g_l[m-2];
	g_u[m-1] = g_l[m-1];

  #if DEBUG_LEVEL > 1
  std::cout << "VSDNLP::get_bounds_info complete" << std::endl;
  m_xl.print("-> m_xl");
  m_xu.print("-> m_xu");
  m_gl.print("-> m_gl");
  m_gu.print("-> m_gus");
  #endif 

	return true;
}

// returns the initial point for the problem
bool VSDNLP::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  #if DEBUG_LEVEL > 1
    std::cout << "VSDNLP::get_starting_point start" << std::endl;
  #endif 
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

	// set starting point
  for (int i = 0; i < n; i++){
    x[i] = m_startingPoint(i);
  }		


	
  #if DEBUG_LEVEL > 1
    std::cout << "VSDNLP::get_starting_point complete" << std::endl;
    m_startingPoint.print("-> m_startingPoint");
  #endif 
  	return true;
}

// returns the value of the objective function
bool VSDNLP::eval_f(Index n, const Number* x, bool new_x, 
		Number& obj_value)
{
  #if DEBUG_LEVEL > 2
    std::cout << "VSDNLP::eval_f start" << std::endl;
  #endif 
  	//assert(n == 4);
	arma::vec paramsShort(n);
	for (int i = 0; i < n; i++){
		paramsShort(i) = x[i];
	}
	
	obj_value = VarSpeedDubins::costFunction(paramsShort, 
		m_vsdcand.get_pathClass(), m_vsdcand.get_pathType(), m_vsdcand.get_orientation(), 
    m_vsdprob.get_R());

  #if DEBUG_LEVEL > 2
    std::cout << "VSDNLP::eval_f complete" << std::endl;
    paramsShort.print("-> paramsShort");
    std::cout << "-> obj_value :" << obj_value << std::endl;
  #endif 
  	return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool VSDNLP::eval_grad_f(Index n, const Number* x, bool new_x, 
		Number* grad_f)
{
	// If the objective functions if f(x), where x is n-dimensional
	// The grad_f is given as a n-dimensional row vector
	// grad_f[0] = df/dx0
	// grad_f[1] = df/dx1
	// ...
	// grad_f[n-1] = df/dx(n-1)
  #if DEBUG_LEVEL > 2
    std::cout << "VSDNLP::eval_grad_f start" << std::endl;
  #endif 
  	// assert(n == 4);
	arma::vec paramsShort(n);
	for (int i = 0; i < n; i++){
		paramsShort(i) = x[i];
	}

	arma::vec costFuncGradient = VarSpeedDubins::costFunctionGradient(
		m_vsdprob.get_R(), m_vsdcand.get_pathClass(), m_vsdcand.get_pathType(), 
    m_vsdcand.get_orientation());

	for (int i = 0; i < n; i++){
		grad_f[i] = costFuncGradient(i); 
	}
  #if DEBUG_LEVEL > 2
    std::cout << "VSDNLP::eval_grad_f complete" << std::endl;
    paramsShort.print("-> paramsShort");
    costFuncGradient.print("-> costFuncGradient");
  #endif 
  	return true;
}

// return the value of the constraints: g(x)
bool VSDNLP::eval_g(Index n, const Number* x, bool new_x, 
		Index m, Number* g)
{
    #if DEBUG_LEVEL > 2
      std::cout << "VSDNLP::eval_g start" << std::endl;
    #endif 
	arma::vec paramsShort(n);
	for (int i = 0; i < n; i++){
		paramsShort(i) = x[i];
	}

	arma::vec endpoint = VarSpeedDubins::pathEndpoint(paramsShort, 
		m_vsdcand.get_pathClass(), m_vsdcand.get_pathType(), 
		m_vsdcand.get_orientation(), m_vsdprob.get_R(), m_vsdprob.get_r());

	arma::vec linConstrEval = m_A*paramsShort;
	for (int i = 0; i < m - 3; i++){
		g[i] = linConstrEval(i);
	}
	g[m-3] = endpoint(0); 
	g[m-2] = endpoint(1); 
	g[m-1] = endpoint(2); 
    #if DEBUG_LEVEL > 2
      std::cout << "VSDNLP::eval_g complete" << std::endl;
    #endif 
  	return true;
}

// return the structure or values of the jacobian
bool VSDNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
// The jacobian is a matrix where the i-th row and j-th col corresponds 
// to te derivative dgi/dxj
// so the matrix itself will have dimensions m x n

//  if (values == NULL) {
//    // return the structure of the jacobian

//    // this particular jacobian is dense
//    iRow[0] = 0;
//    jCol[0] = 0;
//    iRow[1] = 0;
//    jCol[1] = 1;
//    iRow[2] = 0;
//    jCol[2] = 2;
//    iRow[3] = 0;
//    jCol[3] = 3;
//    iRow[4] = 1;
//    jCol[4] = 0;
//    iRow[5] = 1;
//    jCol[5] = 1;
//    iRow[6] = 1;
//    jCol[6] = 2;
//    iRow[7] = 1;
//    jCol[7] = 3;
//  }
//  else {
//    // return the values of the jacobian of the constraints

//    values[0] = x[1]*x[2]*x[3]; // 0,0
//    values[1] = x[0]*x[2]*x[3]; // 0,1
//    values[2] = x[0]*x[1]*x[3]; // 0,2
//    values[3] = x[0]*x[1]*x[2]; // 0,3

//    values[4] = 2*x[0]; // 1,0
//    values[5] = 2*x[1]; // 1,1
//    values[6] = 2*x[2]; // 1,2
//    values[7] = 2*x[3]; // 1,3
//  }
    #if DEBUG_LEVEL > 2
      std::cout << "VSDNLP::eval_jac_g start" << std::endl;
    #endif 
	if (values == NULL) {
		// return the structure of the jacobian as full, dense
		int k = 0;
  		for (Index i = 0; i < m; i++) {
			for (Index j = 0; j < n; j++) {
    			iRow[k] = i; 
				jCol[k] = j; 
				// j iterates first so the values are given col-wise
				// i.e. for a sparse matrix
				// values[0] = jac_g(0,0) = dg0/dx0
				// values[1] = jac_g(0,1) = dg0/dx1, etc...
				k++; 
			}
  		}
	}
  	else {
		// return the values of the jacobian of the constraints
		arma::vec paramsShort(n);
		for (int i = 0; i < n; i++){
			paramsShort(i) = x[i];
		}
		int k = 0;
	  	for (Index i = 0; i < m-3; i++) { // row
			for (Index j = 0; j < n; j++) { // first 4 cols
				values[k] = m_A(i,j);
				k++; 
			}
	  	}

		// 4 x 3 jacobian of the nonlinear constraints
		// i.e.
		// grad_out(0,0) = dgx_da1 = dg0/dx0
		// grad_out(1,0) = dgx_dg1 = dg0/dx1
		// grad_out(2,0) = dgx_dL  = dg0/dx2
		// grad_out(3,0) = dgx_dg2 = dg0/dx3

		arma::mat nlconJacob = VarSpeedDubins::nonlinearEqualityConstraintJacobian(
			paramsShort, m_vsdcand.get_pathClass(), m_vsdcand.get_pathType(), 
			m_vsdcand.get_orientation(), m_vsdprob.get_r(), m_vsdprob.get_R());

		// take transpose 
		for (Index i = 0; i < 3; i++) { // for each row, (each constraint, gi)
			for (Index j = 0; j < n; j++) { // for each indep. variable
				values[k] = nlconJacob(j,i);
				k++; 
			}
    }

		arma::mat valuesMat(m,n);
		int l = 0;
	  	for (Index i = 0; i < m; i++) {
			  for (Index j = 0; j < n; j++) {
				  valuesMat(i,j) = values[l];
				  l++; 
			  }
	  	}
  	}
    #if DEBUG_LEVEL > 2
      std::cout << "VSDNLP::eval_jac_g complete" << std::endl;
    #endif 
  	return true;
}



void VSDNLP::finalize_solution(SolverReturn status,
                Index n, const Number* x, const Number* z_L, const Number* z_U,
                Index m, const Number* g, const Number* lambda,
                Number obj_value,
                const IpoptData* ip_data,
                IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  // For this example, we write the solution to the console
  // std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
	m_optSoln.set_size(4);
  for (Index i=0; i<n; i++) {
	  m_optSoln(i) = x[i];
  }


//  std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
//  for (Index i=0; i<n; i++) {
//    std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
//  }
//  for (Index i=0; i<n; i++) {
//    std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
//  }

  #if DEBUG_LEVEL > 1
  std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
  for (Index i=0; i < n; i++) {
    std::cout << "x[" << i << "] = " << x[i] << std::endl;
  }
  std::cout << std::endl << std::endl << "Objective value" << std::endl;
  std::cout << "f(x*) = " << obj_value << std::endl;
  std::cout << std::endl << "Final value of the constraints:" << std::endl;
  for (Index i=0; i<m ;i++) {
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  }
  m_vsdcand.print();
	arma::vec endpoint = VarSpeedDubins::pathEndpoint(m_optSoln, 
		m_vsdcand.get_pathClass(), m_vsdcand.get_pathType(), 
		m_vsdcand.get_orientation(), m_vsdprob.get_R(), m_vsdprob.get_r());
  #endif

}


