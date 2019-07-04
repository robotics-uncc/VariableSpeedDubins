// ProblemStatement.cpp
// Last Modified:  09-Dec-2015, Artur Wolek

// custom headers
#include<VSDProblem.h>

// constructor
VarSpeedDubins::ProblemStatement::ProblemStatement(){
	m_xInitial = 0;
	m_yInitial = 0;
	m_hInitial = 0;
	m_R = 1.0; // default values
	m_r = 0.3;
}

// destructor
VarSpeedDubins::ProblemStatement::~ProblemStatement(){
}

// main functions

void VarSpeedDubins::ProblemStatement::print(){
	std::cout << "--------------Problem Statement----------------" << std::endl;
	std::cout << "x initial : \t" << m_xInitial << std::endl;
	std::cout << "y initial : \t" << m_yInitial << std::endl;
	std::cout << "h initial : \t" << m_hInitial << std::endl;
	std::cout << "x final : \t" << m_xFinal << std::endl;
	std::cout << "y final : \t" << m_yFinal << std::endl;
	std::cout << "h final : \t" << m_hFinal << std::endl;
	std::cout << "R : \t" << m_R << std::endl;
	std::cout << "r : \t" << m_r << std::endl;
	std::cout << "-----------------------------------------------" << std::endl;
}

// set functions

void VarSpeedDubins::ProblemStatement::set_stateInitial(arma::colvec 
		stateInitial){
	set_xInitial(stateInitial(0));
	set_yInitial(stateInitial(1));
	set_hInitial(stateInitial(2));
}

void VarSpeedDubins::ProblemStatement::set_stateInitial(double xInitial, 
		double yInitial, double hInitial){
	set_xInitial(xInitial);
	set_yInitial(yInitial);
	set_hInitial(hInitial);
}

void VarSpeedDubins::ProblemStatement::set_stateFinal(arma::colvec stateFinal){
	set_xFinal(stateFinal(0));
	set_yFinal(stateFinal(1));
	set_hFinal(stateFinal(2));
}

void VarSpeedDubins::ProblemStatement::set_stateFinal(double xFinal, 
		double yFinal, double hFinal){
	set_xFinal(xFinal);
	set_yFinal(yFinal);
	set_hFinal(hFinal);
}

// get functions
arma::colvec VarSpeedDubins::ProblemStatement::get_stateFinal(){
	arma::colvec stateFinal = {m_xFinal, m_yFinal, m_hFinal};
	return stateFinal;
}














