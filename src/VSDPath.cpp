// VSDPath.cpp
// Last Modefied: 18-Jan-2016, Artur Wolek

#include<VSDPath.h>
#include<VSDUtils.h>

// standard headers
#include<limits> 

// constructor
VarSpeedDubins::Path::Path(){
	// intialize variables
	m_nomSpacing = 0.05;
  m_paramsShort = arma::zeros<arma::vec>(4);
  m_cost = std::numeric_limits<double>::max();
}

// set_paramsLong
void VarSpeedDubins::Path::set_paramsLong(arma::vec paramsLong){
	m_paramsLong = paramsLong;
	VarSpeedDubins::Path::expandParamsLong();
}

// set_paramsShort
void VarSpeedDubins::Path::set_paramsShort(arma::vec paramsShort){
	// store input parameters
	m_paramsShort = paramsShort;
	// expand paramsShort, assumes m_pathClass, m_pathType are already set
	m_paramsLong = VarSpeedDubins::expandParamsShort(paramsShort, m_pathClass,      
                                                   m_pathType);
	VarSpeedDubins::Path::expandParamsLong();
}

// computePathHistory
void VarSpeedDubins::Path::computePathHistory(){	
	VarSpeedDubins::pathHistory(m_pathHistory, m_paramsShort, m_pathClass, 
		                          m_pathType, m_orientation, m_R, m_r, 
                              m_nomSpacing);
	return;
}

// computeEndpoint
void VarSpeedDubins::Path::computeEndpoint(){
	m_endpoint = VarSpeedDubins::pathEndpoint(m_paramsShort, m_pathClass, 
		                                        m_pathType, m_orientation, 
                                            m_R, m_r); 
	return;
}

// print
void VarSpeedDubins::Path::print(){
	std::cout << "---------------------- Path -------------------" << std::endl;
	std::cout << "Path Class \t : " << m_pathClass << std::endl;
	std::cout << "Path Type \t : " << m_pathType << std::endl;
	std::cout << "Path Orientation : " << m_orientation << std::endl;
	std::cout << "x final \t : " << m_endpoint(0) << std::endl;
	std::cout << "y final \t : " << m_endpoint(1) << std::endl;
	std::cout << "h final \t : " << m_endpoint(2) << std::endl;
	std::cout << "(a1,b1,g1) \t : (" << m_a1 << "," << m_b1 << "," << m_g1
		        << ")" << std::endl;
	std::cout << "L \t \t: " << m_L << std::endl;
	std::cout << "(a2,b2,g2) \t : (" << m_a2 << "," << m_b2 << "," << m_g2
		        << ")" << std::endl;
	std::cout << "(a3,b3,g3) \t : (" << m_a3 << "," << m_b3 << "," << m_g3
		        << ")" << std::endl;
	std::cout << "(a4,b4,g4) \t : (" << m_a4 << "," << m_b4 << "," << m_g4
		        << ")" << std::endl;
	std::cout << "cost \t : " << m_cost << std::endl;
	std::cout << "-----------------------------------------------" << std::endl;
	return;
}

// print
void VarSpeedDubins::Path::printCandidate(){
	std::cout << "---------------------- Path -------------------" << std::endl;
	std::cout << "Path Class \t : " << m_pathClass << std::endl;
	std::cout << "Path Type \t : " << m_pathType << std::endl;
	std::cout << "Path Orientation : " << m_orientation << std::endl;
	std::cout << "-----------------------------------------------" << std::endl;
	return;
}



// writePathPlotCommands
void VarSpeedDubins::Path::writePathPlotCommands(std::string fileName, 
		                                             int figureNumber, 
                                                 std::string xVarName, 
                                                 std::string yVarName, 
		                                             std::string style){	
	arma::colvec xData = m_pathHistory.col(0);	
	arma::colvec yData = m_pathHistory.col(1);
	MathTools::writeVectorData(fileName, xData, xVarName);
	MathTools::writeVectorData(fileName, yData, yVarName);
	MathTools::plot1DCurve(fileName, figureNumber, xVarName, yVarName, style);
  MathTools::axisProperty(fileName, "equal");		
	return;			
}

// writeEndpointPlotCommands
void VarSpeedDubins::Path::writeEndpointPlotCommands(std::string fileName, 
		                                                 int figureNumber, 
                                                     std::string xVarName, 
                                                     std::string yVarName, 
		                                                 std::string style){
		MathTools::plotPoint(fileName, figureNumber, m_endpoint(0), m_endpoint(1), 
                         xVarName, yVarName, style);		
	return;			
}

// writeSwitchingPtPlotCommands
void VarSpeedDubins::Path::writeSwitchingPtPlotCommands(std::string fileName, 
		                                                    int figureNumber,         
                                                        std::string xVarName, 
                                                        std::string yVarName, 
		                                                    std::string style){
    arma::mat switchingPts = extremalSwitches(m_paramsShort, m_pathClass, 
		                                          m_pathType, m_orientation, 
                                              m_R, m_r);
    arma::vec xData = switchingPts.row(0).t();
    arma::vec yData = switchingPts.row(1).t();
  	MathTools::writeVectorData(fileName, xData, xVarName);
	  MathTools::writeVectorData(fileName, yData, yVarName);  
	  MathTools::plot1DCurve(fileName, figureNumber, xVarName, yVarName, style);		
	return;			
}

// get_cost()
double VarSpeedDubins::Path::get_cost(){\
  if ( m_solutionStatus.compare("feasible") || 
       m_solutionStatus.compare("optimal") ){
	  m_cost = VarSpeedDubins::costFunction(m_paramsShort, m_pathClass, m_pathType, 
				                                m_orientation, m_R);
  }
	return m_cost;
}

// expandParamsLong
void VarSpeedDubins::Path::expandParamsLong(){
	m_a1 = m_paramsLong(0);
	m_b1 = m_paramsLong(1);
	m_g1 = m_paramsLong(2);
	m_L = m_paramsLong(3);
	m_a2 = m_paramsLong(4);
	m_b2 = m_paramsLong(5);
	m_g2 = m_paramsLong(6);
	m_a3 = m_paramsLong(7);
	m_b3 = m_paramsLong(8);
	m_g3 = m_paramsLong(9);
	m_a4 = m_paramsLong(10);	
	m_b4 = m_paramsLong(11);
	m_g4 = m_paramsLong(12);
}
