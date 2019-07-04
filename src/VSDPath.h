// VSDPath.h
// Provides functions for defining and displaying a Variable Speed Dubins path.
// Last Modefied: 18-Jan-2016, Artur Wolek

#ifndef VSDPATH_H
#define VSDPATH_H

// C++ headers
#include<vector>
#include<string>

// external
#include<armadillo>

// custom headers
#include<VSDProblem.h>

// custom
namespace VarSpeedDubins {

class Path {

private: 

	// Consider the following model of a planar mobile robot with position (x,y) 
	// and heading, h.
	// 
	// 	dx/dt = v cos h
	// 	dy/dt = v sin h
	// 	dh/dt = u		where  |u| <= umax is the turn rate control input
	//					and     0 < vmin <= v <= vmax is the speed control input
	//
	//	Recall that along a circular arc : (turn rate) =  (speed)/(turn radius)
	//	then at |u| = umax, and v = vmin, the turn radius is R = umax/vmin
	// 	and  at |u| = umax, and v = vmax, the turn radius is r = umax/vmax < R
	//					
	// A Variable Speed Dubins Path (VSDP) is the fastest (minimum-time) path 
	// that joins two states (x0, y0, h0) and (xFinal, yFinal, hFinal), 
	// while satisfying the above motion model. 


	// turn radii, as described above
	double m_R;
	double m_r;
	
	// a variable speed Dubins path is a member of the following set:
	// TST: turn-straight-turn
	// TTTT: a sequence of at most four turns
	// UMAX: a fixed turn rate path
	// Dubins: a Dubins-like turn
	std::string m_pathClass; // a string identifying the path class
	
	// each turn "T" as mentioned above is generically composed of a "BCB" 
	// segment where
	// B: is a "bang" turn with |u| = umax, v = vmax, and turn radius R
	// C: is a "cornering" turn with |u| = umax, v = vmin, and turn radius r < R
	// if we consider all possible forms of T, and all possible orientations of 
	// the turns corresponding to the different pathClasses described above, it 
	// turns out there are 72 path types to consider
	std::string m_pathType; // a string identifying the path type
	std::string m_orientation;
	
	// we describe the "B" "C" and "B" segments in a "BCB" turn using the 
	// arc-length parameters alpha (a), beta (b), and gamma (g) corresponding to 
	// each segment respectively the longest sequence is of the form 
	// TTTT = BCB-BCB-BCB-BCB, thus we have four sets of parameters:
	double m_a1, m_a2, m_a3, m_a4;
	double m_b1, m_b2, m_b3, m_b4;
	double m_g1, m_g2, m_g3, m_g4;
	
	// In the TST case the straight segment is described by the length-parameter
	double m_L;

	// all of the above arc-length parameters are grouped into one large vector, 
	// even though only a subset of parameters apply to a particular pathType
	// (e.g. for a "BCB-BCB-BCB" path the parameters a4,b4,g4 are not used)
	arma::vec m_paramsLong;

	// although there may be many parameters describing a variable speed Dubins 
	// path, they do not all vary indepndentalty, and only (at most) four unique 
	// arc-length parameters are required.
	arma::vec m_paramsShort;
	
	// curvature signs
	// for a "rigid" turn it is sufficient to only define the sign of the 
	// initial turn since by definition the following turns are of the opposite 
	// sense
	double m_sgnK1;

	// for a "flexible" turn (of the form TST), the initial and final turn can 
	// vary independently thus a second parameter describes the sign of the 
	// final turn
	double m_sgnK2;

	// the cost of a VSDP is the time it takes to traverse the path
	double m_cost;

	// for the purposes of displaying a VSDP, we must select how many points to 
	// use to represent the path
	int m_numPtsTotal; // total number of points
	double m_nomSpacing;
	
	// further, the number of points for each segment must be determined
	int m_numPts_a1, m_numPts_a2, m_numPts_a3, m_numPts_a4;
	int m_numPts_b1, m_numPts_b2, m_numPts_b3, m_numPts_b4;
	int m_numPts_g1, m_numPts_g2, m_numPts_g3, m_numPts_g4;
	int m_numPts_L;
	
	bool debugFlag = false; // flag used for toggling debug outputs
	
	// private functions
	void expandParamsLong();
	void computePathSpacing(); 
	
	// outputs
	std::string m_solutionStatus;

	// the time and state histories along the path are output
	arma::vec m_endpoint = {0, 0, 0};
	arma::vec m_timeHistory;
	arma::vec m_xHistory;
	arma::vec m_yHistory;
	arma::vec m_hHistory; 
	arma::mat m_pathHistory;
	

public: 
  // constructor and destructor
	Path();
	virtual ~Path(){};

	// Optional?

	
	// required set functions
	void set_pathType(std::string pathType){m_pathType = pathType;};
	void set_pathClass(std::string pathClass){m_pathClass = pathClass;};
	void set_pathOrientation(std::string orientation){
			m_orientation = orientation;};
	void set_cornerTurnRadius(double r){m_r = r;};
	void set_bangTurnRadius(double R){m_R = R;};
	void set_turnRadii(double r, double R){m_r = r; m_R = R;};
	void set_nominalSpacing(double nomSpacing){m_nomSpacing = nomSpacing;};
	// only of the following needs to be set:
	void set_paramsLong(arma::vec paramsLong);
	void set_paramsShort(arma::vec paramsShort);

	// these are not always used
	void set_numPtsTotal(int numPtsTotal){m_numPtsTotal = numPtsTotal;};
	void set_solutionStatus(std::string solutionStatus){
			m_solutionStatus = solutionStatus;};

	// main functions
	void computePathHistory(); 
	void print(); 
	void printCandidate(); 
	void computeEndpoint();
	void writePathPlotCommands(std::string fileName, int figureNumber, 
			std::string xVarName, std::string yVarName, std::string style);
	void writeEndpointPlotCommands(std::string fileName, int figureNumber, 
			std::string xVarName, std::string yVarName, std::string style);
  void writeSwitchingPtPlotCommands(std::string fileName, int figureNumber, 
			std::string xVarName, std::string yVarName, std::string style);

	// get functions
	double get_r(){return m_r;};
	double get_R(){return m_R;};

	std::string get_pathType(){return m_pathType;};
	std::string get_pathClass(){return m_pathClass;};
	std::string get_orientation(){return m_orientation;};
	std::string get_solutionStatus(){return m_solutionStatus;};

	arma::vec get_paramsLong(){return m_paramsLong;};
	arma::vec get_paramsShort(){return m_paramsShort;}; 

	double get_cost();
	double get_xFinal(){return m_endpoint(0);};
	double get_yFinal(){return m_endpoint(1);};
	double get_hFinal(){return m_endpoint(2);};

	arma::vec get_timeHistory(){return m_timeHistory;};
	arma::vec get_xHistory(){return m_xHistory;};
	arma::vec get_yHistory(){return m_yHistory;};
	arma::vec get_hHistory(){return m_hHistory;}; 	

}; // class

} // namespace
#endif
