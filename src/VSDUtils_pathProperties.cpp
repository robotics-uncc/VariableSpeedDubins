// VCSDUtils_pathProperties.cpp
// Last Modefied: 25-Jan-2016, Artur Wolek

#include<VSDUtils.h>


// returns curvature parameters, given pathClass and orientation
void VarSpeedDubins::NLPdimensions(std::string pathClass, std::string pathType,
		int &numParams, int &numConstraints, int &numElementsJacob){
	if ( pathClass.compare("TST") == 0 ){
		numParams = 4;
		// 3 boundary conditions + parameter constraints
		if (pathType.compare("BCB-S-BCB") == 0){ // 1-4
			numConstraints = 5; 
		}
		else if (pathType.compare("BCB-S-BC") == 0){ // 5-8
			numConstraints = 4;
		}
		else if (pathType.compare("BCB-S-B") == 0){ // 9-12
			numConstraints = 5;
		}
		else if (pathType.compare("CB-S-BCB") == 0){ // 13-16
			numConstraints = 4;
		}
		else if (pathType.compare("CB-S-BC") == 0){ // 17-20
			numConstraints = 5;
		}
		else if (pathType.compare("CB-S-B") == 0){ // 21-24
			numConstraints = 5;
		}
		else if (pathType.compare("B-S-BCB") == 0){ // 25-28
			numConstraints = 5;
		}
		else if (pathType.compare("B-S-BC") == 0){ // 29-32
			numConstraints = 5;
		}
	}
	else if ( pathClass.compare("TT") == 0 ){
		numParams = 4;
		if (pathType.compare("BCB-BCB") == 0){ // 33-36
			numConstraints = 6; 
		}
		else if (pathType.compare("BCB-BC") == 0){ // 37-40
			numConstraints = 6;
		}
		else if (pathType.compare("BCB-B") == 0){ // 41-44
			numConstraints = 6;
		}
		else if (pathType.compare("CB-BCB") == 0){ // 45-48
			numConstraints = 6;
		}
		else if (pathType.compare("B-BCB") == 0){ // 49-52
			numConstraints = 6;
		}
	}
	else if ( (pathClass.compare("TTTT") == 0)){
			numParams = 4;
		if (pathType.compare("CB-TT-BC") == 0){ // 53-54
			numConstraints = 6; /// 
		}
		else if (pathType.compare("CB-TT-B") == 0){ // 55-56
			numConstraints = 6; /// 
		}
		else if (pathType.compare("B-TT-BC") == 0){ // 63-64
			numConstraints = 6; /// 
		}
		else if (pathType.compare("B-TT-B") == 0){ // 65-66
			numConstraints = 6; 
		}
	}
	else if ( (pathClass.compare("TTT") == 0)){
			numParams = 4;
		if (pathType.compare("CB-T-BCB") == 0){ // 57-58
			numConstraints = 6; ///
		}
		else if (pathType.compare("CB-T-BC") == 0){ // 59-60
			numConstraints = 6; ///
		}
		else if (pathType.compare("CB-T-B") == 0){ // 61-62
			numConstraints = 6;
		}
		else if (pathType.compare("B-T-BCB") == 0){ // 67-68
			numConstraints = 6;
		}
		else if (pathType.compare("B-T-BC") == 0){ // 69-70
			numConstraints = 6;
		}
		else if (pathType.compare("B-T-B") == 0){ // 71-72
			numConstraints = 6;
		}
		else if (pathType.compare("BCB-T-BC") == 0){ // 73-74
			numConstraints = 6;
		}
		else if (pathType.compare("BCB-T-B") == 0){ // 75-76
			numConstraints = 6;
		}
	}
	numElementsJacob = numParams*numConstraints;
}

// returns curvature parameters, given pathClass and orientation
void VarSpeedDubins::determineCurvatureParams(std::string 
	pathClass, std::string orientation, double &k1, double &k2){
	if ( pathClass.compare("TST") == 0 ){
			if (orientation.compare("LSL") == 0){
					k1 = 1.0; 
					k2 = 1.0;
			}
			else if (orientation.compare("LSR") == 0){
					k1 = 1.0; 
					k2 = -1.0;
			}
			else if (orientation.compare("RSL") == 0){
					k1 = -1.0; 
					k2 = 1.0;
			}
			else if (orientation.compare("RSR") == 0){
					k1 = -1.0; 
					k2 = -1.0;
			}
	}
	else if ( pathClass.compare("TT") == 0 ){
			if (orientation.compare("LL") == 0){	
					k1 = 1.0; 
					k2 = 1.0;		
			}
			else if (orientation.compare("LR") == 0){	
					k1 = 1.0; 
					k2 = -1.0;		
			}
			else if (orientation.compare("RL") == 0){	
					k1 = -1.0; 
					k2 = 1.0;		
			}
			else if (orientation.compare("RR") == 0){	
					k1 = -1.0; 
					k2 = -1.0;		
			}
	}
	else if ( pathClass.compare("TTT") == 0 ){
			if (orientation.compare("LRL") == 0){	
					k1 = 1.0; 
					k2 = -1.0;	// undefined		
			}
			else if (orientation.compare("RLR") == 0){	
					k1 = -1.0; 
					k2 = 1.0;	// undefined		
			}
	}
	else if ( pathClass.compare("TTT") == 0 ){
			if (orientation.compare("LRL") == 0){	
					k1 = 1.0; 
					k2 = -1.0;	// undefined		
			}
			else if (orientation.compare("RLR") == 0){	
					k1 = -1.0; 
					k2 = 1.0;	// undefined		
			}
	}
	else if ( pathClass.compare("TTTT") == 0 ){
			if (orientation.compare("LRLR") == 0){	
					k1 = 1.0; 
					k2 = -1.0;	// undefined		
			}
			else if (orientation.compare("RLRL") == 0){	
					k1 = -1.0; 
					k2 = 1.0;	// undefined		
			}
	}
}

// returns a string of the orientation given the initial and final 
// curvature (k1, k2), and the pathClass
std::string VarSpeedDubins::determineOrientation(double k1, double k2, std::string pathClass){
	std::string orientation;
	if ( pathClass.compare("TST") == 0 ){
		if (k1 == 1 && k2 == 1){
		    orientation = "LSL";
		}
		else if (k1 == 1 && k2 == -1){
		    orientation = "LSR";
		}
		else if (k1 == -1 && k2 == 1){
		    orientation = "RSL";
		}
		else if (k1 == -1 && k2 == -1){
		    orientation = "RSR";
		}
	}
	else if ( pathClass.compare("TT") == 0 ){
		if (k1 == 1 && k2 == 1){
		    orientation = "LL";
		}
		else if (k1 == 1 && k2 == -1){
		    orientation = "LR";
		}
		else if (k1 == -1 && k2 == 1){
		    orientation = "RL";
		}
		else if (k1 == -1 && k2 == -1){
		    orientation = "RR";
		}
	}       
	else if ( pathClass.compare("TTT") == 0 ){
		if (k1 == 1){
		    orientation = "LRL";
		}
		else if (k1 == -1){
		    orientation = "RLR";
		}
	}
	else if ( pathClass.compare("TTTT") == 0 ){
		if (k1 == 1){
		    orientation = "LRLR";
		    }
		else if (k1 == -1){
		    orientation = "RLRL";
		    }
		}	
	return orientation;	
}

std::vector< std::vector<std::string> > VarSpeedDubins::candidateList(){
	std::vector< std::string > candidate;
	std::vector< std::vector<std::string> > candidateList;
	// Total number of types:
	// TST 		: 32
	// TT 		: 20
	// TTT		: 10
	// TTTT		: 2
	// Dubins	: 8
	// ------------------
	// TOTAL 	: 72 path types
	
	// TST, 32 types
	candidateList.resize(84);
	// BCB-S-BCB: 0-3
	candidate = {"TST","BCB-S-BCB","LSL"};
	candidateList.at(0) = candidate;
	candidate = {"TST","BCB-S-BCB","LSR"};
	candidateList.at(1) = candidate;
	candidate = {"TST","BCB-S-BCB","RSL"};
	candidateList.at(2) = candidate;
	candidate = {"TST","BCB-S-BCB","RSR"}; 
	candidateList.at(3) = candidate;	
	// BCB-S-BC: 4-7
	candidate = {"TST","BCB-S-BC","LSL"};
	candidateList.at(4) = candidate;
	candidate = {"TST","BCB-S-BC","LSR"};
	candidateList.at(5) = candidate;
	candidate = {"TST","BCB-S-BC","RSL"};
	candidateList.at(6) = candidate;
	candidate = {"TST","BCB-S-BC","RSR"};
	candidateList.at(7) = candidate;	
	// BCB-S-B: 8-11
	candidate = {"TST","BCB-S-B","LSL"};
	candidateList.at(8) = candidate;
	candidate = {"TST","BCB-S-B","LSR"};
	candidateList.at(9) = candidate;
	candidate = {"TST","BCB-S-B","RSL"};
	candidateList.at(10) = candidate;
	candidate = {"TST","BCB-S-B","RSR"};
	candidateList.at(11) = candidate;	
	// CB-S-BCB: 12-15
	candidate = {"TST","CB-S-BCB","LSL"};
	candidateList.at(12) = candidate;
	candidate = {"TST","CB-S-BCB","LSR"};
	candidateList.at(13) = candidate;
	candidate = {"TST","CB-S-BCB","RSL"};
	candidateList.at(14) = candidate;
	candidate = {"TST","CB-S-BCB","RSR"};
	candidateList.at(15) = candidate;
	// CB-S-BC: 16-19	
	candidate = {"TST","CB-S-BC","LSL"};
	candidateList.at(16) = candidate;
	candidate = {"TST","CB-S-BC","LSR"};
	candidateList.at(17) = candidate;
	candidate = {"TST","CB-S-BC","RSL"};
	candidateList.at(18) = candidate;
	candidate = {"TST","CB-S-BC","RSR"};
	candidateList.at(19) = candidate;	
	// CB-S-B: 20-23
	candidate = {"TST","CB-S-B","LSL"};
	candidateList.at(20) = candidate;
	candidate = {"TST","CB-S-B","LSR"};
	candidateList.at(21) = candidate;
	candidate = {"TST","CB-S-B","RSL"};
	candidateList.at(22) = candidate;
	candidate = {"TST","CB-S-B","RSR"};
	candidateList.at(23) = candidate;	
	// B-S-BCB: 24-27
	candidate = {"TST","B-S-BCB","LSL"};
	candidateList.at(24) = candidate;
	candidate = {"TST","B-S-BCB","LSR"};
	candidateList.at(25) = candidate;
	candidate = {"TST","B-S-BCB","RSL"};
	candidateList.at(26) = candidate;
	candidate = {"TST","B-S-BCB","RSR"};
	candidateList.at(27) = candidate;	
	// B-S-BC: 28-31
	candidate = {"TST","B-S-BC","LSL"};
	candidateList.at(28) = candidate;
	candidate = {"TST","B-S-BC","LSR"};
	candidateList.at(29) = candidate;
	candidate = {"TST","B-S-BC","RSL"};
	candidateList.at(30) = candidate;
	candidate = {"TST","B-S-BC","RSR"};
	candidateList.at(31) = candidate;	
	// TT 
	// BCB-BCB: 32-35
	candidate = {"TT","BCB-BCB","LL"};
	candidateList.at(32) = candidate;
	candidate = {"TT","BCB-BCB","LR"};
	candidateList.at(33) = candidate;
	candidate = {"TT","BCB-BCB","RL"};
	candidateList.at(34) = candidate;	
	candidate = {"TT","BCB-BCB","RR"};	
	candidateList.at(35) = candidate;	
	// BCB-BC: 36-39
	candidate = {"TT","BCB-BC","LL"};
	candidateList.at(36) = candidate;
	candidate = {"TT","BCB-BC","LR"};
	candidateList.at(37) = candidate;	
	candidate = {"TT","BCB-BC","RL"};
	candidateList.at(38) = candidate;
	candidate = {"TT","BCB-BC","RR"};	
	candidateList.at(39) = candidate;
	// BCB-B: 40-43
	candidate = {"TT","BCB-B","LL"};
	candidateList.at(40) = candidate;
	candidate = {"TT","BCB-B","LR"};
	candidateList.at(41) = candidate;
	candidate = {"TT","BCB-B","RL"};
	candidateList.at(42) = candidate;
	candidate = {"TT","BCB-B","RR"};
	candidateList.at(43) = candidate;	
	// CB-BCB: 44-47
	candidate = {"TT","CB-BCB","LL"};
	candidateList.at(44) = candidate;
	candidate = {"TT","CB-BCB","LR"};	
	candidateList.at(45) = candidate;	
	candidate = {"TT","CB-BCB","RL"};
	candidateList.at(46) = candidate;
	candidate = {"TT","CB-BCB","RR"};	
	candidateList.at(47) = candidate;	
	// B-BCB: 48-51
	candidate = {"TT","B-BCB","LL"};
	candidateList.at(48) = candidate;
	candidate = {"TT","B-BCB","LR"};	
	candidateList.at(49) = candidate;	
	candidate = {"TT","B-BCB","RL"};
	candidateList.at(50) = candidate;
	candidate = {"TT","B-BCB","RR"};	
	candidateList.at(51) = candidate;

	// TTTT and TTT
	// CB-TT-BC: 52-53
	candidate = {"TTTT","CB-TT-BC","LRLR"};
	candidateList.at(52) = candidate;
	candidate = {"TTTT","CB-TT-BC","RLRL"};
	candidateList.at(53) = candidate;
	// CB-TT-B: 54-55
	candidate = {"TTTT","CB-TT-B","LRLR"};
	candidateList.at(54) = candidate;
	candidate = {"TTTT","CB-TT-B","RLRL"};
	candidateList.at(55) = candidate;
	// CB-T-BCB: 56-57
	candidate = {"TTT","CB-T-BCB","LRL"};
	candidateList.at(56) = candidate;
	candidate = {"TTT","CB-T-BCB","RLR"};
	candidateList.at(57) = candidate;
	// CB-T-BC: 58-59
	candidate = {"TTT","CB-T-BC","LRL"};
	candidateList.at(58) = candidate;
	candidate = {"TTT","CB-T-BC","RLR"};
	candidateList.at(59) = candidate;
	// CB-T-B: 60-61
	candidate = {"TTT","CB-T-B","LRL"};
	candidateList.at(60) = candidate;
	candidate = {"TTT","CB-T-B","RLR"};
	candidateList.at(61) = candidate;
	// B-TT-BC: 62-63
	candidate = {"TTTT","B-TT-BC","LRLR"};
	candidateList.at(62) = candidate;
	candidate = {"TTTT","B-TT-BC","RLRL"};
	candidateList.at(63) = candidate;
	// B-TT-B: 64-65
	candidate = {"TTTT","B-TT-B","LRLR"};
	candidateList.at(64) = candidate;
	candidate = {"TTTT","B-TT-B","RLRL"};
	candidateList.at(65) = candidate;
	// B-T-BCB: 66-67
	candidate = {"TTT","B-T-BCB","LRL"};
	candidateList.at(66) = candidate;
	candidate = {"TTT","B-T-BCB","RLR"};
	candidateList.at(67) = candidate;	
  // B-T-BC: 68-69
	candidate = {"TTT","B-T-BC","LRL"};
	candidateList.at(68) = candidate;
	candidate = {"TTT","B-T-BC","RLR"};
	candidateList.at(69) = candidate;	
	// B-T-B: 70-71
	candidate = {"TTT","B-T-B","LRL"};
	candidateList.at(70) = candidate;
	candidate = {"TTT","B-T-B","RLR"};
	candidateList.at(71) = candidate;	
	// BCB-T-BC: 72-73
	candidate = {"TTT","BCB-T-BC","LRL"};
	candidateList.at(72) = candidate;
	candidate = {"TTT","BCB-T-BC","RLR"};
	candidateList.at(73) = candidate;	
	// BCB-T-B: 74-75
	candidate = {"TTT","BCB-T-B","LRL"};
	candidateList.at(74) = candidate;
	candidate = {"TTT","BCB-T-B","RLR"};
	candidateList.at(75) = candidate;	

	// Dubins-like : 64-71
	candidate = {"Dubins","B-S-B","LSL"};
	candidateList.at(76) = candidate;
	candidate = {"Dubins","B-S-B","LSR"};
	candidateList.at(77) = candidate;			
	candidate = {"Dubins","B-S-B","RSL"};
	candidateList.at(78) = candidate;
	candidate = {"Dubins","B-S-B","RSR"};
	candidateList.at(79) = candidate;	
	candidate = {"Dubins","B-B-B","LRL"};
	candidateList.at(80) = candidate;
	candidate = {"Dubins","B-B-B","LRL"};
	candidateList.at(81) = candidate;	
	candidate = {"Dubins","C-S-C","LRL"};
	candidateList.at(82) = candidate;
	candidate = {"Dubins","C-S-C","RLR"};
	candidateList.at(83) = candidate;							
	return candidateList;
}


