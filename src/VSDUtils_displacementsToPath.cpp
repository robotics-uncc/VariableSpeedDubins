// VCSDUtils_displacementsToPath.cpp
// Last Modefied: 25-Jan-2016, Artur Wolek

#include<VSDUtils.h>

// returns continuous state history from several vector displacement segments
void VarSpeedDubins::displacementsToPath(arma::vec &segJoined, arma::vec &seg1,
		                                     arma::vec &seg2){
  // check if the inputs are empty
	bool seg1empty = (seg1.n_elem == 0);
	bool seg2empty = (seg2.n_elem == 0);
	// if all segments are non-empty (non-empty, non-empty, non-empty)~(n,n)
	if (!seg1empty && !seg2empty){
		// it is assumed that the first element of each segment is zero
		// so when joining several segments we remove these zeros when 
		// necessary to avoid duplicate entries
		//
		// in this case, seg1 will start with zero
		// when computing seg2shifted, we add the last element of seg1 to all
		// the elements of seg2, excluding the first element of seg2 which is 
		// a zero
		arma::vec seg2shifted = seg1(seg1.n_elem-1) + 
			                      seg2(arma::span(1,seg2.n_elem-1));
		segJoined = arma::join_vert(seg1, seg2shifted);	
	}
	// (n,e)
	else if (!seg1empty && seg2empty){
		segJoined = seg1;
	}	
	// (e,e)
	else if (seg1empty && seg2empty){
    // all inputs are empty
    return;
	}	
	// (e,n)
	else if (seg1empty && !seg2empty){
		segJoined = seg2;
	}	
  return;
}

void VarSpeedDubins::displacementsToPath(arma::vec &segJoined, arma::vec &seg1,
		                                     arma::vec &seg2, arma::vec &seg3){
	displacementsToPath(segJoined, seg1, seg2);
	displacementsToPath(segJoined, segJoined, seg3);
	return;
}

// returns continuous state history from several matrix displacement segments
void VarSpeedDubins::displacementsToPath(arma::mat &segJoined, arma::mat &seg1, 
		                                     arma::mat &seg2){
	// check if segments are empty
  // since the segments always start with zero, we still call a matrix
  // with only 1 row empty
	bool seg1empty = (seg1.n_rows <= 1);
	bool seg2empty = (seg2.n_rows <= 1);
  //std::cout << "seg1empty: " << seg1empty << std::endl;
  //std::cout << "seg2empty: " << seg2empty << std::endl;
	// (n,n) both non-empty
	if (!seg1empty && !seg2empty){	
    //std::cout << "both not empty" << std::endl;
		// delete first row of zeros and shift matrix	
    // Note: this assumes there are at least two rows!
		arma::rowvec a = seg1.row(seg1.n_rows-1);
		arma::mat B = seg2.rows(1,seg2.n_rows-1);
		arma::mat seg2shifted = MathTools::addVectorToMatrixColoumnWise(a, B);
    //std::cout << "join_vert" << std::endl;
		// join to initial segment
		segJoined = arma::join_vert(seg1, seg2shifted);
	}
	// (n,e)
	else if (!seg1empty && seg2empty){	
		segJoined = seg1;
	}
	// (e,n)
	else if (seg1empty && !seg2empty){	
		segJoined = seg2;
	}	
	// (e,e)
	else if (seg1empty && seg2empty){
		std::cout << "displacementsToPath(), Error: all inputs are empty" 
			        << std::endl; 	
	}	
	return;
}

// returns continuous state history from several matrix displacement segments
void VarSpeedDubins::displacementsToPath(arma::mat &segJoined, arma::mat &seg1,
	                                       arma::mat &seg2, arma::mat &seg3){
  //std::cout << "join seg1 and seg2" << std::endl;
  //seg1.print("seg1");
  //seg2.print("Seg2");
	displacementsToPath(segJoined, seg1, seg2);

  ///segJoined.print("segJoined");
  //std::cout << "join seg2 and seg3" << std::endl;
  //segJoined.print("segJoined");
  //seg3.print("Seg3");
	displacementsToPath(segJoined, segJoined, seg3);
	return;
}

void VarSpeedDubins::displacementsToPath(arma::mat &segJoined, arma::mat &seg1,
	                                       arma::mat &seg2, arma::mat &seg3, 
                                         arma::mat &seg4){
	displacementsToPath(segJoined, seg1, seg2);
	displacementsToPath(segJoined, segJoined, seg3);
	displacementsToPath(segJoined, segJoined, seg4);
	return;
}

