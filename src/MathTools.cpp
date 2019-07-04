#include<MathTools.h>

std::vector<int> MathTools::minIndicesWithTolerance(
		                                            std::vector<double> testVector, 
                                                double tolerance){
	// initialize vector of min indices
	std::vector<int> minIndices;
	// store index of smallest value in testVector
	minIndices.push_back(minElement(testVector)); 
	// store smallest value for reference
	double minValueRef = testVector[minIndices[0]]; 
	// set this element to largest possible value since it has been recorded
	testVector[minIndices[0]] = std::numeric_limits<double>::max(); 
	int i = 1;
	// check if there are additional minimum values within the specified tol
	while ( std::abs(minValue(testVector) - minValueRef) <= tolerance ){
		minIndices.push_back(minElement(testVector)); // store additional index
		// set to largest possible value	
		testVector[minIndices[i]] = std::numeric_limits<double>::max(); 
		i++;
	}
	return minIndices;
}

double MathTools::polarDistance(double a, double b){
  // make both angles positive
  a = MathTools::mod(a, 2.0*M_PI);
  b = MathTools::mod(b, 2.0*M_PI);
  // calculate angular distance without wrapping around 2pi
  double delta = std::abs(b-a);
  // check if wrapping around is shorter
  if (delta > 2.0*M_PI - delta){
    return 2.0*M_PI - delta;
  }
  else {
    return delta;
  }
}

double MathTools::distance(std::vector<double> referenceVector, 
		                       std::vector<double> testVector){
	std::vector<double>::size_type vectorSize = referenceVector.size();
	double distanceSquared = (referenceVector[0]-testVector[0]) * 
							             (referenceVector[0]-testVector[0]);
	for (std::vector<double>::size_type i = 1; i < vectorSize ; i++){
		distanceSquared = distanceSquared + (referenceVector[i]-testVector[i]) 
						                          * (referenceVector[i]-testVector[i]);
	}	
	return std::sqrt(distanceSquared);
}

double MathTools::mod(double number, double base){
  return number - std::floor(number/base)*base;
}

void MathTools::append(std::vector<double> &baseVector, 
					             std::vector<double> &appendVector){
	// Modify a std::vector by appending a second vector to its end
	baseVector.insert(std::end(baseVector), std::begin(appendVector), 
	                  std::end(appendVector));
}

int MathTools::minElement(std::vector<double> testVector){
	// Returns an iterator pointing to the element with the smallest value 
	std::vector<double>::iterator it  = std::min_element(testVector.begin(), 
		testVector.end());
	// Calculates the number of elements between first and last.
	return std::distance( testVector.begin(), it); // returns index not a value
}

double MathTools::minValue(std::vector<double> testVector){
	return testVector[minElement(testVector)];
}

std::vector<double> MathTools::linspace(double xmin, double xmax, 
		                                    double spacing, int nSkip){
	// Generate a squence of evenly spaced points within prescribed bounds
	// if desired skip nSkip several points from the begining
	// note: spacing cannot always be exactly enforced
	std::vector<double> outputVector;	
	double delx = xmax - xmin;
	double dir = MathTools::sign(delx); //arma direction of interval
	// check if request spacing is feasible
	if (std::abs(spacing) >= std::abs(delx)){ // not feasible, manually set 
		outputVector.push_back(xmin);
		outputVector.push_back(xmax);
		std::cout << "MathTools::linspace warning: requested spacing is "
				      << "larger than interval" << std::endl;
	}
	else { // spacing is feasible 
		double nptsDecimal = std::abs(delx/spacing)+1.0;
		int nptsInt = (int)nptsDecimal;
		// check if number of skip points is feasible
		if (nSkip >= nptsInt){
			std::cout << "MathTools::linspace error: number of skip points is "
					      << "not feasible, reduce spacing" << std::endl;
		}
		else {
			// check if spacing divides interval into a whole number
			if (MathTools::isInteger(nptsDecimal)){
				for (int i = nSkip; i < nptsInt ; i++){
					outputVector.push_back(xmin + dir*spacing*i);
				}
			}
			else { // set the last point to xmax, overriding spacing request
				for (int i = 0; i < nptsInt-1; i++){
					outputVector.push_back(xmin + dir*spacing*i);
				}
				outputVector.push_back(xmax);
			}
		}
	}
	return outputVector;
}

std::vector<double> MathTools::linspace(double xmin, double xmax, 
		                                    double spacing){
	return MathTools::linspace(xmin, xmax, spacing, 0);
}

bool MathTools::isInteger(double value){
	if (value == floor(value))
		return true;
	else 
		return false;
}

void MathTools::writeVectorData(std::string fileName, std::vector<double> data, 
		                            std::string varName){
  // open for writing
	std::ofstream outfile(fileName, std::ofstream::app); 
  // write coordinates
  if ( data.size() == 0 ){
    throw std::runtime_error("MathTools::writeVectorData: vector to write is empty.");
  }
  outfile << std::setprecision(17) << varName << "= [ " << data[0]; 
  for (int i = 1; i < data.size(); i++){
      outfile << " , " << data[i];
  }	
  outfile << "]; \n";
  outfile.close();	   	// close file
}

void MathTools::plot1DCurve(std::string fileName, int figureNumber, 
		                        std::string xVarName, std::string yVarName, 
                            std::string style){	
	std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << "figure " << figureNumber << " \n ";
	outfile << "plot ("<< xVarName << " , " << yVarName <<",'" << style << 
		"', 'linewidth',3); hold on; \n";
	outfile << "grid on;" << "\n ";
	outfile << "set(gca, 'linewidth', 3, 'fontsize', 16) \n";
    	// close file
    	outfile.close();
	return;
}

void MathTools::axisProperty(std::string fileName, std::string property){
	// open file for writing
	std::ofstream outfile(fileName, std::ofstream::app);
	outfile << "axis(\"" << property << "\") \n";
	// close file
    	outfile.close();
}

void MathTools::writeMatrixData(std::string fileName, 
		                            std::vector< std::vector<double> > &data, 
                                std::string varName){
	std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << varName << "= [ ";
	// assume that data is given as a vector of row vectors
	std::vector< double > currentRow;
    for (int i = 0; i < data.size(); i++){ // for each row
		currentRow = data.at(i);
		  for (int j = 0; j < currentRow.size(); j++){
			  outfile << " , " << currentRow.at(j);
		  }
      outfile << " ; "; // end of row
    }	
	outfile << "]; \n";
    outfile.close();	   	// close file
}

void MathTools::writeMatrixData(std::string fileName, 
		                            std::vector< std::vector<int> > &data, 
                                std::string varName){
	std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << varName << "= [ ";
	// assume that data is given as a vector of row vectors
	std::vector< int > currentRow;
    for (int i = 0; i < data.size(); i++){ // for each row
		currentRow = data.at(i);
		for (int j = 0; j < currentRow.size(); j++){
			outfile << " , " << currentRow.at(j);
		}
        outfile << " ; "; // end of row
    }	
	outfile << "]; \n";
    outfile.close();	   	// close file
}

void MathTools::runOctaveScript(std::string fileName){
	std::string runCommand = "octave --force-gui --persist " + fileName;
	system(runCommand.c_str());
}

void MathTools::addCommand(std::string fileName, std::string command){
  std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << command << "\n"; // write command
  outfile.close(); // close file
}

bool MathTools::checkIntervalIntersect(double a, double b, double c, double d){
  // suppose interval 1 is [a,b]
  // suppose interval 2 is [c,d]
  // first check if a is contained in [c,d]
  bool aIntersects = MathTools::checkBounds(a, c, d);
  // then check if b is contained in [c,d]
  bool bIntersects = MathTools::checkBounds(b, c, d);
  // the only other possibility is that [c,d] is contained entirely inside [a,b]
  bool cd_in_ab =  ( (a <= c) && ( b >= d) );
  // if any of the above conditions hold then the intervals intersect
  return (aIntersects || bIntersects || cd_in_ab);
}

void MathTools::plotPoint(std::string fileName, int figureNumber, double x, 
						  double y, std::string xVarName, std::string yVarName, 
						  std::string style){
	MathTools::writeVariable(fileName, x, xVarName);
	MathTools::writeVariable(fileName, y, yVarName);
	MathTools::plot1DCurve(fileName, figureNumber, xVarName, yVarName, style);
	return;
}

int MathTools::safeCurveSpacing(double arcLength, double nominalSpacing){
	if ( (arcLength != 0) && (arcLength <= nominalSpacing) ){
		return 2; // at least two points to define a curve
	}	
	else if (arcLength == 0){
		return 0; // curve does not exist
	}
	else {
		return (int)floor(arcLength/nominalSpacing) + 1; // will be at least 
	}
}

int MathTools::pascalsTriangle(int n, int k){
  int entry = 1;
  for (int i = 1; i < k; i++){
    entry = entry*(n+1-i)/(i);
  }
  return entry;
}

int MathTools::binomial(int n, int k){
  // binomial(n, k) = n!/[k!(n-k)!]
  return MathTools::pascalsTriangle(n,k+1);
}

bool MathTools::checkBounds(double var, double lower_bound, double upper_bound){
	return ( (var >= lower_bound) && (var <= upper_bound) );
}

void MathTools::writeVariable(std::string fileName, double x, 
		                          std::string varName){
  std::ofstream outfile(fileName, std::ofstream::app); // open for writing
	outfile << varName << "= " << x << "; \n"; // write variable
  outfile.close(); // close file
}

int MathTools::factorial(int n){
  int nfact = 1;
  for (int i = 0; i < n; i++){
    nfact = nfact*(i+1);
  }
  return nfact;
}

bool MathTools::isEven(int x){
  // (x % 2) outputs remainder after dividing x by 2  
  return (x % 2 == 1);
}

void MathTools::debugPrintVector(std::vector<double> varVector, bool debugFlag){
	if (debugFlag){		
		// check size of vector
		std::vector<double>::size_type vectorSize = varVector.size();
		// iterate through each element and print
		std::cout << " Element  :  Value " << std::endl;
		for (std::vector<double>::size_type i = 0; i <= vectorSize - 1; i ++){ 
			std::cout << i << " : " << varVector[i] << std::endl;
		}
	}
}

void MathTools::debugPrintVector(std::vector<double> varVector){
	MathTools::debugPrintVector(varVector, true);
}

void MathTools::debugPrintVector(std::vector<int> varVector){
	MathTools::debugPrintVector(varVector, true);
}

void MathTools::debugPrintVector(std::vector<int> varVector, bool debugFlag){
	if (debugFlag){		
		// check size of vector
		std::vector<int>::size_type vectorSize = varVector.size();
		// iterate through each element and print
		std::cout << " Element  :  Value " << std::endl;
		for (std::vector<int>::size_type i = 0; i <= vectorSize - 1; i ++){ 
			std::cout << i << " : " << varVector[i] << std::endl;
		}
	}
}

void MathTools::runOctaveScript(std::string fileName, std::string options){
	std::string runCommand = "octave " + options + " " + fileName;
	system(runCommand.c_str());
}
