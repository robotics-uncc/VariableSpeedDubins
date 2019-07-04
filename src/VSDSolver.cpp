// VSDSolver.cpp
#include<VSDSolver.h>
#include<limits> 

// constructor
// -----------------------------------------------------------------------------
VarSpeedDubins::Solver::Solver(){
  m_numPaths = 76; // max 76
  m_numPathsDubins = 8;
  m_solveStatus = -1; // default -1 unsolved
  m_solveStatusAll = arma::zeros<arma::vec>(84);
  // initialize to infinite cost
  m_candCosts = arma::ones<arma::vec>(84) * std::numeric_limits<double>::max();
  m_candList = VarSpeedDubins::candidateList();
  m_paramGuessSuppliedFlag = false;
  m_Lmax = 10.0; // default
  m_numSamplesMax = 250;
  m_samplingAlgorithm = "barycentric";
  m_equalCostTolerance = 0.0001;
  m_multipleSolnFlag = 0;
}

// set_problemStatement
// -----------------------------------------------------------------------------
void VarSpeedDubins::Solver::set_problemStatement(VarSpeedDubins::ProblemStatement vsdprob){
  m_vsdprob = vsdprob;
  // set dubins problem and solver with "bang" turns
  m_rdprob_bang.set_stateFinal(m_vsdprob.get_xFinal(), m_vsdprob.get_yFinal(), 
                               m_vsdprob.get_hFinal() );
  m_rdprob_bang.set_minTurningRadius(m_vsdprob.get_R());
  m_rdsolver_bang.set_problemStatement(m_rdprob_bang);
  // set dubins problem and solver with "cornering" turns
  m_rdprob_corn.set_stateFinal(m_vsdprob.get_xFinal(), m_vsdprob.get_yFinal(), 
                               m_vsdprob.get_hFinal() );
  m_rdprob_corn.set_minTurningRadius(m_vsdprob.get_r());
  m_rdsolver_corn.set_problemStatement(m_rdprob_corn);
};

// set_NLP_options
// -----------------------------------------------------------------------------
void VarSpeedDubins::Solver::set_NLP_options(SmartPtr<IpoptApplication> &app){
    app->RethrowNonIpoptException(true);
  	// Change some options
  	// Note: The following choices are only examples, they might not be
  	//       suitable for your optimization problem.
  	app->Options()->SetNumericValue("tol", 1e-7);
  	app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    //app->Options()->SetStringValue("derivative_test","first-order");
  	//app->Options()->SetStringValue("output_file", "ipopt.out");
	  app->Options()->SetIntegerValue("print_level", 0);
    app->Options()->SetStringValue("print_user_options", "no");
  	// The following overwrites the default name (ipopt.opt) of the
  	// options file
  	// app->Options()->SetStringValue("option_file_name", "hs071.opt");
}

// solveAll
// -----------------------------------------------------------------------------
void VarSpeedDubins::Solver::solveAll(){
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  set_NLP_options(app);
  // Initialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    std::cout << std::endl << std::endl 
              << "*** Error during initialization!" << std::endl;
  }
  // iterate through all candidates
	std::vector<std::string> cand;
  std::string pathClass, pathType, pathOrientation;
  for (int i = 0; i < m_numPaths; i++){
    //std::cout << "Computing candidate : " << i << std::endl;  
    // create a new NLP
    m_vsdnlp = new VSDNLP();
    m_vsdnlp->set_problemStatement(m_vsdprob);
    // create a new candidate path
    VarSpeedDubins::Path vsdcand;
		cand = m_candList[i];
		pathClass = cand[0];
		pathType = cand[1];
		pathOrientation = cand[2];
		vsdcand.set_turnRadii(0.3, 1.0);
		vsdcand.set_pathClass(pathClass);
		vsdcand.set_pathType(pathType);
		vsdcand.set_pathOrientation(pathOrientation);
    // pass candidate to nlp
    m_vsdnlp->set_candidate(vsdcand);    
    m_vsdnlp->set_Lmax(m_Lmax);
    // set starting point method
    if (m_paramGuessSuppliedFlag){
      //std::cout << "starting point manually set" << std::endl;
      m_vsdnlp->set_startingPoint(m_paramGuess);
    }    
    else {
      //std::cout << "starting point will be computed" << std::endl;
      m_vsdnlp->set_numSamplesMax(m_numSamplesMax);
      m_vsdnlp->set_samplingAlgorithm(m_samplingAlgorithm);
    }
    m_vsdnlp->initializeCandidate();
    // ask ipopt to solve the problem
    //std::cout << "optimize" << std::endl;
    status = app->OptimizeTNLP(m_vsdnlp);
    // process results
    //std::cout << " save results " << std::endl;
    if (status == Solve_Succeeded) {
      //std::cout << "sovle_succeeded" << std::endl;
      arma::vec optSoln = m_vsdnlp->get_optSoln();
  		vsdcand.set_paramsShort(optSoln);
  		vsdcand.computeEndpoint();	
      //std::cout << " candCosts" << std::endl;
      m_candCosts(i) = vsdcand.get_cost();
      // check additional suboptimality conditions
      bool suboptFlag = VarSpeedDubins::checkSuboptimality(
                                                      vsdcand.get_paramsShort(), 
                                                      vsdcand.get_pathClass(),
                                                      vsdcand.get_pathType(),
                                                      vsdcand.get_R(),
                                                      vsdcand.get_r());
      //std::cout << " subopt flag" << std::endl;
      if (suboptFlag == true){
		    m_solveStatus = 2; // feasible but suboptimal
      }
      else {
		    m_solveStatus = 1; // feasible
      }
      //std::cout << "*** The problem solved!" << std::endl;
      //vsdcand.print();
    }
    else {
      m_solveStatus = 0;
      //std::cout << "*** The problem FAILED!" << std::endl;
    }
    //std::cout << " m_candPaths.push_back(vsdcand); " << std::endl;
    m_candPaths.push_back(vsdcand);
    m_solveStatusAll(i) = m_solveStatus;
  }
  // solve Dubins cases
  //std::cout << "dubins solve()" << std::endl;
  m_rdsolver_bang.solve();
  m_rdsolver_corn.solve();
  // update candidate dubins paths
  m_candPathsDubins.push_back(m_rdsolver_bang.get_LSL());
  m_candPathsDubins.push_back(m_rdsolver_bang.get_LSR());
  m_candPathsDubins.push_back(m_rdsolver_bang.get_RSL());
  m_candPathsDubins.push_back(m_rdsolver_bang.get_RSR());
  m_candPathsDubins.push_back(m_rdsolver_bang.get_LRL());
  m_candPathsDubins.push_back(m_rdsolver_bang.get_RLR());
  m_candPathsDubins.push_back(m_rdsolver_corn.get_LRL());
  m_candPathsDubins.push_back(m_rdsolver_corn.get_RLR());
  //std::cout << "dubins" << std::endl;
  int j = 0;
  for (int i = 76; i < 84; i++){
    // check status
    std::string statusString = (m_candPathsDubins.at(j)).get_solutionStatus();
    if ( (statusString.compare("feasible") == 0) || 
         (statusString.compare("optimal") == 0) ){
      m_solveStatusAll(i) = 1;
    }
    else {
      m_solveStatusAll(i) = 0;
    }
    m_candCosts(i) = (m_candPathsDubins.at(j)).get_cost();
    j++;
  }
  // sort results
  //std::cout << " sortFeasiblePaths() " << std::endl;
  sortFeasiblePaths();
}

// solveGivenNLP
// -----------------------------------------------------------------------------
int VarSpeedDubins::Solver::solveGivenNLP(){
	// create an ipopt application
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  set_NLP_options(app);
  // Initialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();

  if (status != Solve_Succeeded) {
    //std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
    return (int) status;
  }
  if (m_paramGuessSuppliedFlag){
    //std::cout << "starting point manually set" << std::endl;
    m_vsdnlp->set_startingPoint(m_paramGuess);
  }    
  else {
    //std::cout << "starting point will be computed" << std::endl;
    m_vsdnlp->set_numSamplesMax(m_numSamplesMax);
    m_vsdnlp->set_samplingAlgorithm(m_samplingAlgorithm);
  }
  //std::cout << "initializeCandidate" << std::endl;
  m_vsdnlp->initializeCandidate();
  //std::cout << "print" << std::endl;
  //m_vsdnlp->print();
  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(m_vsdnlp);

  if (status == Solve_Succeeded) {
    //std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
    // check additional suboptimality conditions
//    bool suboptFlag = VarSpeedDubins::checkSuboptimality(
//                                                      vsdcand.get_paramsShort(), 
//                                                      vsdcand.get_pathClass(),
//                                                      vsdcand.get_pathType(),
//                                                      vsdcand.get_R(),
//                                                      vsdcand.get_r());
    bool suboptFlag = false;
    if (suboptFlag == true){
      m_solveStatus = 2; // feasible but suboptimal
    }
    else {
      m_solveStatus = 1; // feasible
    }
  }
  else {
    //std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
		m_solveStatus = 0;
  }

  	// As the SmartPtrs go out of scope, the reference count
  	// will be decremented and the objects will automatically
  	// be deleted.
  	return (int) status;
}

// print all
// -----------------------------------------------------------------------------
void VarSpeedDubins::Solver::printAll(){
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "          VarSpeedDubins::Solver Output          " << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "No., Class, Type, Dir., Status, Cost, p1,p2,p3,p4" << std::endl;
  arma::vec params; 
  int numFeasiblePaths = 0;
  for (int i = 0; i < m_numPaths; i++){
    // print if a feasible solution was found
    if (m_solveStatusAll(i) != 0){
      VarSpeedDubins::Path curPath = m_candPaths.at(i);  
      params = curPath.get_paramsShort();
      std::cout << i << ")  " << curPath.get_pathClass() << "," 
                << curPath.get_pathType() << "," << curPath.get_orientation() 
                << "," << m_solveStatusAll(i) << ",\t" << curPath.get_cost() 
                << " \t," << params(0) << "," << params(1) << " , " 
                << params(2) << "," << params(3) << std::endl;
      numFeasiblePaths++;
    }
  }
  // show robust dubins cases
  int j = 0;
  for (int i = m_numPaths; i < m_numPaths + m_numPathsDubins; i++){
    // print if a feasible solution was found
    if (m_solveStatusAll(i) == 1){

      std::string pathClass;
      if (j < 6){
        pathClass = "Dubins-R";    
      }
      else {
        pathClass = "Dubins-r";
      }
      RobustDubins::Path curPath = m_candPathsDubins.at(j); 
      std::cout << i << ")  " << pathClass 
                << "," << curPath.get_pathType() << "," 
                << m_solveStatusAll(i)  << ",\t" << curPath.get_cost() << " \t," 
                << curPath.get_aParamUnsigned() << ","  
                << curPath.get_bParamUnsigned() << ","  
                << curPath.get_cParamUnsigned() << std::endl;
      numFeasiblePaths++;
    }
    j++;
  }
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << " No. of Candidates Found: " <<  numFeasiblePaths << std::endl;
}


// sortFeasiblePaths
// -----------------------------------------------------------------------------
void VarSpeedDubins::Solver::sortFeasiblePaths(){
  //std::cout << "sortFeasiblePaths() start" << std::endl;
  // sort for the lowest cost
  m_optPathTypes = MathTools::minIndicesWithTolerance(m_candCosts, 
                                                      m_equalCostTolerance);
  // set size of output matrices to accomodate solns
  m_optPathParams.set_size(4,m_optPathTypes.size());
  m_optPathCosts.set_size(m_optPathTypes.size());
  if (m_optPathTypes.size() > 1){
    m_multipleSolnFlag = 1;
  }
  //std::cout << "m_optPathTypes" << std::endl;
  //MathTools::debugPrintVector(m_optPathTypes, true);
  //m_candCosts.print("m_candCost");

  // save outputs
  for (int i = 0; i < m_optPathTypes.size(); i++){
    if (m_optPathTypes[i] < 76){  // processing a var speed dubins path
      m_optPathParams.col(i)=m_candPaths.at(m_optPathTypes[i]).get_paramsShort();
      m_optPathCosts(i) = m_candPaths.at(m_optPathTypes[i]).get_cost();
    }
    else { // process a dubins path
      // get dubins path params
      double a = m_candPathsDubins.at(m_optPathTypes[i] - 76).get_aParamUnsigned();
      double b = m_candPathsDubins.at(m_optPathTypes[i] - 76).get_bParamUnsigned();
      double c = m_candPathsDubins.at(m_optPathTypes[i] - 76).get_cParamUnsigned();
      //std::cout << "m_optPathParams " << m_optPathParams.size() << std::endl;   
      //std::cout << " a " << a << std::endl;
      m_optPathParams(0,i) = a;
      m_optPathParams(1,i) = b;
      m_optPathParams(2,i) = c;
      m_optPathParams(3,i) = -1; // undefined
      //std::cout << "m_optPathCosts " << m_optPathCosts.size() << std::endl;   
      m_optPathCosts(i) = m_candPathsDubins.at(m_optPathTypes[i] - 76).get_cost();
    }
  }
  //std::cout << "sortFeasiblePaths() end" << std::endl;
  return;
}

// plotAll
// -----------------------------------------------------------------------------
void VarSpeedDubins::Solver::plotAll(std::string fileName){
	std::string xVarName = "x";
  std::string yVarName = "y";
  // plot all the variable speed dubins paths
  for (int i = 0; i < m_numPaths; i++){
    if (m_solveStatusAll(i) != 0){
      VarSpeedDubins::Path curPath = m_candPaths.at(i);    
      curPath.computePathHistory();
      curPath.computeEndpoint();
      curPath.writeEndpointPlotCommands(fileName,1,xVarName,yVarName,"r");
      curPath.writeSwitchingPtPlotCommands(fileName,1,xVarName,yVarName,"ro");
      if (m_solveStatusAll(i) == 2){
        curPath.writePathPlotCommands(fileName,1,xVarName,yVarName,"m");	
      }
      else if (m_solveStatusAll(i) == 1){
        curPath.writePathPlotCommands(fileName,1,xVarName,yVarName,"b");	
      }
    }
  }
  MathTools::axisProperty(fileName, "equal");
  // plot the dubins paths
  m_rdsolver_bang.writeOctaveCommandsToPlotSolution(fileName, 2);
  m_rdsolver_corn.writeOctaveCommandsToPlotSolution(fileName, 2);
  MathTools::axisProperty(fileName, "equal");
}

//// plotOptimal
//// -----------------------------------------------------------------------------
//void VarSpeedDubins::Solver::plotOptimal(std::string fileName){
//	std::string xVarName = "x_optimal";
//  std::string yVarName = "y_optimal";
//  // plot all the variable speed dubins paths
//  for (int i = 0; i < m_numPaths; i++){
//    if (m_solveStatusAll(i) != 0){
//      VarSpeedDubins::Path curPath = m_candPaths.at(i);    
//      curPath.computePathHistory();
//      curPath.computeEndpoint();
//      curPath.writeEndpointPlotCommands(fileName,1,xVarName,yVarName,"r");
//      curPath.writeSwitchingPtPlotCommands(fileName,1,xVarName,yVarName,"ro");
//      if (m_solveStatusAll(i) == 2){
//        curPath.writePathPlotCommands(fileName,1,xVarName,yVarName,"m");	
//      }
//      else if (m_solveStatusAll(i) == 1){
//        curPath.writePathPlotCommands(fileName,1,xVarName,yVarName,"b");	
//      }
//    }
//  }
//  MathTools::axisProperty(fileName, "equal");
//  // plot the dubins paths
//  m_rdsolver_bang.writeOctaveCommandsToPlotSolution(fileName, 2);
//  m_rdsolver_corn.writeOctaveCommandsToPlotSolution(fileName, 2);
//  MathTools::axisProperty(fileName, "equal");
//}




