#include<MathToolsArma_Combinatronics.h>
#include<MathToolsArma_VectorOperators.h>
#include<MathTools.h>

// permutations
arma::mat MathTools::permute(arma::vec x){
  // get size of x
  int n = x.n_elem;
  int j = 0;
  int nfact = factorial(n);
  arma::mat xperm(n,factorial(n));
  MathTools::heapPermute(n, x, j, xperm);
  return xperm;
}



// Heap's permutation algorithm
void MathTools::heapPermute(int n, arma::vec &x, int &j, arma::mat &xperm){
  if (n == 1){
    (x.t()).print();
    xperm.col(j) = x;
    j++;
  }
  else {
    for (int i = 0; i < n ; i++){
      heapPermute(n - 1, x, j, xperm);
      if (MathTools::isEven(n)){
        MathTools::swap(x, i, n - 1);   
      }
      else { // n is odd
        MathTools::swap(x, 0, n - 1); 
      }
    }
  }
}


// returns matrix with rows of all permutations of integer pairs that sum to n 
arma::mat MathTools::constantSumPairs(int n){
  arma::mat out(n+1,2); // initialize
  // generate sequence: {(n,0), (n-1,1) , ... , (0,n)}
  for (int i = 0; i < n+1; i++){ 
    out(i,0) = i; 
    out(i,1) = n-i;
  }
  return out;
}

// returns a matrix with rows that contain all possible sequences of m integers 
// (including zero) that sum to n 
arma::mat MathTools::constantSumSubsets(int n, int m){
  int numPerms = MathTools::binomial(n+m-1, m-1); // get number of permutations
  arma::mat output = arma::ones<arma::mat>(numPerms, m); // initialize output 
  int curDepth = 0; // col of output matrix being populated
  int curRow = 0; // row of output matrix being populated
  constantSumSubsets(n, m, curDepth, curRow, output);
  return output;
}

// main algorithm: recursively populates the matrix with the desired sequences
void MathTools::constantSumSubsets(int n, int &m, int &curDepth, int &curRow, 
                                   arma::mat &seqMat){
  if (curDepth == m - 2){ // stopping condition     
    // last two cols of seqMat will now be populated
    arma::mat lastCols = MathTools::constantSumPairs(n); 
    // iterate through the list of constant sum pairs
    for (int i = 0; i < lastCols.n_rows; i++){
      // at each iteration insert a row into seqMat
      if (i != 0){ // once curRow is incremented need to copy previous values
        for (int j = 0; j < curDepth; j++){ // populate cols 0 to curDepth 
          seqMat(curRow, j) = seqMat(curRow - 1, j);         
        }
      }
      // insert new values
      seqMat(curRow, curDepth) = lastCols(i,0); // second-to-last col
      seqMat(curRow, curDepth + 1) = lastCols(i,1); // last col
      // increment row counter
      curRow++;
    }
    return;
  }
  else {
    // generate a list of constant sum pairs from n
    arma::mat listPairs = MathTools::constantSumPairs(n);
    // iterate through this list
    for (int i = 0; i < listPairs.n_rows; i++){
      if (i != 0){ // for all but the first pair
        // copy values between cols 0 and curDepth from previous row to curRow
        for (int j = 0; j < curDepth; j++){ 
          seqMat(curRow, j) = seqMat(curRow - 1, j);       
        }
      }
      // insert first value of listPair at the current depth 
      seqMat(curRow, curDepth) = listPairs(i,0); 
      curDepth++; // prepare to insert the second value
      // if the second value is not the last col, then we need to continue 
      // expanding the curDepth col recurisvely 
      MathTools::constantSumSubsets(listPairs(i,1), m, curDepth, curRow, seqMat);
      curDepth--; // return to previos col
    } 
  }
}



