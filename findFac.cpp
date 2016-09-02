//findFac.cpp 
#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <set>
#include "findFac.h"

double findFac(int nInput, int mInput, char * wavePt) { 

  int n,m, isCos, isSin;
  double pi = 3.141592553589;

  // Translating from the program indexing to the function/actual equation indexing
  n = (nInput+1)/2;
  m = (mInput+1)/2;

  if (nInput%2 == 1 || mInput%2 == 1 || nInput == 0 || mInput == 0) { 
    isCos = 1; 
  } else { 
    isCos = 0; 
  } 

  if ((nInput != 0 && nInput%2 == 0) || (mInput != 0 && mInput%2 == 0)) { 
    isSin = 1; 
  } else { 
    isSin = 0; 
  } 

// Result of the coupling integrals, where isSin denotes the prescence of a
// sine in the integral, and isCos denotes the prescence of a cosine in the
// integral - the form of the integral is determined by the equations
  if (strcmp(wavePt, "1") == 0) {

    if (nInput == 0 && mInput == 0) return 2.*pi; 
    if (nInput == mInput) return pi;

  } else if (strcmp(wavePt, "sin") == 0) {
    if ((nInput == 0 && mInput == 2) || (nInput == 2 && mInput == 0))return pi; 

    if(isCos == 1 && isSin == 1) {

      if (abs(n - m) == 1 && abs(nInput - mInput) == 1) return -pi/2; 
      if(abs(n - m) == 1 && abs(nInput - mInput) == 3) return pi/2; 

    } 
  } else if (strcmp(wavePt, "cos") == 0) {

    if (nInput + mInput == 1) return pi; 
    if (((isCos == 0 && isSin == 1) || (isCos == 1 && isSin == 0)) && (abs(m-n) == 1))  return pi/2; 
  }

  return 0.; 
}
