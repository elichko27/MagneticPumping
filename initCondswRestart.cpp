//initCondswRestart.cpp
#include <iostream>
#include <cstring>
#include <cmath>
#include <set>
#include "initCondswRestart.h"
#include "initConds.h"


void initCondswRestart(int &fMatFinalCounter, double * fMatFinalOld, 
		       double * fMatFinalNew, double * fMatFinal, 
		       int &Ntsteps, bool isRestart, int Nvsteps, int n, int expansionLevel, double vthe, 
		       double delV, double vMax, char * distChoice, char * rChoice, char * fileTag, 
		       double tMax, double tMin, double delT, int fMatFinalTMax, 
		       int downSampleT, int downSampleV){
  if (isRestart == 0) {

    for(int i=0; i<Nvsteps*n*expansionLevel;i++) {
      fMatFinalOld[i] = 0.0; 
      fMatFinalNew[i] = 0.0; 
    }
    initConds(fMatFinalOld, Nvsteps, n, expansionLevel, vthe, delV, vMax, distChoice );

    for(int i=0; i<Nvsteps/downSampleV*n*expansionLevel;i++) {
      fMatFinal[i] = 0.0; 
    }

    fMatFinalCounter = 1; 

  } else {
    std::cout << "Should check isRestart - not set up to restart from files yet\n";
  }

}
