#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <set>
#include "initConds.h"

void initConds(double * fMatFinalNew, int Nvsteps, int n, int expansionLevel, 
	       double vthe, double delV, double vMax, char * distChoice ) {

  if (strcmp(distChoice, "Maxwellian") == 0) {
    for(int i=0; i<Nvsteps;i++) 
      for(int j=0; j<n; j++) 
	for(int k=0; k<expansionLevel; k++) {
	  if (j == 0 && k ==0) {
	    fMatFinalNew[k*n*Nvsteps + j*Nvsteps + i] = exp(-1*pow(double(i)*delV,2)/pow(vthe,2)); 
	    //if (i == 0) fMatFinalNew[k*n*Nvsteps + j*Nvsteps + i] = 1.0; 
	  } else {
	    fMatFinalNew[k*n*Nvsteps + j*Nvsteps + i] = 0.0; 
	  }
	}
  } else if (strcmp(distChoice, "Kappa") == 0) {
    double kappa = 3.0; 
    double gammaFac = tgamma(kappa+1.0)/tgamma(kappa-1./2.)/tgamma(3./2.); 
    double pi = 3.1415926535;
    for(int i=0; i<Nvsteps;i++) 
      for(int j=0; j<n; j++) 
	for(int k=0; k<expansionLevel; k++) {
	  if (j == 0 && k ==0) {
	    fMatFinalNew[k*n*Nvsteps + j*Nvsteps + i] = pow(kappa*vthe,3./2.)/2./pi*gammaFac*
	      pow(1.+pow(double(i)*delV,2.)/kappa/pow(vthe, 2.), -1.*kappa-1.); 
	  } else {
	    fMatFinalNew[k*n*Nvsteps + j*Nvsteps + i] = 0.0; 
	  }
	}

  } else {
    std::cout << "Incorrect choice for distChoice - must be Maxwellian or Kappa\n";
  }

}
