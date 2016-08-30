//nuF.cpp
#include <iostream>
#include <cstring>
#include <cmath>
#include <set>
#include "nuF.h"

void nuF(double * nuFac, double nu, int n, double vthe, 
	 double delV, int Nvsteps, char * scattType) {

  double fac1; 

  for(int j = 0; j < n; j++){
    fac1 = -1*double(j)*(1+double(j)); 
    if (j == 0) { 
      for(int i=0;i<Nvsteps;i++) nuFac[j*Nvsteps + i] = 0.0;
    } else if (strcmp(scattType, "Constant") == 0) { 
      for(int i=0;i<Nvsteps;i++) nuFac[j*Nvsteps + i] = nu*fac1; 
    } else if (strcmp(scattType, "Lorentz") == 0) {
      for(int i=1;i<Nvsteps;i++) nuFac[j*Nvsteps + i] = fac1/pow(double(i)*delV/vthe,3.0); 
      nuFac[0] = nuFac[1]; 
    } else {
      std::cout << "The code does not support this choice for scattType at this time. \n"; 
    }
  }

}
