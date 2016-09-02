//appendfMatNew.cpp
#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include "appendfMatNew.h"
#include "gradOpts.h"

void appendfMatNew(int iN, int jN, int kN, double * fMatTempNew, int iO, int jO, int kO, double * fMatTempOld, 
		   double overallFactor, double * facMat, int * expanMat, char * isGrad, int isOneCosSin, 
		   int expansionLevel, int Nvsteps, int n, double delV) { 


  kO = 0; 
  while (expanMat[isOneCosSin*expansionLevel*expansionLevel + kN*expansionLevel + kO] != -1) { 
     // Taking care of the switching of cos and sin for kpar terms
    if (abs(jN - jO) == 1 && kN != 0 && kO != 0 && jN != 0 && jO != 0) { 
      if (jN%2 == 1) { 
	fMatTempNew[kN*Nvsteps*n + jN*Nvsteps + iN] += double((kN+1)/2)*overallFactor
	  *facMat[isOneCosSin*expansionLevel*expansionLevel + (kN+1)*expansionLevel + kO]
	  *fMatTempOld[kO*Nvsteps*n + jO*Nvsteps + iO];
      } else { 
	fMatTempNew[kN*Nvsteps*n + jN*Nvsteps + iN] += -1.0*double((kN+1)/2)*overallFactor
	  *facMat[isOneCosSin*expansionLevel*expansionLevel + (kN-1)*expansionLevel + kO]
	  *fMatTempOld[kO*Nvsteps*n + jO*Nvsteps + iO];
      }

    // Terms that involve gradients
    } else if (strcmp(isGrad, "Grad") == 0) { 
      fMatTempNew[kN*Nvsteps*n + jN*Nvsteps + iN] += overallFactor
	*facMat[isOneCosSin*expansionLevel*expansionLevel + kN*expansionLevel + kO]
	*gradOpts(fMatTempOld, delV, iO, 0, Nvsteps, kO*Nvsteps + jO*Nvsteps); 

    // Terms that don't involve gradients
    } else { 
      fMatTempNew[kN*Nvsteps*n + jN*Nvsteps + iN] += overallFactor
	*facMat[isOneCosSin*expansionLevel*expansionLevel + kN*expansionLevel + kO]
	*fMatTempOld[kO*Nvsteps*n + jO*Nvsteps + iO];
    }
    kO++; 
  } 

}
