// assignFactors.cpp
#include <cstdlib>
#include <iostream>
#include <cmath> 
#include "findFac.h"
#include "assignFactors.h" 

void assignFactors(double * facMat, int * expanMat, int expansionLevel) { 

  char oneC[] = "1"; 
  char sinC[] = "sin"; 
  char cosC[] = "cos"; 

  for(int i = 0; i < expansionLevel; i++) { 
    for(int j = 0; j < expansionLevel; j++) { 
      facMat[0*expansionLevel*expansionLevel + i*expansionLevel + j] = findFac(i,j, oneC)/findFac(i,i,oneC);
      facMat[1*expansionLevel*expansionLevel + i*expansionLevel + j] = findFac(i,j, cosC)/findFac(i,i,oneC);
      facMat[2*expansionLevel*expansionLevel + i*expansionLevel + j] = findFac(i,j, sinC)/findFac(i,i,oneC);

      expanMat[0*expansionLevel*expansionLevel + i*expansionLevel + j] = -1; 
      expanMat[1*expansionLevel*expansionLevel + i*expansionLevel + j] = -1;
      expanMat[2*expansionLevel*expansionLevel + i*expansionLevel + j] = -1;
    }
  }

  int lastSpot = 0; 
  for(int k = 0; k < 3; k++) {
    for(int i = 0; i < expansionLevel; i++) { 
      for(int j = 0; j < expansionLevel; j++) { 
	
	if(facMat[k*expansionLevel*expansionLevel + i*expansionLevel + j] != 0.0) { 
	  expanMat[lastSpot] = j; 
	  lastSpot++;
	}
	
	if (j == (expansionLevel - 1)) lastSpot = double(k)*pow(double(expansionLevel),2.0) + double(i+1)*double(expansionLevel);
	
      }
    }
  }

}
