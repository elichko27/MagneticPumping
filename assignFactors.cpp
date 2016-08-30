// assignFactors.cpp
#include <cstdlib>
#include <iostream>
//#include <findFac.h>
#include "assignFactors.h" 

void assignFactors(double * facMat, int * expanMat, int expansionLevel) { 

  for(int i = 0; i < expansionLevel; i++) { 
    for(int j = 0; j < expansionLevel; j++) { 
      facMat[i*expansionLevel + j] = 1.0;//findFac(i,j);

      expanMat[i*expansionLevel + j] = -1; 
    }
  }

  int lastSpot = 0; 
  for(int i = 0; i < expansionLevel; i++) { 
    for(int j = 0; j < expansionLevel; j++) { 

      if(facMat[i*expansionLevel + j] != 0.0) { 
	facMat[lastSpot] = j; 
	lastSpot++;
      }

      if (j == (expansionLevel - 1)) lastSpot = double(i+1)*double(expansionLevel);

    }
  }

}
