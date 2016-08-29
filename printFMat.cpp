#include <iostream>
#include "printFMat.h"

void printFMat(double * fMatFinalOld, int expansionLevel, int n, int Nvsteps) {
// Printing the last fMatFinalOld just to see what happens
  for (int j = 0; j<Nvsteps*n*expansionLevel; j++) { 
    if(j%(Nvsteps) == 0 && j%n == 0 && j != 0)  std::cout << std::endl; 
    if(j%(Nvsteps) == 0 && j%expansionLevel == 0 && j != 0) std::cout << std::endl;
    if(j%(Nvsteps) == 0 && j%expansionLevel == 0 && j%n == 0) std::cout << std::endl <<"expanLevel = " 
									<< j/Nvsteps/n << std::endl;
    std::cout << " " << fMatFinalOld[j] << " "; 
  }

}
