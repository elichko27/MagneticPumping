//writeResults.cpp
#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <set>
#include <ctime>
#include "writeResults.h"

void writeResults(int i, int downSampleT, int downSampleV, double * fMatFinalOld, double *fMatFinal, 
		  int Nvsteps, int n, int expansionLevel, std::ofstream & outputFile, FILE * ptr_fp, 
		  time_t begin, time_t end) {
  if (i%downSampleT == 0) {
    if (downSampleV == 1) {
      time(&end);
      fwrite(fMatFinalOld, Nvsteps*n*expansionLevel*sizeof(double), 1, ptr_fp);
      outputFile << i << " tSteps in : " << difftime(end, begin) << " seconds" << std::endl;
    } else {
      
      int jcount = 0; 
      for(int j = 0; j<Nvsteps; j+=downSampleV) {
	for(int k = 0; k<n; k++) 
	  for(int m = 0; m<expansionLevel; m++) { 
	    fMatFinal[m*n*Nvsteps/downSampleV + k*Nvsteps/downSampleV + jcount] 
	      = fMatFinalOld[m*n*Nvsteps+k*Nvsteps+j]; 
	    //cout << j << endl; 
	  }
	jcount++; 
	//std::cout << "jcount is: " << jcount << "Nvsteps is: " << Nvsteps << std::endl; 
      }
      std::cout << "Passed assigning fMatFinalOld to fMatFinal" << std::endl; 
      time(&end);
      std::cout << "Passed the time step" << std::endl; 
      fwrite(fMatFinal, Nvsteps/downSampleV*n*expansionLevel*sizeof(double), 1, ptr_fp);
      outputFile << i << " tSteps in : " << difftime(end, begin) << " seconds" << std::endl;
    }
  }
}
