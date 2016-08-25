//gradOpts.cpp
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <set>
#include "gradOpts.h"

double gradOpts(double * vec1, double delT, int j, int jMin, int jMax, int rest) {

  /*if ((abs((vec1[rest + j+1] - vec1[rest + j])/delT) > 10000. && j == jMin) 
      || (abs((vec1[rest + j] - vec1[rest + j-1])/delT) > 100. && j == jMax) 
      || (abs((vec1[rest + j+1] - vec1[rest + j-1])/delT) > 100.)){
    std::cout << std::endl << "delT = " << delT << std::endl 
	      << "vec1[rest + j+1] = " << vec1[rest + j+1] << std::endl 
	      << "vec1[rest + j] = " << vec1[rest + j] << std::endl
	      << "vec1[rest + j-1] = " << vec1[rest + j-1] << std::endl
	      << "rest + j + 1 = " << rest + j + 1 << std::endl
	      << "rest + j = " << rest + j << std::endl
	      << "rest + j - 1 = " << rest + j - 1 << std::endl;
	exit(1); 
	}*/

  if (j == jMin) {
    return (vec1[rest + j+1] - vec1[rest + j])/delT; 
  } else if (j == (jMax-1)) {
    return (vec1[rest + j] - vec1[rest + j-1])/delT; 
  } else {
    return (vec1[rest + j+1] - vec1[rest + j-1])/2./delT; 
  }

}
