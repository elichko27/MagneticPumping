// printMat.cpp
#include <cstdlib>
#include <iostream>
#include "printIntMat.h"

void printIntMat(int * matrix, int n, int x, int y) { 

  for(int i = 0; i < n*x*y; i++) { 
    if (i%(x*y) == 0 && i!=0) std::cout << std::endl; 
    if (i%x == 0 && i != 0) std::cout << std::endl;
    std::cout << matrix[i] << " ";
    
  }

}
