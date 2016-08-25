//R0.cpp 
#include <iostream>
#include <cstring>
#include <cmath>
#include <set>
#include "R0.h"

double R0( double omega, double delR, double t) {
  return omega*delR*sin(omega*t); 
}
