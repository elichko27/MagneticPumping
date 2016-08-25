//simParams.cpp
#include <iostream>
#include <cstring>
#include <cmath>
#include <set>
#include "simParams.h"

void simParams(char * rChoice, char * regionChoice, double &vthe, double &memi, 
	     double &B, double &wpe, double &wce, double &wci, double &beta_e, 
	     double &vA, double &c1, double &vMin, double &vMax, double &omega, 
	     double &delB, double &deln, double &delR, double &kpar, int &Lx, int &l0) 
{ 
  // Overall Constants
  vthe = 0.0707; 
  memi = 1/25.0; 
  B = 0.5; 
  wpe = 1.0; 
  wce = B*wpe; 
  wci = wce*memi; 
  beta_e = 2*pow(wpe,2)*pow(wce,2)*pow(vthe,2); 
  vA = sqrt(memi/beta_e)*vthe; 
  Lx = 320; 
  c1 = 0.0; 
  vMin = 0.0; 
  vMax = 1.4*10.0; //1.4*20.0;

  double pi = 3.1415926535897;

  if (strcmp(regionChoice, "Region1") == 0) {
      double xMin = 403./534.;
      double xMax = 449./534.;
      int Nxsteps = 449 - 403;
      // xVec = linspace(xMin, xMax, Nxsteps);
      //cosFac = mean(cos(xVec*2*pi));
      //sinFac = mean(sin(xVec*2*pi));
      //cos2Fac = mean(cos(2*xVec*2*pi));
      //sin2Fac = mean(sin(2*xVec*2*pi));
    
      double x = 0.75; 
      //spatialFacVec = [1 cos(x*2*pi) sin(x*2*pi) cos(x*2*pi) sin(x*2*pi)];
    }

  if (strcmp(rChoice, "59703") == 0) {
    l0 = 0;
    omega = 0.095786937159141*wci;
    
    delB = 0.1;
    deln = 0.0; 
    delR = std::abs(2/3*deln - delB);
  } else if (strcmp(rChoice, "60703") == 0) {
    l0 = 1;
    omega = 0.136534197049846*wci;   
    
    delB = 0.1; 
    deln = 0.0; 
    delR = std::abs(2/3*deln - delB);
  } else if (strcmp(rChoice, "1DMagneticPumping") == 0) {
    memi =1/100.0; 
    l0 = 0;
    omega = 0.1*wpe;
    Lx = 4;
    
    delB = 0.25;
    deln = 3.*delB; 
    delR = std::abs(2./3.*deln - delB);
  }
  kpar = double(l0)*2.0*pi/double(Lx);

  
}
