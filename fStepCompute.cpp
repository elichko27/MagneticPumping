//fStepCompute.cpp
#include <iostream>
#include <cstring>
#include <cmath>
#include <set>
#include "R0.h"
#include "nuF.h"
#include "gradOpts.h"
#include "Ais.h"
#include "fStepCompute.h"

void fStepCompute( double * fMatTempNew, double * fMatTempOld, double delT, double delV, 
		   double c1, double nu, double omega, double vthe, double delR, 
		   double deln, int Nvsteps, int n, int expansionLevel, double kpar, double t, 
		   char * scattType, char * gradientOption ) {

  //Error catching 
  if (expansionLevel%2 != 1) {
    std::cout << "Invalid Expansion Level - expansion level must be odd!\n";
    //return 1;
  }

  double nuFac[Nvsteps];
  double df0dv, df0cdv, df0sdv, df02cdv, df02sdv, df1dv, df1cdv, df1sdv, df12cdv, df12sdv, 
    df2dv, df2cdv, df2sdv, df22cdv, df22sdv, df3dv, df3cdv, df3sdv, df32cdv, df32sdv, f0, 
    f0c, f0s, f02c, f02s, f1, f1c, f1s, f12c, f12s, f2, f2c, f2s, f22c, f22s, f3, f3c, f3s, 
    f32c, f32s, dfnm2dv, dfnm2cdv, dfnm2sdv, dfnm22cdv, dfnm22sdv, dfndv, dfncdv, dfnsdv, 
    dfn2cdv, dfn2sdv, fnm2, fnm2c, fnm2s, fnm22c, fnm22s, fnm1c, fnm1s, fnm12c, fnm12s, fn, 
    fnc, fns, fn2c, fn2s, fnp1c, fnp1s, fnp12c, fnp12s, fnp2, fnp2c, fnp2s, fnp22c, fnp22s, 
    dfnp2dv, dfnp2cdv, dfnp2sdv, dfnp22cdv, dfnp22sdv;  
  double jd, fac1, fac2, fac3; 

  double dR0 = R0(omega, delR, t);
  double dn0 = R0(omega, deln, t);

  nuF(nuFac, nu, 0, vthe, delV, Nvsteps, scattType); 

  if (expansionLevel == 5) {

    for (int i=0; i<Nvsteps; i++) 
      for(int j=0; j<n; j++) 
	for(int k=0; k<expansionLevel; k++) {

	  if (j == 0) { 
	    if (k == 0) { 

	      f0 = fMatTempOld[0*n*Nvsteps + 0*Nvsteps + i];
	      f2c = fMatTempOld[1*n*Nvsteps + 2*Nvsteps + i];

	      df2cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 2*Nvsteps);
	      df0cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 0*Nvsteps);

	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i] 
		- delT*(c1*f0 + nuFac[i]*f0 + dR0/10.*double(i)*delV*df2cdv
		      + dR0/2.*21./10.*f2c
		      + dn0/2.*double(i)*delV*df0cdv);

	    } else if (k == 1) { 

	      df2dv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + 2*Nvsteps);
	      df22cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + 2*Nvsteps); 
	      df0dv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + 0*Nvsteps); 
	      df02cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + 0*Nvsteps);

	      f2 = fMatTempOld[0*n*Nvsteps + 2*Nvsteps + i]; 
	      f22c = fMatTempOld[3*n*Nvsteps + 2*Nvsteps + i]; 
	      f1s = fMatTempOld[2*n*Nvsteps + 1*Nvsteps + i]; 
	      f0c = fMatTempOld[1*n*Nvsteps + 0*Nvsteps + i]; 

	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i] 
		- delT*(c1*f0c + nuFac[i]*f0c + dR0/5.*double(i)*delV*(df2dv + 1./2.*df22cdv)
		      + dR0*21./10.*(f2 + 1/2*f22c) + kpar/3.*double(i)*delV*f1s
		      + dn0*double(i)*delV*(df0dv + 1./2.*df02cdv));

	    } else if (k == 2) {

	      f0s = fMatTempOld[2*n*Nvsteps + 0*Nvsteps + i]; 
	      f22s = fMatTempOld[4*n*Nvsteps + 2*Nvsteps + i]; 
	      f1c = fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i]; 

	      df22sdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + 2*Nvsteps); 
	      df02sdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + 0*Nvsteps); 

	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i]
		- delT*(c1*f0s + nuFac[i]*f0s + dR0/10.*double(i)*delV*df22sdv + dR0/2.*21./10.*f22s
			- kpar/3.*double(i)*delV*f1c + dn0/2.*double(i)*delV*df02sdv);  

	    } else if (k == 3) { 
	      
	      df2cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps +2*Nvsteps);  
	      df0cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps +0*Nvsteps); 

	      f2c = fMatTempOld[1*n*Nvsteps + 2*Nvsteps + i]; 
	      f12s = fMatTempOld[4*n*Nvsteps + 1*Nvsteps + i];
	      f02c = fMatTempOld[3*n*Nvsteps + 0*Nvsteps + i]; 
	      
	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i] 
		- delT*(c1*f02c + nuFac[i]*f02c + dR0/10.*double(i)*delV*df2cdv + dR0/2.*21./10.*f2c 
			+ 2.*kpar/3.*double(i)*delV*f12s + dn0/2.*double(i)*delV*df0cdv);

	    } else {

	      f02s = fMatTempOld[4*n*Nvsteps + 0*Nvsteps + i]; 
	      f2s = fMatTempOld[2*n*Nvsteps + 2*Nvsteps + i]; 
	      f12c = fMatTempOld[3*n*Nvsteps + 1*Nvsteps + i]; 

	      df2sdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps +2*Nvsteps); 
	      df0sdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps +0*Nvsteps); 
	      
	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i] 
		- delT*(c1*f02s + nuFac[i]*f02s + dR0/10.*double(i)*delV*df2sdv + dR0/2.*21./10.*f2s
			- 2.*kpar/3.*double(i)*delV*f12c + dn0/2.*double(i)*delV*df0sdv); 
	    }
	  } else if (j == 1) { 
	    if (k == 0) { 

	      f1 = fMatTempOld[0*n*Nvsteps + 1*Nvsteps + i];
	      f1c = fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i]; 
	      f3c = fMatTempOld[1*n*Nvsteps + 3*Nvsteps + i];

	      df1cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 1*Nvsteps);
 
	      if (n > 3) {
		df3cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 3*Nvsteps);  
	      } else { 
		df3cdv = 0.0; 
	      }

	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i]
		- delT*(c1*f1 + nuFac[i]*f1 + dR0/2.*double(i)*delV*(2./5.*df1cdv + 9./35.*df3cdv) 
			+ dR0/2*(3./5.*f1c + 177./70.*f3c) + dn0/2.*double(i)*delV*df1cdv);

	    } else if (k == 1) { 
	      f0s = fMatTempOld[2*n*Nvsteps + 0*Nvsteps + i];
	      f1 = fMatTempOld[0*n*Nvsteps + 1*Nvsteps + i]; 
	      f1c = fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i]; 
	      f12c = fMatTempOld[3*n*Nvsteps + 1*Nvsteps + i];
	      f2s = fMatTempOld[2*n*Nvsteps + 2*Nvsteps + i];

	      df1dv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + 1*Nvsteps);  
	      df12cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + 1*Nvsteps);  
	      
	      if (n > 3) { 
		f3  = fMatTempOld[0*n*Nvsteps + 3*Nvsteps + i]; 
		f32c  = fMatTempOld[3*n*Nvsteps + 3*Nvsteps + i];

		df3dv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + 3*Nvsteps);
		df32cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + 3*Nvsteps);
	      } else {
		f3  = 0.0; 
		f32c  = 0.0;

		df3dv = 0.0;
		df32cdv = 0.0;
	      }

	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i]
		- delT*(c1*f1c + nuFac[i]*f1c + dR0*double(i)*delV*((2./5.*df1dv + 9./35.*df3dv) 
								    + 1./2.*(2./5.*df12cdv + 9./35.*df32cdv))
			+ dR0*(3./5.*(f1 + 1./2.*f12c) + 177./70.*(f3 + 1./2.*f32c)) 
			+ kpar*double(i)*delV*(f0s + 2./5.*f2s)
			+ dn0*double(i)*delV*(df1dv + 1./2.*df12cdv)); 

	    } else if (k == 2) {

	      f0c = fMatTempOld[1*n*Nvsteps + 0*Nvsteps + i];
	      f1s = fMatTempOld[2*n*Nvsteps + 1*Nvsteps + i];   
	      f12s = fMatTempOld[4*n*Nvsteps + 1*Nvsteps + i]; 
	      f2c = fMatTempOld[1*n*Nvsteps + 2*Nvsteps + i]; 

	      df12sdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + 1*Nvsteps); 

	      if (n > 3) {
		f32s = fMatTempOld[4*n*Nvsteps + 3*Nvsteps + i];
		df32sdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + 3*Nvsteps);
	      } else {
		f32s = 0.0;
		df32sdv = 0.0;
	      }

	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i]
		- delT*(c1*f1s + nuFac[i]*f1s + dR0/2.*double(i)*delV*(2./5.*df12sdv + 9./35.*df32sdv)
			+ dR0/2*(3./5.*f12s + 177./70.*f32s) - kpar*double(i)*delV*(f0c + 2./5.*f2c)
			+ dn0/2.*double(i)*delV*df12sdv); 

	    } else if (k == 3) { 
	      
	      f02s = fMatTempOld[4*n*Nvsteps + 0*Nvsteps + i];
	      f1c = fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i]; 
	      f12c = fMatTempOld[3*n*Nvsteps + 1*Nvsteps + i];
	      f22s = fMatTempOld[4*n*Nvsteps + 2*Nvsteps + i]; 

	      df1cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 1*Nvsteps); 
	      
	      if (n > 3) {
		f3c = fMatTempOld[1*n*Nvsteps + 3*Nvsteps + i]; 
		df3cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 3*Nvsteps);
	      } else {
		f3c = 0.0; 
		df3cdv = 0.0;
	      }
	      
	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i]
		- delT*(c1*f12c + nuFac[i]*f12c + dR0/2.*double(i)*delV*(2./5.*df1cdv + 9./35.*df3cdv)
			+ dR0/2.*(3./5.*f1c + 177./70.*f3c) + 2.*kpar*double(i)*delV*(f02s + 2./5.*f22s)
			+ dn0/2.*double(i)*delV*df1cdv);

	    } else {

	      f02c = fMatTempOld[3*n*Nvsteps + 0*Nvsteps + i];
	      f1s = fMatTempOld[2*n*Nvsteps + 1*Nvsteps + i];
	      f12s = fMatTempOld[4*n*Nvsteps + 1*Nvsteps + i]; 
	      f22c = fMatTempOld[3*n*Nvsteps + 2*Nvsteps + i];

	      df1sdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + 1*Nvsteps); 

	      if (n > 3) {
		f3s = fMatTempOld[2*n*Nvsteps + 3*Nvsteps + i];
		df3sdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + 3*Nvsteps);
	      } else {
		f3s = 0.0;
		df3sdv = 0.0;
	      }

	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i]
		- delT*(c1*f12s + nuFac[i]*f12s + dR0/2.*double(i)*delV*(2./5.*df1sdv + 9./35.*df3sdv)
			+ dR0/2.*(3./5.*f1s + 177./70.*f3s) - 2.*kpar*double(i)*delV*(f02c + 2./5.*f22c)
			+ dn0/2.*double(i)*delV*df1sdv); 
	    }
	  } else {
	    
	    jd = double(j); 
	    fac1 = 3./2.*jd*(jd-1)/(2.*jd-3.)/(2.*jd-1.);
	    fac2 = 3./2.*(pow(jd+1.,2.)/(2.*jd+1.)/(2.*jd+3.) 
			  + pow(jd,2.)/(2*jd+1)/(2.*jd-1.)) - 1./2.;
	    fac3 = 3./2.*(jd+1.)*(jd+2.)/(2.*jd+3.)/(2.*jd+5.);

	    if (k == 0) { 

	      fn = fMatTempOld[0*n*Nvsteps + j*Nvsteps + i]; 
	      fnm2c = fMatTempOld[1*n*Nvsteps + (j-2)*Nvsteps + i]; 
	      fnc = fMatTempOld[1*n*Nvsteps + j*Nvsteps + i]; 

	      dfnm2cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j-2)*Nvsteps); 
	      dfncdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j)*Nvsteps); 

	      if (n > (j+2)) { 
		fnp2c = fMatTempOld[1*n*Nvsteps + (j+2)*Nvsteps + i];
		dfnp2cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j+2)*Nvsteps); 
	      } else {
		fnp2c = 0.0;
		dfnp2cdv = 0.0;
	      }
	      
	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i]
		- delT*(c1*fn + nuFac[i]*fn + dR0/2.*double(i)*delV*(fac1*dfnm2cdv + fac2*dfncdv + fac3*dfnp2cdv)
			+ dR0/2.*(Ais(j,0)*fnm2c + Ais(j,1)*fnc + Ais(j,2)*fnp2c)
			+ dn0/2.*double(i)*delV*dfncdv); 
 
	    } else if (k == 1) { 

	      fnm2 = fMatTempOld[0*n*Nvsteps + (j-2)*Nvsteps + i]; 
	      fnm22c = fMatTempOld[3*n*Nvsteps + (j-2)*Nvsteps + i]; 
	      fnm1s = fMatTempOld[2*n*Nvsteps + (j-1)*Nvsteps + i]; 
	      fn = fMatTempOld[0*n*Nvsteps + j*Nvsteps + i]; 
	      fnc = fMatTempOld[1*n*Nvsteps + j*Nvsteps + i]; 
	      fn2c = fMatTempOld[3*n*Nvsteps + j*Nvsteps + i];  

	      dfnm2dv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + (j-2)*Nvsteps);
	      dfnm22cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + (j-2)*Nvsteps); 
	      dfndv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + j*Nvsteps); 
	      dfn2cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + j*Nvsteps); 

	      if (n > (j+1)) {
		fnp1s = fMatTempOld[2*n*Nvsteps + (j+1)*Nvsteps + i]; 
	      } else {
		fnp1s = 0.0;
	      }

	      if (n > (j+2)) {
		fnp2 = fMatTempOld[0*n*Nvsteps + (j+2)*Nvsteps + i];
		fnp22c = fMatTempOld[3*n*Nvsteps + (j+2)*Nvsteps + i];

		dfnp2dv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + (j+2)*Nvsteps);
		dfnp22cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + (j+2)*Nvsteps); 
	      } else {

		fnp2 = 0.0;
		fnp22c = 0.0;

		dfnp2dv = 0.0;
		dfnp22cdv= 0.0;
	      }


	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i]
		- delT*(c1*fnc + nuFac[i]*fnc + dR0*double(i)*delV*((fac1*dfnm2dv + fac2*dfndv + fac3*dfnp2dv) 
								    + 1./2.*(fac1*dfnm22cdv + fac2*dfn2cdv + fac3*dfnp22cdv))
			+ dR0*(Ais(j,0)*(fnm2 + 1./2.*fnm22c) + Ais(j,1)*(fn + 1/2*fn2c) + Ais(j,2)*(fnp2 + 1./2.*fnp22c))
			+ kpar*double(i)*delV*(jd/(2.*jd-1.)*fnm1s + (jd+1.)/(2.*jd+3.)*fnp1s)
			+ dn0*double(i)*delV*(dfndv + 1./2.*dfn2cdv));

	    } else if (k == 2) {

	      fnm22s = fMatTempOld[4*n*Nvsteps + (j-2)*Nvsteps + i]; 
	      fnm1c = fMatTempOld[1*n*Nvsteps + (j-1)*Nvsteps + i]; 
	      fns = fMatTempOld[2*n*Nvsteps + j*Nvsteps + i]; 
	      fn2s = fMatTempOld[4*n*Nvsteps + j*Nvsteps + i];

	      dfnm22sdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + (j-2)*Nvsteps);
	      dfn2sdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + j*Nvsteps); 

	      if (n > (j+1)) {
		fnp1c = fMatTempOld[1*n*Nvsteps + (j+1)*Nvsteps + i]; 
	      } else {
		fnp1c = 0.0;
	      }

	      if (n > (j+2)) {
		fnp22s = fMatTempOld[4*n*Nvsteps + (j+2)*Nvsteps + i]; 
		dfnp22sdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + (j+2)*Nvsteps); 
	      } else {
		fnp22s = 0.0; 
		dfnp22sdv = 0.0; 
	      }


	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i]
		- delT*(c1*fns + nuFac[i]*fns
			+ dR0/2.*double(i)*delV*(fac1*dfnm22sdv + fac2*dfn2sdv + fac3*dfnp22sdv)
			+ dR0/2.*(Ais(j,0)*fnm22s + Ais(j,1)*fn2s + Ais(j,2)*fnp22s)
			- kpar*double(i)*delV*(jd/(2.*jd-1.)*fnm1c + (jd+1)/(2.*jd+3.)*fnp1c)
			+ dn0/2.*double(i)*delV*dfn2sdv); 

	    } else if (k == 3) { 
 
	      fnm2c = fMatTempOld[1*n*Nvsteps + (j-2)*Nvsteps + i]; 
	      fnm12s = fMatTempOld[4*n*Nvsteps + (j-1)*Nvsteps + i]; 
	      fnc = fMatTempOld[1*n*Nvsteps + j*Nvsteps + i]; 
	      fn2c = fMatTempOld[3*n*Nvsteps + j*Nvsteps + i]; 

	      dfnm2cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j-2)*Nvsteps);  
	      dfncdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + j*Nvsteps);

	      if (n > (j+1)) {
		fnp12s = fMatTempOld[4*n*Nvsteps + (j+1)*Nvsteps + i]; 
	      } else {
		fnp12s = 0.0; 
	      }

	      if (n > (j+2)) {
		fnp2c = fMatTempOld[1*n*Nvsteps + (j+2)*Nvsteps + i];
		dfnp2cdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j+2)*Nvsteps);
	      } else {
		fnp2c = 0.0;
		dfnp2cdv = 0.0;
	      }

	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i]
		- delT*(c1*fn2c + nuFac[i]*fn2c 
			+ dR0/2.*double(i)*delV*(fac1*dfnm2cdv + fac2*dfncdv + fac3*dfnp2cdv)
			+ dR0/2.*(Ais(j,0)*fnm2c + Ais(j,1)*fnc + Ais(j,2)*fnp2c)
			+ 2.*kpar*double(i)*delV*(jd/(2.*jd-1)*fnm12s + (jd+1.)/(2.*jd+3.)*fnp12s)
			+ dn0/2.*double(i)*delV*dfncdv); 

	    } else {

	      fnm2s = fMatTempOld[2*n*Nvsteps + (j-2)*Nvsteps + i]; 
	      fnm12c = fMatTempOld[3*n*Nvsteps + (j-1)*Nvsteps + i]; 
	      fns = fMatTempOld[2*n*Nvsteps + j*Nvsteps + i];
	      fn2s = fMatTempOld[4*n*Nvsteps + j*Nvsteps + i]; 

	      dfnm2sdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + (j-2)*Nvsteps); 
	      dfnsdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + j*Nvsteps); 

	      if (n > (j+1)) {
		fnp12c = fMatTempOld[3*n*Nvsteps + (j+1)*Nvsteps + i]; 
	      } else {
		fnp12c = 0.0; 
	      }

	      if (n > (j+2)) {
		fnp2s = fMatTempOld[2*n*Nvsteps + (j+2)*Nvsteps + i]; 
		dfnp2sdv = gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + (j+2)*Nvsteps); 
	      } else {
		fnp2s = 0.0; 
		dfnp2sdv = 0.0;
	      }

	      fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i]
		- delT*(c1*fn2s + nuFac[i]*fn2s 
			+ dR0/2.*double(i)*delV*(fac1*dfnm2sdv + fac2*dfnsdv + fac3*dfnp2sdv)
			+ dR0/2.*(Ais(j,0)*fnm2s + Ais(j,1)*fns + Ais(j,2)*fnp2s)
			- 2*kpar*double(i)*delV*(jd/(2.*jd-1.)*fnm12c + (jd+1.)/(2.*jd+3.)*fnp12c)
			+ dn0/2.*double(i)*delV*dfnsdv);

	      /*if (n >= (j+4)) { 
		fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempNew[k*n*Nvsteps + j*Nvsteps + i]
		  + delT*dR0/2*(Ais(j,3)*fMatTempOld[k*n*Nvsteps + j*Nvsteps + i]);
		  }*/
	    } // if (k == 0)

	    // Adding on the infinite terms
	    /*if (j%2 == 0 && j > 2) {
	      fMatTempNew[0*n*Nvsteps + 0*Nvsteps + i] =  fMatTempNew[0*n*Nvsteps + 0*Nvsteps + i] 
		+ delT*dR0/2.*(3./2.)*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i]; 
	      fMatTempNew[1*n*Nvsteps + 0*Nvsteps + i] =  fMatTempNew[1*n*Nvsteps + 0*Nvsteps + i] 
		+ delT*dR0*(3./2.)*(fMatTempOld[0*n*Nvsteps + j*Nvsteps + i] + 1./2.*fMatTempOld[3*n*Nvsteps + j*Nvsteps + i]);
	      fMatTempNew[2*n*Nvsteps + 0*Nvsteps + i] =  fMatTempNew[2*n*Nvsteps + 0*Nvsteps + i] 
		+ delT*dR0/2.*(3./2.)*fMatTempOld[4*n*Nvsteps + j*Nvsteps + i];
	      fMatTempNew[3*n*Nvsteps + 0*Nvsteps + i] =  fMatTempNew[3*n*Nvsteps + 0*Nvsteps + i] 
		+ delT*dR0/2.*(3./2.)*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i];
	      fMatTempNew[4*n*Nvsteps + 0*Nvsteps + i] =  fMatTempNew[4*n*Nvsteps + 0*Nvsteps + i] 
		+ delT*dR0/2.*(3./2.)*fMatTempOld[2*n*Nvsteps + j*Nvsteps + i];


	      for (int m = 2; m<n; m+=2) { //looping through the equations that are to be modified
		if (j >= m+4) {
		  fMatTempNew[0*n*Nvsteps + m*Nvsteps + i] =  fMatTempNew[0*n*Nvsteps + m*Nvsteps + i] 
		    + delT*dR0/2.*(Ais(j,3))*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i]; 
		  fMatTempNew[1*n*Nvsteps + m*Nvsteps + i] =  fMatTempNew[1*n*Nvsteps + m*Nvsteps + i] 
		    + delT*dR0*(Ais(j,3))*(fMatTempOld[0*n*Nvsteps + j*Nvsteps + i] 
					   + 1./2.*fMatTempOld[3*n*Nvsteps + j*Nvsteps + i]);
		  fMatTempNew[2*n*Nvsteps + m*Nvsteps + i] =  fMatTempNew[2*n*Nvsteps + m*Nvsteps + i] 
		    + delT*dR0/2.*(Ais(j,3))*fMatTempOld[4*n*Nvsteps + j*Nvsteps + i];
		  fMatTempNew[3*n*Nvsteps + m*Nvsteps + i] =  fMatTempNew[3*n*Nvsteps + m*Nvsteps + i] 
		    + delT*dR0/2.*(Ais(j,3))*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i];
		  fMatTempNew[4*n*Nvsteps + m*Nvsteps + i] =  fMatTempNew[4*n*Nvsteps + m*Nvsteps + i] 
		    + delT*dR0/2.*(Ais(j,3))*fMatTempOld[2*n*Nvsteps + j*Nvsteps + i];
		}
	      }

	    } else if (j%2 == 1 && j > 3) { 
	      fMatTempNew[0*n*Nvsteps + 1*Nvsteps + i] =  fMatTempNew[0*n*Nvsteps + 1*Nvsteps + i] 
		+ delT*dR0/2.*(15./2.)*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i]; 
	      fMatTempNew[1*n*Nvsteps + 1*Nvsteps + i] =  fMatTempNew[1*n*Nvsteps + 1*Nvsteps + i] 
		+ delT*dR0*(15./2.)*(fMatTempOld[0*n*Nvsteps + j*Nvsteps + i] + 1./2.*fMatTempOld[3*n*Nvsteps + j*Nvsteps + i]);
	      fMatTempNew[2*n*Nvsteps + 1*Nvsteps + i] =  fMatTempNew[2*n*Nvsteps + 1*Nvsteps + i] 
		+ delT*dR0/2.*(15./2.)*fMatTempOld[4*n*Nvsteps + j*Nvsteps + i];
	      fMatTempNew[3*n*Nvsteps + 1*Nvsteps + i] =  fMatTempNew[3*n*Nvsteps + 1*Nvsteps + i] 
		+ delT*dR0/2.*(15./2.)*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i];
	      fMatTempNew[4*n*Nvsteps + 1*Nvsteps + i] =  fMatTempNew[4*n*Nvsteps + 1*Nvsteps + i] 
		+ delT*dR0/2.*(15./2.)*fMatTempOld[2*n*Nvsteps + j*Nvsteps + i];
	      
	      for (int m = 3; m<n; m+=2) { //looping through the equations that are to be modified
		if (j >= m+4) {
		  fMatTempNew[0*n*Nvsteps + m*Nvsteps + i] =  fMatTempNew[0*n*Nvsteps + m*Nvsteps + i] 
		    + delT*dR0/2.*(Ais(j,3))*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i]; 
		  fMatTempNew[1*n*Nvsteps + m*Nvsteps + i] =  fMatTempNew[1*n*Nvsteps + m*Nvsteps + i] 
		    + delT*dR0*(Ais(j,3))*(fMatTempOld[0*n*Nvsteps + j*Nvsteps + i] 
					   + 1./2.*fMatTempOld[3*n*Nvsteps + j*Nvsteps + i]);
		  fMatTempNew[2*n*Nvsteps + m*Nvsteps + i] =  fMatTempNew[2*n*Nvsteps + m*Nvsteps + i] 
		    + delT*dR0/2.*(Ais(j,3))*fMatTempOld[4*n*Nvsteps + j*Nvsteps + i];
		  fMatTempNew[3*n*Nvsteps + m*Nvsteps + i] =  fMatTempNew[3*n*Nvsteps + m*Nvsteps + i] 
		    + delT*dR0/2.*(Ais(j,3))*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i];
		  fMatTempNew[4*n*Nvsteps + m*Nvsteps + i] =  fMatTempNew[4*n*Nvsteps + m*Nvsteps + i] 
		    + delT*dR0/2.*(Ais(j,3))*fMatTempOld[2*n*Nvsteps + j*Nvsteps + i];
		}
	      } // for (int m = 3; m<n; m+=2)
	    } // if (j%2 == 0 && j > 2) */

	  } //if (j == 0)
	} // for (int k=0; k<expansionLevel; k++)
  } // if (expansionLevel == 5)
} // void fStepCompute
