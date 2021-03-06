//fStepCompute2.cpp
#include <iostream>
#include <cstring>
#include <cmath>
#include <set>
#include "R0.h"
#include "nuF.h"
#include "gradOpts.h"
#include "Ais.h"
#include "fStepCompute2.h"

void fStepCompute2( double * fMatTempNew, double * fMatTempOld, double delT, double delV, 
		   double c1, double nu, double omega, double vthe, double delR, 
		   double deln, int Nvsteps, int n, int expansionLevel, double kpar, double t, 
		   char * scattType, char * gradientOption ) {

  //Error catching 
  if (expansionLevel%2 != 1) {
    std::cout << "Invalid Expansion Level - expansion level must be odd!\n";
    //return 1;
  }

  double nuFac[Nvsteps];
  double jd, fac1, fac2, fac3; 
  int j; 

  double dR0 = R0(omega, delR, t);
  double dn0 = R0(omega, deln, t);

  nuF(nuFac, nu, 0, vthe, delV, Nvsteps, scattType); 

  if (expansionLevel == 5) {

    for (int i=0; i<Nvsteps; i++) {

      fMatTempNew[0*n*Nvsteps + 0*Nvsteps + i] = fMatTempOld[0*n*Nvsteps + 0*Nvsteps + i] 
	- delT*(c1*fMatTempOld[0*n*Nvsteps + 0*Nvsteps + i] 
		+ nuFac[i]*fMatTempOld[0*n*Nvsteps + 0*Nvsteps + i] 
		+ dR0/10.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 2*Nvsteps)
		+ dR0/2.*21./10.*fMatTempOld[1*n*Nvsteps + 2*Nvsteps + i]
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 0*Nvsteps));
      
      fMatTempNew[1*n*Nvsteps + 0*Nvsteps + i] = fMatTempOld[1*n*Nvsteps + 0*Nvsteps + i] 
	- delT*(c1*fMatTempOld[1*n*Nvsteps + 0*Nvsteps + i] 
		+ nuFac[i]*fMatTempOld[1*n*Nvsteps + 0*Nvsteps + i] 
		+ dR0/5.*double(i)*delV*(gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + 2*Nvsteps) 
					 + 1./2.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + 2*Nvsteps))
		+ dR0*21./10.*(fMatTempOld[0*n*Nvsteps + 2*Nvsteps + i] 
			       + 1./2.*fMatTempOld[3*n*Nvsteps + 2*Nvsteps + i]) 
		+ kpar/3.*double(i)*delV*fMatTempOld[2*n*Nvsteps + 1*Nvsteps + i]
		+ dn0*double(i)*delV*(gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + 0*Nvsteps)
				      + 1./2.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + 0*Nvsteps)));

      fMatTempNew[2*n*Nvsteps + 0*Nvsteps + i] = fMatTempOld[2*n*Nvsteps + 0*Nvsteps + i]
	- delT*(c1*fMatTempOld[2*n*Nvsteps + 0*Nvsteps + i] 
		+ nuFac[i]*fMatTempOld[2*n*Nvsteps + 0*Nvsteps + i] 
		+ dR0/10.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + 2*Nvsteps)
		+ dR0/2.*21./10.*fMatTempOld[4*n*Nvsteps + 2*Nvsteps + i]
		- kpar/3.*double(i)*delV*fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i] 
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + 0*Nvsteps));  
      
      fMatTempNew[3*n*Nvsteps + 0*Nvsteps + i] = fMatTempOld[3*n*Nvsteps + 0*Nvsteps + i] 
	- delT*(c1*fMatTempOld[3*n*Nvsteps + 0*Nvsteps + i] 
		+ nuFac[i]*fMatTempOld[3*n*Nvsteps + 0*Nvsteps + i] 
		+ dR0/10.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 2*Nvsteps)
		+ dR0/2.*21./10.*fMatTempOld[1*n*Nvsteps + 2*Nvsteps + i] 
		+ 2.*kpar/3.*double(i)*delV*fMatTempOld[4*n*Nvsteps + 1*Nvsteps + i] 
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 0*Nvsteps));
      
      fMatTempNew[4*n*Nvsteps + 0*Nvsteps + i] = fMatTempOld[4*n*Nvsteps + 0*Nvsteps + i] 
	- delT*(c1*fMatTempOld[4*n*Nvsteps + 0*Nvsteps + i] 
		+ nuFac[i]*fMatTempOld[4*n*Nvsteps + 0*Nvsteps + i] 
		+ dR0/10.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + 2*Nvsteps) 
		+ dR0/2.*21./10.*fMatTempOld[2*n*Nvsteps + 2*Nvsteps + i]
		- 2.*kpar/3.*double(i)*delV*fMatTempOld[3*n*Nvsteps + 1*Nvsteps + i] 
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + 0*Nvsteps)); 

      fMatTempNew[0*n*Nvsteps + 1*Nvsteps + i] = fMatTempOld[0*n*Nvsteps + 1*Nvsteps + i]
	- delT*(c1*fMatTempOld[0*n*Nvsteps + 1*Nvsteps + i] 
		+ nuFac[i]*fMatTempOld[0*n*Nvsteps + 1*Nvsteps + i] 
		+ dR0/2.*double(i)*delV*2./5.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 1*Nvsteps) 
		+ dR0/2*3./5.*fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i]  
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 1*Nvsteps)); 
      
      fMatTempNew[1*n*Nvsteps + 1*Nvsteps + i] = fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i]
	- delT*(c1*fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i] 
		+ nuFac[i]*fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i] 
		+ dR0*double(i)*delV*(2./5.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + 1*Nvsteps) 
				      + 1./2.*2./5.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + 1*Nvsteps))
		+ dR0*3./5.*(fMatTempOld[0*n*Nvsteps + 1*Nvsteps + i] + 1./2.*fMatTempOld[3*n*Nvsteps + 1*Nvsteps + i]) 
		+ kpar*double(i)*delV*(fMatTempOld[2*n*Nvsteps + 0*Nvsteps + i] 
				       + 2./5.*fMatTempOld[2*n*Nvsteps + 2*Nvsteps + i])
		+ dn0*double(i)*delV*(gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + 1*Nvsteps) 
				      + 1./2.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + 1*Nvsteps))); 
      
      fMatTempNew[2*n*Nvsteps + 1*Nvsteps + i] = fMatTempOld[2*n*Nvsteps + 1*Nvsteps + i]
	- delT*(c1*fMatTempOld[2*n*Nvsteps + 1*Nvsteps + i] 
		+ nuFac[i]*fMatTempOld[2*n*Nvsteps + 1*Nvsteps + i] 
		+ dR0/2.*double(i)*delV*2./5.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + 1*Nvsteps)
		+ dR0/2*3./5.*fMatTempOld[4*n*Nvsteps + 1*Nvsteps + i] 
		- kpar*double(i)*delV*(fMatTempOld[1*n*Nvsteps + 0*Nvsteps + i] 
				       + 2./5.*fMatTempOld[1*n*Nvsteps + 2*Nvsteps + i])
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + 1*Nvsteps)); 
      
      fMatTempNew[3*n*Nvsteps + 1*Nvsteps + i] = fMatTempOld[3*n*Nvsteps + 1*Nvsteps + i]
	- delT*(c1*fMatTempOld[3*n*Nvsteps + 1*Nvsteps + i] 
		+ nuFac[i]*fMatTempOld[3*n*Nvsteps + 1*Nvsteps + i] 
		+ dR0/2.*double(i)*delV*2./5.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 1*Nvsteps)
		+ dR0/2.*3./5.*fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i] 
		+ 2.*kpar*double(i)*delV*(fMatTempOld[4*n*Nvsteps + 0*Nvsteps + i] 
					  + 2./5.*fMatTempOld[4*n*Nvsteps + 2*Nvsteps + i])
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 1*Nvsteps));

      fMatTempNew[4*n*Nvsteps + 1*Nvsteps + i] = fMatTempOld[4*n*Nvsteps + 1*Nvsteps + i]
	- delT*(c1*fMatTempOld[4*n*Nvsteps + 1*Nvsteps + i] 
		+ nuFac[i]*fMatTempOld[4*n*Nvsteps + 1*Nvsteps + i] 
		+ dR0/2.*double(i)*delV*2./5.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + 1*Nvsteps)
		+ dR0/2.*3./5.*fMatTempOld[2*n*Nvsteps + 1*Nvsteps + i] 
		- 2.*kpar*double(i)*delV*(fMatTempOld[3*n*Nvsteps + 0*Nvsteps + i] 
					  + 2./5.*fMatTempOld[3*n*Nvsteps + 2*Nvsteps + i])
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + 1*Nvsteps));
      
      if (n > 3) {
	fMatTempNew[0*n*Nvsteps + 1*Nvsteps + i] += 
	  -1.*delT*(dR0/2.*double(i)*delV*9./35.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 3*Nvsteps) 
		   + dR0/2*177./70.*fMatTempOld[1*n*Nvsteps + 3*Nvsteps + i]); 
      
	fMatTempNew[1*n*Nvsteps + 1*Nvsteps + i] += 
	  -1.*delT*(dR0*double(i)*delV*(9./35.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + 3*Nvsteps)
				       + 1./2.*9./35.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + 3*Nvsteps))
		   + dR0*177./70.*(fMatTempOld[0*n*Nvsteps + 3*Nvsteps + i] 
				   + 1./2.*fMatTempOld[3*n*Nvsteps + 3*Nvsteps + i])); 
      
	fMatTempNew[2*n*Nvsteps + 1*Nvsteps + i] += 
	  -1.*delT*(dR0/2.*double(i)*delV*9./35.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + 3*Nvsteps)
		   + dR0/2*177./70.*fMatTempOld[4*n*Nvsteps + 3*Nvsteps + i]); 
      
	fMatTempNew[3*n*Nvsteps + 1*Nvsteps + i] += 
	  -1.*delT*(dR0/2.*double(i)*delV*9./35.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 3*Nvsteps)
		   + dR0/2.*177./70.*fMatTempOld[1*n*Nvsteps + 3*Nvsteps + i]);

	fMatTempNew[4*n*Nvsteps + 1*Nvsteps + i] += 
	  -1.*delT*(dR0/2.*double(i)*delV*9./35.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + 3*Nvsteps)
		   + dR0/2.*177./70.*fMatTempOld[2*n*Nvsteps + 3*Nvsteps + i]); 
      }

      j = 2; 
      while (j < n) {
	    
	jd = double(j); 
	fac1 = 3./2.*jd*(jd-1)/(2.*jd-3.)/(2.*jd-1.);
	fac2 = 3./2.*(pow(jd+1.,2.)/(2.*jd+1.)/(2.*jd+3.) 
		      + pow(jd,2.)/(2*jd+1)/(2.*jd-1.)) - 1./2.;
	fac3 = 3./2.*(jd+1.)*(jd+2.)/(2.*jd+3.)/(2.*jd+5.);
	      
	fMatTempNew[0*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[0*n*Nvsteps + j*Nvsteps + i]
	  - delT*(c1*fMatTempOld[0*n*Nvsteps + j*Nvsteps + i] 
		  + nuFac[i]*fMatTempOld[0*n*Nvsteps + j*Nvsteps + i] 
		  + dR0/2.*double(i)*delV*(fac1*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j-2)*Nvsteps) 
					   + fac2*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j)*Nvsteps))
		  + dR0/2.*(Ais(j,0)*fMatTempOld[1*n*Nvsteps + (j-2)*Nvsteps + i] 
			    + Ais(j,1)*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i])
		  + dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j)*Nvsteps)); 

	fMatTempNew[1*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[1*n*Nvsteps + j*Nvsteps + i]
	  - delT*(c1*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i] 
		  + nuFac[i]*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i] 
		  + dR0*double(i)*delV*((fac1*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + (j-2)*Nvsteps) 
					 + fac2*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + (j)*Nvsteps)) 
					+ 1./2.*(fac1*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + (j-2)*Nvsteps) 
						 + fac2*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + (j)*Nvsteps)))
		  + dR0*(Ais(j,0)*(fMatTempOld[0*n*Nvsteps + (j-2)*Nvsteps + i] 
				   + 1./2.*fMatTempOld[3*n*Nvsteps + (j-2)*Nvsteps + i]) 
			 + Ais(j,1)*(fMatTempOld[0*n*Nvsteps + j*Nvsteps + i] 
				     + 1/2*fMatTempOld[3*n*Nvsteps + j*Nvsteps + i]))
		  + kpar*double(i)*delV*jd/(2.*jd-1.)*fMatTempOld[2*n*Nvsteps + (j-1)*Nvsteps + i]
		  + dn0*double(i)*delV*(gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + (j)*Nvsteps) 
					+ 1./2.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + (j)*Nvsteps)));

	fMatTempNew[2*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[2*n*Nvsteps + j*Nvsteps + i]
	  - delT*(c1*fMatTempOld[2*n*Nvsteps + j*Nvsteps + i] 
		  + nuFac[i]*fMatTempOld[2*n*Nvsteps + j*Nvsteps + i]
		  + dR0/2.*double(i)*delV*(fac1*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + (j-2)*Nvsteps) 
					   + fac2*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + (j)*Nvsteps))
		  + dR0/2.*(Ais(j,0)*fMatTempOld[4*n*Nvsteps + (j-2)*Nvsteps + i] 
			    + Ais(j,1)*fMatTempOld[4*n*Nvsteps + j*Nvsteps + i])
		  - kpar*double(i)*delV*jd/(2.*jd-1.)*fMatTempOld[1*n*Nvsteps + (j-1)*Nvsteps + i]
		  + dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + (j)*Nvsteps)); 
	
	fMatTempNew[3*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[3*n*Nvsteps + j*Nvsteps + i]
	  - delT*(c1*fMatTempOld[3*n*Nvsteps + j*Nvsteps + i] 
		  + nuFac[i]*fMatTempOld[3*n*Nvsteps + j*Nvsteps + i] 
		  + dR0/2.*double(i)*delV*(fac1*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j-2)*Nvsteps) 
					   + fac2*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j)*Nvsteps))
		  + dR0/2.*(Ais(j,0)*fMatTempOld[1*n*Nvsteps + (j-2)*Nvsteps + i] 
			    + Ais(j,1)*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i])
		  + 2.*kpar*double(i)*delV*jd/(2.*jd-1)*fMatTempOld[4*n*Nvsteps + (j-1)*Nvsteps + i]
		  + dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j)*Nvsteps)); 

	fMatTempNew[4*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[4*n*Nvsteps + j*Nvsteps + i]
	  - delT*(c1*fMatTempOld[4*n*Nvsteps + j*Nvsteps + i] 
		  + nuFac[i]*fMatTempOld[4*n*Nvsteps + j*Nvsteps + i] 
		  + dR0/2.*double(i)*delV*(fac1*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + (j-2)*Nvsteps) 
					   + fac2*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + (j)*Nvsteps))
		  + dR0/2.*(Ais(j,0)*fMatTempOld[2*n*Nvsteps + (j-2)*Nvsteps + i] 
			    + Ais(j,1)*fMatTempOld[2*n*Nvsteps + j*Nvsteps + i])
		  - 2*kpar*double(i)*delV*jd/(2.*jd-1.)*fMatTempOld[3*n*Nvsteps + (j-1)*Nvsteps + i]
		  + dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + (j)*Nvsteps));

	if (n > (j+1)) {
	  fMatTempNew[1*n*Nvsteps + j*Nvsteps + i] += 
	    -1.*delT*(kpar*double(i)*delV*(jd+1.)/(2.*jd+3.)*fMatTempOld[2*n*Nvsteps + (j+1)*Nvsteps + i]);

	  fMatTempNew[2*n*Nvsteps + j*Nvsteps + i] += 
	    -1.*delT*(-1.*kpar*double(i)*delV*(jd+1)/(2.*jd+3.)*fMatTempOld[1*n*Nvsteps + (j+1)*Nvsteps + i]); 
	
	  fMatTempNew[3*n*Nvsteps + j*Nvsteps + i] += 
	    -1.*delT*(2.*kpar*double(i)*delV*(jd+1.)/(2.*jd+3.)*fMatTempOld[4*n*Nvsteps + (j+1)*Nvsteps + i]); 

	  fMatTempNew[4*n*Nvsteps + j*Nvsteps + i] += 
	    -1.*delT*(-2.*kpar*double(i)*delV*(jd+1.)/(2.*jd+3.)*fMatTempOld[3*n*Nvsteps + (j+1)*Nvsteps + i]);
	}

	if (n > (j+2)) {
	  fMatTempNew[0*n*Nvsteps + j*Nvsteps + i] += 
	  -1.*delT*(dR0/2.*double(i)*delV*fac3*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j+2)*Nvsteps)
		  + dR0/2.*Ais(j,2)*fMatTempOld[1*n*Nvsteps + (j+2)*Nvsteps + i]); 

	fMatTempNew[1*n*Nvsteps + j*Nvsteps + i] += 
	  -1.*delT*(dR0*double(i)*delV*fac3*(gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + (j+2)*Nvsteps) 
					+ 1./2.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + (j+2)*Nvsteps))
		  + dR0*Ais(j,2)*(fMatTempOld[0*n*Nvsteps + (j+2)*Nvsteps + i] 
				  + 1./2.*fMatTempOld[3*n*Nvsteps + (j+2)*Nvsteps + i]));

	fMatTempNew[2*n*Nvsteps + j*Nvsteps + i] += 
	  -1.*delT*(dR0/2.*double(i)*delV*fac3*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + (j+2)*Nvsteps)
		  + dR0/2.*Ais(j,2)*fMatTempOld[4*n*Nvsteps + (j+2)*Nvsteps + i]); 
	
	fMatTempNew[3*n*Nvsteps + j*Nvsteps + i] += 
	  -1.*delT*(dR0/2.*double(i)*delV*fac3*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j+2)*Nvsteps)
		  + dR0/2.*Ais(j,2)*fMatTempOld[1*n*Nvsteps + (j+2)*Nvsteps + i]); 

	fMatTempNew[4*n*Nvsteps + j*Nvsteps + i] += 
	  -1.*delT*(dR0/2.*double(i)*delV*fac3*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + (j+2)*Nvsteps)
		  + dR0/2.*Ais(j,2)*fMatTempOld[2*n*Nvsteps + (j+2)*Nvsteps + i]);
	}

	j++;
      } // while (j < n)

    } //for (int i=0; i<Nvsteps; i++)
  } // if (expansionLevel == 5)

} // void fStepCompute
