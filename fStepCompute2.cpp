//fStepCompute2.cpp
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <cmath>
#include <set>
#include "R0.h"
#include "nuF.h"
#include "gradOpts.h"
#include "Ais.h"
#include "fStepCompute2.h"
#include "appendfMatNew.h"

void fStepCompute2( double * fMatTempNew, double * fMatTempOld, double * nuFac, double delT, double delV, 
		   double c1, double nu, double omega, double vthe, double delR, 
		   double deln, int Nvsteps, int n, int expansionLevel, double kpar, double t, 
		    char * scattType, char * gradientOption, double *facMat, int * expanMat ) {

  //Error catching 
  if (expansionLevel%2 != 1) {
    std::cout << "Invalid Expansion Level - expansion level must be odd!\n";
    //return 1;
  }

  //double nuFac0[Nvsteps];
  //double nuFac1[Nvsteps]; 
  //double nuFacn[Nvsteps]; 
  double jd, fac1, fac2, fac3; 
  int i, j, k; 

  double dR0 = R0(omega, delR, t);
  double dn0 = R0(omega, deln, t); 

  if (expansionLevel == 5) {

    for (int i=0; i<Nvsteps; i++) {

      //nuF(nuFac, nu, 0, vthe, delV, Nvsteps, scattType);

      fMatTempNew[0*n*Nvsteps + 0*Nvsteps + i] = fMatTempOld[0*n*Nvsteps + 0*Nvsteps + i] 
	- delT*(c1*fMatTempOld[0*n*Nvsteps + 0*Nvsteps + i] 
		+ nuFac[0*Nvsteps + i]*fMatTempOld[0*n*Nvsteps + 0*Nvsteps + i] 
		+ dR0/10.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 2*Nvsteps)
		+ dR0/2.*21./10.*fMatTempOld[1*n*Nvsteps + 2*Nvsteps + i]
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 0*Nvsteps));

      if (fMatTempNew[i] > 100){
	std::cout << std::endl << "i = " << i << "\n" 
		  << ", fMatTempOld[0*n*Nvsteps + 0*Nvsteps + i] = " << fMatTempOld[0*n*Nvsteps + 0*Nvsteps + i] << std::endl
		  << ", nuFac[i] = " << nuFac[0*Nvsteps + i] << std::endl 
		  << ", gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 2*Nvsteps) = " 
		  << gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 2*Nvsteps) << std::endl 
		  << ", fMatTempOld[1*n*Nvsteps + 2*Nvsteps + i] = " << fMatTempOld[1*n*Nvsteps + 2*Nvsteps + i] << std::endl 
		  << ", gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 0*Nvsteps) = " 
		  << gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 0*Nvsteps) << std::endl 
		  << ", dR0 = " << dR0 << std::endl 
		  << ", dn0 = " << dn0 << std::endl << std::endl; 
	exit(1); 
      }

      
      fMatTempNew[1*n*Nvsteps + 0*Nvsteps + i] = fMatTempOld[1*n*Nvsteps + 0*Nvsteps + i] 
	- delT*(c1*fMatTempOld[1*n*Nvsteps + 0*Nvsteps + i] 
		+ nuFac[0*Nvsteps + i]*fMatTempOld[1*n*Nvsteps + 0*Nvsteps + i] 
		+ dR0/5.*double(i)*delV*(gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + 2*Nvsteps) 
					 + 1./2.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + 2*Nvsteps))
		+ dR0*21./10.*(fMatTempOld[0*n*Nvsteps + 2*Nvsteps + i] 
			       + 1./2.*fMatTempOld[3*n*Nvsteps + 2*Nvsteps + i]) 
		+ kpar/3.*double(i)*delV*fMatTempOld[2*n*Nvsteps + 1*Nvsteps + i]
		+ dn0*double(i)*delV*(gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + 0*Nvsteps)
				      + 1./2.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + 0*Nvsteps)));

      fMatTempNew[2*n*Nvsteps + 0*Nvsteps + i] = fMatTempOld[2*n*Nvsteps + 0*Nvsteps + i]
	- delT*(c1*fMatTempOld[2*n*Nvsteps + 0*Nvsteps + i] 
		+ nuFac[0*Nvsteps + i]*fMatTempOld[2*n*Nvsteps + 0*Nvsteps + i] 
		+ dR0/10.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + 2*Nvsteps)
		+ dR0/2.*21./10.*fMatTempOld[4*n*Nvsteps + 2*Nvsteps + i]
		- kpar/3.*double(i)*delV*fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i] 
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + 0*Nvsteps));  
      
      fMatTempNew[3*n*Nvsteps + 0*Nvsteps + i] = fMatTempOld[3*n*Nvsteps + 0*Nvsteps + i] 
	- delT*(c1*fMatTempOld[3*n*Nvsteps + 0*Nvsteps + i] 
		+ nuFac[0*Nvsteps + i]*fMatTempOld[3*n*Nvsteps + 0*Nvsteps + i] 
		+ dR0/10.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 2*Nvsteps)
		+ dR0/2.*21./10.*fMatTempOld[1*n*Nvsteps + 2*Nvsteps + i] 
		+ 2.*kpar/3.*double(i)*delV*fMatTempOld[4*n*Nvsteps + 1*Nvsteps + i] 
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 0*Nvsteps));
      
      fMatTempNew[4*n*Nvsteps + 0*Nvsteps + i] = fMatTempOld[4*n*Nvsteps + 0*Nvsteps + i] 
	- delT*(c1*fMatTempOld[4*n*Nvsteps + 0*Nvsteps + i] 
		+ nuFac[0*Nvsteps + i]*fMatTempOld[4*n*Nvsteps + 0*Nvsteps + i] 
		+ dR0/10.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + 2*Nvsteps) 
		+ dR0/2.*21./10.*fMatTempOld[2*n*Nvsteps + 2*Nvsteps + i]
		- 2.*kpar/3.*double(i)*delV*fMatTempOld[3*n*Nvsteps + 1*Nvsteps + i] 
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 2*n*Nvsteps + 0*Nvsteps)); 

      // The f1 terms
      //nuF(nuFac, nu, 1, vthe, delV, Nvsteps, scattType);

      fMatTempNew[0*n*Nvsteps + 1*Nvsteps + i] = fMatTempOld[0*n*Nvsteps + 1*Nvsteps + i]
	- delT*(c1*fMatTempOld[0*n*Nvsteps + 1*Nvsteps + i] 
		+ nuFac[1*Nvsteps + i]*fMatTempOld[0*n*Nvsteps + 1*Nvsteps + i] 
		+ dR0/2.*double(i)*delV*2./5.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 1*Nvsteps) 
		+ dR0/2*3./5.*fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i]  
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 1*Nvsteps)); 
      
      fMatTempNew[1*n*Nvsteps + 1*Nvsteps + i] = fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i]
	- delT*(c1*fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i] 
		+ nuFac[1*Nvsteps + i]*fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i] 
		+ dR0*double(i)*delV*(2./5.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + 1*Nvsteps) 
				      + 1./2.*2./5.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + 1*Nvsteps))
		+ dR0*3./5.*(fMatTempOld[0*n*Nvsteps + 1*Nvsteps + i] + 1./2.*fMatTempOld[3*n*Nvsteps + 1*Nvsteps + i]) 
		+ kpar*double(i)*delV*(fMatTempOld[2*n*Nvsteps + 0*Nvsteps + i] 
				       + 2./5.*fMatTempOld[2*n*Nvsteps + 2*Nvsteps + i])
		+ dn0*double(i)*delV*(gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 0*n*Nvsteps + 1*Nvsteps) 
				      + 1./2.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 3*n*Nvsteps + 1*Nvsteps))); 
      
      fMatTempNew[2*n*Nvsteps + 1*Nvsteps + i] = fMatTempOld[2*n*Nvsteps + 1*Nvsteps + i]
	- delT*(c1*fMatTempOld[2*n*Nvsteps + 1*Nvsteps + i] 
		+ nuFac[1*Nvsteps + i]*fMatTempOld[2*n*Nvsteps + 1*Nvsteps + i] 
		+ dR0/2.*double(i)*delV*2./5.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + 1*Nvsteps)
		+ dR0/2*3./5.*fMatTempOld[4*n*Nvsteps + 1*Nvsteps + i] 
		- kpar*double(i)*delV*(fMatTempOld[1*n*Nvsteps + 0*Nvsteps + i] 
				       + 2./5.*fMatTempOld[1*n*Nvsteps + 2*Nvsteps + i])
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + 1*Nvsteps)); 
      
      fMatTempNew[3*n*Nvsteps + 1*Nvsteps + i] = fMatTempOld[3*n*Nvsteps + 1*Nvsteps + i]
	- delT*(c1*fMatTempOld[3*n*Nvsteps + 1*Nvsteps + i] 
		+ nuFac[1*Nvsteps + i]*fMatTempOld[3*n*Nvsteps + 1*Nvsteps + i] 
		+ dR0/2.*double(i)*delV*2./5.*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 1*Nvsteps)
		+ dR0/2.*3./5.*fMatTempOld[1*n*Nvsteps + 1*Nvsteps + i] 
		+ 2.*kpar*double(i)*delV*(fMatTempOld[4*n*Nvsteps + 0*Nvsteps + i] 
					  + 2./5.*fMatTempOld[4*n*Nvsteps + 2*Nvsteps + i])
		+ dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + 1*Nvsteps));

      fMatTempNew[4*n*Nvsteps + 1*Nvsteps + i] = fMatTempOld[4*n*Nvsteps + 1*Nvsteps + i]
	- delT*(c1*fMatTempOld[4*n*Nvsteps + 1*Nvsteps + i] 
		+ nuFac[1*Nvsteps + i]*fMatTempOld[4*n*Nvsteps + 1*Nvsteps + i] 
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

      // Higher order terms/equations
      j = 2; 
      while (j < n) {
	   
	//nuF(nuFac, nu, j, vthe, delV, Nvsteps, scattType);
	jd = double(j); 
	fac1 = 3./2.*jd*(jd-1)/(2.*jd-3.)/(2.*jd-1.);
	fac2 = 3./2.*(pow(jd+1.,2.)/(2.*jd+1.)/(2.*jd+3.) 
		      + pow(jd,2.)/(2*jd+1)/(2.*jd-1.)) - 1./2.;
	fac3 = 3./2.*(jd+1.)*(jd+2.)/(2.*jd+3.)/(2.*jd+5.);
	      
	fMatTempNew[0*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[0*n*Nvsteps + j*Nvsteps + i]
	  - delT*(c1*fMatTempOld[0*n*Nvsteps + j*Nvsteps + i] 
		  + nuFac[j*Nvsteps + i]*fMatTempOld[0*n*Nvsteps + j*Nvsteps + i] 
		  + dR0/2.*double(i)*delV*(fac1*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j-2)*Nvsteps) 
					   + fac2*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j)*Nvsteps))
		  + dR0/2.*(Ais(j,0)*fMatTempOld[1*n*Nvsteps + (j-2)*Nvsteps + i] 
			    + Ais(j,1)*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i])
		  + dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j)*Nvsteps)); 

	fMatTempNew[1*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[1*n*Nvsteps + j*Nvsteps + i]
	  - delT*(c1*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i] 
		  + nuFac[j*Nvsteps + i]*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i] 
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
		  + nuFac[j*Nvsteps + i]*fMatTempOld[2*n*Nvsteps + j*Nvsteps + i]
		  + dR0/2.*double(i)*delV*(fac1*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + (j-2)*Nvsteps) 
					   + fac2*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + (j)*Nvsteps))
		  + dR0/2.*(Ais(j,0)*fMatTempOld[4*n*Nvsteps + (j-2)*Nvsteps + i] 
			    + Ais(j,1)*fMatTempOld[4*n*Nvsteps + j*Nvsteps + i])
		  - kpar*double(i)*delV*jd/(2.*jd-1.)*fMatTempOld[1*n*Nvsteps + (j-1)*Nvsteps + i]
		  + dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 4*n*Nvsteps + (j)*Nvsteps)); 
	
	fMatTempNew[3*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[3*n*Nvsteps + j*Nvsteps + i]
	  - delT*(c1*fMatTempOld[3*n*Nvsteps + j*Nvsteps + i] 
		  + nuFac[j*Nvsteps + i]*fMatTempOld[3*n*Nvsteps + j*Nvsteps + i] 
		  + dR0/2.*double(i)*delV*(fac1*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j-2)*Nvsteps) 
					   + fac2*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j)*Nvsteps))
		  + dR0/2.*(Ais(j,0)*fMatTempOld[1*n*Nvsteps + (j-2)*Nvsteps + i] 
			    + Ais(j,1)*fMatTempOld[1*n*Nvsteps + j*Nvsteps + i])
		  + 2.*kpar*double(i)*delV*jd/(2.*jd-1)*fMatTempOld[4*n*Nvsteps + (j-1)*Nvsteps + i]
		  + dn0/2.*double(i)*delV*gradOpts(fMatTempOld, delT, i, 0, Nvsteps, 1*n*Nvsteps + (j)*Nvsteps)); 

	fMatTempNew[4*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[4*n*Nvsteps + j*Nvsteps + i]
	  - delT*(c1*fMatTempOld[4*n*Nvsteps + j*Nvsteps + i] 
		  + nuFac[j*Nvsteps + i]*fMatTempOld[4*n*Nvsteps + j*Nvsteps + i] 
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
  } else { 
    
    char isGrad[] = "Grad"; 
    char isntGrad[] = "NoGrad"; 
    
    for(k = 0; k < expansionLevel; k++) { 
      for(j = 0; j < n; j++) { 
	jd = double(j); 
	fac1 = 3./2.*jd*(jd-1)/(2.*jd-3.)/(2.*jd-1.);
	fac2 = 3./2.*(pow(jd+1.,2.)/(2.*jd+1.)/(2.*jd+3.) 
			+ pow(jd,2.)/(2*jd+1)/(2.*jd-1.)) - 1./2.;
	fac3 = 3./2.*(jd+1.)*(jd+2.)/(2.*jd+3.)/(2.*jd+5.);
	
	for(i = 0; i < Nvsteps; i++) { 

	  fMatTempNew[k*n*Nvsteps + j*Nvsteps + i] = fMatTempOld[k*n*Nvsteps + j*Nvsteps + i];
	  
	  appendfMatNew(i, j, k, fMatTempNew, i, j, k, fMatTempOld, 
			-1.0*delT*dR0*fac2*double(i)*delV, facMat, expanMat, isGrad, 
			1, expansionLevel, Nvsteps, n, delV); 
	  appendfMatNew(i, j, k, fMatTempNew, i, j, k, fMatTempOld, 
			-1.0*delT*dR0*Ais(j,1), facMat, expanMat, isntGrad, 
			1, expansionLevel, Nvsteps, n, delV);
	  appendfMatNew(i, j, k, fMatTempNew, i, j, k, fMatTempOld, 
			-1.0*delT*dn0*double(i)*delV, facMat, expanMat, isGrad, 
			1, expansionLevel, Nvsteps, n, delV);

	  if (j-2 >= 0) {
	    appendfMatNew(i, j, k, fMatTempNew, i, j-2, k, fMatTempOld, 
			  -1.0*delT*dR0*fac1*double(i)*delV, facMat, expanMat, isGrad, 
			  1, expansionLevel, Nvsteps, n, delV); 
	    appendfMatNew(i, j, k, fMatTempNew, i, j-2, k, fMatTempOld, 
			  -1.0*delT*dR0*Ais(j,0), facMat, expanMat, isntGrad, 
			  1, expansionLevel, Nvsteps, n, delV);
	  }

	  if (j-1 >= 0 && k != 0 && kpar != 0.) { 
	    appendfMatNew(i, j, k, fMatTempNew, i, j-1, k, fMatTempOld, 
			  -1.0*delT*kpar*double(i)*delV*jd/(2.*jd -1.), facMat, expanMat, isntGrad, 
			  0, expansionLevel, Nvsteps, n, delV);
	  }

	  if (n > j+1 && k != 0 && kpar != 0.) { 
	    appendfMatNew(i, j, k, fMatTempNew, i, j+1, k, fMatTempOld, 
			  -1.0*delT*kpar*double(i)*delV*(jd + 1.)/(2.*jd + 3.), facMat, expanMat, isntGrad, 
			  0, expansionLevel, Nvsteps, n, delV);
	  }
	  
	  if (n > j+2) {
	    appendfMatNew(i, j, k, fMatTempNew, i, j+2, k, fMatTempOld, 
			  -1.0*delT*dR0*fac3*double(i)*delV, facMat, expanMat, isGrad, 
			  1, expansionLevel, Nvsteps, n, delV); 
	    appendfMatNew(i, j, k, fMatTempNew, i, j+2, k, fMatTempOld, 
			  -1.0*delT*dR0*Ais(j,2), facMat, expanMat, isntGrad, 
			  1, expansionLevel, Nvsteps, n, delV);
	  }

	  // Adding in infinite terms
	  if (n > j+4 && n >= 5) { 
	    // count down using for loop and II -= 2 to 0 && append appropriately
	  }

	}
      }
    }
    
  }// if (expansionLevel == 5)

} // void fStepCompute
