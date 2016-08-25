#ifndef FSTEPCOMPUTE2_H
#define FSTEPCOMPUTE2_H

void fStepCompute2( double * fMatTempNew, double * fMatTempOld, double delT, double delV, 
		   double c1, double nu, double omega, double vthe, double delR, 
		   double deln, int Nvsteps, int n, int expansionLevel, double kpar, double t, 
		   char * scattType, char * gradientOption );

#endif
