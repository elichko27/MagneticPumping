#ifndef FSTEPCOMPUTE_H
#define FSTEPCOMPUTE_H

void fStepCompute( double * fMatTempNew, double * fMatTempOld, double delT, double delV, 
		   double c1, double nu, double omega, double vthe, double delR, 
		   double deln, int Nvsteps, int n, int expansionLevel, double kpar, double t, 
		   char * scattType, char * gradientOption );

#endif

