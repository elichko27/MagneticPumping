#ifndef INITCONDSWRESTART_H
#define INITCONDSWRESTART_H

void initCondswRestart(int &fMatFinalCounter, double * fMatFinalOld, 
		       double * fMatFinalNew, double * fMatFinal, double * nuFac,
		       int &Ntsteps, bool isRestart, int Nvsteps, int n, int expansionLevel, double vthe, double nu, 
		       double delV, double vMax, char * distChoice, char * rChoice, char * fileTag, char * scattType, 
		       double tMax, double tMin, double delT, int downSampleT, int downSampleV);
#endif
