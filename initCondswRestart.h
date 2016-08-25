#ifndef INITCONDSWRESTART_H
#define INITCONDSWRESTART_H

void initCondswRestart(int &fMatFinalCounter, double * fMatFinalOld, 
		       double * fMatFinalNew, double * fMatFinal, 
		       int &Ntsteps, bool isRestart, int Nvsteps, int n, int expansionLevel, double vthe, 
		       double delV, double vMax, char * distChoice, char * rChoice, char * fileTag, 
		       double tMax, double tMin, double delT, int fMatFinalTMax,
		       int downSampleT, int downSampleV);
#endif
