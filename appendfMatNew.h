#ifndef APPENDFMATNEW_H
#define APPENDFMATNEW_H

void appendfMatNew(int iN, int jN, int kN, double * fMatTempNew, int iO, int jO, int kO, double * fMatTempOld, 
		   double overallFactor, double * facMat, int * expanMat, char * isGrad, int isOneCosSin, 
		   int expansionLevel, int Nvsteps, int n, double delV);

#endif
