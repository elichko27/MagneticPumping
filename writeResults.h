#ifndef WRITERESULTS_H
#define WRITERESULTS_H

void writeResults(int i, int downSampleT, int downSampleV, double * fMatFinalOld, double *fMatFinal, 
		  int Nvsteps, int n, int expansionLevel, std::ofstream & outputFile, FILE * ptr_fp, 
		  time_t begin, time_t end);

#endif
