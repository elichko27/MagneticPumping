#ifndef CREATEOUTPUTFILE_H
#define CREATEOUTPUTFILE_H

void createOutputFile (std::ofstream & outputFile, bool isRestart, double delT, double delV, 
		       double tMin, double tMax, double vMin, double vMax, int Nvsteps, int Ntsteps, 
		       int downSampleT, int downSampleV, char * rChoice, char * regionChoice, 
		       char * distChoice, char * scattType, char * gradientOption, double vthe, 
		       double memi, double B, double wpe, double wce, double wci, double beta_e, 
		       double vA, double c1, double nu, double omega, double delB, double deln, 
		       double delR, double kpar, int Lx, int l0);
//void createOutputFile (std::ofstream & outputFile, bool isRestart);

#endif
