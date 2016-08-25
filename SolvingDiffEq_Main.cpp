//SolvingDiffEq_Main.cpp
#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <set>
#include <ctime>
#include "initConds.h"
#include "simParams.h"
#include "initCondswRestart.h"
#include "R0.h"
#include "nuF.h"
#include "gradOpts.h"
#include "Ais.h"
#include "fStepCompute2.h"
using namespace std; 


int main () {
  time_t begin, end;
  time(&begin); 

  //////////////////////
  // Input Parameters //
  //////////////////////
  //In the future switch these over to an input file
  int n = 3;                //number of equations - goes to n-1 order
  int expansionLevel = 5;   //number of cos terms in expansion

  double delT = 0.005; 
  double delV = 0.0001; 

  double tMin = 0; 
  double tMax = 600; 

  bool isRestart = 0; 
  int downSampleT = 1000; 
  int downSampleV = 2; 

  int fMatFinalTMax = 5; 

  // Run comparison parameters
  char rChoice[] = "1DMagneticPumping";
  char regionChoice[] = "Region1";
  char distChoice[] = "Maxwellian"; 
  char scattType[] = "Constant"; 
  char gradientOption[] = "Standard"; 
  double nu = 30.0*pow(10.0,-4.0); 

  double vthe, memi, B, wpe, wce, wci, beta_e, vA, c1, vMin, vMax, omega, delB,
    deln, delR, kpar; 
  int Lx, l0; 
  simParams(rChoice, regionChoice, vthe, memi, B, wpe, wce, wci, beta_e, vA, c1, 
	  vMin, vMax, omega, delB, deln, delR, kpar, Lx, l0);

  /////////////////////////
  // Computed Parameters //
  /////////////////////////
  int Ntsteps = int((tMax - tMin)/delT);
  int Nvsteps = int((vMax - vMin)/delV);
  
  int ft; 
  char fileTag [1000];
  ft = sprintf(fileTag, "_PnExpan%d_CosExpan%d_nu%d_dist_dt%d_dv%d_tMax%d_vMax%d___downSampleT%d_downSampleV%d", 
	       /*rChoice,*/ n, expansionLevel, int(nu*pow(1.,5.)),/* distChoice,*/ int(1./delT), int(1./delV), 
	       int(tMax*1.0), int(vMax*10.)/*, scattType, gradientOption*/, downSampleT, downSampleV); 

  ////////////////////////////////////////////////////////
  // Initial Conditions and Saving Output Initilization //
  ////////////////////////////////////////////////////////

  int fMatFinalCounter; 
  //double * fMatFinalOld, fMatFinalNew, fMatFinal; 
  double * fMatFinalOld = new double[Nvsteps*n*expansionLevel];
  double * fMatFinalNew = new double[Nvsteps*n*expansionLevel];
  double * fMatFinal = new double[Nvsteps/downSampleV*n*expansionLevel];

  initCondswRestart(fMatFinalCounter, fMatFinalOld, fMatFinalNew, fMatFinal, 
		    Ntsteps, isRestart, Nvsteps, n, expansionLevel, vthe, delV, vMax, 
		    distChoice, rChoice, fileTag, tMax, tMin, delT, fMatFinalTMax, 
		    downSampleT, downSampleV);

  /*ofstream outfile("cppTestOutput.bin", ios::out | ios::binary);
  if(!outfile) {
    cout << "Cannot open file.";
    return 1;
    }*/ 

  FILE *ptr_fp; 
  if((ptr_fp = fopen("cppTestOutput.bin", "ab")) == NULL)
    {
      printf("Unable to open file!\n");
      exit(1);
    }else printf("Opened file successfully for writing.\n");
  cout << "Ntsteps = " << Ntsteps << ", Nvsteps = " << Nvsteps << endl; 

  /*ptr_fp = fopen("test.bin", "ab");
  if (ptr_fp == NULL) { 
    printf("Unable to open file!\n"); 
    exit(1); 
  } else printf("Opened file successfully for writing.\n");
  fclose(ptr_fp); */

  ////////////////////////////
  // Running the Simulation //
  ////////////////////////////

  for (int i = 0; i<Ntsteps; i++) {
    if (i%1000 == 0) {
      time(&end); 
      cout << i << " tSteps in : " << difftime(end, begin) << " seconds" << endl; 
    }
 
    // Writing the result to a binary file
    if (i%downSampleT == 0) {
      if (downSampleV == 1) {
	//outfile.write((char *) &fMatFinalOld, sizeof fMatFinalOld);
	fwrite(fMatFinalOld, Nvsteps*n*expansionLevel*sizeof(double), 1, ptr_fp);
      } else {

	int jcount = 0; 
	for(int j = 0; j<Nvsteps; j+=downSampleV) {
	  for(int k = 0; k<n; k++) 
	    for(int m = 0; m<expansionLevel; m++) { 
	      fMatFinal[m*n*Nvsteps/downSampleV + k*Nvsteps/downSampleV + jcount] 
		= fMatFinalOld[m*n*Nvsteps+k*Nvsteps+j]; 
	    }
	  jcount++; 
	}

	//outfile.write((char *) &fMatFinal, sizeof fMatFinal);
	fwrite(fMatFinal, Nvsteps/downSampleV*n*expansionLevel*sizeof(double), 1, ptr_fp);
      }
      //cout << "Tstep: " << i << " out of: " << Ntsteps << endl; 
    }

    // Computing the next time step
    fStepCompute2( fMatFinalNew, fMatFinalOld, delT, delV, c1, nu, omega, vthe, delR, 
		   deln, Nvsteps, n, expansionLevel, kpar, double(i)*delT, scattType, gradientOption );

    // Making the new matrix the old matrix
    copy(fMatFinalNew, fMatFinalNew+Nvsteps*n*expansionLevel, fMatFinalOld);  

  } //for (int i = 0; i<Ntsteps; i ++)


  //outfile.close();
  fclose(ptr_fp);
  // Printing the last fMatFinalOld just to see what happens
  for (int j = 0; j<Nvsteps*n*expansionLevel; j++) { 
    if(j%(Nvsteps) == 0 && j%n == 0 && j != 0)  cout << endl; 
    if(j%(Nvsteps) == 0 && j%expansionLevel == 0 && j != 0) cout << endl;
    if(j%(Nvsteps) == 0 && j%expansionLevel == 0 && j%n == 0) cout << endl <<"expanLevel = " << j/Nvsteps/n << endl;
    cout << " " << fMatFinalOld[j] << " "; 
  }

  // Calculating how long the code ran for
  delete[] fMatFinalOld;
  delete[] fMatFinalNew;
  delete[] fMatFinal;
  time(&end);
  double elapsed_secs = difftime(end,begin);
  cout << endl; 
  cout << "Ntsteps = " << Ntsteps << ", Nvsteps = " << Nvsteps << endl; 
  cout << "Program is finished in: " << elapsed_secs << "secs" << endl << endl; 
  return 0; 
}

  // Unused loading code
  /*  int size_of_vec = Nvsteps*n*expansionLevel*2; 
  double fread[size_of_vec]; 
  ifstream infile("cppTestOutput.bin", ios::out | ios::binary);
  infile.read((char *) &fread, sizeof fread); 
  cout << "This is the file I saved: " << endl; 
  for(int i = 0; i<size_of_vec; i++) cout << fread[i] << " "; 
  cout << endl << endl; 
  infile.close(); 
  cout << "This is the last fMatFinal:" << endl; 
  for(int i = 0; i<Nvsteps*n*expansionLevel; i++) cout << fMatFinalNew[i] << " ";
  cout << endl << endl; */
