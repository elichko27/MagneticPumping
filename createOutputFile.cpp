//createOutputFile.cpp
#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <set>
#include <ctime>
#include "createOutputFile.h"

void createOutputFile (std::ofstream & outputFile, bool isRestart, double delT, double delV, 
		       double tMin, double tMax, double vMin, double vMax, int Nvsteps, int Ntsteps, 
		       int downSampleT, int downSampleV, char * rChoice, char * regionChoice, 
		       char * distChoice, char * scattType, char * gradientOption, double vthe, 
		       double memi, double B, double wpe, double wce, double wci, double beta_e, 
		       double vA, double c1, double nu, double omega, double delB, double deln, 
		       double delR, double kpar, int Lx, int l0) {
  
  if (isRestart == 1) {
    outputFile.open("outputFile.txt", std::ios::out | std::ios::app); 
    std::cout << "Restart -> appending already existing outfile" << std::endl; 
  } else {

    outputFile.open("outputFile.txt"); 
    
    // Adding important information to the output file
    outputFile << "Simulation Variables:" << std::endl; 
    outputFile << "delT = " << delT << std::endl; 
    outputFile << "delV = " << delV << std::endl;
    outputFile << "tMin = " << tMin << std::endl; 
    outputFile << "tMax = " << tMax << std::endl; 
    outputFile << "vMin = " << vMin << std::endl; 
    outputFile << "vMax = " << vMax << std::endl;
    outputFile << "Nvsteps = " << Nvsteps << std::endl;
    outputFile << "Ntsteps = " << Ntsteps << std::endl; 
    outputFile << "downSampleT = " << downSampleT << std::endl; 
    outputFile << "downSampleV = " << downSampleV << std::endl << std::endl; 
  
    outputFile << "Physical Variables:" << std::endl;
    outputFile << "rChoice = " << rChoice << std::endl; 
    outputFile << "regionChoice = " << regionChoice << std::endl; 
    outputFile << "distChoice = " << distChoice << std::endl; 
    outputFile << "scattType = " << scattType << std::endl; 
    outputFile << "gradientOption = " << gradientOption << std::endl << std::endl;

    outputFile << "Physical Parameters:" << std::endl; 
    outputFile << "vthe = " << vthe << std::endl; 
    outputFile << "memi = " << memi << std::endl; 
    outputFile << "B = " << B << std::endl;
    outputFile << "wpe = " << wpe << std::endl;
    outputFile << "wce = " << wce << std::endl;
    outputFile << "wci = " << wci << std::endl;
    outputFile << "beta_e = " << beta_e << std::endl;
    outputFile << "vA = " << vA << std::endl;
    outputFile << "c1 = " << c1 << std::endl;
    outputFile << "nu = " << nu << std::endl;
    outputFile << "omega = " << omega << std::endl;
    outputFile << "delB = " << delB << std::endl;
    outputFile << "deln = " << deln << std::endl;
    outputFile << "delR = " << delR << std::endl;
    outputFile << "kpar = " << kpar << std::endl;
    outputFile << "Lx = " << Lx << std::endl;
    outputFile << "l0 = " << l0 << std::endl << std::endl;

    outputFile << "Simulation Output:" << std::endl;
  }
}
