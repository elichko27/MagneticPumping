# MagneticPumping

In this repository is the set of c++ code that I use to compute the solution for the differential equations and also plot the results. 

To compile and run: 

$ g++ -O2 -o SolvingDiffEq_Main ./SolvingDiffEq_Main.cpp ./initConds.cpp ./simParams.cpp ./initCondswRestart.cpp ./R0.cpp ./nuF.cpp ./gradOpts.cpp ./Ais.cpp ./fStepCompute2.cpp
$ ./SolvingDiffEq_Main
