# MagneticPumping

In this repository is the set of c++ code that I use to compute the solution for the differential equations and also plot the results. 

To compile and run: 

$ g++ -O2 -o SolvingDiffEq_Main ./SolvingDiffEq_Main.cpp ./initConds.cpp ./simParams.cpp ./initCondswRestart.cpp ./R0.cpp ./nuF.cpp ./gradOpts.cpp ./Ais.cpp ./fStepCompute2.cpp
$ ./SolvingDiffEq_Main

Reminders for when using git: 

$ git add FILENAME

adds the file to the staging area, then 

$ git commit 

adds it to the branch. You can also do: 

$ git commit -m "Description of changes here"

to have an in-line description of the changes. 

Then when you're sure you want to add it to the main repository type: 

$ git push origin master