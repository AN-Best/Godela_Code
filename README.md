# Code for Godela Challenge
This repository contains the code for the take home challenge. <br>
The main code of interests is the OFOAM_Functions.py, there are three main functions within it: create_mesh, run_solver, view_results as well as some helper functions. The remaining codes are as follows:
- MultiRun.py - calls OFOAM_Functions in a loop with the different files being called
  - runs create_mesh and run_solver
- SingleRun_Clean.py - Just runs a single variant calling the functions
- ViewResults.py - Given the case_dir this plots the results
- SingleRun_Dirty.py - this was my initial attempt at interfaciing with OpenFOAM, it is messy but I thought it would be good to keep as it shows my process/progress 
