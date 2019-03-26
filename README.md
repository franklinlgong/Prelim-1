# Prelim-1

The .pdf file contains the written solutions to all the problems in the prelim.

prelim1problem2system.m contains the differential equations and parameters for part a of problem 2
Prelim_1_Problem_2.m contains the solver for the differential equations in prelim1problem2system.m for part a of problem 2

Prelim1Problem2b.m solves prelim1problem2system.m and stores the averaged concentration values for each protein over several intervals in a matrix. This is combined with manually shifting the parameter values and matrix name for use in Prelim1Problem2bStore.m which solves the finite difference between the forward perturbed and reverse perturbed matrices against an unperturbed value to give a manually stored set of matrices containing scaled sensitivity coefficients.

Prelim1Problem3b.m contains the linear programming functions for problem 3b and returns a 6 point graph of the protein concentration vs inducer concentration. 
