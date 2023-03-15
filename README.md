# Thesis
This repository will house all image analysis programs related to my thesis, INVESTIGATING THE ROLE OF LIPID DROPLETS AS INTRACELLULAR MECHANICAL STRESSORS IN NON-ALCOHOLIC FATTY LIVER DISEASE, Abigail Erin Loneker, 2023. 


calcNuclearIrregulairty is a function to quantify nuclear deformation, globally with the nuclear irregularity measure, and locally by estimating curvature along the linearized membrane boundary. Additional details on the rationale for this can be found in a preprint available on bioRxiv https://doi.org/10.1101/2022.08.27.505524, or in the text of my thesis (to be made available through UPenn Scholarly Commons https://repository.upenn.edu/etd.html)

The irregularity program was modified to estimate nuclear membrane fluctuations (see folder NuclearMembraneFluctuations). The program linearizes the membrane boundary at each time point and then computes a mean squared displacement for each point along the membrane boundary. Additional details can be found in the text of my thesis. 
Programs for this analysis are found in the NuclearMembraneFluctuations folder. NuclearFlux_Main is the main program, which calls calcNuclearFluctuations to linearize the membrane boundaries. Both programs are needed for this analysis.

lipidDropletStiffness is a program to estimate the stiffness of lipid droplets embedded in polyacrylamide gels and subjected to specified stresses or strains. The program requires image stacks of lipid droplet cross-sections before (first image) and after compression (all subsequent images). The analysis is based on the Eshelby inclusion problem for spherical inclusions. Additional details can be found in the text of my thesis.
