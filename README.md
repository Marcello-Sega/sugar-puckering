





# sugar-puckering

<p align="center">
  <img src="http://www.science.unitn.it/~sega/Sugars_files/shapeimage_1.png" />
</p>

This repository contains 

1. a variant of the gromos 45a4 force field for carbohydrates, 
which improves puckering free energies [J. Chem. Phys. 133, 095104 (2010)](http://dx.doi.org/10.1063/1.3476466). 
The files include Gromacs topologies for all hexopyranoses 
(alpha and beta: Glc, Gal, Man, All, Tal, Gul, Alt, Ido) 
including the corrected parametrization reported in [J. Chem. Phys. 134, 149901 (2011)](http://dx.doi.org/10.1063/1.3578611)

2. the code (CP_reconstruct.c) used to reconstruct the atomic positions of 
6-membered rings, given the Cremer-Pople coordinates as input. The program computes also Strauss-Pickett 
out-of-plane dihedrals.


3. Some notes (CPinversion.pdf) on how to recover atomic coordinates from Cremer-Pople ones
