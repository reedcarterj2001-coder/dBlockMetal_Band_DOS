# dBlockMetal_Band_DOS
Scripts for solid state Spring 2026 term project on electronic energy band and density of states calculations applied to the d-block transition metals.

1. The user runs controller.m and selects which metal they want by typing its atomic symbol into the input (not case sensitive).
2. controller.m then calls metalDatabase. for information on the metal matching the symbol the user input and uses convertMetalUnits.m
   to put energies in eV and lattice parameters in Angstroms.
3. controller.m takes the metal (now in eV and A) and gives it to buildCrystal.m which reads the structure tagged to the metal (BCC, FCC, or HCP) and generates
   the appropriate geometric properties:
   - Primitive lattice vectors and basis
   - High symmetry k-points
   - k-path
   - Nearest neighbor list
   - Nearest neighbor distance
   - Number of atoms per primitive cell
   - Reciprocal lattice vectors using reciprocalLattice.m
4. controller.m feeds the crystal and a partiion value to generateKPath.m to get the k-points
   along the path and the positions of the symmetry labels for the band plot.
5. The k-path points are then sent to computeSpectum.m to obtain various information like the bands.
   computeSpectrum.m is used so that the Hamiltonian only has to be computed once reducing the timecost
   of the program. The Hamiltonian is found using the aptly named buildHamiltonian.m which relies on
   SlaterKosterMatrix.m for the onsite and neghbor-hopping physics. The Slater-Koster parameters were
   taken from Papaconstantopoulos's "Handbook of the Band Structure of Elemental Solids" second edition
   published by Springer.
6. After the spectrum is returned, the DOS is computed with user-set binning and Gaussian smearing
   parameters. The parameters are in controller.m directly above the call for computeDOS.m. The electron
   spin degeneracy and DOS normalization are handled here.
8. From the DOS, the Fermi energy is found using findFermi.m
9. Just about all of the information is then funneled to plotBandDOS.m for plotting. This is where the
   plots are shifted by the Fermi energy to prevent shift related bugs. This also calls
   computeOrbitalWeights.m to get the orbital contributions to each band for a tritone heatmap in the
   band plot. The script ends with figure exporting functionality using MATLAB's exportgraphics function
   and set to export as a .tif. The user needs to add their figure directory as the value of outdir.
