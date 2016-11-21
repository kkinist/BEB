# BEB
Calculate electron-ionization total ionization cross sections for molecules.

1.  Starting from simple coordinates, input files are generated for the Gaussian09 quantum chemistry program (G09). 
You must have a license for G09 to get any results. 
2.  There is a shell script to run all the G09 calculations. 
3.  The resulting output files are parsed to create a 'BUN' file of molecular data. 
4.  The BUN file is processed to run the BEB calculation of total ionization cross section (electron impact). 
