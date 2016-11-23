#!/bin/bash
# Run G09 jobs for BEB calculation (input files prepared by beb_g09build.py)
# Abort if the geometry optimization leads to a saddle point. 
# KKI
# 
if [ $# -lt 1 ]; then
	echo 'Please supply a filename root (i.e., specify molecule)'
	exit
fi
for suff in "_opt" "_bu" "_bupp" "_ept1" "_ept2" "_cc" "_cc1hi" "_cc1lo" "_cc2hi" "_cc2lo" ; do
	fin="$1$suff.gjf"
	fout="$1$suff.out"
	if [ -e "$fout" ]; then
		echo $fout already exists
	elif [ -e "$fin" ]; then
		echo Running $fin
		g09 < $fin > $fout
	fi
    if [ "$suff" = "_opt" ]; then
        # check number of imaginary vib. freqs
        nim=`nimag.pl $fout`
        if [[ $nim =~ [1-9]$ ]]; then
            echo "$nim"
            echo "----Adjust the geometry and try again."
            exit
        fi
    fi
done
echo "Calculations complete for: " $1
