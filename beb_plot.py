#!/usr/bin/python3
#
# Given a BUN file and the script beb_tbl.pl, create a plot of 
#   TICS vs. incident electron energy
# K. K. Irikura, NIST 2017
# 
# save plot to a PNG file, 4/27/2018
#
import sys
import re
import os
import matplotlib.pyplot as plt
#
# process command-line arguments
if len(sys.argv) < 2:
    sys.exit( 'Usage:  beb_plot.py <BUNfile> <graph title>' )
try:
    finp = open( sys.argv[1] )
    fname = sys.argv[1]
except:
    # assume user forgot file suffix, probably .bun
    finp = open( sys.argv[1] + '.bun' )
    fname = sys.argv[1] + '.bun'
froot = os.path.splitext( sys.argv[1] )[0]
# look for molecular formula in BUN file, using format from beb_g09parse.py
regex = re.compile(r'# Molecule is .* (\S+)\)')
for line in finp.readlines():
    m = regex.match(line)
    if m:
        # label and stoichiometry should be on this line
        #   label should be the same as 'froot'
        stoich = m.group(1)
    else:
        # this BUN file not generated by beb_g09parse.py; use filename root
        m = re.match(r'(.+)\.?', sys.argv[1])
        stoich = m.group(1)
try:
    figtitle = sys.argv[2]
except:
    # it was probably left blank; use the stoichiometry
    figtitle = stoich
#
# run beb_tbl.pl
result = os.popen('perl beb_tbl.pl {:s}'.format(fname)).read()
regtbl = re.compile(r'Energy .* BEB')
regblank = re.compile(r'^\s*$')
intbl = False   # flag
for line in result.split('\n'):
    if intbl:
        if regblank.match(line):
            # done reading data
            intbl = False
            break
        # read a line of data
        data = line.split()
        energy.append(data[0])
        xsec.append(data[1])
    if regtbl.search(line):
        # header for data table
        intbl = True
        # initialize data lists
        energy = []
        xsec = []
#
# create the plot
plt.semilogx(energy, xsec, 'o-')
plt.title(figtitle)
plt.xlabel('Incident energy / eV')
plt.ylabel('Total ionization cross section / 10$^{-20}$ m$^2$')
# save the plot as PNG
plotfile = stoich + '_BEB_plot.png'
plt.savefig(plotfile)
# display the plot
plt.show()
