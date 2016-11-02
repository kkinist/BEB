#!/usr/bin/python3
# Create Gaussian09 input files needed for BEB calculation on a neutral target.
# Basic method:  G. E. Scott and K. K. Irikura, Surf. Interf. Anal. 37, 973-977 (2005)
#
# Required input file format (whitespace-delimited): 
#   line 1: charge, spin-multiplicity [integers]
#   next: element symbol, cartesian triple [one atom per line]
#   final blank line
#
#   Input files as examples and for testing:
#       h2o.inp (singlet, non-degenerate orbitals)
#       ch4.inp (singlet, degenerate orbitals)
#       no2.inp (doublet, non-degen)
#       ch3.inp (doublet, degen)
#       ch2.inp (triplet, non-degen)
#       o2.inp  (triplet, degen)
# Command-line arguments can be used to specify memory and number of processors:
#   -m<mwords> -n<nprocs>
#
# BEB is defined for neutral and singly-positive targets but only neutrals are
#   handled by this software.  
# Latest change: 18 August 2016
#
# DISCLAIMER:  
#   Certain commercial software is identified in order to specify procedures
#   completely.  In no case does such identification imply recommendation or
#   endorsement by the National Institute of Standards and Technology, nor does
#   it imply that the material or equipment identified is necessarily the best
#   available for the purpose.
#
# Software citation: 
#   BEB automation software, Karl K. Irikura, National Institute of Standards
#   and Technology (2016). 
#
import sys
import os
import re
from beb_subs import *
##
def get_elements( coords ):
    # parse the list of lines from the input file
    # return a list of element symbols found
    amongatoms = 0
    atoms = []
    for line in coords:
        m = re.search( r'^\s*([A-Z][a-z]?)\s*', line )
        if m:
            amongatoms = 1
            atoms.append( m.group(1) )
        if amongatoms and line.isspace():
            # all atoms have been read
            break
    return atoms
##
def sort_atoms( atoms, name='z' ):
    # return three lists of Z: light, medium, and heavy atoms (they need different
    #   basis sets)
    # parameter 'name' specifies whether to return atomic numbers or symbols
    atoms = set( atoms )
    atno = list( map( elz, atoms ) )
    lowz = []
    midz = []
    highz = []
    for z in atno:
        if z < 11:
            lowz.append( z )
        elif z < 37:
            midz.append( z )
        else:
            highz.append( z )
    if name == 'symbol':
        # return element symbols instead of atomic numbers
        lowz = list( map( elz, lowz ) )
        midz = list( map( elz, midz ) )
        highz = list( map( elz, highz ) )
    return lowz, midz, highz
##
def specify_basis( atoms, step=1 ):
    # determine appropriate basis set in BEB calculation
    # return:
    #   formatted string suitable for G09 input file,
    #   if step==4: also return dict of number of orbitals replaced by ECPs: nopp{symbol}
    # first arg is coordinates, to be parsed for element symbols
    # second arg is BEB step number:
    #   1 = geometry optimization
    #   2 = binding and kinetic energies -- all-electron
    #   3 = as above, but with ECP on most atoms (if needed)
    #   4 = correlated valence binding energies
    #   5 = high-level 1st vertical ionization energy (IE)
    #   6 = correlated 2nd vertical IE (two alternatives)   !!! NOT YET DONE !!!
    light = [ '6-31G(d)', '6-311G(d,p)', '6-311G(d,p)', '6-311+G(d,p)', 'cc-pVTZ' ]
    middle = [ '6-31G(d)', '6-311G(d,p)', 'SDDall', '6-311+G(d,p)', 'cc-pVTZ' ]
    heavy = [ '3-21G*', '3-21G*', 'SDD', 'SDD', 'SDD' ]
    # main-group pseudopotentials need d-polarization functions added (from 6-311G**)
    dpol = { 'Na': 0.175, 'Mg': 0.175,
        'Al': 0.325, 'Si': 0.45, 'P': 0.55, 'S': 0.65, 'Cl': 0.75, 'Ar': 0.85,
        'K': 0.229, 'Ca': 0.26, 'Ga': 0.169, 'Ge': 0.228, 'As': 0.264,
        'Se': 0.305, 'Br': 0.451, 'Kr': 0.395, 'I': 0.302 }
    plusd = [ 'Al', 'Si', 'P', 'S', 'Cl', 'Ar' ]
    basisstring = ''
    ecpstring = ''
    lowz, midz, highz = sort_atoms( atoms, 'symbol' )
    if len( lowz ) > 0:
        basisstring += " ".join( lowz ) + ' 0\n'
        basisstring += light[ step - 1 ]
        basisstring += '\n****\n'
    if len( midz ) > 0:
        if re.search( 'cc-pVTZ', middle[step-1] ):
            # do we need cc-pVTZ, cc-pV(T+d)Z, or both?
            # add tight d-function to old cc basis set
            needd = list( set(plusd) & set(midz) )
            noneed = list( set(midz) - set(plusd) )
            if len( needd ):
                basisstring += '@/home/irikura/basis/t_plusd.gbs\n'
            if len( noneed ):
                basisstring += " ".join( noneed ) + ' 0 \n'
                basisstring += middle[step-1] + '\n****\n'
        else:
            basisstring += " ".join( midz ) + ' 0\n'
            basisstring += middle[ step - 1 ]
            if re.search( 'SDD', middle[step-1] ):
                # need to specify pseudopotential
                ecpstring += " ".join( midz ) + ' 0\n'
                ecpstring += 'SDDall\n'
                # may need to add d-polarization function to valence basis
                for el in midz:
                    if el in dpol:
                        # add d-polarization
                        basisstring += '\n D   1 1.00       0.000000000000\n'
                        basisstring += '      %16.10e  0.1000000000D+01' % dpol[el]
                        break
            basisstring += '\n****\n'
    if len( highz ) > 0:
        # polarization functions available only for some elements; put one per line
        isecp = re.search( 'SDD', heavy[step-1] )
        for el in highz:
            basisstring += el + ' 0\n'
            basisstring += heavy[ step - 1 ]
            if el in dpol and isecp:
                # d-polarization function is available for this element--use it
                basisstring += '\n D   1 1.00       0.000000000000\n'
                basisstring += '      %16.10e  0.1000000000D+01' % dpol[el]
            basisstring += '\n****\n'
        if isecp:
            # need to specify pseudopotential
            ecpstring += " ".join( highz ) + ' 0\n'
            ecpstring += 'SDD\n'
    basisstring += '\n'
    if ( len(ecpstring) > 0 ):
        basisstring += ecpstring + '\n'
    if (step == 4):
        # EPT calculation needs to know how many orbitals are replaced by ECPs; create dict
        nopp = {}
        for el in atoms:
            nopp[el] = 0
        for el in midz:
            if re.search('SDD', middle[3]):
                nopp[el] = n_sdd(el) / 2
        for el in highz:
            if re.search('SDD', heavy[3]):
                nopp[el] = n_sdd(el) / 2
        return basisstring, nopp
    return basisstring
##
def need_ecp( basisstring ):
    # return 1 if basisstring specifies ECP-basis, else 0
    regx = re.compile( 'SDD' )
    mch = regx.search( basisstring )
    if mch:
        return 1
    return 0
##
def need_6D(bstring):
    # return True if only 6-31G* basis sets are in use, else False
    # want to match: 6-31g* 6-31g(d 6-31+g* 6-31++g*
    d6 = True
    pat = '6-31[+]{0,2}[Gg](\*|\(d)'
    lines = bstring.split('\n')
    for i in range(len(lines)):
        if re.search('0\s*$', lines[i]):
            # atom specification line; look at next line
            i += 1
            d6 = d6 and re.search(pat, lines[i])
    return d6
#
# MAIN
#
# For the Gaussian09 program, create the following input files. Some molecules
# will not require all these computations. 
#
import re
from beb_subs import *
##
fn_cc = froot + '_cc.gjf'
fgjf = open( fn_cc, 'w' )
print( 'created file ', fn_cc )
basis = specify_basis( elem, 5 )
directive = '# CCSD(T)/GEN geom=check scf=xqc'
if need_ecp( basis ):
    directive += ' pseudo=read\n'
else:
    directive += '\n'
comment = '\nBEB step ' + str(step) + ': CCSD(T) for neutral ' + froot + '\n\n'
fgjf.write( header + directive + comment )
fgjf.write( coords[0] + '\n' + basis )
#
# sub-step B: vertical single ionization
c2 = comment.replace( 'neutral', 'singly ionized' )
dir2 = directive.replace( 'GEN', 'CHKBAS' )
dir2 = dir2.replace( 'pseudo=read', '' )
fgjf.write( '--Link1--\n' + header + dir2 + c2 )
charge += 1
mult += 1       # consider two spin states if the neutral is a radical
fgjf.write( "%d %d\n\n" % (charge, mult) )
if mult > 2:
    # the neutral was not a singlet; consider lower-spin ion
    c3 = c2.replace( 'singly', 'low-spin singly' )
    fgjf.write( '--Link1--\n' + header + dir2 + c3 )
    fgjf.write( "%d %d\n\n" % (charge, mult-2) )
fgjf.close()
#
# create input files for double ionization threshold (step 6)
# provide two alternatives:
#   EPT/big-basis of cation (fairly cheap)
#   CCSD(T)/big-basis of low- and high-spin dications (most expensive step in whole series)
#
step = '6A'
fn_ept2 = froot + '_ept2.gjf'
fgjf = open(fn_ept2, 'w')
print( 'created file ', fn_ept2 )
c2 = '\nBEB step ' + str(step) + ': EPT for IE of ' + froot + ' cation \n\n'
fgjf.write( header + directive_ept + c2 )
fgjf.write( "%d %d\n\n" % (charge, mult) )
fgjf.write( basis )
# 'mult' is already one higher than the input neutral value; adjust 'nbet' to match
nbet -= 1;
fgjf.write( "%d %d %d %d\n\n" % (nalp, nalp, nbet, nbet) )
if mult > 2:
    # the neutral was not a singlet; also consider low-spin cation reference
    c2 = c2.replace( 'cation ', 'low-spin cation ' )
    fgjf.write( '--Link1--\n' + header + directive_ept + c2 )
    fgjf.write( "%d %d\n\n" % (charge, mult-2) )
    fgjf.write( basis )
    # adjust 'nalp' and 'nbet' for low spin
    nbet += 1
    nalp -= 1
    fgjf.write( "%d %d %d %d\n\n" % (nalp, nalp, nbet, nbet) )
fgjf.close()
#
step = '6B'
fn_dbl = froot + '_dbl.gjf'
fgjf = open(fn_dbl, 'w')
print( 'created file ', fn_dbl )
c2 = '\nBEB step ' + str(step) + ': CCSD(T) for doubly ionized ' + froot + ' (high-spin)\n\n'
charge += 1
mult += 1
dir2 = dir2.replace( 'CCSD(T)', 'CCSD(T,maxcyc=100)' )  # dications can be hard to converge
fgjf.write( header + dir2 + c2 )
fgjf.write( "%d %d\n\n" % (charge, mult) )
# consider lower-spin ion
c3 = c2.replace( 'high-spin', 'low-spin' )
fgjf.write( '--Link1--\n' + header + dir2 + c3 )
fgjf.write( "%d %d\n\n" % (charge, mult-2) )
fgjf.close()

#!/usr/bin/python3
# Create Gaussian09 input files needed for BEB calculation on a neutral target.
# Basic method:  G. E. Scott and K. K. Irikura, Surf. Interf. Anal. 37, 973-977 (2005)
#
# Required input file format (whitespace-delimited): 
#   line 1: charge, spin-multiplicity [integers]
#   next: element symbol, cartesian triple [one atom per line]
#   final blank line
#
#   Input files as examples and for testing:
#       h2o.inp (singlet, non-degenerate orbitals)
#       ch4.inp (singlet, degenerate orbitals)
#       no2.inp (doublet, non-degen)
#       ch3.inp (doublet, degen)
#       ch2.inp (triplet, non-degen)
#       o2.inp  (triplet, degen)
# Command-line arguments can be used to specify memory and number of processors:
#   -m<mwords> -n<nprocs>
#
# BEB is defined for neutral and singly-positive targets but only neutrals are
#   handled by this software.  
# Latest change: 18 August 2016
#
# DISCLAIMER:  
#   Certain commercial software is identified in order to specify procedures
#   completely.  In no case does such identification imply recommendation or
#   endorsement by the National Institute of Standards and Technology, nor does
#   it imply that the material or equipment identified is necessarily the best
#   available for the purpose.
#
# Software citation: 
#   BEB automation software, Karl K. Irikura, National Institute of Standards
#   and Technology (2016). 
#
import sys
import os
fn_cc = froot + '_cc.gjf'
fgjf = open( fn_cc, 'w' )
print( 'created file ', fn_cc )
basis = specify_basis( elem, 5 )
directive = '# CCSD(T)/GEN geom=check scf=xqc'
if need_ecp( basis ):
    directive += ' pseudo=read\n'
else:
    directive += '\n'
comment = '\nBEB step ' + str(step) + ': CCSD(T) for neutral ' + froot + '\n\n'
fgjf.write( header + directive + comment )
fgjf.write( coords[0] + '\n' + basis )
#
# sub-step B: vertical single ionization
c2 = comment.replace( 'neutral', 'singly ionized' )
dir2 = directive.replace( 'GEN', 'CHKBAS' )
dir2 = dir2.replace( 'pseudo=read', '' )
fgjf.write( '--Link1--\n' + header + dir2 + c2 )
charge += 1
mult += 1       # consider two spin states if the neutral is a radical
fgjf.write( "%d %d\n\n" % (charge, mult) )
if mult > 2:
    # the neutral was not a singlet; consider lower-spin ion
    c3 = c2.replace( 'singly', 'low-spin singly' )
    fgjf.write( '--Link1--\n' + header + dir2 + c3 )
    fgjf.write( "%d %d\n\n" % (charge, mult-2) )
fgjf.close()
#
# create input files for double ionization threshold (step 6)
# provide two alternatives:
#   EPT/big-basis of cation (fairly cheap)
#   CCSD(T)/big-basis of low- and high-spin dications (most expensive step in whole series)
#
step = '6A'
fn_ept2 = froot + '_ept2.gjf'
fgjf = open(fn_ept2, 'w')
print( 'created file ', fn_ept2 )
c2 = '\nBEB step ' + str(step) + ': EPT for IE of ' + froot + ' cation \n\n'
fgjf.write( header + directive_ept + c2 )
fgjf.write( "%d %d\n\n" % (charge, mult) )
fgjf.write( basis )
# 'mult' is already one higher than the input neutral value; adjust 'nbet' to match
nbet -= 1;
fgjf.write( "%d %d %d %d\n\n" % (nalp, nalp, nbet, nbet) )
if mult > 2:
    # the neutral was not a singlet; also consider low-spin cation reference
    c2 = c2.replace( 'cation ', 'low-spin cation ' )
    fgjf.write( '--Link1--\n' + header + directive_ept + c2 )
    fgjf.write( "%d %d\n\n" % (charge, mult-2) )
    fgjf.write( basis )
    # adjust 'nalp' and 'nbet' for low spin
    nbet += 1
    nalp -= 1
    fgjf.write( "%d %d %d %d\n\n" % (nalp, nalp, nbet, nbet) )
fgjf.close()
#
step = '6B'
fn_dbl = froot + '_dbl.gjf'
fgjf = open(fn_dbl, 'w')
print( 'created file ', fn_dbl )
c2 = '\nBEB step ' + str(step) + ': CCSD(T) for doubly ionized ' + froot + ' (high-spin)\n\n'
charge += 1
mult += 1
dir2 = dir2.replace( 'CCSD(T)', 'CCSD(T,maxcyc=100)' )  # dications can be hard to converge
fgjf.write( header + dir2 + c2 )
fgjf.write( "%d %d\n\n" % (charge, mult) )
# consider lower-spin ion
c3 = c2.replace( 'high-spin', 'low-spin' )
fgjf.write( '--Link1--\n' + header + dir2 + c3 )
fgjf.write( "%d %d\n\n" % (charge, mult-2) )
fgjf.close()

