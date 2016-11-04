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
# NEED MORE TESTS: SULFUR, ZINC, ETC.
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
#from beb_subs import *
sys.path.append('/media/sf_bin3')
from g09_subs3 import *
from chem_subs import *
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
    # first arg is coordinates, to be parsed for element symbols
    # second arg is the computational job number:
    #   1 = geometry optimization and vibrational frequencies
    #   2 = binding and kinetic energies -- all-electron
    #   3 = as above, but with ECP on most atoms (if needed)
    #   4 = correlated (EPT/OVGF) valence binding energies
    #       also return dict of number of orbitals replaced by ECPs: nopp{symbol}
    #   5 = correlated (EPT/OVGF) IE of cation for VIE2 value
    #       also return dict of number of orbitals replaced by ECPs: nopp{symbol}
    #   6 = CCSD(T) for neutral (target)
    #   7 = CCSD(T) for ion for VIE1 -- high-spin
    #   8 = CCSD(T) for ion for VIE1 -- low-spin
    #   9 = CCSD(T) for dication for VIE2 -- high-spin
    #  10 = CCSD(T) for dication for VIE2 -- low-spin
    #
    # the next three lists show the basis set to be used for each step number
    light = ['', '6-31G(d)', '6-311G(d,p)', '6-311G(d,p)', '6-311+G(d,p)', '6-311+G(d,p)', 
        'cc-pVTZ', 'cc-pVTZ', 'cc-pVTZ', 'cc-pVTZ', 'cc-pVTZ']
    middle = ['', '6-31G(d)', '6-311G(d,p)', 'SDDall', '6-311+G(d,p)', '6-311+G(d,p)',
        'cc-pVTZ', 'cc-pVTZ', 'cc-pVTZ', 'cc-pVTZ', 'cc-pVTZ']
    heavy = ['', '3-21G*', '3-21G*', 'SDD', 'SDD', 'SDD',
        'SDD', 'SDD', 'SDD', 'SDD', 'SDD']
    # main-group pseudopotentials need d-polarization functions added (from 6-311G**)
    dpol = { 'Na': 0.175, 'Mg': 0.175,
        'Al': 0.325, 'Si': 0.45, 'P': 0.55, 'S': 0.65, 'Cl': 0.75, 'Ar': 0.85,
        'K': 0.229, 'Ca': 0.26, 'Ga': 0.169, 'Ge': 0.228, 'As': 0.264,
        'Se': 0.305, 'Br': 0.451, 'Kr': 0.395, 'I': 0.302 }
    # the 3p-block elements need +d added to the Dunning basis sets
    plusd = [ 'Al', 'Si', 'P', 'S', 'Cl', 'Ar' ]
    basisstring = ''
    ecpstring = ''
    lowz, midz, highz = sort_atoms( atoms, 'symbol' )
    if len( lowz ) > 0:
        basisstring += " ".join( lowz ) + ' 0\n'
        basisstring += light[step]
        basisstring += '\n****\n'
    if len( midz ) > 0:
        if re.search( 'cc-pVTZ', middle[step] ):
            # do we need cc-pVTZ, cc-pV(T+d)Z, or both?
            # add tight d-function to old cc basis set
            needd = list( set(plusd) & set(midz) )
            noneed = list( set(midz) - set(plusd) )
            if len( needd ):
                basisstring += '@/home/irikura/basis/t_plusd.gbs\n'
            if len( noneed ):
                basisstring += " ".join( noneed ) + ' 0 \n'
                basisstring += middle[step] + '\n****\n'
        else:
            basisstring += " ".join( midz ) + ' 0\n'
            basisstring += middle[step]
            if re.search( 'SDD', middle[step] ):
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
        isecp = re.search( 'SDD', heavy[step] )
        for el in highz:
            basisstring += el + ' 0\n'
            basisstring += heavy[step]
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
    if (step in [4,5]):
        # EPT calculation needs to know how many orbitals are replaced
        #   by ECPs (aka pseudopotentials, PPs); create dict
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
##
def newfile(froot, filext):
    # open a file named using the two arguments
    # return the filehandle
    fname = '{:s}_{:s}.gjf'.format(froot, filext)
    fgjf = open( fname, 'w' )
    print( 'creating file', fname )
    return fgjf
#
# MAIN
#
# For the Gaussian09 program, create the following input files. Some molecules
# will not require all these computations. 
#
#
# read command-line input
if len(sys.argv) < 2:
    sys.exit( 'Usage:  beb_g09build.py <inputfile.inp> [-n<cpus>] [-m<Mwords>]' )
try:
    finp = open( sys.argv[1] )
except:
    # assume user forgot file suffix, probably .inp
    finp = open( sys.argv[1] + '.inp' )
froot = os.path.splitext( sys.argv[1] )[0]
nprocs = 8  # default value
mem = 200   # default value (Mwords)
for arg in sys.argv:
    # nprocs specification leads with '-n'
    m = re.match(r'-n(\d+)', arg)
    if m:
        nprocs = int(m.group(1))
    # memory specification leads with '-m' for Mwords
    m = re.match(r'-m(\d+)', arg)
    if m:
        mem = int(m.group(1))
# read input file
coords = finp.readlines()
finp.close()
elem = get_elements( coords )
(charge, mult) = ( int(x) for x in coords[0].split() )
if ( charge != 0):
    print( 'charge =', charge )
    sys.exit( 'The target molecule must be electrically neutral.' )
#
# filename extensions for each computational job step
filext = ['', 'opt', 'bu', 'bupp', 'ept1', 'ept2',
    'cc', 'cc1hi', 'cc1lo', 'cc2hi', 'cc2lo']
#
# create input file for geometry optimization and frequencies (step 1)
#
step = 1
fgjf = newfile(froot, filext[step])
header =  '%chk={:s}.chk\n%nprocs={:d}\n%mem={:d}mw\n'.format(froot.replace(",", "-"),
    nprocs, mem)
directive = '# B3LYP/GEN OPT FREQ scf=xqc'
comment = '\nBEB step {:d}: geom and freqs for {:s}\n\n'.format(step, froot)
basis = specify_basis( elem, step )
if need_6D(basis):
    directive += ' 6D\n'
else:
    directive += '\n'
fgjf.write( header + directive + comment )
fgjf.write( "".join( coords ) )
fgjf.write( basis )
fgjf.close()
#
# create all-electron input file for Hartree-Fock B + U (step 2)
# [include Mulliken population analysis]
#
step = 2
fgjf = newfile(froot, filext[step])
directive = '# HF/GEN geom=check IOP(6/81=3) scf=xqc pop=reg\n'
comment = '\nBEB step {:d}: Hartree-Fock B and U for {:s}\n\n'.format(step, froot)
fgjf.write( header + directive + comment )
basis = specify_basis( elem, step )
fgjf.write( coords[0] + '\n' + basis )  # coords[0] is the "charge mult" line
fgjf.close()
#
# create pseudopotential input file for Hartree-Fock B + U (step 3), if needed
#
step = 3
lowz, midz, hiz = sort_atoms( elem )
if ( len(midz) + len(hiz) ) > 0:
    # there are atoms with n>3: create PP files for B+U calc.
    fgjf = newfile(froot, filext[step])
    directive = '# HF/GEN geom=check IOP(6/81=3) pseudo=read scf=xqc\n'
    comment = '\nBEB step {:d}: pseudopotential B and U for {:s}\n\n'.format(step, froot)
    fgjf.write( header + directive + comment )
    basis = specify_basis( elem, step )
    fgjf.write( coords[0] + '\n' + basis )
    fgjf.close()
#
# create input file for correlated (EPT) binding energies B (step 4)
#
step = 4
fgjf = newfile(froot, filext[step])
basis, nopp = specify_basis( elem, step )
directive = '# EPT(OVGF+P3,ReadOrbitals)/GEN geom=check scf=xqc'
if need_ecp( basis ):
    directive += ' pseudo=read\n'
else:
    directive += '\n'
comment = '\nBEB step {:d}: EPT B for {:s}\n\n'.format(step, froot)
fgjf.write( header + directive + comment )
fgjf.write( coords[0] + '\n' + basis )
# figure out the orbital window to get all valence orbs
nel = sum( map( elz, elem ) )
nalp = (nel + mult - 1) // 2     # number of alpha electrons
nbet = nel - nalp
ncore = 0
for el in elem:
    ncore += max( n_core(el), nopp[el] )
ncore //= 2
nalp -= ncore   # number of valence alpha electrons
nbet -= ncore
if mult > 1:
    fgjf.write('1 {:d} 1 {:d}\n\n'.format(nalp, nbet))
    #fgjf.write( "1 %d 1 %d\n\n" % (nalp, nbet) )
else:
    fgjf.write('1 {:d} {:d} {:d}\n\n'.format(nalp, nbet, nbet))
    #fgjf.write( "1 %d %d %d\n\n" % (nalp, nbet, nbet) )
fgjf.close()
directive_ept = directive
#
# For ions, assume minimal values for 'high' and 'low' spin.
# [0] is for even spin multiplicities (i.e., odd number of electrons)
# [1] is for odd spin multiplicities (i.e., even number of electrons)
ion_mult = [{}] * 2
ion_mult[0] = { 'hi' : 4, 'lo' : 2 }
ion_mult[1] = { 'hi' : 3, 'lo' : 1 }
#
# create input file for EPT VIE2 (step 5)
#
i = (mult + 2) % 2  # for +2 ion
# cation multiplicity must be between dication options
ref_mult = (ion_mult[i]['hi'] + ion_mult[i]['lo']) // 2
step = 5
fgjf = newfile(froot, filext[step])
basis, nopp = specify_basis( elem, step )
comment = '\nBEB step {:d}: cation EPT for VIE2 for {:s}\n\n'.format(step, froot)
fgjf.write( header + directive_ept + comment )
fgjf.write('1 {:d}\n'.format(ref_mult) + '\n' + basis )
#fgjf.write( coords[0] + '\n' + basis )
# figure out the orbital window to get the HOMO(s)
nalp = (nel-1 + ref_mult-1) // 2     # number of alpha electrons in cation
nbet = nel-1 - nalp
nalp -= ncore   # number of valence alpha electrons
nbet -= ncore
fgjf.write('{:d} {:d} {:d} {:d}\n\n'.format(nalp, nalp, nbet, nbet))
fgjf.close()
#
# create input file for neutral CCSD(T) (step 6)
#
step = 6
fgjf = newfile(froot, filext[step])
basis = specify_basis( elem, step )
directive = '# CCSD(T)/GEN geom=check scf=xqc'
if need_ecp( basis ):
    directive += ' pseudo=read\n'
else:
    directive += '\n'
comment = '\nBEB step {:d}: CCSD(T) for neutral {:s}\n\n'.format(step, froot)
fgjf.write( header + directive + comment )
fgjf.write( coords[0] + '\n' + basis )
fgjf.close()
#
# create input file for high-spin cation CCSD(T) (step 7)
#
step = 7
fgjf = newfile(froot, filext[step])
basis = specify_basis( elem, step )
comment = '\nBEB step {:d}: high-spin CCSD(T) for cation of {:s}\n\n'.format(step, froot)
fgjf.write( header + directive + comment )
i = (mult + 1) % 2
fgjf.write('1 {:d}\n\n'.format(ion_mult[i]['hi']) + basis)
fgjf.close()
#
# create input file for low-spin cation CCSD(T) (step 8)
#
step = 8
fgjf = newfile(froot, filext[step])
basis = specify_basis( elem, step )
comment = '\nBEB step {:d}: low-spin CCSD(T) for cation of {:s}\n\n'.format(step, froot)
fgjf.write( header + directive + comment )
i = (mult + 1) % 2
fgjf.write('1 {:d}\n\n'.format(ion_mult[i]['lo']) + basis)
fgjf.close()
#
# create input file for high-spin dication CCSD(T) (step 9)
#
step = 9
fgjf = newfile(froot, filext[step])
basis = specify_basis( elem, step )
comment = '\nBEB step {:d}: high-spin CCSD(T) for dication of {:s}\n\n'.format(step, froot)
fgjf.write( header + directive + comment )
i = (mult + 2) % 2
fgjf.write('2 {:d}\n\n'.format(ion_mult[i]['hi']) + basis)
fgjf.close()
#
# create input file for low-spin dication CCSD(T) (step 10)
#
step = 10
fgjf = newfile(froot, filext[step])
basis = specify_basis( elem, step )
comment = '\nBEB step {:d}: low-spin CCSD(T) for dication of {:s}\n\n'.format(step, froot)
fgjf.write( header + directive + comment )
i = (mult + 2) % 2
fgjf.write('2 {:d}\n\n'.format(ion_mult[i]['lo']) + basis)
fgjf.close()

