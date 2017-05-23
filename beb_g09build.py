#!/usr/bin/python3
#
# Create Gaussian09 input files needed for BEB calculation on a neutral target.
# Basic method:  G. E. Scott and K. K. Irikura, Surf. Interf. Anal. 37, 973-977 (2005)
#     (DOI: 10.1002/sia.2091)
#
# BEB is defined for neutral and singly-positive targets but only neutrals are
#   handled by this software.  
#
############
# Required input file format (whitespace-delimited): 
#   line 1: charge, spin-multiplicity [integers]
#   next: element symbol, cartesian triple [one atom per line]
# 
# Command-line arguments can be used to specify memory and number of processors:
#   beb_g09build.py <input file name> -m<Mwords> -n<nprocs>
#   (1 Mword = 8 Mbyte)
#
# The resulting Gaussian input files (*.gjf) can be run using the script
#   "run_beb.sh".  This requires the additional script, "nimag.pl", to verify
#   that the geometry optimization has succeeded. 
#
#   Only the job "*_bu.gjf" is required; the others provide more refined
#   molecular-orbital parameters. 
#
# The Gaussian output files are then processed using the script "beb_g09parse.py". 
#   That produces a *.bun file of molecular-orbital parameters.  The BUN file is
#   input to the script "beb_tbl.pl", which produces total ionization cross sections
#   using the BEB model. 
############
#
# Software citation: 
#   BEB automation software, Karl K. Irikura, National Institute of Standards
#   and Technology (2017). 
#
# DISCLAIMER:  
#   Certain commercial software is identified in order to specify procedures
#   completely.  In no case does such identification imply recommendation or
#   endorsement by the National Institute of Standards and Technology, nor does
#   it imply that the material or equipment identified is necessarily the best
#   available for the purpose.
#
# Input files as examples and for testing:
#   ch2.inp (triplet, non-degen)
#   ch3.inp (doublet, degen)
#   ch4.inp (singlet, degenerate orbitals)
#   gacl.inp (singlet, mixed n>2 and/or ECP)
#   h2o.inp (singlet, non-degenerate orbitals)
#   h2s.inp (singlet, ECP or special n>2)
#   no2.inp (doublet, non-degen)
#   o2.inp  (triplet, degen)
#   ZnEt2.inp (singlet, ECP and special n>2)
#   BiFClBrIAt.inp (singlet; ECP; special n>2; heavy atoms)
#   cef4.inp (singlet; f-element; symmetry)
#   thf3.inp (doublet; 5f; symmetry)
#   Li.inp (atom with one valence electron--pathological case with manual workaround)
#       --suggested by Prof. K. N. Joshipura (Sardar Patel Univ.)
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
        # this is where 'light', 'medium', and 'heavy' atoms are defined
        if z < 11:
            lowz.append( z )
        elif z < 37:
            midz.append( z )
        elif z < 104:
            highz.append( z )
        else:
            sys.exit('I cannot cope with elements heavier than Lr because of ' +
                'a lack of all-electron basis sets. Sorry.')
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
    # second arg is the computational step number:
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
    # also return:
    #   list of strings describing the basis set ([El, basis] for each element)
    #
    # the next three lists show the basis set to be used for each step number
    light = ['', '6-31G(d)', '6-311G(d,p)', '6-311G(d,p)', '6-311+G(d,p)', '6-311+G(d,p)', 
        'cc-pVTZ', 'cc-pVTZ', 'cc-pVTZ', 'cc-pVTZ', 'cc-pVTZ']
    middle = ['', '6-31G(d)', '6-311G(d,p)', 'SDDall', '6-311+G(d,p)', '6-311+G(d,p)', 
        'cc-pVTZ', 'cc-pVTZ', 'cc-pVTZ', 'cc-pVTZ', 'cc-pVTZ']
    heavy = ['', 'SDD', 'assorted', 'SDD', 'SDD', 'SDD',
        'SDD', 'SDD', 'SDD', 'SDD', 'SDD']
    # main-group pseudopotentials need d-polarization functions added (from 6-311G**)
    dpol = { 'Na': 0.175, 'Mg': 0.175,
        'Al': 0.325, 'Si': 0.45, 'P': 0.55, 'S': 0.65, 'Cl': 0.75, 'Ar': 0.85,
        'K': 0.229, 'Ca': 0.26, 'Ga': 0.169, 'Ge': 0.228, 'As': 0.264,
        'Se': 0.305, 'Br': 0.451, 'Kr': 0.395, 'I': 0.302 }
    # the 3p-block elements need +d added to the Dunning basis sets
    plusd = [ 'Al', 'Si', 'P', 'S', 'Cl', 'Ar' ]
    # 'basdir' is the location of the special basis set files, "*.gbs"
    basdir = ''
    basisstring = ''
    ecpstring = ''
    basisdescr = {}
    lowz, midz, highz = sort_atoms( atoms, 'symbol' )
    if len( lowz ) > 0:
        basisstring += " ".join( lowz ) + ' 0\n'
        basisstring += light[step]
        basisstring += '\n****\n'
        for el in lowz:
            basisdescr[el] = light[step]
    if len( midz ) > 0:
        for el in midz:
            basisdescr[el] = middle[step]
        if re.search( 'cc-pVTZ', middle[step] ):
            # do we need cc-pVTZ, cc-pV(T+d)Z, or both?
            needd = list( set(plusd) & set(midz) )
            noneed = list( set(midz) - set(needd) )
            if len(needd):
                # cc-pV(T+d)Z basis sets will be read from a file at run time
                basfile = basdir + 't_plusd.gbs'
                basisstring += '@{:s}\n'.format(basfile)
                if not os.path.exists(basfile):
                    # basis-set file is not on this machine; OK if it's found at run time
                    print('*** Warning: basis-set file not found: {:s}\n'.format(basfile))
                for el in needd:
                    basisdescr[el] = basisdescr[el].replace('T', '(T+d)')
            if len(noneed):
                # plain cc-pVTZ basis sets
                basisstring += " ".join(noneed) + ' 0\n'
                basisstring += middle[step]
                basisstring += '\n****\n'
        else:
            # not a cc-pVTZ type of basis
            basisstring += " ".join( midz ) + ' 0\n'
            basisstring += middle[step]
            if re.search( 'SDD', middle[step] ):
                # need to specify pseudopotential
                ecpstring += " ".join( midz ) + ' 0\n'
                ecpstring += 'SDDall\n'
                # may need to add d-polarization function to valence basis
                for el in midz:
                    if el in dpol:
                        # add d-polarization for each element separately
                        basisstring += '\n****\n'
                        basisstring += '{:s} 0'.format(el)
                        basisstring += '\n D   1 1.00       0.000000000000\n'
                        basisstring += '      {:16.10e}  0.1000000000D+01'.format(dpol[el])
                        basisdescr[el] += '+d'
            basisstring += '\n****\n'
    if len( highz ) > 0:
        # polarization functions available only for some elements; put one per line
        isecp = re.search( 'SDD', heavy[step] )
        isspecial = re.search( 'assorted', heavy[step] )
        if isspecial:
            # the basis sets will be read at run time
            basfile = basdir + 'dzp_dkh.gbs'
            for basfile in ['dzp_dkh.gbs', 'sarc_dkh.gbs', 'ano_rcc_dz.gbs']:
                basisstring += '@{:s}{:s}\n'.format(basdir, basfile)
                if not os.path.exists(basfile):
                    # basis-set file is not on this machine; OK if it's found at run time
                    print('*** Warning: basis-set file not found: {:s}\n'.format(basfile))
        for el in highz:
            if isspecial:
                # check atomic number to learn basis set
                z = elz(el)
                if z in [87, 88]:
                    # this is my DZ truncation of the ANO-RCC sets
                    basisdescr[el] = 'ANO-RCC-DZ'
                elif z in list(range(37,58)) + list(range(72, 87)):
                    basisdescr[el] = 'DZP-DKH'
                else:
                    # hopefully this covers all other cases encountered here
                    basisdescr[el] = 'SARC-DKH'
            else:
                # not an 'assorted' basis set
                basisdescr[el] = heavy[step]
                basisstring += el + ' 0\n'
                basisstring += heavy[step]
                if el in dpol and isecp:
                    # d-polarization function is available for this element--use it
                    basisstring += '\n D   1 1.00       0.000000000000\n'
                    basisstring += '      %16.10e  0.1000000000D+01' % dpol[el]
                    basisdescr[el] += '+d'
                basisstring += '\n****\n'
        if isecp:
            # need to specify pseudopotential
            ecpstring += " ".join( highz ) + ' 0\n'
            ecpstring += 'SDD\n'
    basisstring += '\n'
    if ( len(ecpstring) > 0 ):
        basisstring += ecpstring + '\n'
    if (step in [4,5]):
        # EPT calculation needs to know if any core orbitals are
        #   replaced by ECPs (aka pseudopotentials, PPs); create dict
        nopp = {}
        for el in atoms:
            nopp[el] = 0
        for el in midz:
            if re.search('SDD', middle[step]):
                nopp[el] = n_sdd(el) // 2
        for el in highz:
            if re.search('SDD', heavy[step]):
                nopp[el] = n_sdd(el) // 2
        return basisstring, nopp, basisdescr
    return basisstring, basisdescr
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
    return fgjf, fname
##
def update_log():
    # update log file (descriptions of calculations) using global variables
    nparen = 1  # maximum number of parenthetical arguments to display
    task = ['', 'Optimize geometry to a local minimum (no imaginary vibrational frequencies)',
        'All-electron calculation of binding energies and kinetic energies',
        'Orbital kinetic energies with effective core potentials (pseudopotentials)',
        'Correlated binding energies for outer orbitals',
        'Doubly vertical ionization energy of cation',
        'Correlated total ground-state energy',
        'Correlated total energy of high-spin, vertical cation',
        'Correlated total energy of low-spin, vertical cation',
        'Correlated total energy of high-spin, doubly vertical dication',
        'Correlated total energy of low-spin, doubly vertical dication' ]
    s = '{:s} (step {:d})'.format(fname, step)
    flog.write('\n{:s}\n'.format(s))
    s = '-' * len(s)
    flog.write('{:s}\n'.format(s))
    s = task[step]
    flog.write('{:s}\n'.format(s))
    # extract theoretical method from 'directive'
    regtheory = re.compile(r'#\s*(.+?)/')
    reg6D = re.compile(r'\s6D')
    regparen = re.compile(r'\(.*\)')
    regdk = re.compile(r'int\S+(DKH\S*)', re.IGNORECASE)
    m = regtheory.match(directive)
    if m:
        theory = m.group(1)
        # if there are parenthetical arguments, only display the first 'nparen' of them
        m = regparen.search(theory)
        if m:
            # may not display all these sub-commands
            scmd = m.group(0).strip('()').split(',')
            n = min(nparen, len(scmd))
            theory = regparen.sub('(' + ','.join(scmd[0:n]) + ')', theory)
        # check for Douglas-Kroll-Hess hamiltonian
        m = regdk.search(directive)
        if m:
            theory += '-{:s}'.format(m.group(1).upper())
        s = 'Theoretical method: {:s}'.format(theory)
    else:
        # failed to read theoretical method; print directive instead
        s = '{:s}\n'.format(directive)
    flog.write('{:s}\n'.format(s))
    # report basis sets 
    s = 'Elem\tBasis'
    flog.write('{:s}\n'.format(s))
    for b in basdescr:
        s = '{:s}\t{:s}'.format(b, basdescr[b])
        flog.write('{:s}\n'.format(s))
    if reg6D.search(directive):
        s = 'd-functions used in sets of 6'
        flog.write('{:s}\n'.format(s))
    return
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
nprocs = 1  # default value
mem = 1000   # default value (Mwords)
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
# remove trailing blank lines, then add a single blank line
regblank = re.compile(r'^\s*$')
while regblank.match(coords[-1]):
    coords.pop()
coords.append('\n')
elem = get_elements( coords )
lowz, midz, hiz = sort_atoms( elem )
(charge, mult) = ( int(x) for x in coords[0].split() )
if ( charge != 0):
    print( 'charge =', charge )
    sys.exit( 'The target molecule must be electrically neutral.' )
#
# create log file describing each computation
fname_log = '{:s}_build.log'.format(froot)
flog = open(fname_log, 'w')
flog.write('Descriptions of the Gaussian calculations are below.\n')
flog.write('See corresponding *.txt files for descriptions of ' \
    'non-standard basis sets.\n')
#
# filename extensions for each computational job step
filext = ['', 'opt', 'bu', 'bupp', 'ept1', 'ept2',
    'cc', 'cc1hi', 'cc1lo', 'cc2hi', 'cc2lo']
#
# create input file for geometry optimization and frequencies (step 1)
#
step = 1
fgjf, fname = newfile(froot, filext[step])
header =  '%chk={:s}.chk\n%nprocs={:d}\n%mem={:d}mw\n'.format(froot.replace(",", "-"),
    nprocs, mem)
directive = '# B3LYP/GEN gfinput OPT FREQ scf=xqc'
comment = '\nBEB step {:d}: geom and freqs for {:s}\n\n'.format(step, froot)
basis, basdescr = specify_basis( elem, step )
if re.search('SDD', basis):
    # must read pseudopotential
    directive += ' pseudo=read'
if need_6D(basis):
    directive += ' 6D\n'
else:
    directive += '\n'
fgjf.write( header + directive + comment )
fgjf.write( "".join( coords ) )
fgjf.write( basis )
fgjf.close()
update_log()    # update logfile (description of calculations) 
#
# create all-electron input file for Hartree-Fock B + U (step 2)
# [include Mulliken population analysis for identifying orbitals with n>2]
# Only this step is absolutely required.
#
step = 2
fgjf, fname = newfile(froot, filext[step])
directive = '# HF/GEN gfinput geom=check IOP(6/81=3) scf=xqc pop(all,thresh=1)'
# if there are heavy elements, include scalar relativistic effects via DKH2
if len(hiz) > 0:
    # there are heavy elements; use DKH2 to get "better" BU values
    directive += ' int=DKH2\n'
else:
    directive += '\n'
comment = '\nBEB step {:d}: Hartree-Fock B and U for {:s}\n\n'.format(step, froot)
fgjf.write( header + directive + comment )
basis, basdescr = specify_basis( elem, step )
fgjf.write( coords[0] + '\n' + basis )  # coords[0] is the "charge mult" line
fgjf.close()
update_log()
#
# create pseudopotential input file for Hartree-Fock B + U (step 3), if needed
# [include Mulliken population analysis for identifying orbitals with n>2]
#
step = 3
if ( len(midz) + len(hiz) ) > 0:
    # there are atoms with n>3: create PP files for B+U calc.
    fgjf, fname = newfile(froot, filext[step])
    directive = '# HF/GEN gfinput geom=check IOP(6/81=3) pseudo=read scf=xqc pop(all,thresh=1)\n'
    comment = '\nBEB step {:d}: pseudopotential B and U for {:s}\n\n'.format(step, froot)
    fgjf.write( header + directive + comment )
    basis, basdescr = specify_basis( elem, step )
    fgjf.write( coords[0] + '\n' + basis )
    fgjf.close()
    update_log()
#
# create input file for correlated (EPT) binding energies B (step 4)
# [include Mulliken population analysis to sort orbitals consistently
#  with previous calculations]
#
step = 4
fgjf, fname = newfile(froot, filext[step])
basis, nopp, basdescr = specify_basis( elem, step )
directive = '# EPT(OVGF+P3,ReadOrbitals)/GEN gfinput geom=check guess=check scf=xqc pop(all,thresh=1)'
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
    ncore += max( n_core(el, 'g09'), 2 * nopp[el] )
ncore //= 2
nalp -= ncore   # number of valence alpha electrons
nbet -= ncore
if mult > 1:
    fgjf.write('1 {:d} 1 {:d}\n\n'.format(nalp, nbet))
else:
    fgjf.write('1 {:d} {:d} {:d}\n\n'.format(nalp, nbet, nbet))
fgjf.close()
directive_ept = directive
update_log()
#
# For ions, assume 'high' and 'low' spin value close to neutral.
# [0] is for even spin multiplicities (i.e., odd number of electrons)
# [1] is for odd spin multiplicities (i.e., even number of electrons)
ion_mult = [{}] * 2
if mult > 1:
    ion_mult[0] = { 'hi' : mult+1, 'lo' : mult-1 }  # for cation
else:
    ion_mult[0] = { 'hi' : 4, 'lo' : 2 }
if mult > 2:
    ion_mult[1] = { 'hi' : mult, 'lo' : mult-2 }  # for dication
else:
    ion_mult[1] = { 'hi' : mult+2, 'lo' : mult }
#
# create input file for EPT VIE2 (step 5)
#
# reference cation multiplicity must be between dication options
ref_mult = (ion_mult[1]['hi'] + ion_mult[1]['lo']) // 2
step = 5
fgjf, fname = newfile(froot, filext[step])
basis, nopp, basdescr = specify_basis( elem, step )
comment = '\nBEB step {:d}: cation EPT for VIE2 for {:s}\n\n'.format(step, froot)
fgjf.write( header + directive_ept + comment )
fgjf.write('1 {:d}\n'.format(ref_mult) + '\n' + basis )
# figure out the orbital window to get the HOMO(s)
nalp = (nel-1 + ref_mult-1) // 2     # number of alpha electrons in cation
nbet = nel-1 - nalp
nalp -= ncore   # number of valence alpha electrons
nbet -= ncore
fgjf.write('{:d} {:d} {:d} {:d}\n\n'.format(nalp, nalp, nbet, nbet))
fgjf.close()
update_log()
#
# create input file for neutral CCSD(T) (step 6)
#
step = 6
fgjf, fname = newfile(froot, filext[step])
basis, basdescr = specify_basis( elem, step )
directive = '# CCSD(T,maxcyc=100)/GEN gfinput geom=check scf=xqc'
if need_ecp( basis ):
    directive += ' pseudo=read\n'
else:
    directive += '\n'
comment = '\nBEB step {:d}: CCSD(T) for neutral {:s}\n\n'.format(step, froot)
fgjf.write( header + directive + comment )
fgjf.write( coords[0] + '\n' + basis )
fgjf.close()
update_log()
#
# create input file for high-spin cation CCSD(T) (step 7)
#
step = 7
fgjf, fname = newfile(froot, filext[step])
basis, basdescr = specify_basis( elem, step )
comment = '\nBEB step {:d}: high-spin CCSD(T) for cation of {:s}\n\n'.format(step, froot)
directive_ch = directive.rstrip() + ' guess=check\n'  # avoid excited states
fgjf.write( header + directive_ch + comment )
fgjf.write('1 {:d}\n\n'.format(ion_mult[0]['hi']) + basis)
fgjf.close()
update_log()
#
# create input file for low-spin cation CCSD(T) (step 8)
#
step = 8
fgjf, fname = newfile(froot, filext[step])
basis, basdescr = specify_basis( elem, step )
comment = '\nBEB step {:d}: low-spin CCSD(T) for cation of {:s}\n\n'.format(step, froot)
# not using 'guess=check' here:  it gave a bad result for CS molecule
fgjf.write( header + directive + comment )
fgjf.write('1 {:d}\n\n'.format(ion_mult[0]['lo']) + basis)
fgjf.close()
update_log()
#
# create input file for high-spin dication CCSD(T) (step 9)
#
step = 9
fgjf, fname = newfile(froot, filext[step])
basis, basdescr = specify_basis( elem, step )
comment = '\nBEB step {:d}: high-spin CCSD(T) for dication of {:s}\n\n'.format(step, froot)
fgjf.write( header + directive_ch + comment )
fgjf.write('2 {:d}\n\n'.format(ion_mult[1]['hi']) + basis)
fgjf.close()
update_log()
#
# create input file for low-spin dication CCSD(T) (step 10)
#
step = 10
fgjf, fname = newfile(froot, filext[step])
basis, basdescr = specify_basis( elem, step )
comment = '\nBEB step {:d}: low-spin CCSD(T) for dication of {:s}\n\n'.format(step, froot)
fgjf.write( header + directive + comment )  # no 'guess=check' for low-spin
fgjf.write('2 {:d}\n\n'.format(ion_mult[1]['lo']) + basis)
fgjf.close()
update_log()
#
flog.close()
print('Descriptions of the calculations are in the file "{:s}".'.format(fname_log))
if nprocs == 1:
    print('Your calculations will use 1 processor and {:d} Mw of memory.'.format(mem))
else:
    print('Your calculations will use {:d} processors and {:d} Mw of memory.'.format(nprocs, mem))
