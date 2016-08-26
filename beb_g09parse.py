#!/usr/bin/python3
# Read G09 output files and create BUN file for BEB calculation
# Expected files are those created by beb_g09build.py
# Karl Irikura, last change 18 August 2016 
#
import sys
import os
import re
from beb_subs import *
##
def read_cc( fhandl ):
    # read multi-step CCSD(T) output for BEB calculation
    # return list of lists with lowest energies, [ ionlevel, spin, energy ]
    fgau.seek(0)  # rewind file
    regx = re.compile( '^ BEB step (5|6B): CCSD\(T\) for (\w+) ' )
    regxel = re.compile( 'NAE=\s+(\d+) NBE=\s+(\d+) ' )
    regxcc = re.compile( r'^ CCSD\(T\)=\s*([-]?\d\.\d+D[-+]\d+)\b' )
    ecc = {}  # key is ionlevel, value is [ spin, energy ]
    for line in fhandl:
        mch = regxel.search( line )
        if mch:
            # numbers of alpha and beta electrons
            nalp = int( mch.group(1) )
            nbet = int( mch.group(2) )
            spin = nalp - nbet + 1
        mch = regx.search( line )
        if mch:
            # which level of ionization
            ionlevel = mch.group(2)
            try:
                ecc[ionlevel]
            except:
                ecc[ionlevel] = []
        mch = regxcc.search( line )
        if mch:
            # CCSD(T) energy
            energy = float( mch.group(1).replace('D','E') )
            ecc[ionlevel].append( [ spin, energy ] )
    # find the lowest energy spin state for ions
    bestcc = {}
    for ioniz in ecc:
        # there may be more than one spin state for this ionization level
        # pick the spin state with lowest energy
        nstate = len( ecc[ioniz] )
        if nstate > 0:
            ilow = 0
            for istate in range(nstate):
                if ecc[ioniz][istate][1] < ecc[ioniz][ilow][1]:
                    ilow = istate
            bestcc[ioniz] = ecc[ioniz][ilow]
        else:
            # number of states calculated for this ionization level is < 1
            print( 'I was not expecting %d spin states for calculation on %s' % (nstate,ioniz) )
    return bestcc
##
def select_ept( eptlist ):
    # select highest-level EPT results with acceptable pole strength
    # return data structure of similar type as 'bu', 'buae', and 'bupp' in main program
    okpole = 0.75
    rank = { 'HF': 0, 'OVGF': 1, 'P2': 2, 'P3': 3 }
    na = nb = -1
    # count alpha and beta orbitals; 'eptlist' might not be sorted
    for x in eptlist:
        if x[0] == 'alpha' and x[1] > na:
            na = x[1]
        if x[0] == 'beta' and x[1] > nb:
            nb = x[1]
    na += 1
    nb += 1
    sel = { 'alpha': [ [] ]*na, 'beta': [ [] ]*nb }
    for line in eptlist:
        (spin, iorb, meth, b, ps) = line
        if len( sel[spin][iorb] ) == 0:
            # create new entry
            if ps >= okpole:
                sel[spin][iorb] = [ b, meth ]
            else:
                # PS too low--just record a zero
                sel[spin][iorb] = [ 0, meth ]
        else:
            # possibly replace value with higher-level one
            if ( rank[meth] > rank[ sel[spin][iorb][1] ] ):
                # better method
                if ps >= okpole:
                    # replace value
                    sel[spin][iorb] = [ b, meth ]
                elif sel[spin][iorb][0] == 0:
                    # method failed, but inferior method did too; replace method
                    sel[spin][iorb][1] = meth
    return sel
##
## MAIN
##
if len(sys.argv) < 2:
    sys.exit( 'Usage:  beb_g09parse.py <file rootname>\n\tE.g., beb_g09parse.py c2h5oh' )
mol = sys.argv[1]
ppcore = 0
# concatenate both ccsd(t) output files into one file
fncc = mol + '_cccat.out'
with open(fncc, 'w') as outfile:
    for suff in ['_cc', '_dbl']:
        fn = mol + suff + '.out'
        try:
            with open(fn) as infile:
                outfile.write(infile.read())
        except:
            # probably user gave bad command-line argument
            os.unlink(fncc)  # don't leave debris
            raise
# parse each type of output file in turn
for suff in [ '_geom', '_bu', '_bupp', '_ept', '_cccat', '_ept2' ]:
    fn = mol + suff + '.out'
    try:
        fgau = open( fn, 'r' )
    except:
        # hopefully this file was not needed
        continue
    print( 'Reading file ', fn )
    if suff == '_geom':
        # check that geometry optimization converged
        if opt_success( fgau ):
            print( '\tgeometry optimization converged' )
        else:
            sys.exit( '\tgeometry optimization FAILED\n' )
        # check that nimag = 0
        arch = get_arch( fgau, 2 )      # freq calc is 2nd step in this file
        nimag = arch2nimag( arch )
        if nimag > 0:
            msg = '\tgeometry is NOT a minimum: nimag = %d' % nimag
            sys.exit( msg )
        if nimag < 0:
            sys.exit( '\tno valid value of NIMAG was found\n' )
        print( '\tnimag = 0 is good' )
        # extract stoichiometry, charge, and spin multiplicity
        (stoich, charge, mult) = read_stoich( fgau )[0]
        elems = stoich2dict( stoich )
        ncore = n_core( elems ) / 2  # number of conventional core orbitals
        nalp = int( read_regex( r'NAE=\s+(\d+)', fgau )[0] )  # total alpha electrons
    if suff == '_bu' or suff == '_bupp':
        # extract HF orbital and kinetic energies
        bu = read_oeke( fgau )
        uhf = len( bu['beta'] )  # flag to indicate separate beta orbitals
        if uhf:
            print( '\t%d alpha and %d beta orbital energies ' % (len(bu['alpha']), len(bu['beta'])) )
        else:
            print( '\t%d closed-shell orbital energies ' % len(bu['alpha']) )
        # orbital energies are negative binding energies
        if suff == '_bu':
            # save all-electron BU in separate variable
            buae = bu
        else:
            # save ECP BU in separate variable
            bupp = bu
            ppcore = len(buae['alpha']) - len(bupp['alpha'])
            print( '\tPP replaced %d orbitals' % ppcore )
            # note any light-atom core orbitals remaining
            ppskip = ncore - ppcore
        for x in bu:
            # print B and U values in eV; key 'x' is 'alpha' or 'beta'
            if len( bu[x] ) < 1 :
                continue
            print( '\t\t', x, '[-B, U]/eV' )
            for y in bu[x]:
                # 'y' is the list [B, U] for one orbital; print in eV and rounded
                yev = [ round(t*ev, 2) for t in y ]
                print( '\t\t', yev )
    if suff == '_ept':
        # extract correlated orbital energies (negative of binding energies)
        nfc = int( read_regex( r'NFC=\s+(\d+)', fgau )[0] )
        print( '\t%d core orbitals frozen in EPT calculation' % nfc )
        nae = int( read_regex( r'NAE=\s+(\d+)', fgau )[0] )
        nfcept = nalp - nae + nfc # number of frozen/ECP orbitals not correlated in EPT calc.
        ept = read_ept( fgau )
        bcorr = select_ept( ept )
        for x in ept:
            print( '\t\t', x )
        print( '\tSelected values' )
        # print EPT-correlated B values in eV
        iemin = -999
        for x in bcorr:
            if len( bcorr[x] ) < 1 :
                continue
            print( '\t\t', x, 'B/eV' )
            for i, y in enumerate( bcorr[x] ):
                print( '\t\t %d\t%.2f (%s)' % (i, -ev*y[0], y[1]) )
                if y[0] > iemin:
                    iemin = y[0]
                    ie1_ept = [ y[0], x, y[1] ] # orbital energy, 'alpha' or 'beta', calc type
    if suff == '_ept2':
        # extract ionization energy of ion (add to IE1 to get IE2)
        noa = int( read_regex( r'NOA=\s+(\d+)', fgau )[0] )
        nob = int( read_regex( r'NOB=\s+(\d+)', fgau )[0] )
        nval = {
            'alpha': noa-1,
            'beta':  nob-1,
        }
        print( '\tNOA =', noa, ' and NOB =', nob )
        ept2 = read_ept( fgau )
        bcorr2 = select_ept( ept2 )
        for x in ept2:
            print( '\t\t', x )
        print( '\tSelected values' )
        # print EPT-correlated B values in eV
        print( '\t\t\t\teV' )
        iemin = 999
        for x in bcorr2:
            if len( bcorr2[x] ) < 1 :
                continue
            i = nval[x]
            y = bcorr2[x][i]
            ie = -ev * y[0]
            print( '\t\t%-5s\t%d\t%.2f (%s)' % (x, i, ie, y[1]) )
            if ie < iemin:
                iemin = ie
                ie2_ept = [ y[0], x, y[1] ] # orbital energy, 'alpha' or 'beta', calc. type
    if suff == '_cccat':
        # get CCSD(T) thresholds for single- and double-ionization
        cc = read_cc( fgau )
        for x in cc:
            print( '\t', x, cc[x] )
    fgau.close()
os.unlink(fncc)
# ionization thresholds from EPT calculations 
print( 'Vertical ionization thresholds from EPT calculations' )
print( '\tNeutral is', spinname(mult) )
mult1 = mult + 1
if mult > 1 and ie1_ept[1] == 'alpha':
    mult1 = mult - 1
iev = -ev * ie1_ept[0]
print( '\t+1 ion is %s, higher by %.2f eV (%s)' % (spinname(mult1), iev, ie1_ept[2]) )
mult2 = mult1 + 1
if mult1 > 1 and ie2_ept[1] == 'alpha':
    mult2 = mult1 - 1
iev = -ev * ie2_ept[0]
print( '\t+2 ion is %s, higher by %.2f eV (%s)' % (spinname(mult2), iev, ie2_ept[2]) )
try:
    ie1 = cc['singly'][1] - cc['neutral'][1]
    ie1spin = spinname( cc['singly'][0] )
    print( 'Vertical IE from CCSD(T) = %.2f eV to %s' % (ev*ie1, ie1spin) )
    ie1_calc = 'CCSD(T)'
except:
    print( 'No ionization threshold available from CCSD(T). Using EPT value instead.' )
    ie1 = -ie1_ept[0]
    ie1spin = spinname(mult1)
    print( '\t(%.2f eV to the %s)' % (ev*ie1, ie1spin) )
    ie1_calc = ie1_ept[2]
try:
    ie2 = cc['doubly'][1] - cc['neutral'][1]
    ie2spin = spinname( cc['doubly'][0] )
    print( 'Vertical IE2 from CCSD(T) = %.2f eV to %s' % (ev*ie2, ie2spin) )
    ie2_calc = 'CCSD(T)'
except:
    print( 'No double-ionization threshold available from CCSD(T).  Using EPT value.' )
    ie2 = ie1 - ie2_ept[0]
    ie2spin = spinname(mult2)
    print( '\t(%.2f eV to the %s)' % (ev*ie2, ie2spin) )
    ie2_calc = ie2_ept[2]
# now prepare the BUN file
fn = mol + '.dat'
fbun = open( fn, 'w' )
fbun.write( '# Results from automated beb_g09build.py and beb_g09parse.py\n' )
if charge:
    msg = '# Molecule is \'%s\' (%s %s%+d)\n' % (mol, spinname(mult), stoich, charge)
else:
    # don't specify zero charge
    msg = '# Molecule is \'%s\' (%s %s)\n' % (mol, spinname(mult), stoich)
fbun.write( msg )
msg = '# Double-ionization threshold = %.2f eV from %s (to %s)\n' % (ie2*ev, ie2_calc, ie2spin)
fbun.write( msg )
fbun.write( '#MO\tB/eV\tU/eV\tN\tQ\tDblIon\tSpecial\tRemarks\n' )
q = 1  # for BEB approximation to BED
special = 'none'  # for ECP-based computation on neutral molecule
nalp = len( buae['alpha'] )
nbet = len( buae['beta'] )
# check for degnerate orbitals--needed if using CCSD(T) for HOMO
# relies upon the orbitals being ordered
degen = { 'alpha': [], 'beta': [] }
# each element of degen[spin] is a list of degenerate orbitals (by index)
print( 'checking for degenerate orbitals' )
for spin in ['alpha', 'beta']:
    print( 'BUAE[%s]: ' % (spin) )
    for imo in range(len(buae[spin])):
        [b, u] = buae[spin][imo]
        print( '\t%f\t%f' % (b, u) )
        # check for energy match with previous orbital
        if ( imo > 0 and b == buae[spin][imo-1][0] and u == buae[spin][imo-1][1] ):
            # this is degenerate with the previous orbital
            degen[spin][-1].append( imo )
        else:
            # start a new list of degenerate orbitals
            degen[spin].append( [imo] )
    print( 'degeneracy lists:' )
    for t in degen[spin]:
        print( '\t', t )
# done checking for degen
#
# Prepare the list of orbital data
orblist = []
#
class Orbital:
    # one line of orbital data, to be printed to the BUN file
    def __init__(self, nmo=0, b=0, u=0, n=0, q=0, dblion=0, special=0, remark=0):
        self.mo = nmo   # orbital number (starting with 1), maybe followed by a or b
        self.b = b  # orbital binding energy in eV
        self.u = u  # orbital kinetic energy in eV
        self.n = n  # number of electrons in orbital
        self.q = q  # value of dipole constant (=1 for BEB approximation to BED)
        self.dblion = dblion    # Yes or No
        self.special = special  # none or ion
        self.remark = remark    # source of B value
    def outstring(self):
        line = '%s\t%.2f\t%.2f\t%d\t%d\t%s\t%s\t%s\n' % (self.mo, self.b,
            self.u, self.n, self.q, self.dblion, self.special, self.remark)
        return line
#
n = 2  # electrons per orbital
otyp = ''
if uhf:
    n = 1
for spin in ['alpha', 'beta']:
    nel = ( nalp if spin == 'alpha' else nbet )
    if uhf:
        otyp = ( 'a' if spin == 'alpha' else 'b' )
    for iorb in range( nel ):
        mo = str( iorb + 1 ) + otyp
        dblion = 'No'
        remark = ''
        if ppcore and iorb - ppcore >= ppskip:
            # use pseudopotential U for valence electrons
            u = bupp[spin][iorb-ppcore][1]
        else:
            # use all-electron B and U
            b = buae[spin][iorb][0]
            u = buae[spin][iorb][1]
        if iorb >= nfcept:
            # use valence binding energies from EPT calc, even if they are only Koopmans
            # (because the basis set is a little bigger than above)
            try:
                b = bcorr[spin][iorb-nfcept][0]
            except:
                print( 'Baddd' )
                print( 'spin =', spin, 'iorb =', iorb, 'nfcept =', nfcept )
                print( 'bcorr: ', bcorr )
            meth = bcorr[spin][iorb-nfcept][1]
            if meth != 'HF':
                remark = 'B from ' + meth
        if iorb == nel - 1:
            # use IE1 threshold as binding energy, if it's the right spin
            if spin == 'alpha':
                if cc['neutral'][0] == 1 or cc['singly'][0] < cc['neutral'][0]:
                    b = -ie1
                    remark = 'B from %s' % (ie1_calc)
            else:
                # spin beta
                if cc['singly'][0] > cc['neutral'][0]:
                    b = -ie1
                    remark = 'B from %s' % (ie1_calc)
        # change sign of binding energy for output
        b = -b
        if b > ie2:
            dblion = 'Yes'
        # convert b and u to eV
        b *= ev
        u *= ev
        orblist.append( Orbital(mo, b, u, n, q, dblion, special, remark) )
        #line = '%s\t%.2f\t%.2f\t%d\t%d\t%s\t%s\t%s\n' % (mo, b, u, n, q, dblion, special, remark)
        #fbun.write( line )
# apply degeneracy information
print( 'apply deneracies' )
for spin in ['alpha', 'beta']:
    for group in degen[spin]:
        print( '\tx: ', group )
        u = set( orblist[imo].u for imo in group )
        if len(u) > 1:
            # there should be exactly 1 value of U for each degenerate set or orbitals
            # so there is a problem with this set; don't change anything
            break
        b = set( orblist[imo].b for imo in group )
        if len(b) > 1:
            # B values are not all the same
            # accept only the data for the highest-index orbital, but:
            #   keep only the lowest orbital number
            #   update N to be the sum for all the degenerate orbitals
            iaccept = max( group )
            orbno = min( orblist[imo].mo for imo in group )
            ntot = sum( orblist[imo].n for imo in group )
            print( 'iaccept =', iaccept, 'and orbno =', orbno, ' and ntot =', ntot )
            for imo in group:
                if imo == iaccept:
                    # this orbital will be printed
                    orblist[imo].mo = orbno
                    orblist[imo].n = ntot
                else:
                    # this orbital will not be printed
                    orblist[imo].n = 0
# write the MO data to the file to finish up
for mo in orblist:
    if mo.n > 0:
        fbun.write( mo.outstring() )
fbun.close()

