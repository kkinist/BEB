# Routines for use in running BEB calculations.
# Karl Irikura 
#
import re
import string
#
amu_au = 1 / 5.4857990946E-4    # amu expressed in a.u. (viz., electron masses)
au_wavenumber = 219474.6313708  # hartree expressed in wavenumbers
au_joule = 4.35974434E-18       # hartree expressed in joule
avogadro = 6.02214129e23 # Avogadro constant
planck = 6.62606957e-34  # Planck constant in J.s
clight = 299792458.     # speed of light in m/s
bohr = 0.5291772109     # Bohr radius in angstrom
ev = 27.21138505        # Hartree energy expressed in eV
au_kjmol = au_joule * avogadro / 1000   # hartree expressed in kJ/mol
ev_wavenumber = au_wavenumber / ev      # eV expressed in cm**-1
#
def elz( ar ):
    # return atomic number given an elemental symbol, or
    # return elemental symbol given an atomic number 
    symb = [ 'n',
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
        'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
        'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
        'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
        'Cs', 'Ba',
            'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
                 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
               'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt',
                 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
        'Fr', 'Ra',
            'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
                 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
               'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt' ]
    if type( ar ) == int:
        return symb[ ar ]
    return symb.index( ar )
##
def n_core( atno ):
    # given Z value (or element symbol) return number of core electrons
    # also accept arg of type 'dict', as returned by function stoich2dict() in this file
    ncore = 0
    if type( atno ) == dict:
        for el in atno:
            ncore += atno[el] * n_core(el)
        return ncore
    if type( atno ) != int:
        # convert symbol to Z value
        atno = elz( atno )
    core = {
        3  :  2,
        11 : 10,
        19 : 18,
        31 : 28,
        37 : 36,
        49 : 46,
        55 : 54,
        81 : 78,
        87 : 86
    }
    for ki in sorted( core ):
        if atno >= ki:
            ncore = core[ki]
    return ncore
##
def stoich2dict( stoich ):
    # Starting from a Gaussian-style stoichiometry (e.g., 'C2H6Si'), return a
    #   hash with formula[element] = number
    formula = {}
    uc = string.ascii_uppercase
    lc = string.ascii_lowercase
    num = string.digits
    newel = ''
    newnum = ''
    for i in range( len(stoich) ):
        if stoich[i] in uc:
            # uppercase--new element symbol
            if len( newel ):
                # finish up the previous symbol
                formula[ newel ] = ( 1 if len(newnum) == 0 else int(newnum) )
            newel = stoich[i]
            newnum = ''
        elif stoich[i] in lc:
            # second letter of a symbol
            newel += stoich[i]
        elif stoich[i] in num:
            # one digit of an atom count
            newnum += stoich[i]
    # finish the last
    formula[ newel ] = ( 1 if len(newnum) == 0 else int(newnum) )
    return formula
##
def opt_success( fhandl ):
    # return True for a successful geometry optimization
    regx = re.compile( 'Optimization completed.' )
    for line in fhandl:
        m = regx.search( line )
        if m:
            # successful optimization
            return True
    return False
##
def get_arch( fhandl, stepnum=1 ):
    # read archive block; 2nd arg is which archive (1 = first, 2 = second, etc)
    # return list of major (double-delimited) fields
    fhandl.seek( 0 )	# rewind file
    gtext = fhandl.readlines()	# read whole file
    #archini = re.compile( r'^ 1[\|\\]1[\|\\]' )	# beginning of archive block
    #archend = re.compile( r'[\|\\][\|\\]@$' )
    archini = re.compile( r'^ ((1\|)|(1\\)){2}' )	# beginning of archive block
    archend = re.compile( r'((\|)|(\\)){2}@\s*$' )
    arch = ''
    acount = 0
    for p in gtext:
        mch = archini.match( p )
        if mch:
            acount += 1
            if acount < stepnum:
                continue
            # remove the leading space and trailing newline
            #arch = p.strip()
            arch = p[1:].rstrip('\n').rstrip('\r')
        elif len( arch ) > 0:
            # add to growing archive block
            #arch += p.strip()
            arch += p[1:].rstrip('\n').rstrip('\r');
            if archend.search( p ):
                # that was the last line of the archive block
                break
    if acount < stepnum:
        # request cannot be satisfied; error
        return None
    arch = arch.replace( '\\', '|' )	# replace backslashes with up-bars
    # split it using the delimiter '||'
    arch = arch.split( '||' )
    # 'arch' is now a list of archive fields
    return arch
##
def arch2nimag( arch ):
    # find and return the value of NIMAG in an archive list
    regx = re.compile( 'NImag=(\d+)' )
    for line in arch:
        mch = regx.search( line )
        if mch:
            return int( mch.group(1) )
    # return nonsense value of -1 if no valid value was found
    return -1
##
def read_stoich( fhandl ):
    # return stoichiometry, charge, and spin multiplicity as list of lists
    fhandl.seek(0)
    charge = 0
    mult = 1
    stoichcm = []
    regx = re.compile( r' Stoichiometry\s+([A-Z]\S+)' )
    for line in fhandl:
        mch = regx.search( line )
        if mch:
            stoichstr = mch.group(1)
            t = re.split( '[(]', stoichstr )
            stoich = t[0]
            if len( t ) > 1:
                # look for spin multiplicity (number followed by right-parenthesis)
                mch = re.search( '(\d+)\)', t[1] )
                if mch:
                    mult = int( mch.group(1) )
                else:
                    mult = 1
                # look for charge (number followed by sign)
                mch = re.search( '(\d+)([-+])', t[1] )
                if mch:
                    charge = int( mch.group(1) )
                    if mch.group(2) == '-':
                        charge = -charge
                else:
                    charge = 0
            stoichcm.append( [stoich, charge, mult] )
    return stoichcm
##
def read_regex( regex, fhandl, idx=1 ):
    # Return list of all matching fields #idx for specified regular expression
    fhandl.seek(0)
    matches = []
    regx = re.compile( regex )
    for line in fhandl:
        mch = regx.search( line )
        if mch:
            matches.append( mch.group( idx ) )
    return matches
##
def read_oeke( fhandl, virtual=0 ):
    # read orbital energies and kinetic energies from HF calculation (for BEB application)
    # use virtual=1 to include info for virtual orbitals
    fhandl.seek(0)
    bu = { 'alpha': [], 'beta': [] }
    # for lists 'alpha' and 'beta', index is (orbital no. - 1) and each element is in
    #  turn a list [B U <'O'|'V'>], with the third field only included if virtual==1
    regxstart = re.compile( 'Orbital energies and kinetic energies \((\w+)\):' )
    regxend = re.compile( 'Total kinetic energy from orbitals=' )
    regxocc = re.compile( 'O$' )
    regxdata = re.compile( '\s+(\d+)\s+.*\d\.\d' )
    inorb = 0 # flag
    for line in fhandl:
        if inorb:
            mch = regxend.search( line )
            if mch:
                # all done
                break
            # parse next line of the interesting output block
            # is it a data line?
            mch = regxdata.search( line )
            if mch:
                # yes, it's a data line
                iorb = int(mch.group(1))   # orbital number listed in file
                s = line.split()
                b = float(s[2])  # negative for bound orbitals
                u = float(s[3])
                # occupied or virtual?
                occ = regxocc.search( s[1] )
                lbl = 'O' if occ else 'V'
                if virtual:
                    # virtual orbital data are requested, so include O|V label
                    bu[spin].append( [ b, u, lbl ] )
                elif occ:
                    bu[spin].append( [ b, u ] )
                else:
                    continue
                if iorb != len( bu[spin] ):
                    msg = 'iorb = %d but length of bu[%s] = %d' % (iorb, spin, len(bu[spin]) )
                    sys.exit( msg )
        mch = regxstart.search( line )
        if mch:
            spin = mch.group(1)
            inorb = 1
    return bu
##
def spinname( m ):
    # given a spin multiplity (m = 2S+1), return the text name
    name = [ 'spinless', 'singlet', 'doublet', 'triplet', 'quartet', 'quintet', 'sextet',
        'septet', 'octet', 'nonet', 'decet', 'undecet', 'duodecet' ]
    m = int( m )
    if m in range(12):
        return name[m]
    else:
        return str(m) + '-tet'
##
def read_ept( fhandl ):
    # read OVGF and P3 results from EPT calculation
    # return two lists: labels  and list-of-lists with corresponding values
    fhandl.seek(0)
    searchstring = [ 'Koopmans theorem:',
        'Outer Valence Approximation:', 'Converged 3rd order P3 pole:',
        'Converged second order pole:' ]
    label = [ 'HF', 'OVGF', 'P3', 'P2' ]
    nmethod = len( label )
    energy = []
    regx = []
    regxonum = re.compile( 'Summary of results for (\w+)\s+spin-orbital\s+(\d+)\s' )
    regxdegen = re.compile( 'Orbitals\s+(\d+)\s+and\s+(\d+)\s+are degenerate.  Skipping orbital\s+(\d+)' )
    for i in range( len(searchstring) ):
        sstr = searchstring[i]
        if i > 0:
            # OVGF or P3 or P2: make $1 the energy (in Hartree) and $2 the pole strength
            regx.append( re.compile( sstr + r'\s+(.*)\sau.*(\d\.\d+)\s\(PS\)' ) )
        else:
            # Koopmans: make $1 the energy (in Hartree)
            regx.append( re.compile( sstr + r'\s+(.*)\sau\s' ) )
    for line in fhandl:
        for imeth in range( nmethod ):
            mch = regxonum.search( line )  # look for number and spin of (active) orbital
            if mch:
                spin = mch.group(1)
                iorb = int( mch.group(2) ) - 1  # decrement to mimic index
            mch = regx[imeth].search( line )
            if mch:
                t = [ spin, iorb, label[imeth], float( mch.group(1).replace('D','E') ) ]
                if imeth > 0:
                    # save the pole strength
                    t.append( float( mch.group(2) ) )
                else:
                    # Koopmans--just call 'pole strength' = 1
                    t.append( 1 )
                energy.append( t )
            mch = regxdegen.search( line )  # look for skipped orbitals
            if mch:
                iskip = int( mch.group(3) ) - 1  # index of skipped orbital
                ia = int( mch.group(1) ) - 1
                ib = int( mch.group(2) ) - 1
                if iskip == ia:
                    idup = ib
                else:
                    idup = ia
                # idup is the index of the orbital whose info will be copied
                for t in energy:
                    if t[1] == idup:
                        t2 = list( t )  # copy to new list
                        t2[1] = iskip
                        energy.append( t2 )
    # remove duplicates (inefficiently)
    for i in range( len(energy)-1 ):
        for j in reversed( range(i+1,len(energy)) ):
            if energy[i] == energy[j]:
                # remove the duplicate entry (j)
                energy.pop( j )
    return energy
##
