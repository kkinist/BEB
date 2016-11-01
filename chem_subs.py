# Routines for general quantum chemistry (no particular software package)
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
eV_per_hartree = 27.21138602  # from NIST website 10/25/2016
au_kjmol = au_joule * avogadro / 1000   # hartree expressed in kJ/mol
ev_wavenumber = au_wavenumber / eV_per_hartree      # eV expressed in cm**-1
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
    # if 'atno' is a stoichiometric dict of {'el' : number}, then return the sum for
    #   the whole molecule
    ncore = 0
    if type(atno) == str:
        # convert symbol to Z value
        atno = elz( atno )
    if type(atno) == dict:
        for el, natom in atno.items():
            ncore += n_core(el) * natom
        return ncore
    core = {
		# these are the maximum atomic numbers-1 (Z-1) that have
		#   the given number of core elecrons
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
def read_regex( regex, fhandl, idx=1 ):
    # Return something from a line matchine a regular expression.
    #   First arg is the regular expression; idx is the match-group
    #	to return.  Return a list of values from all matching lines. 
    fhandl.seek(0)
    matches = []
    regx = re.compile( regex )
    for line in fhandl:
        mch = regx.search( line )
        if mch:
            matches.append( mch.group( idx ) )
    return matches
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
def max_not_exceed(bigser, target):
    # args are: (1) a pandas Series
    #           (2) a target value
    # return the largest value in 'bigser' that does not exceed 'target'
    # This is useful for matching up line numbers.
    smaller = bigser[bigser <= target]
    return smaller.max()
##
def hartree_eV(energy, direction='to_eV', multiplier=1):
    # convert from hartree to eV or the reverse (if direction == 'from_eV')
    # use multiplier = -1 to change the sign, too
    if direction == 'to_eV':
        return multiplier * energy * eV_per_hartree 
    elif direction == 'from_eV':
        return multiplier * energy / eV_per_hartree
    else:
        # illegal direction
        return 'unrecognized direction = {:s} in routine hartree_eV'.format(direction)
##
