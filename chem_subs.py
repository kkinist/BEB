# Routines for general quantum chemistry (no particular software package)
# Python3 and pandas
# Karl Irikura 
#
import re
import string
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
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
def elz(ar, choice=None):
    # return atomic number given an elemental symbol, or
    # return elemental symbol given an atomic number 
    # If 'choice' is specified as 'symbol' or 'Z', return that.
    # if 'ar' is a list, then return a corresponding list
    symb = ['n',
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
               'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt']
    if type(ar) == str and not re.match(r'^\d+$', ar):
        # this looks like an element symbol
        ar = ar.title()  # Title Case
        if choice == 'symbol':
            return ar
        else:
            return symb.index(ar)
    if type(ar) == list:
        # process a list of atoms
        vals = []
        for el in ar:
            vals.append(elz(el), choice)
        return vals
    # if we got here, the argument an atomic number
    Z = int(ar)
    if choice == 'Z':
        return Z
    else:
        return symb[Z]
##
def n_core(atno, code=None):
    # given Z value (or element symbol) return number of core electrons
    # if 'atno' is a stoichiometric dict of {'el' : number}, then return the sum for
    #   the whole molecule
    # if the optional argument, 'code', is specified, the number will be the default
    #   for that quantum chemistry code
    ncore = 0
    if type(atno) == str:
        # convert symbol to Z value
        atno = elz(atno)
    if type(atno) == dict:
        # a molecular formula
        for el, natom in atno.items():
            ncore += n_core(el) * natom
        return ncore
    if code == 'g09':
        # default for Gaussian09 frozen-core calculations
        core = {
            # these are the minimum atomic numbers (Z) that have
            #   the given number of core elecrons (Z : ncore)
            3  :  2,
            11 : 10,
            19 : 18,
            37 : 36,
            55 : 54, # this is a guess
            87 : 86  # this is a guess
        }
    else:
        core = {
            # these are the minimum atomic numbers (Z) that have
            #   the given number of core elecrons (Z : ncore)
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
    for ki in sorted(core):
        if atno >= ki:
            ncore = core[ki]
    return ncore
##
def read_regex(regex, fhandl, idx=1):
    # Return something from a line matchine a regular expression.
    #   First arg is the regular expression; idx is the match-group
    #	to return.  Return a list of values from all matching lines. 
    fhandl.seek(0)
    matches = []
    regx = re.compile(regex)
    for line in fhandl:
        mch = regx.search(line)
        if mch:
            matches.append(mch.group(idx))
    return matches
##
def spinname(m):
    # given a spin multiplity (m = 2S+1), return the text name
    name = [ 'spinless', 'singlet', 'doublet', 'triplet', 'quartet', 'quintet', 'sextet',
        'septet', 'octet', 'nonet', 'decet', 'undecet', 'duodecet' ]
    m = int(m)
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
    if direction == 'to_eV':
        return multiplier * energy * eV_per_hartree
    elif direction == 'from_eV':
        return multiplier * energy / eV_per_hartree
    else:
        # illegal direction
        return 'unrecognized direction = {:s} in routine hartree_eV'.format(direction)
##
def starting_n(Ltype, nppe=0):
    # given an orbital-angular momentum type ('s', 'p', etc.), 
    # return the lowest possible principal quantum number (1, 2, etc.)
    # The optional second argument is the number of electrons that have
    #   been replaced by an ECP/pseudopotential
    # This routine only handles the common cases
    nmin = {'s': 1, 'p': 2, 'd': 3, 'f': 4, 'g': 5, 'h': 6}
    cases = [2, 10, 18, 28, 36, 46, 54, 60, 68, 78, 92]
    if nppe > 0:
        # Some electrons have been replaced by ECP; adjust the explicit
        #   shell numbers accordingly
        if (not nppe in cases):
            print('*** Unhandled number of ECP-replaced electrons ***')
            print('\tnppe = {:d} in routine "starting_n"'.format(nppe))
            # But go ahead and apply the algorithm, anyway!
        # determine number of shells replaced
        rcore = {'s': 0, 'p': 0, 'd': 0, 'f':0}
        resid = nppe
        nf = (resid - 28) // 32 # number of f shells replaced
        if nf > 0:
            rcore['f'] = nf
            resid -= nf * 14
        nd = (resid - 10) // 18 # number of d shells replaced
        if nd > 0:
            rcore['d'] = nd
            resid -= nd * 10
        np = (resid - 2) // 8  # number of p shells replaced
        if np > 0:
            rcore['p'] = np
            resid -= np * 6
        ns = resid // 2  # number of s shells replaced
        rcore['s'] = ns
        resid -= ns * 2
        if resid != 0:
            print('*** Unexpected residual electrons in routine "starting_n" ***')
        for L in rcore:
            nmin[L] += rcore[L]
    return nmin[Ltype.lower()]
##
def L_degeneracy(Ltype):
    # given an orbital-angular momentum type ('s', 'p', etc.), 
    # return the degeneracy (1, 3, etc.)
    degen = {'s': 1, 'p': 3, 'd': 5, 'f': 7, 'g': 9, 'h': 11, 'i': 13}
    return degen[Ltype.lower()]
##
def combine_MOspin(df, col1='MO', col2='Spin', colnew='Label'):
    # Given a pandas DataFrame, combine a numeric 'MO' field with
    #   a 'Spin' field ('alpha' or 'beta') to create a new 'MO' field
    #   that is a combination like '1a' or '5b'.
    # Return that new DataFrame.
    dfret = df.copy()
    dfret[colnew] = df.apply(lambda x: str(x[col1])+(x[col2])[0], axis=1)
    return dfret
##
class Atom:
    # Minimal description: element symbol + cartesian coordinates
    def __init__(self, el, xyz):
        # 'el' : Element symbol or atomic number
        # 'xyz': cartesian coordinates as list or numpy array 
        self.el = elz(el, choice='symbol')
        self.xyz = np.array(xyz)
    def copy( self ):
        newatom = Atom(self.el, self.xyz)
        return newatom
    def print( self ):
        # print routine
        print('{:s}\t{:9.5f}\t{:9.5f}\t{:9.5f}'.format(self.el, self.xyz[0], self.xyz[1], self.xyz[2]))
        return
##
def distance(pos1, pos2):
    # return distance between two vectors (numpy)
    # return NaN if the vectors have different dimensionality
    if len(pos1) != len(pos2):
        print('Unequal vector dimensions in "distance": dim1 = {:d}, dim2 = {:d}')
        return np.nan 
    v = pos2 - pos1
    d = np.sqrt(np.dot(v, v))
    return d
##
class Geometry:
    # a list of Atoms
    def __init__(self, *args, intype='1list'):
        # three input types are recognized:
        #   '2lists'    : a list of elements and a list of coordinate triples
        #   '1list'     : a list of [el, x, y, z] quadruples
        #   'atlist'    : a list of Atoms
        #   'DataFrame' : a pandas DataFrame with four columns (Z, x, y, z)
        self.atom = []
        if len(args) == 0:
            # return an empty Geometry
            return
        if intype == 'atlist':
            # argument is already a list of Atoms
            self.atom = list(args[0])
            return
        if intype == '1list':
            # argument is a list of quadruples, [el, x, y, z]
            for quad in args[0]:
                at = Atom(quad[0], quad[1:4])
                self.atom.append(at)
            return
        if intype == '2lists':
            # first argument is a list of elements
            # second argument is a list of triples
            nsymb = len(args[0])
            nxyz = len(args[1])
            if nsymb != nxyz:
                print('*** Inconsistent #symb = {:d} and #xyz = {:d} in Geometry initialization'.format(nsymb, nxyz))
                return  # empty 
            for iat in range(nsymb):
                at = Atom(args[0][iat], args[1][iat])
                self.atom.append(at)
            return
        if intype == 'DataFrame':
            # argument is a four-column pandas DataFrame (Z, x, y, z)
            for iat in range(len(args[0].index)):
                elxyz = args[0].iloc[iat]
                at = Atom(elxyz[0], elxyz[1:].tolist())
                self.atom.append(at)
    def addatom(self, atom):
        self.atom.append(atom)
        return
    def delatom(self, iatom):
        del self.atom[iatom]
        return
    def copy(self, elements=[]):
        # If a list of elements is provided, the new Geometry includes only them
        newgeom = Geometry()
        for a in self.atom:
            if (len(elements) == 0) or (a.el in elements):
                newgeom.addatom(a.copy())
        return newgeom
    def natom(self):
        return len(self.atom)
    def print(self):
        # printing routine
        print('el\t     x\t\t     y\t\t     z')
        for atom in self.atom:
            atom.print()
        return
    def distance(self, i, j):
        # distance between atoms i and j
        try:
            d = distance(self.atom[i].xyz, self.atom[j].xyz)
            return d
        except IndexError:
            s = '*** Illegal atom number in Geometry.distance(): ' + \
                'i = {:d}, j = {:d}'.format(i, j)
            print(s)
            return np.nan
    def distmat(self):
        # 2D array of interatomic distances
        xyz = [a.xyz for a in self.atom]
        dmat = cdist(xyz, xyz, metric='euclidean')
        return dmat
##
