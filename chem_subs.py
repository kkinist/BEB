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
            vals.append(elz(el, choice))
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
def mulpopDiff(df1, df2):
    print('--- OBSOLETE ---')
    # Compare two pandas DataFrames with Mulliken population data,
    #   as returned by routine 'read_AOpop_in_MOs()' in 'g09_subs.py'
    # Return a DataFrame:
    #   (1) MO number
    #   (2) Energy difference (E2 - E1)
    #   (3) Popdiff : numerical difference between AO populations (per MO)
    #   (4) AOpops : DataFrame, columns = AO labels, index = ['P1', 'P2']
    #     (4a) AO contribution in first orbital set
    #     (4b) AO contribution in second orbital set
    #
    # The two input DataFrames may not have the same list of MO numbers
    #   only consider MOs they have in common
    MOnums = sorted(set(list(df1.MO)) & set(list(df2.MO)))  # set intersection
    jointMO = []
    ediff = []
    popdiff = []
    AOpops = []  # a list of DataFrames
    for MO in MOnums:
        orb1 = df1[df1.MO == MO]
        orb2 = df2[df2.MO == MO]
        e1 = orb1.iloc[0]['Energy']
        e2 = orb2.iloc[0]['Energy']
        ediff.append(e2 - e1)
        jointMO.append(MO)
        # compare the Mulliken populations
        mulpop1 = {}
        mulpop2 = {}
        # create a label for each AO that looks like '#5-p' for a p-orbital on atom #5
        for irow in range(len(orb1.index)):
            s = '#{:d}-{:s}'.format(orb1.iloc[irow]['Atom#'], orb1.iloc[irow]['L'])
            c = orb1.iloc[irow]['Contrib']
            mulpop1[s] = c
        for irow in range(len(orb2.index)):
            s = '#{:d}-{:s}'.format(orb2.iloc[irow]['Atom#'], orb2.iloc[irow]['L'])
            c = orb2.iloc[irow]['Contrib']
            mulpop2[s] = c
        # for each AO, subtract the populations
        pdiff = 0
        # make a list of the AOs used by either orbital
        aolist = sorted(set(list(mulpop1.keys()) + list(mulpop2.keys())))
        # create DataFrame for this MO, all zeros at first
        dfpop = pd.DataFrame(data=np.zeros((2,len(aolist))), columns=aolist)
        for ao in aolist:
            try:
                pdiff += abs(mulpop2[ao] - mulpop1[ao])
                dfpop.iloc[0][ao] = mulpop1[ao]
                dfpop.iloc[1][ao] = mulpop2[ao]
            except:
                # probably an AO present in only one of the two orbitals
                try:
                    pdiff += abs(mulpop1[ao])
                    dfpop.iloc[0][ao] = mulpop1[ao]
                    dfpop.iloc[1][ao] = 0
                except:
                    # must be the other orbital
                    pdiff += abs(mulpop2[ao])
                    dfpop.iloc[1][ao] = mulpop2[ao]
                    dfpop.iloc[0][ao] = 0
        popdiff.append(pdiff)
        dfpop.insert(0, 'orbset', ['P1', 'P2'])
        AOpops.append(dfpop.set_index('orbset'))
    data = list(zip(jointMO, ediff, popdiff, AOpops))
    cols = ['MO', 'E2-E1', 'P2-P1', 'AOpops']
    dfdiff = pd.DataFrame(data=data, columns=cols)
    return dfdiff
##
def JSdm(P, Q, base=4):
    # Jensen-Shannon divergence metric; base=4 gives range = [0, 1]
    # P and Q are *discrete* PDFs (with same data type)
    # Allowed data types: tuple; list; dict; numpy array
    # P and Q must be same length, except when dict
    # Return:
    #   (1) metric (float)
    #   (2) messages (list of string)
    #
    message = []
    if type(P) != type(Q):
        print('*** P and Q must be same data type in routine JSdm() ***')
        return (None, None)
    if (type(P) == list) or (type(P) == tuple) or (type(P) == np.ndarray):
        P = np.array(P).astype(float)
        Q = np.array(Q).astype(float)
        allkeys = []   # length will be tested later, to infer input type
    elif type(P) == dict:
        # make a sorted list of all the keys
        allkeys = sorted(set(list(P.keys()) + list(Q.keys())))
        Plist = []
        Qlist = []
        for key in allkeys:
            try:
                Plist.append(P[key])
            except:
                # probably key is not present in this dict
                Plist.append(0)
            try:
                Qlist.append(Q[key])
            except:
                Qlist.append(0)
        if P.keys() != Q.keys():
            message.append('Different key lists merged for P and Q')
        # convert list to numpy array
        P = np.array(Plist).astype(float)
        Q = np.array(Qlist).astype(float)
    else:
        print('*** Unhandled data type in routine JSdm():', type(P))
        return (None, None)
    # No negative values are allowed
    if len(np.where(P < 0)[0]) or len(np.where(Q < 0)[0]):
        print('*** Negative values not allowed in routine JSdm() ***')
        return (None, None)
    # P and Q must have the same number of elements
    if len(P) != len(Q):
        print('*** P and Q must have same length in routine JSdm() ***')
        return (None, None)
    # Normalize both PDFs (L1-normalization)
    Plen = P.sum()
    Qlen = Q.sum()
    if (Plen == 0) or (Qlen == 0):
        print('*** P and Q may not be all zeros in routine JSdm() ***')
        return (None, None)
    P /= Plen
    Q /= Qlen
    pqsum = P + Q
    # find any zeros in (P+Q) and delete corresponding elements in P, Q, and P+Q
    nullidx = np.where(pqsum == 0)[0]
    if len(nullidx > 0):
        # delete the troublesome elements
        if len(allkeys) > 0:
            # input was dict
            message.append('Deleted null elements with indices ' + str([allkeys[i] for i in nullidx]))
        else:
            # input was list-like
            message.append('Deleted null elements with indices ' + str(nullidx))
        P = np.delete(P, nullidx)
        Q = np.delete(Q, nullidx)
        pqsum = np.delete(pqsum, nullidx)
    # compute the JSDM
    # P or Q may still contain zeros, so don't take straight logarithm
    #   instead, use x*ln(y) = ln(y**x) and convention 0**0 = 1
    s1 = 2 * P / pqsum
    s2 = 2 * Q / pqsum
    s1 = s1 ** P
    s2 = s2 ** Q
    s1 = np.log(s1) / np.log(base)
    s2 = np.log(s2) / np.log(base)
    dsq = (s1 + s2).sum()
    return np.sqrt(dsq), message
##
def matchMOsByPop(df1, df2, thresh=0.1):
    print('--- OBSOLETE ---')
    # Compare two pandas DataFrames with Mulliken population data,
    #   as returned by routine 'read_AOpop_in_MOs()' in 'g09_subs.py'
    # Argument 'thresh' is the cutoff for AO contribution differences
    #   to be considered bad.
    # Return a dict of MO number correspondences, based upon Mulliken
    #   AO contributions.  (Weighted) orbital energies are also considered,
    #   as a proxy for principal quantum number. The dict only includes
    #   orbitals that appear to be mismatched. 
    # Keys are MO labels in df2, values are MO labels in df1.
    #
    # Look for better matches: first contruct matrix of AO-pop differences
    Eweight = 0.1  # weighting for energy differences
    momap = {}
    popDiff = mulpopDiff(df1, df2)
    bigDiff = popDiff[(popDiff['P2-P1'] + popDiff['E2-E1'].abs()) > thresh]
    ndim = len(bigDiff.index)
    if (ndim == 0):
        # no problems to correct
        return momap
    diffmat = np.zeros((ndim, ndim))
    idxlist = list(bigDiff.index)  # these are the row/column numbers
    for iorb in bigDiff.index:
        imo = bigDiff.loc[iorb]['MO']
        refdP = bigDiff.loc[iorb]['P2-P1']
        refE = df1[df1.MO == imo].iloc[0]['Energy']
        P1 = bigDiff.loc[iorb]['AOpops'].loc['P1']
        for jorb in bigDiff.index:
            jmo = bigDiff.loc[jorb]['MO']
            jE = df2[df2.MO == jmo].iloc[0]['Energy']
            dE = jE - refE
            # compute the sum|P2-P1| difference in AO contributions
            P2 = bigDiff.loc[jorb]['AOpops'].loc['P2']
            aos = list(bigDiff.loc[iorb]['AOpops'].columns.values) + list(bigDiff.loc[jorb]['AOpops'].columns.values)
            aos = sorted(set(aos))
            # include |E2-E1|*Eweight in the sum
            dP = abs(dE) * Eweight
            for ao in aos:
                try:
                    dP += abs(P2.loc[ao] - P1.loc[ao])
                except:
                    # probably 'ao' is missing from one of the orbitals
                    try:
                        dP += abs(P1.loc[ao])
                    except:
                        # must be the other orbitals
                        dP += abs(P2.loc[ao])
            diffmat[idxlist.index(iorb), idxlist.index(jorb)] = dP
    claimed = []  # list of orbitals in set2 as they are paired
    pairing = list(range(ndim))  # default is no rearrangement
    # start pairing with the closest matches 
    for iorb in diffmat.min(axis=1).argsort():
        # loop over rows, in order by best available match
        for jorb in diffmat[iorb, :].argsort():
            # loop over columns, in order by best match
            if jorb in claimed:
                # this orbital already paired
                continue
            # this is a pairing
            claimed.append(jorb)
            pairing[iorb] = jorb
            break  # done with this first-set MO
    # convert into a mapping of MO numbers
    for i in range(ndim):
        imo = bigDiff.iloc[i]['MO']
        j = pairing[i]
        jmo = bigDiff.iloc[j]['MO']
        if imo != jmo:
            # report only these non-identity mappings
            momap[jmo] = imo  # key is the MO number in the 2nd set
    return momap
##
def AOpopdiffmats(df1, df2):
    # Compare two pandas DataFrames with Mulliken population data,
    #   as returned by routine 'read_AOpop_in_MOs()' in 'g09_subs.py'
    # Return two numpy 2D-arrays:
    #   (1) JSdm() differences in AO populations (Jensen-Shannon divergence metric)
    #   (2) (E2-E1) orbital energy differences
    # Also return two lists of MO numbers:
    #   (3) MO numbers in df1 (rows of matrices)
    #   (4) MO numbers in df2 (columns of matrics)
    MOlist1 = sorted(set(df1.MO))
    MOlist2 = sorted(set(df2.MO))
    nmo1 = len(MOlist1)
    nmo2 = len(MOlist2)
    dPmat = np.zeros((nmo1, nmo2))
    dEmat = np.zeros((nmo1, nmo2))
    for imo in MOlist1:
        # looping over MOs in first set
        idx = MOlist1.index(imo)  # row number in returned matrices
        orb1 = df1[df1.MO == imo]
        E1 = orb1.iloc[0]['Energy']
        # convert AO populations into a dict
        mulpop1 = {}
        # create a label for each AO that looks like '#5-p' for a p-orbital on atom #5
        for ao in orb1.index:
            s = '#{:d}-{:s}'.format(orb1.loc[ao]['Atom#'], orb1.loc[ao]['L'])
            c = orb1.loc[ao]['Contrib']
            if c < 0:
                # this should not happen
                print('*** Warning: negative AO contribution found in MO #{:d} in df1 ***'.format(imo))
                c = abs(c)
            mulpop1[s] = c
        # loop over orbitals in second set
        for jmo in MOlist2:
            jdx = MOlist2.index(jmo)  # column number in returned matrices
            orb2 = df2[df2.MO == jmo]
            E2 = orb2.iloc[0]['Energy']
            dEmat[idx, jdx] = E2 - E1  # signed difference
            # construct dict of AO populations as above
            mulpop2 = {}
            for ao in orb2.index:
                s = '#{:d}-{:s}'.format(orb2.loc[ao]['Atom#'], orb2.loc[ao]['L'])
                c = orb2.loc[ao]['Contrib']
                if c < 0:
                    # this should not happen
                    print('*** Warning: negative AO contribution found in MO #{:d} in df2 ***'.format(jmo))
                    c = abs(c)
                mulpop2[s] = c
            # get JSdm distance between the two AO population vectors
            dist = JSdm(mulpop1, mulpop2)
            dPmat[idx, jdx] = dist[0]
    return dPmat, dEmat, MOlist1, MOlist2
##
def orbitalPopMatch(df1, df2, Eweight=0.1, diagBias=0.001):
    # Find which MOs correspond between two calculations. 
    # Note: Cannot distinguish degenerate orbitals! 
    # Compare two pandas DataFrames with Mulliken population data,
    #   as returned by routine 'read_AOpop_in_MOs()' in 'g09_subs.py'
    # Argument 'Eweight' is the weight to give to energy differences.
    # Argument 'diagBias' is the preference to give to the existing
    #   orbital numbering.
    # Return a dict of MO number correspondences. The dict only includes
    #   orbitals that appear to be mismatched. 
    # Keys are MO labels in df2, values are MO labels in df1.
    #
    dPmat, dEmat, MOs1, MOs2 = AOpopdiffmats(df1, df2)
    # count the MOs in each orbital set
    norb1 = len(MOs1)
    norb2 = len(MOs2)
    nmo = min(norb1, norb2)
    # use unsigned energy differences
    diffmat = dPmat + Eweight * np.fabs(dEmat)
    # install the bias toward perserving the existing numbering
    # Note: Gaussian prints the populations only to 0.01 precision
    for i in range(norb1):
        imo = MOs1[i]
        try:
            j = MOs2.index(imo)
            diffmat[i, j] -= diagBias
        except:
            # probably an orbital missing from set 2
            pass
    # find closest distance for each row
    rowmin = diffmat.min(axis=1)
    # sort by increasing distance (i.e., best matches first)
    rowlist = rowmin.argsort()
    # truncate to smallest dimension
    rowlist = rowlist[0 : nmo]
    claimed = []  # list of orbitals in set2 as they are paired
    pairing = {}  # mapping between orbital indices (not MO numbers/labels)
    for iorb in rowlist:
        # loop over matrix rows, starting with those with best available matches
        for jorb in diffmat[iorb, :].argsort():
            # loop over columns, starting with best match
            if jorb in claimed:
                # this orbital already paired
                continue
            # this is a pairing
            claimed.append(jorb)
            pairing[iorb] = jorb
            break  # done with this first-set MO
    # convert into a mapping of MO numbers
    momap = {}
    for i in pairing.keys():
        imo = MOs1[i]  # MO number from first set
        j = pairing[i]
        jmo = MOs2[j]  # MO number from second set
        if imo != jmo:
            # report only non-identity mappings
            momap[jmo] = imo  # key is the MO number in the 2nd set
    return momap 
##
