# Routines for use in processing quantum chemistry data.
# Karl Irikura, started 8/26/2008 
#
import re
import math
import sys
from numpy import *
import itertools
import string
import datetime
from scipy.spatial.distance import cdist
##
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

def signum( x ):
        # return 1 if x > 0, -1 if x < 0, and 0 if x == 0
        if x > 0:
                return 1
        elif x < 0:
                return -1
        else:
                return 0
##
def vlength( ar ):
        # return length of numpy array (vector) (8/27/08)
        d = math.sqrt( dot( ar, ar ) )
        return d
##
def distance( pos1, pos2 ):
        # return distance between two cartesian triples (8/27/08)
        v = pos2 - pos1
        d = vlength( v )
        return d
##
def nvec( ar, r = 1.0 ):
        # return copy of numpy array 'ar' normalized to length 'r' (8/27/08)
        s = vlength( ar )
        if s == 0.0:
            return zeros( alen(ar) )
        d = r / vlength( ar )
        g = ar * d
        return g
##
def vproject( u, v, along = True ):
        # return the component of array 'u' along 'v' (5/24/11)
        # if third arg is False, then return the component of 'u' perpendicular to 'v'
        b = nvec( v )
        a = dot( u, b ) * b
        perp = u - a
        if along:
                return a
        else:
                return perp
##
def angleabc(a, b, c, units='rad'):
    # return the angle a-b-c in radians, where all are cartesian triples (4/9/15)
    v1 = a - b
    v2 = c - b
    s = dot( v1, v2 )
    s /= vlength( v1 )
    s /= vlength( v2 )
    theta = math.acos( s )
    if units == "deg":
        theta *= 180 / pi
    return theta
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
class Atom:
        # Element symbol (string); cartesian coordinates (numpy array);
        # mass (float)
        # 8/27/08
        def __init__( self, el, x, y, z, m ):
            self.el = el
            self.xyz = array( ( float(x), float(y), float(z) ) )
            self.mass = float( m )
        def copy( self ):
            # 5/24/11
            newatom = Atom( self.el, self.xyz[0], self.xyz[1], self.xyz[2], self.mass )
            return newatom
        def prt( self ):
            # print routine
            print '%s\t%9.5f\t%9.5f\t%9.5f\t%9.5f u' % ( self.el, self.xyz[0], self.xyz[1], self.xyz[2], self.mass )
            return
        def hbonding(self):
            # this is an explicit, fixed list; return True or False
            hbonders = ['N', 'O', 'F', 'S']
            return (self.el in hbonders)
##
class Hbond:
    # Stand-alone hydrogen bond (not attached to its molecule)
    # 4/2/15
    def __init__(self, idonor, ih, iacceptor, rDA=0.0, rHA=0.0, aDHA=0.0):
        # required arguments are indices into some Geometry object
        # optional arguments are DA = donor-acceptor distance, 
        #   HA = H-acceptor distance, DHA = donor-H-acceptor angle (degrees)
        self.donor = idonor
        self.Hatom = ih
        self.acceptor = iacceptor
        self.DA = rDA
        self.HA = rHA
        self.DHA = aDHA
        self.criteria = []  # list of strings that specify H-bonding criteria
    def triple(self):
        # just return tuple (idonor, ih, iacceptor)
        return (self.donor, self.Hatom, self.acceptor)
    def flabel(self, molecule):
        # return list of formatted atom labels in order D, H, A
        labeld = molecule.alabel(self.donor)
        labelh = molecule.alabel(self.Hatom)
        labela = molecule.alabel(self.acceptor)
        return [labeld, labelh, labela]
    def prt(self, molecule):
        # text display
        if molecule == 'header':
            print 'D\tH\tA\tr(DA)\tr(HA)\ta(DHA)\tCriteria'
        else:
            # Look up atoms within 'molecule'
            pbuf = '\t'.join(self.flabel(molecule))
            pbuf += '\t{:.3f}\t{:.3f}\t{:5.1f}\t'.format(self.DA, self.HA, self.DHA)
            pbuf += ', '.join(self.criteria)
            print pbuf
##
def xyz2Atom( atno, xyz ):
	# given Z value (or element symbol) and list [x, y, z], return an Atom
        if type( atno ) == int:
        	el = elz( atno )
        else:
        	# we were probably given an element symbol, not an atomic number
	        el = atno
	        atno = elz( el )
        m = atomic_weight( atno )
        return Atom( el, xyz[0], xyz[1], xyz[2], m )
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
        # starting from a Gaussian-style stoichiometry (e.g., C2H6Si), return a
        #  hash with formula[element] = number
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
class Geometry:
        # a list of Atoms (8/28/08)
        def __init__( self ):
                self.atom = []
                return
        def addatom( self, atom ):
                self.atom.append( atom )
                return
        def delatom( self, iatom ):
                del self.atom[iatom]
                return
        def copy( self, elements=[] ):
                # 5/24/11
                # If a list of elements is provided, the new Geometry includes only them (10/21/14)
                newgeom = Geometry()
                for a in self.atom:
                        if (len(elements) == 0) or (a.el in elements):
                                newgeom.addatom( a.copy() )
                return newgeom
        def natom( self ):
                return len( self.atom )
        def alabel(self, i, iref=1):
            # return formatted label for atom #i (where first atom is numbered 'iref')
            n = i + iref
            symbol = self.atom[i].el
            return '{:s}_{:d}'.format(symbol, n)
        def mass( self ):
                # molecular mass
                m = 0.
                for a in self.atom:
                        m += a.mass
                return m
        def com( self ):
                # molecular mass and coordinates of center-of-mass
                cent = zeros( 3 )       # initialize to zeroes
                for i in range( len( self.atom ) ):
                        cent += self.atom[i].mass * self.atom[i].xyz
                cent /= self.mass()
                return cent
        def center( self ):
                # translate to put center-of-mass at the origin
                cent = self.com()
                for a in self.atom:
                        a.xyz -= cent
                return
        def length( self, i, j ):
                # distance between atoms i and j
                try:
                    d = distance( self.atom[i].xyz, self.atom[j].xyz )
                    return d
                except IndexError:
                    s = '**illegal atom number: i = {}, j = {}, natom = {}\n'.format( i, j, self.natom() )
                    sys.stderr.write( s )
                    return -1
        def stoichiometry(self):
                # stoichiometry string (without charge or spin multiplicity)
                order = ['C', 'H', 'N', 'O', 'F', 'Cl', 'S', 'P']
                # build hash of elements and their atom counts
                acount = {}
                for a in self.atom:
                        try:
                                acount[a.el] += 1
                        except:
                                acount[a.el] = 1
                stoich = ''
                for e in order:
                        if e in acount.keys():
                                stoich += "{:s}{:d}".format(e, acount[e])
                # alphabetical for elements not specified in 'order'
                others = []
                for e in acount.keys():
                        if not e in order:
                                others.append(e)
                if len(others):
                        for e in others.sort():
                                stoich += "{:s}{:d}".format(e, acount[e])
                return stoich
        def distmat( self ):
                # 2D array of interatomic distances (12/18/14)
                xyz = [a.xyz for a in self.atom]
                dmat = cdist(xyz, xyz, metric='euclidean')
                return dmat
        def distmat_looping( self ):
                # 2D array of interatomic distances (5/16/12) SLOW
                n = len( self.atom )
                dmat = zeros( (n, n) )
                for i in range( n ):
                        for j in range( i ):
                                dmat[i][j] = self.length( i, j )
                                dmat[j][i] = dmat[i][j]
                return dmat
        def distmat_eig( self ):
                # return eigenvalues and eigenvectors of distance matrix
                # sort them by eigenvalue
                dmat = self.distmat()
                vals, vecs = linalg.eigh( dmat )
                vecs = vecs[argsort(vals)]
                vals = sort(vals)
                return vals, vecs
        def del_dupl( self ):
                # remove duplicated atoms (same element and coordinates)
                # return the number of atoms deleted
                n = self.natom()
                ndel = 0
                thresh = 0.01   # distance threshold
                for i in range( n - 1 ):
                        for j in reversed( range( i+1, n ) ):
                                if self.atom[i].el == self.atom[j].el:
                                        d = self.length(i,j)
                                        if d < thresh:
                                                # delete atom j
                                                self.delatom(j)
                                                ndel += 1
                        n = self.natom()
                return ndel
        def bonded( self, i, j ):
                # return True if bonded, else False (based on distance only) (3/2/10)
                r0 = r0_ref( self.atom[i].el, self.atom[j].el )
                tol = 1.3       # tolerated amount of bond stretching
                if self.length( i, j ) < r0 * tol:
                        return True
                return False
        def vec( self, i, j ):
                # return vector j-i (3/2/10)
                v = self.atom[j].xyz - self.atom[i].xyz
                return v
        def normvec( self, i, j ):
                # return normalized interatomic vector j-i
                return nvec( self.vec( i, j ) )
        def angle( self, i, j, k, units="rad" ):
                # angle between atoms i-j-k, in radians
                theta = angleabc(self.atom[i].xyz, self.atom[j].xyz, self.atom[k].xyz, units)
                return theta
        def angle_old( self, i, j, k, units="rad" ):
                # angle between atoms i-j-k, in radians
                v1 = self.atom[i].xyz - self.atom[j].xyz
                v2 = self.atom[k].xyz - self.atom[j].xyz
                s = dot( v1, v2 )
                s /= vlength( v1 )
                s /= vlength( v2 )
                theta = math.acos( s )
                if units == "deg":
                        theta *= 180 / pi
                return theta
        def dihedral( self, i, j, k, l, type="linear", units="rad" ):
                # calculate dihedral angle in radians (optionally in degrees)
                # type='linear':  connectivity is i-j-k-l
                # type='branched': connectivity is i-j<kl (both k and l bonded to j)
                a = self.vec( j, i )
                b = self.vec( j, k )
                c = self.vec( l, k )
                if type == "branched":
                        c = self.vec( j, l )
                        a, b = b, a
                b = nvec( b )   # normalize the mid-vector
                x = a - b * dot( a, b ) # component of a normal to b
                z = c - b * dot( c, b )
                x = nvec( x )
                z = nvec( z )
                if ( vlength(x) == 0.0) or ( vlength(z) == 0.0):
                        # something is linear; dihedral is undefined
                        return float('nan')
                s = dot( x, z )
                if s < -1:
                        s = -1
                if s > 1:
                        s = 1
                phi = math.acos( s )    # in range [0, pi)
                s = cross( x, z )       # vector cross-product to get sign
                s = dot( s, b )
                s = signum( s )
                phi *= s        # include sign (right-handed definition)
                if s == 0:
                        # x and z are parallel
                        if dot( x, z ) > 0:
                                phi = 0
                        else: 
                                phi = pi
                if units == "deg":
                        phi *= 180 / pi
                return phi
        def grow_atom(self, el, j, k, l, r, theta, phi):
                # add another atom of specified element at position (r, theta, phi)
                #   relative to atoms (j, k, l)
                #   angles specified in degrees
                a = math.radians( theta )
                d = math.radians( phi )
                loc = self.atom[ j ].xyz.copy() # position of connected atom
                # determine the vector from atom j to the new atom i
                # first component/basis vector (parallel to j-k)
                n1 = self.normvec( j, k )
                displ = n1 * math.cos( a )
                # second component (~parallel to k-l)
                n2 = self.normvec( k, l )
                n2 = n2 - dot( n2, n1 ) * n1    # project out the component along n1
                n2 = nvec( n2 ) # renormalize
                displ += n2 * math.sin( a ) * math.cos( d )
                # third component (perpendicular to other two)
                n3 = cross( n1, n2 )
                displ += n3 * math.sin( a ) * math.sin( d )
                loc += r * displ
                [x, y, z] = list( loc )
                m = atomic_weight(el)
                self.addatom( Atom(el, x, y, z, m) )
                return 
        def inertia_tensor( self ):
                # moment-of-inertia tensor (without any unit conversions) (4/7/11)
                inertia = zeros( (3,3) )
                for a in self.atom:
                        inertia[0][0] += a.mass * ( a.xyz[1] * a.xyz[1] + a.xyz[2] * a.xyz[2] )
                        inertia[1][1] += a.mass * ( a.xyz[0] * a.xyz[0] + a.xyz[2] * a.xyz[2] )
                        inertia[2][2] += a.mass * ( a.xyz[1] * a.xyz[1] + a.xyz[0] * a.xyz[0] )
                        inertia[0][1] -= a.mass * a.xyz[0] * a.xyz[1]
                        inertia[0][2] -= a.mass * a.xyz[0] * a.xyz[2]
                        inertia[1][2] -= a.mass * a.xyz[1] * a.xyz[2]
                inertia[1][0] = inertia[0][1]
                inertia[2][0] = inertia[0][2]
                inertia[2][1] = inertia[1][2]
                inertia = inertia.T
                return inertia
        def nuclear_repulsion( self, unit="bohr" ):
                # nuclear repulsion energy (hartree) (10/20/2014)
                n = len( self.atom )
                repulsion = 0.0
                for i in range( n ):
                        zi = elz( self.atom[i].el )
                        for j in range( i ):
                                rij = self.length( i, j )
                                zj = elz( self.atom[j].el )
                                repulsion += zi * zj / rij
                if unit == "angstrom":
                        repulsion *= bohr
                return repulsion
        def rotation( self ):
                # return rotational constants, moments of inertia, and principal axes
                ### around the center of mass ###
                centered = self.copy()
                centered.center()
                imat = centered.inertia_tensor()
                moment, axes = linalg.eigh( imat )
                # convert moment to kg.m^2, assuming distances in angstrom and masses in u
                moment /= 1.0e20 * avogadro * 1000.0
                rotconst = planck / ( 8 * pi * pi * clight * moment )   # now in units (1/m)
                rotconst *= clight * 1.0e-9      # now in GHZ
                return rotconst, moment, axes
        def rgyr( self ):
                # return radius of gyration (10/20/2014)
                mrr = 0
                cent = self.com()
                for a in self.atom:
                        rvec = a.xyz - cent
                        mrr += a.mass * dot( rvec, rvec )
                mrr /= self.mass()
                return sqrt( mrr )
        def sayvetz( self, fcmw ):
                # THIS ROUTINE DOES NOT YET WORK ! 
                # apply Sayvetz restrictions, following Califano's book
                # fcmw is the matrix of mass-weighted force constants for vibrational analysis
                self.center()   # move center of mass to the origin
                dimen = 3 * self.natom()
                dmat = zeros( (dimen, dimen) )  # this is the array of D vectors
                # define D1, D2, and D3
                for i in range( 3 ):
                        for iat in range( self.natom() ):
                                j = 3 * iat + i
                                dmat[i][j] = sqrt( self.atom[iat].mass )
                # define D4, D5, and D6
                for iat in range( self.natom() ):
                        msqrt = sqrt( self.atom[iat].mass )
                        x0 = self.atom[iat].xyz[0]
                        y0 = self.atom[iat].xyz[1]
                        z0 = self.atom[iat].xyz[2]
                        k = 3 * iat
                        dmat[3][k+1] = -msqrt * z0
                        dmat[3][k+2] = msqrt * y0
                        dmat[4][k] = msqrt * z0
                        dmat[4][k+2] = -msqrt * x0
                        dmat[5][k] = -msqrt * y0
                        dmat[5][k+1] = msqrt * x0
                # normalize
                for i in range( 6 ):
                        n = sqrt( dot( dmat[i], dmat[i] ) )
                        dmat[i] /= n
                pmat = eye( dimen ) - dmat
                # normalize pmat
                for i in range( 6 ):
                        n = sqrt( dot( pmat[i], pmat[i] ) )
                        pmat[i] /= n
                pmat = eye( dimen )
                dmat = dot( fcmw, pmat.T )
                dmat = dot( pmat, dmat )
                return dmat
        def prt( self, start=0 ):
                # printing routine; number atoms beginning with 'start'
                print 'n\tel\t    mass\t     x\t\t     y\t\t     z'
                for i in range( self.natom() ):
                        at = self.atom[i]
                        print '%d\t%s\t%9.5f\t%9.5f\t%9.5f\t%9.5f' % ( i+start, at.el, at.mass, at.xyz[0], at.xyz[1], at.xyz[2] )
                print 'Center of mass:\t%9.5f' % self.mass(),
                print '\t%9.5f\t%9.5f\t%9.5f' % ( tuple( self.com() ) )
                return
        def prtxmol( self ):
                # print in Xmol's XYZ format
                print self.natom()
                print 'molecular mass is', self.mass()
                for i in range( self.natom() ):
                        at = self.atom[i]
                        print '%s\t%9.5f\t%9.5f\t%9.5f' % ( at.el, at.xyz[0], at.xyz[1], at.xyz[2] )
                return
        def write_xyz( self, filename='mol.xyz', comment='XYZ file (U. Minnesota XMOL)' ):
                # write coordinates to a file in XYZ format
                fp = open( filename, 'w' )
                s = '%d\n%s\n' % ( self.natom(), comment )
                fp.write( s )
                for i in range( self.natom() ):
                        s = '%-3s %.6f %.6f %.6f\n' % ( self.atom[i].el, self.atom[i].xyz[0], self.atom[i].xyz[1], self.atom[i].xyz[2] )
                        fp.write( s )
                fp.close()
                return
        def write_mol(self, molfile='mol.mol', comment='MOL file (MDL format)'):
                # write to file in Accelrys/MDL's MOL V2000 format
                # no charge information is included
                fmol = open(molfile, 'w')
                s = 'file generated by \'mol_subs.py\' for ' + self.stoichiometry() + '\n'
                fmol.write(s)
                s = '\ton ' + datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y") + '\n'
                fmol.write(s)
                # next is 'counts' line
                fmol.write(comment + '\n')
                natom = self.natom()
                connex = self.connection_table()
                nbond = int(connex.sum()/2)
                chiral = 1
                s = "{:3d}{:3d}{:3d}{:3d}{:3d}{:3d}{:3d}{:3d}{:3d}{:3d}{:3d}{:>6s}\n".format(natom, 
                        nbond, 0, 0, chiral, 0, 0, 0, 0, 0, 999, 'V2000')
                fmol.write(s)
                # next are atoms
                for a in self.atom:
                        s = "{:10.4f}{:10.4f}{:10.4f} ".format(a.xyz[0], a.xyz[1], a.xyz[2])
                        mdiff = 0
                        charge = 0
                        hadd = 0
                        s += "{:3s}{:2d}{:3d}{:3d}{:3d}".format(a.el, mdiff, charge, 0, hadd)
                        h0 = 0
                        s += "{:3d}{:3d}{:3d}{:3d}{:3d}".format(0, 0, h0, 0, 0)
                        s += "{:3d}{:3d}{:3d}\n".format(0, 0, 0)
                        fmol.write(s)
                # next are bonds
                for i in range(natom):
                        for j in range(i+1,natom):
                                # call everything a single bond!
                                if connex[i][j]:
                                        bondtype = 1
                                        s = "{:3d}{:3d}{:3d}".format(i+1, j+1, bondtype)
                                        s += "{:3d}{:3d}{:3d}{:3d}\n".format(0, 0, 0, 0)
                                        fmol.write(s)
                # that's enough
                fmol.write('M  END\n')
                fmol.close()
                return
        def xyzvec( self, mweight=False ):
                # return numpy 1D array of all cartesian coordinates with optional mass weighting
                v = self.atom[0].xyz
                if mweight:
                        v *= sqrt( self.atom[0].mass )
                for i in range( len( self.atom ) ):
                        if i == 0:
                                continue
                        if mweight:
                                v = hstack( [ v, self.atom[i].xyz * sqrt(self.atom[i].mass) ] )
                        else:
                                v = hstack( [ v, self.atom[i].xyz ] )
                return v
        def symbols(self):
                # return python list of element symbols
                sym = []
                for i in range(self.natom()):
                        sym.append(self.atom[i].el)
                return sym
        def connection_table(self):
                # return a connection table:  a 2D array indicating bonded distances (= 0 or 1)
                # 12/18/2014
                tol = 1.3  # bonding criterion: less than 130% of standard interatomic distance
                dmat = self.distmat() / tol
                connex = zeros_like(dmat)
                for i in range(self.natom()):
                        for j in range(i):
                                # j < i
                                if dmat[i][j] < r0_ref(self.atom[i].el, self.atom[j].el):
                                        connex[i][j] = 1
                                        connex[j][i] = 1
                return connex
        def connection_table_looping( self ):
                # return a connection table:  a 2D array indicating bonded distances (= 0 or 1)
                # 3/21/2012
                # SLOW
                dimen = self.natom()
                connex = zeros( (dimen, dimen) )
                for i in range( dimen ):
                        for j in range( i ):
                                # j < i
                                if ( self.bonded( i, j ) ):
                                        connex[i][j] = 1
                                        connex[j][i] = 1
                return connex
        def bonded_list(self, iatom):
            # return a list of indices to atoms that are bonded to iatom
            connex = self.connection_table()
            blist = where(connex[iatom] > 0)[0]   # 'where' returns a tuple
            blist = blist.tolist()  # convert to regular python list
            return blist
        def zmpair(self, j, npair=1):
            # return a list of pairs (k, l) of atoms that can be used
            #   to define the position of a new atom bonded to atom j
            # If npair > 1, return that many pairs if possible.
            # If npair < 1, return all possible pairs
            connex = self.connection_table()
            klis = where(connex[j] > 0)[0]
            llis = []
            klpair = []
            for ik in range(len(klis)):
                llis.append( [l for l in where(connex[klis[ik]] > 0)[0] if l != j] )
                for l in llis[ik]:
                    klpair.append( (klis[ik], l) )
            np = len(klpair)
            if (npair < 1) or (npair >= np):
                # return all pairs found
                return klpair
            # a limited number of pairs was requested: must eliminate some 
            # find any near-linear geometries
            angl = []
            adev = []  # deviation of angle from 0 or 180 degrees
            for ip in range(np):
                a = self.angle(j, klpair[ip][0], klpair[ip][1], units="deg")
                angl.append(a)
                b = min(180. - a, a)
                adev.append(b)
            # sort by decreasing adev
            ii = sorted(range(np), key=lambda x: adev[x], reverse=True)
            # return the requested number of pairs
            if npair > 1:
                return [klpair[i] for i in ii[:npair]]
            else:
                # special case for 1 pair: return bare tuple
                return klpair[0]
        def findElement(self, el, valent=-1):
            # return a list of indices for specified chemical element and specified
            #   number of bonded atoms (-1 to ignore number of bonded atoms)
            #   4/9/15
            sel = []
            for iat in range(self.natom()):
                if self.atom[iat].el == el:
                    sel.append(iat)
            if valent >= 0:
                # restrict to specified number of bonded atoms
                connex = self.connection_table()
                nconn = sum(connex, axis=1)
                for i in reversed(range(len(sel))):
                    # loop backwards because elements will be deleted
                    iat = sel[i]
                    if nconn[iat] != valent:
                        del sel[i]
            return sel
        def protonated(self, Arg=True, His=True):
            # return a list of atoms that appear to be protonated, with description
            #   (list of pairs [atom index, description] )
            # Assumptions are appropriate (I hope) for closed-shell peptides
            # Flags indicate whether to include guanidinium and imidazolium groups
            proted = []
            connex = self.connection_table()
            nconn = sum(connex, axis=1)
            for iat in range(self.natom()):
                irow = (connex[iat,:] == 1)
                # 'jlis' is the array of indices for atoms bonded to atom[iat]
                jlis = arange(self.natom())[irow]
                # 'nlis' is the list of atomic symbols for atoms bonded to atom[iat]
                nlis = map(lambda j: self.atom[j].el, jlis)
                if self.atom[iat].el == 'N':
                    # nitrogen atom
                    if nconn[iat] == 4:
                        # quaternary; call it protonated even if it's a fixed charge
                        # this catches Lysine, N-terminus, N-protonated amides,
                        #   Glutamine, Asparagine
                        proted.append([iat, 'tetravalent N'])
                if self.atom[iat].el == 'C':
                    # carbon atom; check for guanidine and imidazole
                    if nconn[iat] == 3:
                        if nlis == ['N', 'N', 'N'] and Arg:
                            # guanidine or guanidinium:  Arginine
                            if array_equal(nconn[jlis], array([3, 3, 3])):
                                # protonated gaunidinium
                                # all three N atoms are considered protonated
                                for jat in jlis:
                                    proted.append([jat, 'gaunidinium N'])
                        if sorted(nlis) == ['H', 'N', 'N'] and His:
                            # neutral or protonated imidazole:  Histidine
                            nbeta = sum(nconn[jlis])
                            if nbeta == 7:
                                # protonated imidazole; add both N atoms
                                for jat in jlis:
                                    if self.atom[jat].el == 'N':
                                        proted.append([jat, 'imidazolium N'])
                if self.atom[iat].el == 'O':
                    # oxygen atom
                    if nconn[iat] == 3:
                        # protonated OH group (carboxylic acid or alcohol)
                        proted.append([iat, 'trivalent O'])
                    if nconn[iat] == 2:
                        # could be neutral OH (there are no ethers or esters)
                        #   or protonated carbonyl (amide, carboxylic acid;
                        #   there are no esters, ketones or aldehydes)
                        if sorted(nlis) == ['C', 'H']:
                            # O is bonded to H and to C; check valency of the C
                            ic = jlis[0]
                            if self.atom[ic].el != 'C':
                                # the other one must be the C atom
                                ic = jlis[1]
                            if nconn[ic] == 3:
                                # C atoms is bonded to three atoms; protonated carbonyl
                                #   (amide or carboxylic acid) or neutral carboxylic acid
                                crow = (connex[ic,:] == 1)
                                clis = arange(self.natom())[crow]  # array of atoms bonded to C
                                cbonded = map(lambda j: self.atom[j].el, clis) # element symbols
                                cbonded = sorted(cbonded)
                                if cbonded == ['C', 'N', 'O']:
                                    # protonated amide; find the N atom and also add it to the list
                                    for jat in clis:
                                        if self.atom[jat].el == 'N':
                                            proted.append([iat, 'amidium OH'])  # add the O atom
                                            proted.append([jat, 'amidium NH'])  # add the N atom
                                if cbonded == ['C', 'O', 'O']:
                                    # neutral or protonated carboxylic acid
                                    for jat in clis:
                                        if jat != iat and self.atom[jat].el == 'O':
                                            if nconn[jat] == 2:
                                                # the other O atom also has two bonds
                                                # protonated carboxylic acid
                                                proted.append([iat, 'carboxylium O'])
                                                proted.append([jat, 'carboxylium O'])
            return proted
        def hydrogen_bonds(self, rHA=2.5, aDHA=90.):
            # Analyze H-bonding by multiple methods:
            #   "torshin" : r(H-A) <= 2.5 ang, a(D-H-A) >= 90. deg
            #         from Torshin et al., Protein Eng. 15(5), 359-363 (2002)
            #   "oommachen" : r(D-A) <= 4.5 ang, a(D-H-A) >=90. deg
            #         from Oommachen et al., JPCB 112, 5702 (2008)
            #   "patriksson" : r(D-A) <= 3.5 ang, a(D-H-A) >= 150. deg
            #         from Patrikkson et al., IJMS 248, 124 (2006)
            #   "custom" :  r(H-A) and a(D-H-A) values passed to this function
            #         default is identical with 'torshin'
            #
            # This version uses the Hbond object and returns a list of them.  KKI 4/2/2015
            #
            crit = {'torshin' : ['HA', 2.5, 90.],
                'oommachen' : ['DA', 4.5, 90.],
                'patriksson' : ['DA', 3.5, 150.],
                'custom' : ['HA', rHA, aDHA]
                }
            hbonds = []
            n = self.natom()
            dmat = self.distmat()
            connex = self.connection_table()
            haccept = [ [] ] * n   # each element is a list of acceptors for a D or H ("custom")
            newhb = Hbond(-1, -1, -1)  # nonsense value
            for i in range(n):
                # loop over all atoms in the Molecule, as H-bond donors
                if not self.atom[i].hbonding():
                    continue
                # look for H atoms attached to atom i
                hi = []  # list of H atoms on potential donors
                for ih in range(n): 
                    if connex[ih, i] and self.atom[ih].el == 'H':
                        hi.append(ih)
                for j in range(n):
                    # look for potential acceptors
                    if j == i:
                        # acceptor and donor can't be the same atom
                        continue
                    if not self.atom[j].hbonding():
                        continue
                    for ih in hi:
                        # check distance criteria
                        rh = dmat[ih, j]
                        r = dmat[i, j]
                        for meth in crit.keys():
                            if crit[meth][0] == 'HA' and rh > crit[meth][1]:
                                # H and Acceptor are too distant
                                continue
                            if crit[meth][0] == 'DA' and r > crit[meth][1]:
                                # Donor and Acceptor are too distant
                                continue
                            # distance is OK for this method; now check the angle
                            ang = self.angle(i, ih, j, units='deg')
                            if ang > crit[meth][2]:
                                # angle is also OK; save this hydrogen bond
                                if i != newhb.donor or ih != newhb.Hatom or j != newhb.acceptor:
                                    # it's a new hydrogen bond
                                    newhb = Hbond(i, ih, j, r, rh, ang)
                                # add this criterion's name to list (whether new or existing H-bond)
                                newhb.criteria.append(meth)
                        # done looping over criteria; add this H-bond to the list
                        if newhb.donor >= 0:
                            hbonds.append(newhb)
                            newhb = Hbond(-1, -1, -1)
            return hbonds
        def hydrogen_bonds_old2(self, rHA=2.5, aDHA=90., score=False, r0=1.6):
            # Analyze H-bonding by multiple methods:
            #   "torshin" : r(H-A) <= 2.5 ang, a(D-H-A) >= 90. deg
            #         from Torshin et al., Protein Eng. 15(5), 359-363 (2002)
            #   "oommachen" : r(D-A) <= 4.5 ang, a(D-H-A) >=90. deg
            #         from Oommachen et al., JPCB 112, 5702 (2008)
            #   "patriksson" : r(D-A) <= 3.5 ang, a(D-H-A) >= 150. deg
            #         from Patrikkson et al., IJMS 248, 124 (2006)
            #   "custom" :  r(H-A) and a(D-H-A) values passed to this function
            #         default is identical with 'torshin'
            # return [donor, h-atom, acceptor, distance(H-A)/ang, angle/deg] for each 
            # rHA is the DH-A distance, angstrom
            # aDHA is the D-H-A angle, degree
            #
            # If score==True, also return H-bonding score I made up:
            #   "Hscore" : exp(1-r/r0)*[1-sin(theta)]**2 where theta = DHA angle
            #         This was inspired partly by Fig. 2 in 
            #         Kortemme et al., J. Mol. Biol. 326, 1239-59 (2003)
            #   Scores are summed over all the "custom" H-bonds
            # KKI 3/13/2015
            #
            crit = {'torshin' : ['HA', 2.5, 90.],
                'oommachen' : ['DA', 4.5, 90.],
                'patriksson' : ['DA', 3.5, 150.],
                'custom' : ['HA', rHA, aDHA]
                }
            scoretype = ['Hscore2'] # only one so far
            hpair = {}
            for meth in crit.keys():
                hpair[meth] = []
            n = self.natom()
            dmat = self.distmat()
            connex = self.connection_table()
            sc = {}
            for st in scoretype:
                sc[st] = 0.0
            haccept = [ [] ] * n   # each element is a list of acceptors for a D or H ("custom")
            for i in range(n):
                # loop over all atoms in the Molecule, as H-bond donors
                if not self.atom[i].hbonding():
                    continue
                # look for H atoms attached to atom i
                hi = []  # list of H atoms on potential donors
                for ih in range(n): 
                    if connex[ih, i] and self.atom[ih].el == 'H':
                        hi.append(ih)
                for j in range(n):
                    # look for potential acceptors
                    if j == i:
                        # acceptor and donor can't be the same atom
                        continue
                    if not self.atom[j].hbonding():
                        continue
                    for ih in hi:
                        # check distance criteria
                        rh = dmat[ih, j]
                        r = dmat[i, j]
                        for meth in crit.keys():
                            if crit[meth][0] == 'HA' and rh > crit[meth][1]:
                                # H and Acceptor are too distant
                                continue
                            if crit[meth][0] == 'DA' and r > crit[meth][1]:
                                # Donor and Acceptor are too distant
                                continue
                            # distance is OK for this method; now check the angle
                            ang = self.angle(i, ih, j, units='deg')
                            if ang > crit[meth][2]:
                                # angle is also OK; save this pair
                                hpair[meth].append([i, ih, j, rh, ang])
                                if meth == 'custom':
                                    # keep track of acceptors
                                    haccept[i].append(j)
                                    haccept[ih].append(j)
                                    # calculate 'Hscore' score
                                    sc['Hscore2'] += exp(1-(rh/r0)) * (1-sin(ang*pi/180))**2
                                    # if other score types are tested, put them here
            if score:
                return hpair, sc
            else:
                return hpair
        def hydrogen_bonds_old(self, rHA=2.5, aDHA=90., score=False, r0=1.6):
            # return a list of [donor, acceptor, distance(H-A)/ang, angle/deg]
            # default criteria from Torshin et al., Protein Eng. 15(5), 359-363 (2002)
            # rHA is the DH-A distance criterion, angstrom
            # aDHA is the D-H-A angle criterion, degree
            # if score==True, also return an H-bonding score I made up inspired
            # partly by Fig. 2 in Kortemme et al., J. Mol. Biol. 326, 1239-59 (2003)
            #   sc = SUM{ exp(1-r/r0)*[1-sin(theta)]**2 } where theta = DHA angle
            hpair = []
            n = self.natom()
            dmat = self.distmat()
            connex = self.connection_table()
            sc = 0.0
            for i in range(n):
                if not self.atom[i].hbonding():
                    continue
                # look for H atoms attached to atom i
                hi = []  # list of H atoms on potential donors
                for ih in range(n): 
                    if self.atom[ih].el == 'H' and connex[ih, i]:
                        hi.append(ih)
                for j in range(i):
                    if not self.atom[j].hbonding():
                        continue
                    # is atom i a donor and atom j an acceptor?
                    for ih in hi:
                        r = dmat[ih, j]
                        if r < rHA:
                            # distance is OK
                            ang = self.angle(i, ih, j, units='deg')
                            if ang > aDHA:
                                # angle is also OK; save this pair
                                hpair.append([i, j, r, ang])
                                sc += exp(1-(r/r0)) * (1-sin(ang*pi/180))**2
                    # look for H atoms attached to atom j
                    hj = []
                    for ih in range(n):
                        if self.atom[ih].el == 'H' and connex[ih, j]:
                            hj.append(ih)
                    # is atom j a donor and atom i an acceptor?
                    for ih in hj:
                        r = dmat[ih, i]
                        if r < rHA:
                            # distance is OK
                            ang = self.angle(j, ih, i, units='deg')
                            if ang > aDHA:
                                # angle is also OK; save this pair
                                hpair.append([j, i, r, ang])
                                sc += exp(1-(r/r0)) * (1-sin(ang*pi/180))**2
            if score:
                return hpair, sc
            else:
                return hpair
        def list_angles( self ):
                # return a list of bond angles
                # each element of the list is a list: [i, j, k, degrees]
                connex = self.connection_table()
                natom = self.natom()
                anglist = []
                for i in range ( natom ):
                        # find atoms that are bonded to atom i
                        bonds = filter( lambda j: connex[i][j], range( natom ) )
                        if len( bonds ) > 1:
                                # atom i is a vertex
                                for ang in itertools.combinations( bonds, 2 ):
                                        deg = self.angle( ang[0], i, ang[1], units="deg" )
                                        anglist.append( [ (ang[0], i, ang[1]), deg ] )
                return anglist
        def list_dihedrals( self ):
                anglist = self.list_angles()
                dihelist = []
                # find pairs of angles that have two atoms in common
                for i in range( len(anglist) ):
                        if anglist[i][1] > 179.0:
                                # this is too close to linear; don't build a dihedral
                                continue
                        for j in range( i ):
                                if anglist[j][1] > 179.0:
                                        # too close to linear
                                        continue
                                olap = list( set(anglist[i][0]) & set(anglist[j][0]) )  # intersection
                                both = list( set(anglist[i][0]) | set(anglist[j][0]) )  # union
                                if len(olap) == 2:
                                        # this pair of angles forms a dihedral
                                        if anglist[i][0][1] == anglist[j][0][1]:
                                                # the middle atoms are the same, so this is a branched dihedral
                                                b = anglist[i][0][1]
                                                # the first atom is the other one common to both angles
                                                a = ( olap[0] if b == olap[1] else olap[1] )
                                                [c, d] = sort( list( set(both) - set(olap) ) )
                                                val = self.dihedral(a, b, c, d, type='branched', units='deg' )
                                        else:
                                                # this is a linear dihedral
                                                [a, d] = sort( list( set(both) - set(olap) ) ) # unique atoms are terminal
                                                if a in anglist[i][0]:
                                                        if a == anglist[i][0][0]:
                                                                [b, c] = anglist[i][0][1:]
                                                        else:
                                                                [c, b] = anglist[i][0][:2]
                                                else:
                                                        if a == anglist[j][0][0]:
                                                                [b, c] = anglist[j][0][1:]
                                                        else:
                                                                [c, b] = anglist[j][0][:2]
                                                val = self.dihedral(a, b, c, d, type='linear', units='deg' )
                                        #print '\t', a, b, c, d
                                        dihelist.append( [ (a, b, c, d), val ] )
                return dihelist
##
def xyz2Geometry( atnos, xyzs ):
    # args: list of atomic numbers; list of coordinates [x1, y1, z1, x2, y2, z2,...]
    # return a Geometry
    # 9/16/2014
    #
    # check for compatible list lengths
    natom = len( atnos )
    nxyz = len( xyzs )
    if nxyz != 3 * natom:
        print 'Incompatible numbers of atoms and of coordinates:'
        print 'natom =', natom, 'nxyz =', nxyz
        return Null
    # build Geometry one Atom at a time
    molecule = Geometry()
    for i in range( natom ):
        atno = atnos[i]
        xyz = xyzs[3*i:3*i+3]
        atom = xyz2Atom( atno, xyz )
        molecule.addatom( atom )
    return molecule
##
def isblank( line ):
        # return True if line is blank, else False
        mch = re.compile( r'^\s*$' ).match( line )
        return mch
##
def vib_calc( mass, fc ):
        # given list of masses and 2D numpy.array of cartesian force constants,
        # return list of vibrational frequencies (cm^-1) and mode vectors
        # last change 4/6/2011
        mwt = []
        for m in mass:
                mwt.extend( [1/math.sqrt(m)] * 3 )
        wmat = outer( mwt, mwt )
        # apply the mass-weighting matrix to the force constants
        wfc = multiply( fc, wmat )
        wfc /= amu_au
        eigval, eigvec = linalg.eigh( wfc )
        esign = map( signum, eigval )   # save the sign of each eigenvalue
        eigval = multiply( esign, eigval )      # all values are now positive
        eigval = map( math.sqrt, eigval )
        eigval = multiply( esign, eigval )      # imaginary frequencies are "negative"
        eigval *= au_wavenumber;
        # make eigenvectors the rows
        eigvec = eigvec.T
        return eigval, eigvec
##
def mass_unweight( mass, vibmode ):
        # given atomic masses and normal-mode eigenvectors,
        # return list of reduced masses and pure cartesian mode vectors
        # 4/7/2011
        mwt = []
        for m in mass:
                mwt.extend( [1/math.sqrt(m)] * 3 )
        wmat = outer( mwt, mwt )
        vec = multiply( vibmode, wmat )
        # renormalize and get reduced masses
        redm = []
        for i in range( len(vibmode) ):
                s = dot( vec[i], vec[i] )
                redm.extend( [s] )
                s = sqrt( s )
                vec[i] /= s
        redm = array( redm )
        redm = 1 / redm
        print redm
        return redm, vec
##
def read_zmt( fhandl ):
        # read a z-matrix from a file that contains nothing else
        # (but might start with a Gaussian-style charge-multiplicity pair)
        # Can be comma- or space-delimited.
        # Return zmat (list of lists) and zvar (dictionary) separately
        # 8/27/2008; last change 1/9/09
        inmat = invar = False   # flags
        zmat = []
        zvar = {}
        fhandl.seek(0)  # rewind file, just in case
        for line in fhandl.readlines():
                if isblank( line ):
                        if inmat:
                                inmat = False
                        if invar:
                                # all done
                                break
                else:
                        # not a blank line
                        line = line.strip()     # remove leading  and trailing whitespace
                        if re.compile( "=" ).search( line ):
                                # equals sign indicates a variable assignment
                                invar = True
                                row = line.split( "=" )
                                zvar[ row[0] ].append( float( row[1] ) )
                        else:
                                # a z-matrix line
                                line = line.replace( ',', ' ' ) # replace commas with spaces
                                row = line.split()
                                if len( row ) == 2:
                                        # Gaussian-style charge, multiplicity pair; ignore
                                        continue
                                # convert odd fields to integers
                                for i in range( len(row) ):
                                        if i % 2 == 1:
                                                row[i] = int( row[i] )
                                        elif i > 0:
                                                # name of a z-matrix variable
                                                nom = row[i]
                                                if i == 2:
                                                        zvar[ nom ] = [ 'dist' ]        # a distance
                                                elif i == 4:
                                                        zvar[ nom ] = [ 'angl' ]        # a bond angle
                                                elif i == 6:
                                                        # a dihedral angle; remove any negative sign
                                                        nom = nom.replace( '-', '' )    
                                                        zvar[ nom ] = [ 'dihe' ]        # a dihedral angle
                                                else:
                                                        zvar[ nom ] = [ 'unkn' ]        # shouldn't happen
                                zmat.append( row )
        return zmat, zvar
##
class MolFile:
        # basic data from an MDL MOL file, with masses (atomic weights) added to atoms 
        # 3/21/2012
        def __init__( self, header, nbond, geom, ctable ):
                self.header = header    # list of 3 lines of text
                self.nbond = int( nbond )
                self.geom = geom        # a Geometry object
                self.ctable = ctable    # a 2D array indicating connectivity (non-zero for a connection)
        def copy( self ):
                newMolFile = MolFile( self.header, self.natom, self.nbond, self.geom, self.ctable )
                return newMolFile
##
def read_MOL( fhandl ):
        # read the interesting data from an MDL MOL file (not strict about specification)
        # return a MolFile object (5/16/2012)
        iline = 0
        natom = 0
        nbond = 0
        fhandl.seek(0)
        geom = Geometry()
        header = [ '', '', '' ]
        while 1:
                line = fhandl.readline()
                if not line:
                        break
                if iline < 3:
                        header[iline] = line
                elif iline == 3:
                        # '%3d' formatted fields may run together
                        natom = int( line[0:3] )
                        nbond = int( line[3:6] )
                        ctable = zeros( (natom, natom) )
                elif iline < (4 + natom):
                        # an atom's coordinates and elemental symbol
                        sline = line.split()
                        xyz = array( sline[0:3] )
                        el = sline[3]
                        z = elz( el )
                        m = atomic_weight( z )  # this is not from the MOL file
                        atom = Atom( el, xyz[0], xyz[1], xyz[2], m )
                        geom.addatom( atom )
                elif iline < (4 + natom + nbond):
                        # the connection table
                        i = int( line[0:3] ) - 1
                        j = int( line[3:6] ) - 1
                        b = int( line[6:9] )
                        ctable[i][j] = b
                        ctable[j][i] = b
                # ignore any extra lines
                iline += 1
        if natom != geom.natom():
                print 'Warning: inconsistent atom count in MOL file'
                print 'natom =', natom, ' but geom.natom =', geom.natom()
        rval = MolFile( header, nbond, geom, ctable )
        return rval
##
def read_xyz( fhandl ):
    # read an XMol XYZ file, return a Geometry object and the file comment
    natom = 0
    lineno = 0
    fhandl.seek(0)
    geom = Geometry()
    comment = ''
    for line in fhandl.readlines():
        lineno += 1
        if lineno == 1:
            natom = int( line.strip() )
        elif lineno == 2:
            comment = line
        else:
            # read an atom and its cartesians
            sline = line.split()
            if len( sline ) < 4:
                continue
            el = sline[0]
            z = elz( el )
            m = atomic_weight( z )
            xyz = array( sline[1:4] )
            atom = Atom( el, xyz[0], xyz[1], xyz[2], m )
            geom.addatom( atom )
    if natom != geom.natom():
            print 'Warning: inconsistent atom count in XYZ file'
            print 'natom =', natom, ' but geom.natom =', geom.natom()
    return (geom, comment)
##
def xyz2zvar( zmt, geom, units="rad" ):
        # given a z-mat and geometry object, determine
        # the values of z-matrix coordinates.
        # ignore symmetry and modify variable names as needed
        # return a dictionary of var. names and values
        # 8/27/2008; last change 1/7/09
        print '%d atoms in input zmt' % len(zmt)
        print '%d atoms in input geom' % geom.natom()
        #print 'input zmt =', zmt
        zvar = {}       # variable names as keys
        for i in range( 1, len( zmt ) ):
                # row number 0 has no variables so skip it
                zname = zmt[ i ][2]     # name for this distance
                j = zmt[ i ][1] - 1     # index of atom at other end of the distance
                if zvar.has_key( zname ):
                        # this variable has already been assigned; append to the list of values
                        zvar[ zname ].append( geom.length( i, j ) )
                else:
                        zvar[ zname ] = [ 'dist', geom.length( i, j ) ] # make it a list instead of a scalar
                if len( zmt[i] ) > 3:
                        # there is also an angle
                        zname = zmt[ i ][4]
                        k = zmt[ i ][3] - 1     # index of third atom in angle
                        if zvar.has_key( zname ):
                                # add to list of derived values for this angle
                                zvar[ zname ].append( geom.angle( i, j, k, units ) )
                        else:
                                zvar[ zname ] = [ 'angl', geom.angle( i, j, k, units ) ]
                if len( zmt[i] ) > 5:
                        # a dihedral
                        zname = zmt[ i ][6]
                        l = zmt[ i ][5] - 1     # index of fourth atom in dihedral
                        if re.compile( '-' ).match( zname ):
                                # a negative-signed variable
                                signum = -1
                                zname = zname.replace( '-', '' )
                        else:
                                signum = +1
                        if zvar.has_key( zname ):
                                # add to list of values for this dihedral
                                zvar[ zname ].append( geom.dihedral( i, j, k, l, type='linear', units=units ) * signum )
                        else:
                                zvar[ zname ] = [ 'dihe', geom.dihedral( i, j, k, l, type='linear', units=units ) * signum ]
        return zvar
##
def choose_val( kee, zvar ):
        # return the value of the specified z-matrix variable
        # if it's not in the list, expect it to be numeric constant
        # if not, die
        if kee in zvar.keys() and len( zvar[ kee ] ) > 1:
                return( zvar[kee][ 1 ] )
        # check for a negative-signed variable name
        if kee[0] == '-':
                kee = kee[1:]   # strip the leading negative sign
                return( -zvar[kee][ 1 ] )
        try:
                v = float( kee )
                return v
        except:
                print 'Expected a number but found', kee, 'in z-matrix'
                exit(1)
##
def zmt2xyz( zmt, zvar, mass, verbose='default' ):
        # given a z-matrix and values of its variables, and a list of atomic masses,
        # return a (cartesian) geometry object
        # each 'row' of zvar starts with a var type ('dist' or 'angl' or 'dihe')
        # 11/5/08
        if not verbose == 'silent':
                print '%d atoms in input zmt' % len( zmt )
        cart = Geometry()
        for i in range( len( zmt ) ):
                # generate atomic coordinates sequentially
                row = zmt[ i ]  # current line of the z-matrix
                if i == 0:
                        # put first atom at the origin
                        cart.addatom( Atom( row[0], 0, 0, 0, mass[i] ) )
                        continue
                if i == 1:
                        # put second atom along the z-axis
                        z = zvar[ row[2] ][ 1 ] # value of the z-matrix distance variable
                        cart.addatom( Atom( row[0], 0, 0, z, mass[i] ) )
                        continue
                if i == 2:
                        # put third atom in xz plane
                        r = zvar[ row[2] ][ 1 ] # distance to connected atom
                        a = zvar[ row[4] ][ 1 ] # angle to distal atom
                        a = math.radians( a )   # convert angle to radians
                        j = row[ 1 ] - 1        # index of connected atom
                        [x, y, z] = list( cart.atom[ j ].xyz )
                        if row[ 1 ] == 1:
                                # connected to atom at the origin; distal atom is along positive z-axis
                                z += r * math.cos( a )
                        else:
                                # connected to atom on positive z-axis; distal atom is at origin
                                z -= r * math.cos( a )
                        x += r * math.sin( a )
                        cart.addatom( Atom( row[0], x, 0, z, mass[i] ) )
                        continue
                # i >= 3
                j = row[ 1 ] - 1        # index for connected atom
                k = row[ 3 ] - 1        # index for next-neighboring atom
                l = row[ 5 ] - 1        # index for distal atom
                r = zvar[ row[2] ][ 1 ] # distance to atom j
                a = zvar[ row[4] ][ 1 ] # angle to atom k
                d = choose_val( row[6], zvar )  # dihedral to atom l
                a = math.radians( a )
                d = math.radians( d )
                loc = cart.atom[ j ].xyz.copy() # position of connected atom
                # determine the vector from atom j to the new atom i
                # first component/basis vector (parallel to j-k)
                n1 = cart.normvec( j, k )
                displ = n1 * math.cos( a )
                # second component (~parallel to k-l)
                n2 = cart.normvec( k, l )
                n2 = n2 - dot( n2, n1 ) * n1    # project out the component along n1
                n2 = nvec( n2 ) # renormalize
                displ += n2 * math.sin( a ) * math.cos( d )
                # third component (perpendicular to other two)
                n3 = cross( n1, n2 )
                displ += n3 * math.sin( a ) * math.sin( d )
                loc += r * displ
                [x, y, z] = list( loc )
                cart.addatom( Atom( row[0], x, y, z, mass[i] ) )
        if verbose != 'silent':
                cart.prt()
        return cart
##
def dx_dzvar( zmt, zvar, mass ):
        # using finite (double-sided) differences, compute the derivatives of the cartesian
        # coordinates with respect to the z-matrix variables
        # return ordered list of z-matrix variables and corresdponding 2D numpy array of derivatives
        # 1/9/09
        step = { 'dist' : 0.001, 'angl' : 0.1, 'dihe' : 0.1 }   # step sizes in angstrom, degrees
        zv = zvar.copy()        # local copy of dictionary
        zlis = []       # ordered list of z-matrix variables
        dlis = []
        j = 0   # index for this z-matrix variable
        for z in zvar.keys():
                if len( zvar[ z ] ) < 2:
                        # this is probably a numeric constant; skip it
                        continue
                print 'j =', j, 'z =', z, 'zlis =', zlis
                zlis.append( z )
                zv[ z ] = zvar[ z ][ : ]        # local copy of list
                # increment this zvar 
                zv[ z ][ 1 ] += step[ zvar[ z ][ 0 ] ]
                geom_plus = zmt2xyz( zmt, zv, mass, 'silent' )
                # decrement this zvar
                zv[ z ][ 1 ] -= 2 * step[ zvar[ z ][ 0 ] ]
                geom_minus = zmt2xyz( zmt, zv, mass, 'silent' )
                delta = geom_plus.xyzvec() - geom_minus.xyzvec()
                delta /= 2 * step[ zvar[ z ][ 0 ] ]
                dlis.append( delta )
                j += 1
        dxdz = array( dlis )
        return zlis, dxdz
##
def r0_ref( elem1, elem2 ):
        # return single-bonded distances between elements (Angstrom)
        # from b3lyp/6-31g* calculations on molecules specified (3/2/10)
        #   added covalent radii 3/21/2012
        if ( elem1 > elem2 ):
                # put the elements in ascending lexical order
                t = elem1
                elem1 = elem2
                elem2 = t
        if elem1 == 'C':
                if elem2 == 'C':
                        # C-C bond from C2H6
                        return 1.5306
                if elem2 == 'H':
                        # C-H bond from CH4
                        return 1.0936
                if elem2 == 'N':
                        # C-N bond from CH3NH2
                        return 1.4658
                if elem2 == 'O':
                        # C-O bond from CH3OH
                        return 1.4192
        if elem1 == 'H':
                if elem2 == 'H':
                        # H-H bond from H2
                        return 0.743
                if elem2 == 'N':
                        # N-H bond from CH3NH2
                        return 1.0189
                if elem2 == 'O':
                        # O-H bond from CH3OH
                        return 0.9691
        if elem1 == 'N':
                if elem2 == 'N':
                        # N-N bond from N2H4
                        return 1.4374
                if elem2 == 'O':
                        # N-O bond from NH2OH
                        return 1.4481
        if elem1 == 'O':
                if elem2 == 'O':
                        # O-O bond from HOOH
                        return 1.456
        # unknown case; estimate from rough covalent radii
        z1 = elz( elem1 )
        z2 = elz( elem2 )
        r1 = atomic_radius( z1 )
        r2 = atomic_radius( z2 )
        rsum = r1 + r2
        return rsum
##
def elz( ar ):
        # return atomic number given an elemental symbol, or
        # return elemental symbol given an atomic number (3/9/10)
        symb = [ 'n', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
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
def atomic_weight( iz ):
    # return atomic weight given Z (3/21/2012) or elemental symbol (9/16/2014)
    # values are from the NIST 2003 periodic table
    wt = [ 0, 1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
        22.989770, 24.3050, 26.981538, 28.0855, 30.973761, 32.076, 35.453, 39.948,
        39.0983, 40.078, 44.955910, 47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.933200, 58.6934,
        63.546, 65.409, 69.723, 72.64, 74.92160, 78.96, 79.904, 83.798,
        85.4678, 87.62, 88.90585, 91.224, 92.90638, 95.94, 98, 101.07, 102.90550, 106.42,
        107.8682, 112.411, 114.818, 118.710, 121.760, 127.60, 126.90447, 131.293,
        132.90545, 137.327,
        138.9055, 140.116, 140.90765, 144.24, 145, 150.36, 151.964, 157.25, 158.92534,
        162.500, 164.93032, 167.259, 168.93421, 173.04, 174.967,
        178.49, 180.9479, 183.84, 186.207, 190.23, 192.217, 195.078,
        196.96655, 200.59, 204.3833, 207.2, 208.98038, 209, 210, 222,
        223, 226,
        227, 232.0381, 231.03588, 238.02891, 237, 244, 243, 247, 247,
        251, 252, 257, 258, 259, 262,
        261, 262, 266, 264, 277, 268 ]
    if type( iz ) == int:
        return wt[ iz ]
    else:
        # probably an elemental symbol
        z = elz( iz )
        return wt[ z ]
##
def atomic_radius( iz ):
        # return covalent atomic radius given Z (3/21/2012)
        # values are from Wikipedia (attributed to Slater 1964);
        #   I filled blanks with a guess (e.g., Z-1 value)
        r = [ 0, 0.25, 0.25, 1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50, 0.50,
                         1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 1.00,
                         2.20, 1.80, 1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35,
                         1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 1.15,
                         2.35, 2.00, 1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40,
                         1.60, 1.55, 1.55, 1.45, 1.45, 1.40, 1.40, 1.40,
                         2.60, 2.15,
                         1.95, 1.85, 1.85, 1.85, 1.85, 1.85, 1.85, 1.80, 1.75,
                         1.75, 1.75, 1.75, 1.75, 1.75, 1.75,
                         1.55, 1.45, 1.35, 1.35, 1.30, 1.35, 1.35,
                         1.35, 1.50, 1.90, 1.80, 1.60, 1.90, 1.90, 1.90,
                         2.80, 2.15,
                         1.95, 1.80, 1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75,
                         1.75, 1.75, 1.75, 1.75, 1.75, 1.75,
                         1.75, 1.75, 1.75, 1.75, 1.75, 1.75 ]
        if type(iz) == int:
            return r[ iz ]
        else:
            # convert symbol to nuclear charge
            z = elz( iz )
            return r[z]
##
def vdw_radius(iz):
    # return van der Waals radius given Z (3/13/2015)
    # values (angstrom) are from Mantina et al., JPCA 113, 5806 (2009)
    # many values are not provided; here they return 0 
    r = [0] * 118
    r[1:21] = [ 1.10, 1.40, 
        1.81, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54,
        2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88,
        2.75, 2.31 ]
    r[31:39] = [ 1.87, 2.11, 1.85, 1.90, 1.83, 2.02, 3.03, 2.49 ]
    r[49:57] = [ 1.93, 2.17, 2.06, 2.06, 1.98, 2.16, 3.43, 2.68 ]
    r[81:89] = [ 1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 3.48, 2.83 ]
    if type(iz) == int:
        return r[ iz ]
    else:
        # convert symbol to nuclear charge
        z = elz( iz )
        return r[z]
##
def beyer_swinehart(freqs, emax=5000):
    # count vibrational states (harmonic; wavenumbers) (10/23/15)
    # grain size taken as 1 wavenumber so that array index = energy
    # args: freqs is list of harmonic frequencies; emax is energy ceiling
    # B-S algorithm modified to include the ground state
    # return: numpy array of integer counts
    emax = int(emax)
    count = zeros(emax)
    count[0] = 1
    ifreqs = map(int, freqs)
    for freq in ifreqs:
        #count[freq] += 1
        for i in range(freq, emax):
            j = count[i-freq]
            count[i] += j
    return count
##
def z_vib(freqs, temp):
    # return harmonic vibrational partition function at specified T
    # freqs should be in wavenumber and temp in kelvin
    # 10/23/15
    z = 1.0
    kt = temp * 0.695040 # temperature in wavenumbers
    for freq in freqs:
        x = - freq / kt
        q = 1 - exp(x)
        z /= q
    return z
##
def FC_vertical(energy, amplit, thresh=0.95):
    # args are numpy arrays
    # my procedure for finding the "vertical" energy within a FC curve
    # return the vertical energy
    amax = max(amplit)
    thresh *= amax
    nul = zeros_like(amplit)
    amod = where(amplit > thresh, amplit, nul)
    evert = dot(energy, amod) / sum(amod)
    return evert
##
def smoothing(x, y, x2, style='gau', width=-1, normalize=True):
    # I'M NOT SURE THIS IS WORKING
    # return smoothed y values for (x,y) data series (numpy arrays)
	#   ouput is over the smoothed range defined by x2
    # no sorting necessary
    # styles: 'exp' for exponential; 'gau' for gaussian
    # width parameter (sigma) defaults to 1% of x-range
    if len(x) != len(y):
        # bad input data
        return None
    xlo = min(x)
    xhi = max(x)
    if (width <= 0):
        width = (xhi - xlo) * 0.01
    y2 = zeros_like(x2)
    for i in range(len(y)):
        dx = x2 - x[i]
        if style == 'gau':
            dx = (dx/width)**2
            t = exp(-dx)
        if style == 'exp':
            dx = abs(dx/width)
            t = exp(-dx)
        if normalize:
            t = t / t.sum()
        y2 = y2 + t * y[i]
    return y2
##
