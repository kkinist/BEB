# Routines specific to Gaussian09.
# Karl Irikura, started 8/26/2008
#
import re
import math
from mol_subs import *	# one of my files
code_name = 'Gaussian09'
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
def read_std_orientation( fhandl ):
	# return atomic numbers and coordinates (as list of floats) (3/9/2010)
	# returns the last set of coordinates
	fhandl.seek( 0 )	# rewind file
	flg = False
	regx_rot = re.compile( 'Rotational constants ' )
	regx_num = re.compile( '\d+' )
	regx_stdo = re.compile( 'Standard orientation:' )
	for line in fhandl:
		mch = regx_rot.search( line )
		if mch:
			# past the relevant text block
			flg = False
		if flg:
			# inside the relevant text block
			mch = regx_num.search( line )
			if mch:
				# read this line
				line.rstrip()
				row = line.split()
				atno.append( row[1] )
				xyz.extend( map( float, row[3:] ) )
		mch = regx_stdo.search( line )
		if mch:
			# start of relevant text block
			flg = True
			xyz = []
			atno = []
	atno = map( int, atno )
	return atno, xyz
##
def read_orientation( fhandl, ioccur=1, coords='input' ):
	# return atomic numbers and coordinates (as list of floats) (2/17/2015)
	# 'ioccur' is which occurence to return (starting from '1' for the first)
	# set ioccur = -1 to avoid file rewind and take next occurrence
	if ioccur != -1:
            fhandl.seek( 0 )	# rewind file
	flg = False
	iblock = 0
	atno = []
	xyz = []
	regx_stop = re.compile( 'Distance matrix|Rotational constants' )
	regx_num = re.compile( '\d+' )
	regx_inpo = re.compile( 'Input orientation:' )
        if coords == 'standard':
                regx_inpo = re.compile( 'Standard orientation:' )
	for line in fhandl:
		mch = regx_stop.search( line )
		if mch:
			# past a 'standard orientation' text block
			flg = False
			if len(atno) > 0:
			    # all done
			    break
		if flg and (iblock == ioccur or ioccur == -1):
			# inside the requested text block
			mch = regx_num.search( line )
			if mch:
				# read this line
				line = line.rstrip()
				row = line.split()
				atno.append( row[1] )
				xyz.extend( map( float, row[3:] ) )
		mch = regx_inpo.search( line )
		if mch:
			# start of relevant text block
			flg = True
			iblock += 1
			xyz = []
			atno = []
	atno = map( int, atno )
	return atno, xyz
##
def read_mulliken( fhandl, hflag ):
        # read (last) mulliken charges in file
        # return list of element symbols and list of charges
        # 11/27/2012
	fhandl.seek( 0 )	# rewind file
	flg = False
        regx_end = re.compile( 'Sum of Mulliken' )
        if hflag:
                # read the charges that have H atoms folded in
                regx_start = re.compile( 'Mulliken charges with hydrogens summed into heavy atoms:' )
        else:
                regx_start = re.compile( 'Mulliken atomic charges:' )
        for line in fhandl:
                if flg:
                        # in a mulpop section
                        mch = regx_end.search( line )
                        if mch:
                                # end of this mulpop section
                                flg = False
                        else:
                                # continue reading mulpop section
                                field = line.split()
                                if len( field ) == 3:
                                        # read this line
                                        el.append( field[1] )
                                        charge.append( float( field[2] ) )
                else:
                        # not in a mulpop section
                        mch = regx_start.search( line )
                        if mch:
                                # starting a mulpop section
                                el = []
                                charge = []
                                flg = True
        return ( el, charge )
##
def read_masses( fhandl ):
	# read frequency file 
	# return number of atoms and list of atomic masses
	fhandl.seek( 0 )	# rewind file
	m = []
	regx = re.compile( 'has atomic number' )
	for line in fhandl:
		mch = regx.search( line )
		if mch:
			# extract three numbers: atom no., atomic no., mass
			nums = re.compile( r'[\.\d]+' ).findall( line )
			m.append( float( nums[2] ) )
	return ( int(nums[0]), m )
##
def read_link0(fhandl):
    # return the Link0 part of a Gaussian input file as a list
    fhandl.seek(0)  # rewind file
    ln0 = []
    regx = re.compile('\s*(%.*)')
    regstop = re.compile('\s*#')
    for line in fhandl:
        m = regx.match(line)
        if m:
            ln0.append(m.group(1))
        m = regstop.match(line)
        if m:
            # command line; done with link0 statements
            break
    return ln0
##
def read_ghz(fhandl):
    # return the last triple of rotational constants
    fhandl.seek(0)
    regx = re.compile(' Rotational constants \(GHZ\):')
    gline = ''
    for line in fhandl:
        m = regx.search(line)
        if m:
            gline = line
    if gline:
        line = gline.split()
        return map(float, line[-3:])
    else:
        # failure
        return [0, 0, 0]
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
def arch2xyz( arch ):
	# extract cartesian coordinates from archive from get_arch
	# return as a one-dimensional list of floats
	# also return elemental symbols
	xyz = []
	symb = []
	for line in arch[3].split( r'|' ):
		# split each line at commas
		row = line.split( ',' )
		j = len(row)
		if j in [4,5]:
			xyz.extend( map( float, row[j-3:] ) )
			symb.append( row[0] )
	return symb, xyz
##
def arch2chargespin(arch):
        # extract molecular charge and spin multiplicity from archive (from get_arch)
        # return as [charge, mult] 
        charge = 0
        mult = 1
        fields = arch[3].split(r'|')
        row = fields[0].split(',')
        return map(int, row)
##
def arch2theorybasis(arch):
    # return theory and basis set, given archive from get_arch
    fields = arch[0].split(r'|')
    theory = fields[4]
    basis = fields[5]
    return theory, basis
##
def arch2cmd(arch):
    # return command (not including theory/basis) from archive list
    fields = arch[0].split(r'|')
    cmd = fields[3]
    return cmd
##    
def read_Geometry_arch(arch):
        # read archive block, extract coordinates, return Geometry object
        symb, xyz = arch2xyz(arch)
        geom = xyz2Geometry(symb, xyz)
        return geom
##
def arch2fc( arch ):
	# extract cartesian force constants from archive from get_arch
	# return full matrix as one-dimensional list
	row = arch[5].split( ',' )	# split at commas
	row = map( float, row )
	dimen = ( math.sqrt( 8 * len( row ) + 1 ) - 1 ) / 2
	dimen = int( round( dimen ) )	# number of coordinates
	fc = range( dimen * dimen )
	for i in range( dimen ):
		for j in range( dimen ):
			if ( j <= i ):
				n = dimen * i + j
				m = i * (i+1) / 2 + j
				fc[n] = row[m]
				n = dimen * j + i	# other triangle
				fc[n] = row[m]
	return fc
##
def arche( arch ):
    # extract electronic energy from archive from get_arch
    # return hash of values and list of methods (for ordering)
    energy = {}   # hash of energy (name, value) pairs
    method = [] # list of theoretical methods
    #regx = re.compile( '^([a-zA-Z]+)=([-]?\d+\.\d+)$' )
    #regx = re.compile( '^(\w+)=([-]?\d+\.\d+)$' )
    regx = re.compile( '^(\w+(?:\([Tt]\))?)=([-]?\d+\.\d+)$' )
    for line in arch:
        if re.match( 'Version=', line ):
            # look for energy on this line
            line = line.split( r'|' )
            # everything with 'string=' followed by floating-point number
            for pair in line:
                m = regx.match( pair )
                # exclude some 'energies' by name
                if m and m.group(1) not in ['Thermal', 'ZeroPoint']:
                    # report this pair
                    energy[ m.group(1) ] = float( m.group(2) )
                    method.append( m.group(1) )
    return energy, method
##
def write_gjf(fname, link0, cmd, symb, xyz, charge=0, mult=1, comment='comment'):
    # create Gaussian input file
    try:
        fgjf = open(fname, 'w')
    except:
        print 'Problem creating file', fname
        return 'ugh'
    for ln0 in link0:
        # expect a list of subcommands
        fgjf.write(ln0 + '\n')
    fgjf.write('# ' + cmd + '\n')
    fgjf.write('\n' + comment + '\n')
    fgjf.write("\n{:d} {:d}\n".format(charge, mult))
    for i in range(len(symb)):
        j = 3*i
        s = "{:s} {:f} {:f} {:f}\n".format(symb[i], xyz[j], xyz[j+1], xyz[j+2])
        fgjf.write(s)
    fgjf.write('\n')
    fgjf.close()
    return 0
##
def geom_gjf(fname, geom, cmd='PM6', charge=0, mult=1, comment='comment', link0=['']):
    # create Gaussian input file from Geometry object
    symb = geom.symbols()
    xyz = geom.xyzvec().tolist()
    write_gjf(fname, link0, cmd, symb, xyz, charge, mult, comment)
    return 0
##
def read_natoms( fhandl ):
	# read the number of atoms (8/27/12)
	fhandl.seek( 0 )	# rewind file
	natom = 0
	regx = re.compile( ' NAtoms=\s+(\d+)\s' )
	for line in fhandl:
		mch = regx.search( line )
		if mch:
			natom = mch.group(1)
			break
	return int( natom )
##
def read_omega( fhandl, intensities=True ):
    # read harmonic frequencies and intensities from a Gaussian09 output file
    # not for VPT2 outputs; use read_fundamental instead
    # return two lists (8/27/2012) (add intensities 3/19/15)
    fhandl.seek( 0 )	# rewind file
    omega = []
    inten = []
    regx = re.compile( ' Frequencies --' )
    regi = re.compile( ' IR Inten    --' )
    for line in fhandl:
		mch = regx.search( line )
		if mch:
			fields = line.split()
			omega.extend( fields[2:] )
		mch = regi.search(line)
		if mch:
		      fields = line.split()
		      inten.extend(fields[-3:])
    omega = map(float, omega)
    inten = map(float, inten)
    if intensities:
        return omega, inten
    else:
        return omega
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
                        if not cmp( energy[i], energy[j] ):
                                # remove the duplicate entry (j)
                                energy.pop( j )
        return energy
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
def read_zpe( fhandl ):
	# read ZPE from harmonic or VPT2 output file
	fhandl.seek( 0 )
	zpe = 0.0
	regx_har = re.compile( ' Zero-point vibrational energy' )	# in J/mol
	regx_anh = re.compile( ' ZPEtot     =' )	# in cm-1
	for line in fhandl:
		mch = regx_har.search( line )
		if mch:
			line = line.split()
			zpe = float( line[3] ) / 1000.0	# convert to kJ/mol
			zpe *= au_wavenumber / au_kjmol	# convert to cm-1
			continue
		mch = regx_anh.search( line )
		if mch:
			line = line.split()
			zpe = float( line[2] )	# already in cm-1
			continue
	return zpe
##
def read_fundamental( fhandl ):
	# read harmonic and anharmonic frequencies
	fhandl.seek( 0 )	# rewind file
	omega = []
	nu = []
	regx_begin = re.compile( ' Fundamental Bands' )
	regx_end = re.compile( ' Overtones' )
	in_fund = False
	for line in fhandl:
		mch = regx_end.search( line )
		if mch:
			in_fund = False
		if in_fund:
			# parse a line
			line = line.split()
			omega.append( float( line[1] ) )
			nu.append( float( line[2] ) )
			continue
		mch = regx_begin.search( line )
		if mch:
			in_fund = True
	return (omega, nu)
##
def read_overtones( fhandl ):
	# read overtones and combination bands; return them as lists of float
	# they are 2D lists
	fhandl.seek( 0 )
	overt = []
	combi = []
	regx_ov = re.compile( ' Overtones \(DE w.r.t' )
	regx_c = re.compile( ' Combination Bands' )
	regx_end = re.compile( ' ZPE\(harm\) ' )
	in_ov = False
	in_c = False
	for line in fhandl:
		mch = regx_end.search( line )
		if mch:
			break
		if in_ov:
			# read about one (2-0) overtone
			mch = regx_c.search( line )
			if mch:
				# combinations section
				in_ov = False
				in_c = True
				continue
			line = line.replace( '(', ' ' )	# replace left parenthesis by space
			line = line.split()
			if len( line ) != 7:
				# something is wrong with this line
				continue
			a = [ int(line[0]), float(line[3]) ]	# mode number, then anharmonic value
			overt.append( a )
			continue
		if in_c:
			# read about one (1,1') combination band
			line = line.replace( '(', ' ' )
			line = line.split()
			if len( line ) != 9:
				# something is wrong with this line (maybe blank)
				continue
			a = [ int(line[0]), int(line[2]), float(line[5]) ]	# mode numbers, then anharmonic value
			combi.append( a )
			continue
		mch = regx_ov.search( line )
		if mch:
			# overtones section
			in_ov = True
			continue
	return ( overt, combi )
##		
def read_xij( fhandl ):
	# read vibrational anharmonicity constants and ZPE
	# return Xij as 2D list in lower-triangular form
	fhandl.seek( 0 )
	x = []
	x0 = 0.0
	regx_begin = re.compile( ' X matrix of Anharmonic Constant' )
	regx_end = re.compile( ' ZPEharm    =|Max\.Diff\.' )
	regx_float = re.compile( '\d\.\d' )
	regx_x0 = re.compile( ' x0         =' )
	in_x = False
	for line in fhandl:
		mch = regx_end.search( line )
		if mch:
			in_x = False
			continue
		if in_x:
			mch = regx_float.search( line )
			if mch:
				# parse this line
				line = line.split()
				i = int( line[0] )	# 1-adjusted row number
				if i > len( x ):
					x.append( [] )	# add a new, empty row
				x[i-1].extend( line[1:] )	# add to new or existing row
				continue
		mch = regx_begin.search( line )
		if mch:
			in_x = True
			continue
		mch = regx_x0.search( line )
		if mch:
			# read X0 in cm-1
			line = line.split()
			x0 = float( line[2] )
	# convert from string to float
	for i in range( len(x) ):
		x[i] = map( float, x[i] )
	# fill empty elements with zero
	for i in range( len(x) ):
		for j in range( i+1, len(x) ):
			x[i].append( 0.0 )
	return (x, x0)
##
def n_sdd( atno, all=0 ):
        # given Z value (or element symbol) return number of electrons replaced by SDD
        #   pseudopotential (SDDall if all==1)
        if type( atno ) != int:
                # convert symbol to Z value
	        atno = elz( atno )
        if all:
                zmin = 1
        else:
                zmin = 19
        if atno < zmin:
                return 0
        core = {
                3  :  2,
                11 : 10,
                31 : 28,
                49 : 46,
                57 : 28,
                72 : 60,
                81 : 78,
                89 : 60
        }
        npp = 0
        for ki in sorted(core):
                if atno >= ki:
                        npp = core[ki]
        return npp
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
def read_franck_condon(fhandl):
    # Read Franck-Condon curve from Gaussian output file
    # return two lists (x, y)
    fhandl.seek(0)
    rbegin = re.compile('Final Spectrum')
    rend = re.compile('GradGrad')
    flg = False
    energy = []
    amplit = []
    for line in fhandl:
        mch = rend.search(line)
        if mch:
            flg = False
        if flg:
            # Look for data line
            fields = line.replace('D','E').split()
            if len(fields) == 2:
                # could be good
                try:
                    x = float(fields[0])
                    y = float(fields[1])
                    energy.append(x)
                    amplit.append(y)
                except:
                    print 'two bad fields: ', fields
        mch = rbegin.search(line)
        if mch:
            flg = True
    return energy, amplit
