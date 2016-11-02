# Routines specific to Gaussian09.
# Karl Irikura, starting from scratch using python3 and pandas (Sept. 2016)
# Routines are poorly tested, unfit for distribution.
####
import sys
import re
import pandas as pd
import matplotlib.pyplot as plt
####
def read_g09_command(fhandl):
    # read all Gaussian command-lines
    # return a DataFrame: (1) line numbers (base 1),
    #   (2) byte numbers (file positions for seek/tell),
    #   (3) the command-lines of interest
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    cmd = []
    fline = []
    fpos = []
    lineno = 0
    regxdash = re.compile('^\s*-+$')
    regxcmd = re.compile(r'^ ?# ?\w+')
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        m = regxdash.search(line)
        if m:
            # found a dashed line
            # does a command-line follow? 
            line = fhandl.readline()
            lineno += 1
            mcmd = regxcmd.match(line)
            if mcmd:
                # yes, this looks promising
                cmdline = line.strip()
                clineno = lineno
                pos = fhandl.tell()
                # read up to 3 more lines
                for i in range(3):
                    line = fhandl.readline()
                    lineno += 1
                    m = regxdash.search(line)
                    if m:
                        # dashed line indicates end of command-line
                        cmd.append(cmdline)
                        fline.append(clineno)
                        fpos.append(pos)
                        break
                    else:
                        # this is a continuation line
                        cmdline += ' ' + line.strip()
                else:
                    # too many 'continuation lines' -- don't believe it
                    print('Gaussian09 command line too long in "read_g09_command":')
                    print(cmdline)
                    fhandl.seek(byte_start) # restore file pointer to original position
                    return(['confusion'], [0], [0])
    data = list(zip(fline, fpos, cmd))
    cols = ['line', 'byte', 'Command']
    df = pd.DataFrame(data=data, columns=cols)
    fhandl.seek(byte_start) # restore file pointer to original position
    return df
##
def read_g09_rev(fhandl):
    # read and return the first instance of:
    #   major version; revision; early date; late date
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    lineno = 0
    regxdash = re.compile(r'\*{40}')
    regxvers = re.compile(r'^ ?Gaussian (\d+):\s+(\S+)Rev(\S+) (\S+)')
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        m = regxdash.search(line)
        if m:
            # found a line of asterisks
            # does a version line follow?
            line = fhandl.readline()
            lineno += 1
            mvers = regxvers.match(line)
            if mvers:
                # yes
                year = mvers.group(1)
                major = mvers.group(2)
                minor = mvers.group(3)
                date1 = mvers.group(4)
                # read the later date from the next line
                line = fhandl.readline()
                lineno += 1
                date2 = line.strip()
                # check for a closing line of asterisks
                line = fhandl.readline()
                lineno += 1
                m = regxdash.search(line)
                if not m:
                    print('Warning: no closing line of asterisks in "read_g09_rev"')
                break
    fhandl.seek(byte_start) # restore file pointer to original position
    return(year, major, minor, date1, date2)
##
def read_g09_comments(fhandl):
    # read all Gaussian comment lines
    # return a DataFrame: (1) line numbers (base 1),
    #   (2) file positions (bytes),
    #   (3) the user comments as requested
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    comment = []
    fline = []
    fpos = []
    lineno = 0
    pos = 0
    iline = 0
    regxdash = re.compile('^\s*-+$')
    regxtop = re.compile(r'99//99;')
    comm = ''
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        if len(comm) > 0 and iline == 0:
            # the comment block has passed
            comment.append(comm)
            fline.append(clineno)
            fpos.append(pos)
            comm = ''
        if iline == 1:
            # inside a comment block (may occupy multiple lines)
            if len(comm) == 0:
                # first line of comment
                clineno = lineno
                pos = fhandl.tell()
            comm += line.strip()
        if iline > 0:
            # near a comment; look for a long dashed line
            m = regxdash.search(line)
            if m:
                # found a dashed line
                iline -= 1
        m = regxtop.search(line)
        if m:
            # found a route block; comment should follow
            iline = 2
    data = list(zip(fline, fpos, comment))
    cols = ['line', 'byte', 'Comment']
    df = pd.DataFrame(data=data, columns=cols)
    fhandl.seek(byte_start) # restore file pointer to original position
    return df
##
def read_g09_charge_mult(fhandl):
    # read all charge/multiplicity lines, return DataFrame:
    #   (1) line number(s), (2) byte number(s),
    #   (3) charge(s), (4) spin multiplicity(s)
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    charge = []
    mult = []
    fline = []
    fpos = []
    lineno = 0
    regx = re.compile(r'Charge\s+=\s+(-?\d+)\s+Multiplicity\s+=\s+(\d+)')
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        m = regx.search(line)
        if m:
            # found charge and mult
            charge.append(int(m.group(1)))
            mult.append(int(m.group(2)))
            fline.append(lineno)
            fpos.append(fhandl.tell())
    data = list(zip(fline, fpos, charge, mult))
    cols = ['line', 'byte', 'Charge', 'Mult']
    df = pd.DataFrame(data, columns=cols)
    fhandl.seek(byte_start) # restore file pointer to original position
    return df
##
def read_g09_std_orient(fhandl):
    # read all blocks of coordinates labelled "(standard|input) orientation"
    # and return a DataFrame:
    #   (1) line number for the "standard orientation" string,
    #   (2) byte number for same,
    #   (3) unit for the coordinates,
    #   (4) coordinates as a DataFrame (atomic number, x, y, z)
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    coords = []
    units = []
    fline = []
    fpos = []
    xyz = []
    lineno = 0
    regxdash = re.compile('-{60}')
    regx = re.compile('(Standard|Input) orientation:')
    regxunit = re.compile(r'Coordinates\s+\((\w+)\)')
    instd = False
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        if instd:
            # now reading a set of coordinates; count dashed lines
            m = regxdash.search(line)
            if m:
                # found another dashed line
                ndashline += 1
            else:
                # not a dashed line
                if ndashline == 1:
                    # past one dashed line; look for units
                    m = regxunit.search(line)
                    if m:
                        # found units
                        units.append(m.group(1))
                if ndashline == 2:
                    # past two dashed lines; read atomic coordinates: [Z, x, y, z]
                    fields = line.split()
                    xyz.append([int(fields[1]), float(fields[3]),
                        float(fields[4]), float(fields[5])])
                if ndashline == 3:
                    # completely past the coordinates; close out this set
                    df_crd = pd.DataFrame(data=xyz, columns=['Z', 'x', 'y', 'z'])
                    coords.append(df_crd)
                    xyz = []
                    ndashline = 0
                    instd = False
        else:
            # not currently reading a set of coordinates
            m = regx.search(line)
            if m:
                # found a standard orientation
                fline.append(lineno)
                fpos.append(fhandl.tell())
                instd = True    # set flag
                ndashline = 0   # counter
    data = list(zip(fline, fpos, units, coords))
    cols = ['line', 'byte', 'Unit', 'Coordinates']
    df_coord = pd.DataFrame(data=data, columns=cols)
    fhandl.seek(byte_start) # restore file pointer to original position
    return df_coord
##
def read_g09_stoichiometry(fhandl):
    # read stoichiometry string(s) and parse into a dict. Return a DataFrame:
    #   (1) line number, (2) byte number,
    #   (3) stoichiometry string without parenthetical info,
    #   (4) dict of {element: count} pairs
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    raw = []
    stoichlist = []
    stoich = {}
    fline = []
    fpos = []
    lineno = 0
    regx = re.compile(r'Stoichiometry\s+([A-Za-z0-9]+)')
    regxel = re.compile(r'([A-Z][a-z]?)(\d*)')
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        m = regx.search(line)
        if m:
            # found a stoichiometry string
            formula = m.group(1)
            raw.append(formula)
            fline.append(lineno)
            fpos.append(fhandl.tell())
            # parse it into individual elements
            elems = regxel.findall(formula)
            if elems:
                for pair in elems:
                    ct = int('0' + pair[1]) # add string '0' in case it's null
                    if pair[0] in stoich.keys():
                        # this element is already listed; add to its tally
                        stoich[pair[0]] += max(ct, 1)
                    else:
                        # add this element to the dict
                        stoich[pair[0]] = max(ct, 1)
                stoichlist.append(stoich.copy())
                stoich = {}
    data = list(zip(fline, fpos, raw, stoichlist))
    cols = ['line', 'byte', 'Stoich', 'Elements']
    df = pd.DataFrame(data=data, columns=cols)
    fhandl.seek(byte_start) # restore file pointer to original position
    print(df)
    return df
##
def read_g09_rotational(fhandl):
    # read all sets of rotational constants. Return a DataFrame:
    #   (1) line number(s), (2) byte number(s),
    #   (3) measurement unit,
    #   (4) list of rotational constants
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    rotat = []
    unit = []
    fline = []
    fpos = []
    lineno = 0
    regx = re.compile(r'Rotational constants \((\w+)\):')
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        m = regx.search(line)
        if m:
            # found a set of rotational constants
            fline.append(lineno)
            fpos.append(fhandl.tell())
            unit.append(m.group(1))
            # save only non-zero constants
            fields = line.split()
            rset = []
            for c in fields[-3:]:
                cf = float(c)
                if cf > 0:
                    rset.append(cf)
            rotat.append(rset.copy())
    data = list(zip(fline, fpos, unit, rotat))
    cols = ['line', 'byte', 'Unit', 'Rotational Constants']
    df = pd.DataFrame(data=data, columns=cols)
    fhandl.seek(byte_start) # restore file pointer to original position
    return df
##
def read_g09_scfe(fhandl):
    # read all SCF energies and return a DataFrame:
    #   (1) line number(s), (2) byte number(s),
    #   (3) name of SCF method, (4) number of SCF cycles,
    #   (5) SCF energy (hartree)
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    energy = []
    method = []
    ncycle = []
    fline = []
    fpos = []
    lineno = 0
    regx = re.compile(r'SCF Done:\s+E\((\S+)\)\s+=\s+')
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        m = regx.search(line)
        if m:
            # found an SCF energy
            fline.append(lineno)
            fpos.append(fhandl.tell())
            method.append(m.group(1))
            # extract the energy value and the number of scf cycles
            fields = line.split()
            ncycle.append(int(fields[-2]))
            energy.append(float(fields[-5]))
    data = list(zip(fline, fpos, method, ncycle, energy))
    cols = ['line', 'byte', 'Method', 'Cycles', 'Energy']
    df = pd.DataFrame(data=data, columns=cols)
    fhandl.seek(byte_start) # restore file pointer to original position
    return df
##
def read_g09_postHF(fhandl):
    # read all HF and post-HF energies and return a DataFrame:
    #   (1) line number(s), (2) byte number(s),
    #   (3) method name, (4) energy,
    #   repeat #3 and #4 as needed for different post-HF energies
    # Handles (so far): HF, MP2, MP3, CCSD, CCSD(T)
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    energy = []
    method = []
    fline = []
    fpos = []
    lineno = 0
    reghf = re.compile('SCF Done:\s+E\(([RU]?HF)\)\s+=\s+(\S+)')
    regmp2 = re.compile('E2\s+=\s+\S+\s+EUMP2\s+=\s+(\S+)')
    regmp3 = re.compile('E3=\s+\S+\s+EUMP3=\s+(\S+)')
    regccsdt = re.compile('CCSD\(T\)=\s+(\S+)')
    # CCSD is tricky because unconverged energies are printed the same way
    #   as the converged energy, so the context must be checked
    regccsd = re.compile('DE\(Corr\)=\s+\S+\s+E\(CORR\)=\s+(\S+)\s+Delta')
    ccsd_mem = []  # lineno, byte, energy
    regccconv = re.compile('Largest amplitude=')
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        m = reghf.search(line)
        if m:
            # HF energy; save it
            fline.append(lineno)
            fpos.append(fhandl.tell())
            method.append(m.group(1))
            energy.append(float(m.group(2)))
        m = regmp2.search(line)
        if m:
            # MP2 energy; save it
            fline.append(lineno)
            fpos.append(fhandl.tell())
            method.append('MP2')
            energy.append(float(m.group(1).replace('D','E')))
        m = regmp3.search(line)
        if m:
            # MP3 energy; save it
            fline.append(lineno)
            fpos.append(fhandl.tell())
            method.append('MP3')
            energy.append(float(m.group(1).replace('D','E')))
        m = regccsdt.search(line)
        if m:
            # CCSD(T) energy; save it
            fline.append(lineno)
            fpos.append(fhandl.tell())
            method.append('CCSD(T)')
            energy.append(float(m.group(1).replace('D','E')))
        m = regccsd.search(line)
        if m:
            # CCSD energy; remember it but don't save it yet
            ccsd_mem = [lineno, fhandl.tell(), float(m.group(1).replace('D','E'))]
        m = regccconv.search(line)
        if m:
            # CCSD has converged; now save the remembered CCSD info
            fline.append(ccsd_mem[0])
            fpos.append(ccsd_mem[1])
            method.append('CCSD')
            energy.append(ccsd_mem[2])
    data = list(zip(fline, fpos, method, energy))
    cols = ['line', 'byte', 'Method', 'Energy']
    df = pd.DataFrame(data=data, columns=cols)
    fhandl.seek(byte_start) # restore file pointer to original position
    return df
##
def read_g09_identify_scanned(fhandl):
    # identify variables that are being scanned and return a DataFrame:
    #   (1) line number(s), (2) byte number(s),
    #   (3) coordinate label,
    #   (4) coordinate description
    #   (5) initial value of coordinate
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    label = []
    descr = []
    initval = []
    fline = []
    fpos = []
    lineno = 0
    regx = re.compile(r'\s*!\s+(\w+)\s+(\S+)\s+(-?\d*\.\d+)\s+Scan\s+!')
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        m = regx.search(line)
        if m:
            # found a scanned variable
            fline.append(lineno)
            fpos.append(fhandl.tell())
            label.append(m.group(1))
            descr.append(m.group(2))
            initval.append(float(m.group(3)))
    data = list(zip(fline, fpos, label, descr, initval))
    cols = ['line', 'byte', 'Label', 'Description', 'Initial Value']
    df = pd.DataFrame(data=data, columns=cols)
    fhandl.seek(byte_start) # restore file pointer to original position
    return df
##
def read_g09_optim_param(fhandl):
    # read optimized value of coordinates, return a DataFrame:
    #   (1) line number(s) of "optimized parameters" string,
    #   (2) byte number(s),
    #   (3) a DataFrame of the coordinates:
    #       (3a) coordinate label,
    #       (3b) coordinate description,
    #       (3c) optimized value of coordinate
    # Coordinates may be internals or individual cartesians
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    coordset = []
    fline = []
    fpos = []
    lineno = 0
    regx = re.compile(r'\s+!\s+Optimized Parameters\s+!')
    regxdash = re.compile('-{70}')
    regxparen = re.compile(r'[)(]')
    inopt = False
    ndash = 0
    valcols = ['Label', 'Description', 'Value']
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        if inopt:
            # already inside a block of parameters
            if ndash == 2:
                # done with this block
                inopt = False
                crds = list(zip(labelset, descrset, valset))
                coordset.append(pd.DataFrame(data=crds, columns=valcols))
            if ndash == 1:
                # split the line and save fields
                fields = line.split()
                if len(fields) > 3:
                    # a line with real data (not the terminating, dashed line)
                    # check for parentheses
                    if regxparen.search(line):
                        # there are parentheses: internal coordinates
                        labelset.append(fields[1])
                        descrset.append(fields[2])
                        valset.append(float(fields[3]))
                    else:
                        # no parentheses: these are individual cartesians
                        labelset.append(fields[1])
                        valset.append(float(fields[2]))
                        descrset.append('-')
            # check for a delimited dashed line
            m = regxdash.search(line)
            if m:
                # another long dashed line
                ndash += 1
        # check to see if this line contains the main search string
        m = regx.search(line)
        if m:
            # found a block of optimized parameters
            inopt = True
            ndash = 0
            labelset = []
            descrset = []
            valset = []
            fline.append(lineno)
            fpos.append(fhandl.tell())
    data = list(zip(fline, fpos, coordset))
    cols = ['line', 'byte', 'Optimized Coordinates']
    df = pd.DataFrame(data=data, columns=cols)
    fhandl.seek(byte_start) # restore file pointer to original position
    return df
##
def extract_g09_scan(fhandl):
    # assemble scan information, return a sorted DataFrame:
    #   (1) line number(s) of "optimized parameters" string,
    #   (2) byte number(s) of same,
    #   (3) a tuple of names of the scanned variable(s),
    #   (4) a tuple of descriptions of the scanned variable(s),
    #   (5*) column(s) containing the variable value(s),
    #   (6) the energy,
    #   (7) a DataFrame of the cartesian coordinates:
    #       (7a) atomic number,
    #       (7b-d) x, y, z
    #--------------------------------------
    # find which variables were scanned
    byte_start = fhandl.tell()
    df_scancrd = read_g09_identify_scanned(fhandl)
    nvar = len(df_scancrd)
    if nvar < 1:
        sys.exit('No scanned variables were identified!')
    # discard the line and byte numbers before display
    df_scancrd.drop(['line', 'byte'], axis=1, inplace=True)
    # get the values of internal coordinates through the scan
    #   use these line and byte numbers
    df_optcrd = read_g09_optim_param(fhandl)
    npoints = len(df_optcrd)
    # discard the coordinates that were not scanned
    for i in range(npoints):
        point = df_optcrd.iloc[i]['Optimized Coordinates']
        point = point.loc[point.Description.isin(df_scancrd.Description)]
        if len(point) < 1:
            sys.exit('No scanned variables found among coordinates!')
        # now replace the full set of coordinates with only the scanned coordinates
        df_optcrd.set_value(i, 'Optimized Coordinates', point)
    # In each record of df_optcrd, replace the 'Optimized Coordinates' dataframe
    #   with three or more columns: (var labels), (var descriptions), var values [one column per var]
    tuplbl = list(tuple(df_optcrd.iloc[i]['Optimized Coordinates']['Label'].tolist()) for i in range(npoints))
    tupdescr = list(tuple(df_optcrd.iloc[i]['Optimized Coordinates']['Description'].tolist()) for i in range(npoints))
    tupval = list(tuple(df_optcrd.iloc[i]['Optimized Coordinates']['Value'].tolist()) for i in range(npoints))
    del df_optcrd['Optimized Coordinates']
    df_optcrd['Label'] = tuplbl
    df_optcrd['Description'] = tupdescr
    # split up the tuple of variable values and assign each to a new column
    vals = list(zip(*tupval))
    for ivar in range(nvar):
        colname = tupdescr[0][ivar]
        df_optcrd[colname] = vals[ivar]
    # Get the converged SCF energies
    df_scfe = read_g09_scfe(fhandl)
    nenergy = len(df_scfe)
    # Get the converged Cartesian geometries
    df_cart = read_g09_std_orient(fhandl)
    ncart = len(df_cart)
    #
    # For each geometry, choose the preceding energy/cartesians with the largest line number 
    #   this is to match the scanned parameters with the right energies and cartesians
    #
    nrg = []
    cart = []
    ienergy = 0  # index for searching df_scfe
    icart = 0  # index for search df_cart
    for igeom in range(npoints):
        # loop over converged geometries
        linetarget = df_optcrd.iloc[igeom].line
        while (ienergy < nenergy) and (df_scfe.iloc[ienergy].line < linetarget):
            ienergy += 1
        # line number 'ienergy' is now too large by one
        nrg.append(df_scfe.iloc[ienergy-1].Energy)
        while (icart < ncart) and (df_cart.iloc[icart].line < linetarget):
            icart += 1
        # line number 'icart' is now too large by one
        cart.append(df_cart.iloc[icart-1].Coordinates)
    df_optcrd['Energy'] = nrg  # add the energies to the DataFrame
    df_optcrd['Coordinates'] = cart  # add the Cartesians
    fhandl.seek(byte_start) # restore file pointer to original position
    return df_optcrd
##
def combine_g09_scans(filenames, relative=0, outfile=None, plot=False):
    # combine scan information, return a sorted DataFrame:
    #   (1) a tuple of names of the scanned variable(s),
    #   (2) a tuple of descriptions of the scanned variable(s),
    #   (3*) column(s) containing the corresponding variable values,
    #   (4) the energy,
    #   (5) filename
    #   (6) a DataFrame of the cartesian coordinates:
    #       (6a) atomic number,
    #       (6b-d) x, y, z
    #--------------------------------------
    # If a positive (float) value is given for 'relative', then compute
    #   relative energies multiplied by the value given. (So if the Gaussian
    #   energies are in hartree, 'relative=2625.5' yields kJ/mol.)
    # If a filename is specified for 'outfile', write a tab-delimited file.
    #   (Tab-delimited because variable descriptions contains commas.)
    # If a plot is requested, and there are fewer than three
    #   scanned variables, then draw a plot of energy.
    nfiles = len(filenames)
    print('Combining {:d} files'.format(nfiles))
    scandata = []  # list of DataFrames from the files
    for file in filenames:
        # loop through the files to extract scan information
        try:
            fhandl = open(file, 'r')
        except:
            print('*** Failure reading file {:s}'.format(file))
            # proceed to the next file in the list
            continue
        print('<<< {:s} >>>'.format(file))
        df_scan = extract_g09_scan(fhandl)
        fhandl.close()
        # remove the 'line' and 'byte' fields, add the filename
        del df_scan['line']
        del df_scan['byte']
        # put the filename just before the coordinates
        df_scan.insert(len(df_scan.columns)-1, 'FileName', file)
        # add this DataFrame to the list
        scandata.append(df_scan)
        if len(scandata) == 1:
            # save the variables that were scanned
            lbl = df_scan.iloc[0].Label
            nvar = len(lbl)
            descr = df_scan.iloc[0].Description
        else:
            # this is not the first DataFrame in the list
            # Check that the scanned variables are the same as in the first file
            lbl2 = df_scan.iloc[0].Label
            descr2 = df_scan.iloc[0].Description
            if descr2 != descr:
                print('*** change in scanned coordinate! ***')
                print(df_scan[['Label', 'Description']].head(1))
            if lbl2 != lbl:
                print('-- change in name of scanned variable --')
    # combine into a single DataFrame
    alldat = pd.concat(scandata)
    # sort by values of scanned variables
    alldat.sort_values(by=list(descr), inplace=True)
    # Are relative energies requested?
    if relative > 0:
        erel = alldat.Energy  # just the Energy column
        emin = erel.min()
        erel = (erel - emin) * relative
        # add column just before Energy
        icol = alldat.columns.get_loc('Energy')
        alldat.insert(icol, 'E_rel', erel)
    # 
    if outfile:
        # write a tab-separated file of variable descriptions,
        #   their values, and the corresponding energies
        # first make a copy stripped of unwanted columns
        try:
            dcopy = alldat.copy()
            del dcopy['Label']
            del dcopy['Description']
            del dcopy['Coordinates']
            dcopy.set_index(descr[0])
            dcopy.to_csv(outfile, sep='\t', index=False)
            print('>>> Tab-delimited file {:s} written <<<'.format(outfile))
        except:
            print('*** Unable to write tab-delimited file', outfile, '***')
    #
    if plot:
        if nvar > 1:
            print('*** Plotting is coded only for a single scanned variable ***')
        else:
            # create a scatter plot with lines connecting adjacent points
            if 'E_rel' in alldat.columns:
                # plot relative energy if available
                alldat.plot.line(x=descr[0], y='E_rel', style='-o')
            else:
                # plot raw energy
                alldat.plot.line(x=descr[0], y='Energy', style='-o')
            plt.show()
    return alldat
##
def opt_success(fhandl):
    # Return True if the string 'Optimization completed.' is found, else False
    byte_start = fhandl.tell()
    fhandl.seek(0)
    complete = False
    if 'Optimization completed.' in fhandl.read():
        complete = True
    fhandl.seek(byte_start) # restore file pointer to original position
    return complete
##
def read_g09_freqs(fhandl):
    # Return a list, with one row for each frequency block found.  Each
    # row contains three elements:
    #   (1) line number for "Harmonic frequencies" string,
    #   (2) byte number for same,
    #   (3) a DataFrame with seven columns:
    #       (3a) mode number,
    #       (3b) irrep,
    #       (3c) frequency in cm^-1,
    #       (3d) reduced mass,
    #       (3e) force constant,
    #       (3f) IR intensity,
    #       (3g) cartesian normal mode vector as a DataFrame:
    #           (3g1) atomic number (Z),
    #           (3g2-4) displacement coordinates x, y, z
    # 
    byte_start = fhandl.tell()
    fhandl.seek(0)
    freqblock = []
    lineno = 0
    regx = re.compile(r'Harmonic frequencies \(cm\*\*-1\)')
    regblank = re.compile(r'^\s+$')
    regnumbers = re.compile(r'^[\d\s]+$')
    regfreq = re.compile(r'^\s*Frequencies\s')
    regmass = re.compile(r'^\s*Red\. masses\s')
    regforce = re.compile(r'^\s*Frc consts\s')
    reginten = re.compile(r'^\s*IR Inten\s')
    regvec = re.compile(r'\s*Atom\s+AN\s+X\s+Y\s+Z')
    inblock = nmblock = False
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        if inblock:
            # inside a frequency block 
            mblank = regblank.match(line)
            mmodenum = regnumbers.match(line)
            if mblank or mmodenum:
                # close out any preceding vector block
                if len(nmvec):
                    df_nmblock = pd.DataFrame(nmvec, columns=nmhdr)
                    df_nmblock = df_nmblock.apply(pd.to_numeric, errors='ignore')
                    # separate the eigenvectors in this DataFrame
                    for i in range(len(df_nmblock.columns)//3):
                        selhdr = ['{:s}{:d}'.format(x, i+1) for x in ['X','Y','Z']]
                        selhdr.insert(0, 'AN')
                        df = df_nmblock.ix[:,selhdr]
                        # rename the columns
                        df.columns = list('Zxyz')
                        # add to the list
                        vector.append(df)
                    nmvec = []
                if mblank:
                    # blank line indicates end of the block
                    inblock = nmblock = False
                    # pull the data together
                    data = list(zip(mode, irrep, freq, redm, frc, inten, vector))
                    if False:
                        print('mode: ', mode)
                        print('irrep: ', irrep)
                        print('freq: ', freq)
                        print('redm: ', redm)
                        print('frc: ', frc)
                        print('inten: ', inten)
                        print('vector: ', vector)
                        print('data:')
                        print(data)
                    cols = ['Mode', 'Irrep', 'Freq', 'Red. mass', 'Frc const',
                        'IR inten', 'Vector']
                    df = pd.DataFrame(data=data, columns=cols)
                    # convert strings to numbers where possible
                    df = df.apply(pd.to_numeric, errors='ignore')
                    freqblock.append( [blockline, blockbyte, df] )
                else:
                    # normal mode numbers
                    mode.extend(line.split())
                    nmblock = False
                    nmvec = []
                    # also read the following line (irreducible representations, at best)
                    line = fhandl.readline()
                    lineno += 1
                    irrep.extend(line.split())
                continue
            else:
                # parse inside the block
                m = regfreq.match(line)
                if m:
                    # harmonic frequencies in cm^-1
                    freq.extend( line.split()[2:] )
                    continue
                m = regmass.match(line)
                if m:
                    # reduced masses
                    redm.extend( line.split()[3:] )
                    continue
                m = regforce.match(line)
                if m:
                    # force constants
                    frc.extend( line.split()[3:] )
                    continue
                m = reginten.match(line)
                if m:
                    # IR intensities
                    inten.extend( line.split()[3:] )
                    continue
                if nmblock:
                    # parse a line from a normal-mode vectors block
                    nmvec.append( line.split() )
                    continue
                m = regvec.match(line)
                if m:
                    # normal mode vector block
                    nmblock = True
                    nmhdr = line.split()  # column headings
                    # add numerals to headings to remove duplication
                    nmhdr[2:] = ['{:s}{:d}'.format(nmhdr[i], (i+1)//3) for i in range(2,len(nmhdr))]
                    continue
        m = regx.search(line)
        if m:
            # found a frequency block
            inblock = True
            blockline = lineno
            blockbyte = fhandl.tell()
            # lists below are buffers; each element is for one normal mode
            mode = []
            irrep = []
            freq = []
            redm = []
            frc = []
            inten = []
            vector = []
            # 'nmvec' is a buffer for a group of (three) mode vectors
            nmvec = []
    fhandl.seek(byte_start) # restore file pointer to original position
    return freqblock
##
def get_nimag(fhandl):
    # Read the (harmonic) vibrational frequency information
    # Return the number of imaginary (rendered as negative) frequencies
    # If more than one set of frequencies was found in the Gaussian file,
    #   return a list of nimag values.
    freqtbl = read_g09_freqs(fhandl)
    nimag = []
    nfreq = []
    for row in freqtbl:
        tbl = row[2]  # discard first two elements (line# and byte#)
        nfreq.append(len(tbl))
        negfreq = tbl[tbl.Freq < 0]
        nimag.append(len(negfreq))
    if min(nfreq) > 1:
        # there are at least some frequencies in each block
        return max(nimag)
    else:
        # failed to parse frequencies
        return -1
##
def read_g09_electrons(fhandl):
    # read all lines with electron count, return DataFrame:
    #   (1) line number(s), (2) byte number(s),
    #   (3) number of alphas, (4) number of betas,
    #   (5) total number of electrons
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    alpha = []
    beta = []
    total = []
    fline = []
    fpos = []
    lineno = 0
    regx = re.compile(r'(\d+) alpha electrons\s+(\d+) beta electrons')
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        m = regx.search(line)
        if m:
            # found electron counts
            a = int(m.group(1))
            b = int(m.group(2))
            alpha.append(a)
            beta.append(b)
            total.append(a+b)
            fline.append(lineno)
            fpos.append(fhandl.tell())
    data = list(zip(fline, fpos, alpha, beta, total))
    cols = ['line', 'byte', 'Alpha', 'Beta', 'Total']
    df = pd.DataFrame(data, columns=cols)
    fhandl.seek(byte_start) # restore file pointer to original position
    return df
##
def read_g09_oeke(fhandl, virtual=False):
    # read orbital energies and kinetic energies (for BEB application)
    # set virtual=True to include virtual orbitals
    # only read the first set of such data in the target file
    # Return a pandas DataFrame with a row for each orbital like:
    #   (1) line number, (2) byte number,
    #   (3) orbital number, (4) spin (alpha or beta), (5) orbital label,
    #   (6) energy,
    #   (7) kinetic energy
    #
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    orbnum = []
    lbl = []
    energy = []
    ke = []
    fline = []
    fpos = []
    spin = []
    lineno = 0
    regstart = re.compile('Orbital energies and kinetic energies')
    regend = re.compile('Total kinetic energy from orbitals')
    regocc = re.compile(r'O$')
    regbeta = re.compile(r'\(beta\)')
    inblock = False
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        if inblock:
            m = regend.search(line)
            if m:
                # end of block
                data = list(zip(fline, fpos, orbnum, spin, lbl, energy, ke))
                cols = ['line', 'byte', 'Orbital', 'Spin', 'Label', 'Energy', 'KE']
                table = pd.DataFrame(data, columns=cols)
                break   # only read the first such bock of data
            # parse this line
            fields = line.split()
            if len(fields) != 4:
                # not a data line
                m = regbeta.search(line)
                if m:
                    # beta orbitals are coming up
                    curspin = 'beta'
                continue
            # a good line--is it occupied or virtual?
            m = regocc.search(fields[1])
            if m or virtual:
                # report on this orbital
                orbnum.append(int(fields[0]))
                spin.append(curspin)
                lbl.append(fields[1])
                energy.append(float(fields[2]))
                ke.append(float(fields[3]))
                fline.append(lineno)
                fpos.append(fhandl.tell())
        m = regstart.search(line)
        if m:
            # beginning of data block (alpha orbitals)
            inblock = True
            curspin = 'alpha'
    fhandl.seek(byte_start) # restore file pointer to original position
    return table
##
def read_g09_ept(fhandl):
    # Read correlated orbital energies from electron propagator calculation.
    # Only read the first set of such data in the target file.
    # Energies in hartree.
    # NOTE: Gaussian sometimes reports PS = 1.000 when the calculation actually failed.
    #   so change PS = 1.000 to PS = 0.000 if there was a preceding warning.
    # Return a pandas DataFrame with a row for each orbital like:
    #   (1) line number, (2) byte number,
    #   (3) orbital number (starting with 1 for the lowest core orbital),
    #   (4) spin (alpha or beta), 
    #   (5) Koopmans value (uncorrelated),
    #   (6) OVGF second-order energy,
    #   (7) OVGF second-order pole strength,
    #   (8) OVGF third-order energy,
    #   (9) OVGF third-order pole strength,
    #  (10) OVGF outer-valence approx. energy,
    #  (11) OVGF outer-valence approx. pole strength,
    #  (12) P3 second-order energy,
    #  (13) P3 second-order pole strength,
    #  (14) P3 third-order energy,
    #  (15) P3 third-order pole strength.
    #
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    lineno = 0
    oline = []  # for OVGF
    opos = []   # for OVGF
    pline = []  # for P3
    ppos = []   # for P3
    ospin = []
    pspin = []
    oorb = []  # valence orbital numbers as printed by Gaussian
    porb = []  # valence orbital numbers as printed by Gaussian
    ovgf2 = []  # energy, hartree
    ovgf3 = []
    ovgf = []
    p32 = []
    p33 = []
    okoop = []
    pkoop = []
    ovgf2ps = []  # pole strength
    ovgf3ps = []
    ovgfps = []
    p32ps = []
    p33ps = []
    warnPS = 0  # counter/flag
    degenPair = []  # pairs of degenerate orbitals, one skipped by Gaussian
    #   one element of degenPair looks like ('alpha'|'beta', orbnum1, orbnum2)
    degenTo = []  # for each degen. pair, this is the orbital that was skipped
    regx = re.compile('^\s*Summary of results for\s+(alpha|beta)\s+spin-orbital\s+(\d+)\s+(OVGF|P3):')
    regncore = re.compile('^\s*NBasis=\s.*\sNFC=\s+(\d+)\s')
    regkoop = re.compile('^\s*Koopmans theorem:\s+(\S+)')
    reg2nd = re.compile('^\s*Converged second order pole:\s+(\S+).*eV\s+(\S+)')
    reg3rd = re.compile('^\s*Converged (?:third|3rd).*pole:\s+(\S+).*eV\s+(\S+)')
    regova = re.compile('^\s*Outer Valence Approximation:\s+(\S+).*eV\s+(\S+)')
    regdegen = re.compile('^\s*Orbitals\s+(\d+)\s+and\s+(\d+)\s+are degenerate.  Skipping orbital\s+(\d+)')
    regwarn = re.compile('WARNING')
    regblank = re.compile('^\s*$')
    regPSbad = re.compile('1\.000 \(PS\)')
    ncore = -1
    block = ''  # either 'OVGF' or 'P3'
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        if ncore < 0:
            m = regncore.match(line)
            if m:
                # get the number of frozen-core orbitals
                ncore = int(m.group(1))
        else:
            # already found the number of frozen cores; look for results
            m = regdegen.match(line)
            if m:
                # a degeneracy statement
                degenPair.append( ( ospin[-1], int(m.group(1))+ncore, int(m.group(2))+ncore) )
                degenTo.append(int(m.group(3)) + ncore)
            m = regwarn.search(line)
            if warnPS:
                # a warning has been issued recently by the EPT module
                if regblank.match(line):
                    # blank line, decrement warning monitor
                    warnPS -= 1
                else:
                    if regPSbad.search(line):
                        # pole strength of 1.000 is probably a mistake;
                        #   change it to zero
                        line = line.replace('1.000', '0.000')
            if m:
                # the EPT module issued a warning about the following orbital
                # decrement 'warnPS' with each blank line and ignore it when < 1
                warnPS = 2
            if block == '':
                # look for a data block
                m = regx.match(line)
                if m:
                    # this is a data block
                    block = m.group(3)
                    if block == 'OVGF':
                        oline.append(lineno)
                        opos.append(fhandl.tell())
                        ospin.append(m.group(1))
                        oorb.append(int(m.group(2)))
                    if block == 'P3':
                        pline.append(lineno)
                        ppos.append(fhandl.tell())
                        pspin.append(m.group(1))
                        porb.append(int(m.group(2)))
            else:
                # already in an energy data block
                if block == 'OVGF':
                    m = regkoop.match(line)
                    if m:
                        # Koopmans result (no correlation energy)
                        okoop.append(float(m.group(1).replace('D','E')))
                        continue
                    m = reg2nd.match(line)
                    if m:
                        # 2nd-order result
                        ovgf2.append(float(m.group(1).replace('D','E')))
                        ovgf2ps.append(float(m.group(2)))
                        continue
                    m = reg3rd.match(line)
                    if m:
                        # 3rd-order result
                        ovgf3.append(float(m.group(1).replace('D','E')))
                        ovgf3ps.append(float(m.group(2)))
                        continue
                    m = regova.match(line)
                    if m:
                        # outer-valence approximation
                        ovgf.append(float(m.group(1).replace('D','E')))
                        ovgfps.append(float(m.group(2)))
                        continue
                if block == 'P3':
                    m = regkoop.match(line)
                    if m:
                        # Koopmans result (no correlation energy)
                        pkoop.append(float(m.group(1).replace('D','E')))
                        continue
                    m = reg2nd.match(line)
                    if m:
                        # 2nd-order result
                        p32.append(float(m.group(1).replace('D','E')))
                        p32ps.append(float(m.group(2)))
                        continue
                    m = reg3rd.match(line)
                    if m:
                        # 3rd-order result
                        p33.append(float(m.group(1).replace('D','E')))
                        p33ps.append(float(m.group(2)))
                        continue
                # get here if we ran out of data in a data block
                block = ''
    if len(oline) > 0 and len(pline) < 1:
        # we have OVGF data but no P3 data
        data = list(zip(oline, opos, oorb, ospin, okoop, ovgf2, ovgf2ps,
            ovgf3, ovgf3ps, ovgf, ovgfps))
        cols = ['line', 'byte', 'Orbital', 'Spin', 'Koopmans', 'OVGF-2nd', 'OVGF2 PS', 
            'OVGF-3rd', 'OVGF3 PS', 'OVGF', 'OVGF PS']
        df_ovgf = pd.DataFrame(data=data, columns=cols)
        df_ovgf['Orbital'] += ncore # reset the numbering to include the core
        df = df_ovgf
    if len(oline) < 1 and len(pline) > 0:
        # we have P3 data but no OVGF data
        data = list(zip(pline, ppos, porb, pspin, pkoop, 
            p32, p32ps, p33, p33ps))
        cols = ['line', 'byte', 'Orbital', 'Spin', 'Koopmans', 
            'P3-2nd', 'P3-2 PS', 'P3-3rd', 'P3-3 PS']
        df_p3 = pd.DataFrame(data=data, columns=cols)
        df_p3['Orbital'] += ncore
        df = df_p3
    if len(oline) > 0 and len(pline) > 0:
        # we have both OVGF and P3 data
        if len(oline) != len(pline):
            # different numbers of orbitals--should not happen
            fhandl.seek(byte_start) # restore file pointer to original position
            return 'Inconsistency: {:d} orbitals from OVGF but {:d} from P3'.format(len(oline), len(pline))
        # merge OVGF and P3 into one dataframe; use the positions from OVGF
        data = list(zip(oline, opos, oorb, ospin, okoop, ovgf2, ovgf2ps,
            ovgf3, ovgf3ps, ovgf, ovgfps,
            p32, p32ps, p33, p33ps))
        cols = ['line', 'byte', 'Orbital', 'Spin', 'Koopmans', 'OVGF-2nd', 'OVGF2 PS', 
            'OVGF-3rd', 'OVGF3 PS', 'OVGF', 'OVGF PS',
            'P3-2nd', 'P3-2 PS', 'P3-3rd', 'P3-3 PS']
        df = pd.DataFrame(data=data, columns=cols)
        df['Orbital'] += ncore
    fhandl.seek(byte_start) # restore file pointer to original position
    # apply any degeneracy information
    for i in range(len(degenTo)):
        degenFrom = degenPair[i][1]
        if degenFrom == degenTo[i]:
            # choose the other one
            degenFrom = degenPair[i][2]
        # duplicate data for the degenerate orbital
        rowFrom = df[(df['Orbital'] == degenFrom) & (df['Spin'] == degenPair[i][0])]
        df = df.append(rowFrom, ignore_index=True)
        lastrow = len(df.index) - 1
        df.set_value(lastrow, 'Orbital', degenTo[i])
    return df
##
def read_best_ept(fhandl, minPS=0.80):
    # select the "best" values from EPT/OVGF calculations
    # a minimum pole strength of 'minPS' is required
    # Energies in hartree.
    # Return a pandas DataFrame with a row for each orbital like:
    #   (1) orbital number (starting with 1 for the lowest core orbital),
    #   (2) spin (alpha or beta), 
    #   (3) method (e.g., 'OVGF-2nd')
    #   (4) energy
    #   (5) pole strength
    #
    df = read_g09_ept(fhandl)
    # take Koopmans values as default
    df['Method'] = 'Koopmans'
    df['Energy'] = df['Koopmans']
    df['PS'] = 0.0
    # 'trust' is the assumed decreasing ordering of reliability 
    trust = ['P3-3', 'P3-2', 'OVGF', 'OVGF3', 'OVGF2']
    for meth in trust:
        col = '{:s} PS'.format(meth)     # label for columns with pole strengths
        idxps = df.columns.get_loc(col)  # column number for pole strength
        ecol = df.columns[idxps - 1]     # the corresonding energy
        choose = (df[col] > minPS) & (df['Method'] == 'Koopmans')
        df.ix[choose, 'Method'] = ecol
        df.ix[choose, 'Energy'] = df[ecol]
        df.ix[choose, 'PS'] = df[col]
    return df[['Orbital', 'Spin', 'Method', 'Energy', 'PS']]
##
