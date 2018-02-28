# Routines for use in BEB calculations
# --needed by scripts 'beb_g09build.py' and 'beb_g09parse.py'
# Karl Irikura, NIST
#   (file created by beb_subs_extract.pl)
#
import sys
import re
import numpy as np
import pandas as pd
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
    if min(nfreq) > 0:
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
    #   (3) orbital number,
    #   (4) spin ('alpha' or 'beta' or 'both'),
    #   (5) orbital label,
    #   (6) energy,
    #   (7) kinetic energy
    # To get appropriate output from Gaussian, use IOP(6/81=3)
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
    # If there are only alph orbitals, presume this is an RHF calculation
    #   and call the spin 'both'
    if not 'beta' in set(table['Spin']):
        # no 'beta'
        table['Spin'] = 'both'
    return table
##
def read_g09_ept(fhandl):
    # Read correlated orbital energies from electron propagator calculation.
    #   Gaussian keyword: EPT(OVGF+P3)
    # Only read the first set of such data in the target file.
    # Energies in hartree.
    # NOTE: Gaussian sometimes reports PS = 1.000 when the calculation actually failed.
    #   Change PS to 0.000 if there was a preceding warning.
    # Return a pandas DataFrame with a row for each orbital like:
    #   (1) line number, (2) byte number,
    #   (3) orbital number (starting with 1 for the lowest core orbital),
    #   (4) spin (alpha or beta or 'both'), 
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
    regPS = re.compile(r'\d\.\d{3} \(PS\)')
    ncore = -1
    nppe = read_ECP_electrons(fhandl)  # get numbers of ECP-replaced electrons on each center
    nppe = sum(nppe.values())  # the total number of electrons replaced
    pporb = nppe // 2  # the number of core orbitals replaced
    block = ''  # either 'OVGF' or 'P3'
    inEPT = False  # flag to indicate that this is EPT/OVGF, if in a multi-step output file
    regEPT = re.compile(r'\s?#.*(EPT|OVGF)', re.IGNORECASE)
    regcpu = re.compile(r'Job cpu time:')  # signals the end of a computational step
    while True:
        line = fhandl.readline()
        if not line:
            break
        lineno += 1
        if not inEPT:
            # we are not yet reading the results of EPT/OVGF
            m = regEPT.match(line)
            if m:
                # yes, we found the relevant calculation
                inEPT = True
            continue
        # we reach here only if within an EPT/OVGF calculation
        if regcpu.search(line):
            # this is the end of the calculation
            break
        if ncore < 0:
            m = regncore.match(line)
            if m:
                # get the number of frozen-core explicit orbitals
                ncore = int(m.group(1))
                # add it to the number of ECP-replaced orbitals
                pporb += ncore  # this is the total number of missing core orbitals
        else:
            # already found the number of frozen cores; look for results
            m = regdegen.match(line)
            if m:
                # a degeneracy statement
                degenPair.append( ( ospin[-1], int(m.group(1))+pporb, int(m.group(2))+pporb) )
                degenTo.append(int(m.group(3)) + pporb)
            if warnPS > 0:
                # a warning has been issued recently by the EPT module
                if regblank.match(line):
                    # blank line, decrement warning monitor
                    warnPS -= 1
                else:
                    # zero out any pole strength to avoid using the following values
                    if regPS.search(line):
                        line = regPS.sub('0.000 (PS)', line)
            if regwarn.search(line):
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
                        # check for missing P3 results; fill with zeros
                        if len(pkoop) < len(pline):
                            pkoop.append(0.)
                        if len(p32) < len(pline):
                            p32.append(0.)
                        if len(p32ps) < len(pline):
                            p32ps.append(0.)
                        if len(p33) < len(pline):
                            p33.append(0.)
                        if len(p33ps) < len(pline):
                            p33ps.append(0.)
                    if block == 'P3':
                        pline.append(lineno)
                        ppos.append(fhandl.tell())
                        pspin.append(m.group(1))
                        porb.append(int(m.group(2)))
                        # check for missing OVGF results; fill with zeros
                        if len(okoop) < len(oline):
                            okoop.append(0.)
                        if len(ovgf2) < len(oline):
                            ovgf2.append(0.)
                        if len(ovgf2ps) < len(oline):
                            ovgf2ps.append(0.)
                        if len(ovgf3) < len(oline):
                            ovgf3.append(0.)
                        if len(ovgf3ps) < len(oline):
                            ovgf3ps.append(0.)
                        if len(ovgf) < len(oline):
                            ovgf.append(0.)
                        if len(ovgfps) < len(oline):
                            ovgfps.append(0.)
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
    if False:
        print('list lengths:')
        print('\toline: ', len(oline))
        print('\topos: ', len(opos))
        print('\toorb: ', len(oorb))
        print('\tospin: ', len(ospin))
        print('\tokoop: ', len(okoop))
        print('\tovgf2: ', len(ovgf2))
        print('\tovgf2ps: ', len(ovgf2ps))
        print('\tovgf3: ', len(ovgf3))
        print('\tovgf3ps: ', len(ovgf3ps))
        print('\tovgf: ', len(ovgf))
        print('\tovgfps: ', len(ovgfps))
        print('\tpline: ', len(pline))
        print('\tppos: ', len(ppos))
        print('\tporb: ', len(porb))
        print('\tpspin: ', len(pspin))
        print('\tpkoop: ', len(pkoop))
        print('\tp32: ', len(p32))
        print('\tp32ps: ', len(p32ps))
        print('\tp33: ', len(p33))
        print('\tp33ps: ', len(p33ps))
        print('\tovgf:', ovgf)
    if len(oline) > 0 and len(pline) < 1:
        # we have OVGF data but no P3 data
        data = list(zip(oline, opos, oorb, ospin, okoop, ovgf2, ovgf2ps,
            ovgf3, ovgf3ps, ovgf, ovgfps))
        cols = ['line', 'byte', 'Orbital', 'Spin', 'Koopmans', 'OVGF-2nd', 'OVGF2 PS', 
            'OVGF-3rd', 'OVGF3 PS', 'OVGF', 'OVGF PS']
        df_ovgf = pd.DataFrame(data=data, columns=cols)
        df_ovgf['Orbital'] += pporb # reset the numbering to include the core
        df = df_ovgf
    if len(oline) < 1 and len(pline) > 0:
        # we have P3 data but no OVGF data
        data = list(zip(pline, ppos, porb, pspin, pkoop, 
            p32, p32ps, p33, p33ps))
        cols = ['line', 'byte', 'Orbital', 'Spin', 'Koopmans', 
            'P3-2nd', 'P3-2 PS', 'P3-3rd', 'P3-3 PS']
        df_p3 = pd.DataFrame(data=data, columns=cols)
        df_p3['Orbital'] += pporb
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
        df['Orbital'] += pporb
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
        df.iat[lastrow, df.columns.get_loc('Orbital')] = degenTo[i]
    if not (df.Spin == 'beta').any():
        # this is an RHF case; change spin labels to 'both'
        df.Spin = 'both'
    return df
##
def read_best_ept(fhandl, minPS=0.80):
    # select the "best" values from EPT/OVGF calculations
    # a minimum pole strength of 'minPS' is required
    # Energies in hartree.
    # Return a pandas DataFrame with a row for each orbital like:
    #   (1) orbital number (starting with 1 for the lowest core orbital),
    #   (2) spin (alpha or beta or 'both'), 
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
        try:
            idxps = df.columns.get_loc(col)  # column number for pole strength
        except KeyError:
            print('--No data for EPT method {:s}'.format(meth))
            continue
        ecol = df.columns[idxps - 1]     # the corresonding energy
        choose = (df[col] > minPS) & (df['Method'] == 'Koopmans')
        df.ix[choose, 'Method'] = ecol
        df.ix[choose, 'Energy'] = df[ecol]
        df.ix[choose, 'PS'] = df[col]
    return df[['Orbital', 'Spin', 'Method', 'Energy', 'PS']]
##
def n_sdd(atno, all=False):
    # given Z value (or element symbol) return number of electrons replaced by SDD
    #   pseudopotential (SDDall if all==True)
    if type(atno) != int:
        # convert symbol to Z value
        atno = elz(atno)
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
def read_AOpop_in_MOs(fhandl):
    # Read the output from the Gaussian09 keyword:
    #   pop=AllOrbitals
    # Return a pandas DataFrame with the following columns:
    #   (1) orbital number (starting from 1)
    #   (2) 'alpha' or 'beta' or 'both'
    #   (3) 'occ' or 'virt'
    #   (4) orbital energy (hartree)
    #   (5) element symbol of atom
    #   (6) atom number of atom (starting from 1)
    #   (7) AO orbital type ('s', 'p', 'd', etc.)
    #   (8) AO contribution
    # So each MO may be spread across multiple rows of the table. 
    # Long lines may contain embedded line breaks. 
    #
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    regstart = re.compile('Atomic contributions to (Alpha|Beta ) molecular orbitals:')
    regline1 = re.compile(' (occ|vir) ')  # beginning of data line(s)
    regblank = re.compile('^\s*$')
    regAO = re.compile('([A-Z][a-z]?)(\d+)-([spdfghijkl])')
    reginicap = re.compile('^\s*[A-Z]')
    MOnum = []
    spin = []
    occ = []
    OE = []
    symb = []
    atnum = []
    L = []
    contrib = []
    inblock = False
    #endData = False
    prevline = ''
    while True:
        line = fhandl.readline()
        if not line:
            # reached the end of the input file
            break
        if inblock:
            if regblank.match(line):
                # blank line indicates the end of a data block
                inblock = False
                #endData = True
            else:
                # in the data block and not a blank line
                if regline1.search(line):
                    # this is the first of (possibly continuing) data lines
                    # any previous line is complete and should be parsed
                    pass
                else:
                    # this is a continuation line; join it with the previous line
                    # if the continuation line begins with a capital letter, put
                    #   a space in between
                    if reginicap.match(line):
                        prevline = prevline.strip() + ' ' + line.strip()
                    else:
                        prevline = prevline.strip() + line.strip()
                    # don't parse anything yet, just read the next line of input
                    continue
            if len(prevline) == 0:
                # there is no previous line to parse; read next line of input
                prevline = line
                continue
            # parse the (previous) line
            fields = prevline.split()
            # fifth field is word 'is'; remaining fields are contribution strings
            for s in fields[5:]:
                sf = s.split('=')
                # now sf[0] is the AO string and sf[1] is the contribution
                m = regAO.match(sf[0])
                if m:
                    # first field is spin label
                    spin.append(fields[0].lower())
                    # second field is occupation status
                    if fields[1] == 'occ':
                        occ.append('occ')
                    else:
                        occ.append('virt')
                    # third field is MO number
                    MOnum.append(int(fields[2]))
                    # fourth field is OE string
                    s = fields[3].split('=')
                    OE.append(float(s[1]))
                    symb.append(m.group(1))
                    atnum.append(int(m.group(2)))
                    L.append(m.group(3))
                    contrib.append(float(sf[1]))
                else:
                    # something is wrong
                    print('*** Failure parsing AO label:', sf[0])
                    print(line)
            # done parsing the previous, complete line; replace with current line
            prevline = line
        else:
            # not in a data block
            #if endData:
            #    break
            if regstart.search(line):
                # this is the beginning of the data block 
                inblock = True
                continue
    data = list(zip(MOnum, spin, occ, OE, symb, atnum, L, contrib))
    cols = ['Orbital', 'Spin', 'Occ', 'Energy', 'Elem', 'Atom#', 'L', 'Contrib']
    df = pd.DataFrame(data, columns=cols)
    if not (df.Spin == 'beta').any():
        # RHF case; change spin-labels to 'both'
        df.Spin = 'both'
    fhandl.seek(byte_start) # restore file pointer to original position
    return df
##
def read_ECP_electrons(fhandl):
    # Read (one time) the number of electrons, on each atom, that have been replaced
    #   by an effective core potential (ECP; aka pseudopotential, "PP")
    # Return a dict where the key is the atom ordinal number (starting from 1)
    #   and the value is the number of electrons replaced on that atom
    #
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    regstart = re.compile(r'^\s*Pseudopotential Parameters\s*$')
    regdash = re.compile(r'==========================================')
    regx = re.compile(r'^\s{,6}\d+\s+\d+\b')
    regnatom = re.compile(r'\s*NAtoms=\s+(\d+)\s+')
    inblock = False
    centerno = []
    atz = []
    nval = []
    natom = 0
    while True:
        line = fhandl.readline()
        if not line:
            break
        # Look for the number of atoms.  Only needed if no ECP is found.
        m = regnatom.match(line)
        if m:
            natom = int(m.group(1))
        if inblock:
            if regdash.search(line):
                ndash -= 1  # decrement counter
            if ndash < 1:
                # end of ECP data block; finish up
                inblock = False
                # Read only the first ECP data block in the file
                break
            # still within the data block--is it a line of interest?
            m = regx.match(line)
            if m:
                # yes, parse this line
                fields = line.split()
                centerno.append(int(fields[0]))
                atz.append(int(fields[1]))
                if len(fields) == 3:
                    # there is a PP on this center
                    nval.append(int(fields[2]))
                else:
                    # no pseudopotential here; set nval = Z
                    nval.append(int(fields[1]))
        else:
            # check for beginning of data block
            m = regstart.match(line)
            if m:
                inblock = True
                ndash = 3  # counter; down to 0 at end of data block
    # subtract to obtain the number of electrons replaced
    nppe = {}
    for i in range(len(centerno)):
        nppe[centerno[i]] = atz[i] - nval[i]
    if len(nppe) == 0:
        # no ECP was found; just fill with zeros
        for i in range(natom):
            nppe[i+1] = 0
    fhandl.seek(byte_start) # restore file pointer to original position
    return nppe
##
def count_g09_success(fhandl):
    # Return the number of occurrences of "Normal Termination "
    #
    byte_start = fhandl.tell()
    fhandl.seek(0)  # rewind file
    regx = re.compile('Normal termination ')
    nsuccess = 0
    while True:
        line = fhandl.readline()
        if not line:
            break
        if regx.search(line):
            nsuccess += 1
    fhandl.seek(byte_start) # restore file pointer to original position
    return nsuccess
##
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
    try:
        Z = int(ar)
    except:
        print('Error taking int of ar =', ar, type(ar))
        return None
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
def combine_MOspin(df, col1='Orbital', col2='Spin', colnew='MO'):
    # Given a pandas DataFrame, combine a numeric 'Orbital' field with
    #   a 'Spin' field ('alpha' or 'beta') to create a new 'MO' field
    #   that is a combination like '1a' or '5b'.
    # Return that new DataFrame.
    abbrev = {'alpha': 'a', 'beta': 'b', 'both': ''}
    dfret = df.copy()
    dfret[colnew] = df.apply(lambda x: str(x[col1])+abbrev[x[col2]], axis=1)
    return dfret
##
def JSdm(P, Q, base=4):
    # Jensen-Shannon divergence metric; base=4 gives range = [0, 1]
    # P and Q are *discrete* PDFs (with same data type)
    # Allowed data types: tuple; list; dict; 1D numpy array
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
def AOpopdiffmats(df1, df2):
    # Compare two pandas DataFrames with Mulliken population data,
    #   as returned by routine 'read_AOpop_in_MOs()' in 'g09_subs.py'
    # Return two numpy 2D-arrays:
    #   (1) JSdm() differences in AO populations (Jensen-Shannon divergence metric)
    #   (2) (E2-E1) orbital energy differences
    # Also return two lists of MO numbers:
    #   (3) MO labels in df1 (rows of matrices)
    #   (4) MO labels in df2 (columns of matrics)
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
                # treat negative AO pop as a new variable (by changing its label)
                s += '-neg'
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
                    # negative AO pop
                    s += '-neg'
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
    #   Keys are MO labels in df2, values are MO labels in df1.
    # Do not mix alpha with beta orbitals.
    #
    momap = {}
    if (df1['Spin']  == 'alpha').any() & (df1['Spin'] == 'beta').any():
        # this is a UHF case; keep alpha and beta orbitals separate
        for sp in ['alpha', 'beta']:
            set1 = df1[df1['Spin'] == sp]
            set2 = df2[df2['Spin'] == sp]
            momap.update(orbitalPopMatch(set1, set2, Eweight=Eweight, diagBias=diagBias))
        return momap 
    # simple, single-spin case
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
    for i in pairing.keys():
        imo = MOs1[i]  # MO number from first set
        j = pairing[i]
        jmo = MOs2[j]  # MO number from second set
        if imo != jmo:
            # report only non-identity mappings
            momap[jmo] = imo  # key is the MO number in the 2nd set
    return momap 
##
def relabelOrbitals(df, momap):
    # re-label MOs based upon a mapping provided by 'orbitalPopMatch()'
    # Return value: the DataFrame with orbitals re-labeled
    #
    # loop once through the rows, changing MO labels
    for idx in df.index:
        imo = df.loc[idx, 'MO'] 
        if imo in momap.keys():
            # change this MO label
            df.loc[idx, 'MO'] = momap[imo]
    return df
##
