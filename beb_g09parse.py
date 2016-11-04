#!/usr/bin/python3
# Read G09 output files and create BUN data file for BEB calculation
# Expected files are those created by beb_g09build.py
# Karl Irikura, new version started Sept. 2016
#
import sys
import os
import re
#from beb_subs import *
sys.path.append('/media/sf_bin3')
from g09_subs3 import *
from chem_subs import *
import pandas as pd
###
def weakest_binding(df_orbitals):
    # given a DataFrame of orbital energies (etc.), find the most weakly bound
    # return values: (1) binding energy in eV (molecular IE),
    #   (2) orbital spin, (3) theoretical method
    #
    # keep only the row with the weakest binding energy
    IE = max(df_orbitals['Energy']) 
    df = df_orbitals.loc[df_orbitals['Energy'] == IE]
    # change sign and convert to eV
    IE = hartree_eV(-IE, 'to_eV')
    try:
        theory = df.iloc[0]['Method']
    except:
        # no theory; this is probably HF
        theory = 'Koopmans'
    # determine corresponding ground state of dication
    orbspin = df.iloc[0]['Spin']
    return IE, orbspin, theory
##
def ion_multiplicity(nalp_neut, nbet_neut, orb_spin):
    # given N_alpha and N_beta for neutral molecule, and spin of orbital ('alpha'
    #   or 'beta') being ionized, return: 
    #   (1) the spin multiplicity of the resulting ion,
    #   (2) the number of alpha and (3) beta electrons in the ion
    nalp_ion = nalp_neut
    nbet_ion = nbet_neut
    if orb_spin == 'alpha':
        nalp_ion -= 1
    elif orb_spin == 'beta':
        nbet_ion -= 1
    else:
        print('*** Unrecognized orbital spin for ionization: {:s}'.format(orb_spin))
    # require that N_alpha >= N_beta 
    if (nalp_ion < nbet_ion):
        # swap the values
        (nalp_ion, nbet_ion) = (nbet_ion, nalp_ion)
    ionmult = nalp_ion - nbet_ion + 1
    return ionmult, nalp_ion, nbet_ion
##
def make_label(n, a):
    abbrev = {'alpha': 'a', 'beta': 'b'}
    return('{:d}{:s}'.format(n, abbrev[a]))
##
def merge_degenerate(df_orbs):
    # orbitals are considered degenerate if their binding and kinetic energies are identical
    # combine data for degenerate orbitals in DataFrame, merely by increasing
    #   the corresponding electron count
    dupl = df_orbs.duplicated(['Energy', 'KE'])
    uniq = df_orbs.loc[~dupl]
    degen = df_orbs.loc[dupl]
    for i in uniq.index:
        # increase electron count for each matching orbital in 'ndegen' 
        urow = uniq.loc[i].copy()
        for j in degen.index:
            drow = degen.loc[j]
            if (urow.Energy == drow.Energy) and (urow.KE == drow.KE):
                # this is a match
                urow.N += drow.N
        # replace the electron count in 'uniq'
        uniq.set_value(i, 'N', urow.N)
    return uniq
##
#
# Read molecule name from command line
if len(sys.argv) < 2:
    sys.exit('This script is for parsing Gaussian jobs created by script "beb_g09build.py"' +
        '\nUsage:  beb_g09parse.py <file rootname>\n\tExample:  beb_g09parse.py c2h5oh')
mol = sys.argv[1]
# Results will be collected in pandas DataFrame
df_bun = pd.DataFrame()
# initialize CCSD(T) energies to zero (calculations may not be available)
cc_neut = cc_ion_hi = cc_ion_lo = cc_dicat_hi = cc_dicat_lo = 0
ecp = False  # flag
# parse each output job file in turn
for suff in 'opt bu bupp ept1 ept2 cc cc1hi cc1lo cc2hi cc2lo'.split():
    fname = '{:s}_{:s}.out'.format(mol, suff)
    try:
        fgau = open(fname, 'r')
    except:
        # This Gaussian output file is missing
        if suff == 'bu':
            sys.exit('The file {:s} is missing but is required.'.format(fname))
        else:
            print('\n*** Optional file not found: {:s} ***'.format(fname))
            continue
    print('\n>>> Reading file {:s} <<<'.format(fname))
    if suff == 'opt':
        # check that geometry optimization converged
        if opt_success(fgau):
            print('Geometry optimization converged: ', end='')
        else:
            sys.exit('\tGeometry optimization FAILED\n')
        # check that there are no imaginary frequencies
        nimag = get_nimag(fgau)
        if nimag > 0:
            print('NOT a minimum ({:d} imaginary vibrational frequencies)'.format(nimag))
            sys.exit()
        elif nimag < 0:
            print('No vibrational frequencies found')
        else:
            print('local minimum with nimag = 0')
    if suff == 'bu':
        # extract stoichiometry, charge, and spin multiplicity
        df = read_g09_stoichiometry(fgau)
        stoich = df.iloc[0].Stoich
        elems = df.iloc[0].Elements
        df = read_g09_charge_mult(fgau)
        charge = df.iloc[0].Charge
        mult = df.iloc[0].Mult
        ncore = n_core(elems) // 2  # number of conventional core orbitals
        nelec = read_g09_electrons(fgau).iloc[0]
        nalp = nelec.Alpha
        nbet = nelec.Beta
        # extract HF orbital energies and kinetic energies for occupied orbitals
        # orbital energies are negative binding energies
        buae = read_g09_oeke(fgau)[['Orbital', 'Spin', 'Energy', 'KE']]
        print('HF orbital and kinetic energies (hartree):')
        print(buae.to_string(index=False))
        # In case there are no ECP calculations, cope with n>2 possibility
        zvals = elz(list(elems.keys()))
        heavies = max(zvals) > 10
        if heavies:
            # Deal with this!
            buae['Special'] = 'deal'
        else:
            # just create the 'Special' column
            buae['Special'] = 'none'
        uhf = ('beta' in buae.Spin.tolist() )   # flag to indicate separate beta orbitals
        # extract crude molecular ionization threhold
        VIE, iespin, iemeth = weakest_binding(buae)
        mult_ion, nalp_ion, nbet_ion = ion_multiplicity(nalp, nbet, iespin)
        print('Crude VIE = {:.2f} eV to {:s} cation.'.format(VIE, spinname(mult_ion)))
        VIE_method = 'Koopmans'
        VIE2 = 40.0  # guess for vertical double ionization energy (eV)
        VIE2_method = 'guessing'
        mult_dicat = 0  # unknown
    if suff == 'bupp':
        # get valence B/U from HF/ECP calculation
        bupp = read_g09_oeke(fgau)[['Orbital', 'Spin', 'Energy', 'KE']]
        print('HF/ECP orbital (-B) and kinetic (U) energies (hartree):')
        bupp['Method'] = 'ECP'
        # adjust orbital numbers to account for core
        ppcore = len(buae[buae.Spin == 'alpha'].index) - len(bupp[bupp.Spin == 'alpha'].index) 
        print('ECP (pseudopotential) replaced {:d} electron pairs'.format(ppcore))
        bupp['Orbital'] += ppcore
        print(bupp.to_string(index=False))
        # note any light-atom core orbitals remaining
        ppskip = ncore - ppcore
        if ppskip > 0:
            print('({:d} core orbitals remain)'.format(ppskip))
        # extract VIE but don't replace the all-electron value from the previous file
        VIE_pp, ppspin, ppmeth = weakest_binding(bupp)
        mult_pp, nalp_pp, nbet_pp = ion_multiplicity(nalp-ppcore, nbet-ppcore, ppspin)
        print('Crude VIE = {:.2f} eV to {:s} cation.'.format(VIE_pp, spinname(mult_pp)))
        ecp = True
    if suff == 'ept1':
        # extract correlated orbital energies (negative of binding energies)
        ept = read_best_ept(fgau, minPS=0.75)
        # degeneracies will cause a funny orbital ordering--sort before printing
        ept = ept.sort_values(by=['Spin', 'Orbital'])
        # also put columns in the same order as the HF results just displayed
        ept = ept[['Orbital', 'Spin', 'Energy', 'PS', 'Method']]
        print(ept.to_string(index=False))
        # infer cation ground state
        VIE, iespin, VIE_method = weakest_binding(ept)
        mult_ion, nalp_ion, nbet_ion = ion_multiplicity(nalp, nbet, iespin)
        print('VIE = {:.2f} eV to the {:s} cation from {:s} theory.'.format(VIE,
            spinname(mult_ion), VIE_method))
    if suff == 'ept2':
        # get electron counts for cation
        ecation = read_g09_electrons(fgau).iloc[0]
        nalpion = ecation.Alpha
        nbetion = ecation.Beta
        # find ground state of (vertical) dication and (VIE2 - VIE) energy
        ept2 = read_best_ept(fgau, minPS=0.75)
        print(ept2.to_string(index=False))
        # extract molecular ionization energy
        IEcat, ie2spin, ie2meth = weakest_binding(ept2)
        mult_dicat, nalpdbl, nbetdbl = ion_multiplicity(nalpion, nbetion, ie2spin)
        print('Second IE = {:.2f} eV to the {:s} dication from {:s} theory.'.format(IEcat,
            spinname(mult_dicat), ie2meth))
        # infer double-ionization threshold, VIE2
        # was this calculation referenced to the ground state of the ion?
        if mult_ion == nalpion - nbetion + 1:
            # yes
            VIE2 = VIE + IEcat
            VIE2_method = '({:s} + {:s})'.format(VIE_method, ie2meth)
        else:
            # no, it was referenced to an excited state; it should only be too high by S=1
            if mult_ion == nalpion - nbetion - 1:
                # was high-spin; choose the highest beta orbital energy from ept1
                VIE_exc, iespin, exc_method = weakest_binding(ept.loc[ept['Spin'] == 'beta'])
                VIE2 = VIE_exc + IEcat
                VIE2_method = '({:s} + {:s})'.format(exc_method, ie2meth)
            else:
                print('--- I am confused about the EPT spin states to combine for VIE2 ---')
        print('Estimated VIE2 = {:.2f} eV from {:s} theory.'.format(VIE2, VIE2_method))
    if re.match('cc\S*', suff):
        # one of the CCSD(T) calculations
        # use only the first CCSD(T) energy in the file
        df = read_g09_postHF(fgau)  # this includes HF, MP2, etc.
        df = df.loc[df['Method'] == 'CCSD(T)']
        df_mult = read_g09_charge_mult(fgau)
        m = df_mult.iloc[0]['Mult']
        q = df_mult.iloc[0]['Charge']
        # check the charges
        mregx = re.match('cc(\d)\S+', suff)
        if mregx:
            # ion or dication
            qx = int(mregx.group(1))
            if q != qx:
                print('---Warning: charge = {:d} but expected {:d}---'.format(q, qx))
        if suff == 'cc':
            # neutral molecule; check charge and multiplicity
            if q != 0:
                print('---Warning: charge = {:d} for the neutral molecule---'.format(q))
            if m != mult:
                print('---Warning: mult = {:d} but expected {:d}---'.format(m, mult))
            cc_neut = df.iloc[0]['Energy']
            print('CCSD(T) = {:.6f} for neutral molecule'.format(cc_neut))
        elif suff == 'cc1hi':
            # cation, high-spin
            cc_ion_hi = df.iloc[0]['Energy']
            mult_ion_hi = m
            print('CCSD(T) = {:.6f} for {:s} ion'.format(cc_ion_hi, spinname(mult_ion_hi)))
        elif suff == 'cc1lo':
            # cation, low-spin
            cc_ion_lo = df.iloc[0]['Energy']
            mult_ion_lo = m
            print('CCSD(T) = {:.6f} for {:s} ion'.format(cc_ion_lo, spinname(mult_ion_lo)))
        elif suff == 'cc2hi':
            # dication, high-spin
            cc_dicat_hi = df.iloc[0]['Energy']
            mult_dicat_hi = m
            print('CCSD(T) = {:.6f} for {:s} dication'.format(cc_dicat_hi, spinname(mult_dicat_hi)))
        elif suff == 'cc2lo':
            # dication, low-spin
            cc_dicat_lo = df.iloc[0]['Energy']
            mult_dicat_lo = m
            print('CCSD(T) = {:.6f} for {:s} dication'.format(cc_dicat_lo, spinname(mult_dicat_lo)))
        else:
            # unknown and unexpected situation
            print('Unrecognized CC file suffix: ', suff)
# done reading Gaussian output files
fgau.close()
# process available CCSD(T) energies
cc_ion = cc_dicat = 0   # for the lower-energy states
if cc_neut != 0:
    # neutral energy is available; try to compute VIE and VIE2
    if (cc_ion_hi + cc_ion_lo) != 0:
        # there is at least one ion energy available; choose the lower
        if cc_ion_hi < cc_ion_lo:
            cc_ion = cc_ion_hi
            mult_ion = mult_ion_hi
        else:
            cc_ion = cc_ion_lo
            mult_ion = mult_ion_lo
        # calculate and use this VIE 
        VIE = hartree_eV(cc_ion - cc_neut)
        print('\nVIE = {:.2f} eV to {:s} cation from CCSD(T)'.format(VIE, spinname(mult_ion)))
        VIE_method = 'CCSD(T)'
    if (cc_dicat_hi + cc_dicat_lo) != 0:
        # there is at least one dication energy available; choose the lower
        if cc_dicat_hi < cc_dicat_lo:
            cc_dicat = cc_dicat_hi
            mult_dicat = mult_dicat_hi
        else:
            cc_dicat = cc_dicat_lo
            mult_dicat = mult_dicat_lo
        VIE2 = hartree_eV(cc_dicat - cc_neut)
        print('VIE2 = {:.2f} eV to {:s} dication from CCSD(T)'.format(VIE2, spinname(mult_dicat)))
        VIE2_method = 'CCSD(T)'
#
# done parsing Gaussian output files
# assemble the different DataFrames; start with the all-electron HF results
# Make an index, 'Label', that works for UHF and RHF
BUtable = buae.copy()
BUtable['Method'] = 'Koopmans'
# If this is UHF-based, replace orbital numbers with labels like '1a', '1b'
if uhf:
    # re-label orbitals 
    BUtable['Label'] = BUtable.apply(lambda x: make_label(x['Orbital'], x['Spin']), axis=1)
else:
    # RHF case; simple numeric label
    BUtable['Label'] = BUtable['Orbital']
# make 'Label' the index so merging the DataFrames is easier
BUtable.set_index(['Label'], inplace=True)
# Use kinetic energies from ECP calculation, if available
if ecp:
    if uhf:
        bupp['Label'] = bupp.apply(lambda x: make_label(x['Orbital'], x['Spin']), axis=1)
    else:
        bupp['Label'] = bupp['Orbital']
    bupp.set_index(['Label'], inplace=True)
    # replace the kinetic energies
    BUtable.update(bupp[['KE', 'Method']])
    # remove n>2 special cases
    BUtable['Special'] = 'none'
# If EPT calculations were done, use those values where available
try:
    if uhf:
        ept['Label'] = ept.apply(lambda x: make_label(x['Orbital'], x['Spin']), axis=1)
    else:
        ept['Label'] = ept['Orbital']
    ept.set_index(['Label'], inplace=True)
    # Replace with EPT values as available
    BUtable.update(ept['Energy'])
    # If Method == ECP, don't erase that info
    BUtable.loc[BUtable['Method'] != 'ECP'].update(ept['Method'])
    ecprows = BUtable['Method'] == 'ECP'
    BUtable.loc[ecprows, 'Method'] = 'B ' + ept['Method'] + '; U ECP'
except:
    # no EPT calculations are available
    pass
#
# prepare BUtable for output
#
# put Label on the same footing as the other columns
BUtable.reset_index(inplace=True)
# add columns: 'N', 'Q', 'DblIon'
BUtable['N'] = 1 if uhf else 2
BUtable['Q'] = 1
BUtable['DblIon'] = 'No'
# make energies positive and convert to eV
BUtable['Energy'] = BUtable['Energy'].apply(hartree_eV, args=('to_eV', -1))
BUtable['KE'] = BUtable['KE'].apply(hartree_eV)
# mark 'DblIon' column where B > VIE2
BUtable.ix[BUtable['Energy'] > VIE2, 'DblIon'] = 'Yes'
# Consolidate degenerate orbitals
BUtable = merge_degenerate(BUtable)
if VIE_method == 'CCSD(T)':
    # replace binding energy for HOMO
    if uhf and (mult_ion > mult):
        # apply VIE to the highest beta orbital
        homo_spin = 'beta'
    else:
        # apply VIE to the highest alpha orbital
        homo_spin = 'alpha'
    homo = BUtable.loc[BUtable['Spin'] == homo_spin]['Orbital'].idxmax()
    BUtable.loc[homo, 'Energy'] = VIE
    BUtable.loc[homo, 'Method'] = VIE_method
# sort by orbital number (with alphas and betas interleaved)
BUtable = BUtable.sort_values(by='Orbital')
# rename and rearrange columns as needed for input to 'beb_tbl.pl'
BUtable = BUtable.rename(index=str, columns={'Label': '#MO', 'KE': 'U/eV', 'Energy': 'B/eV'})
BUtable = BUtable.rename(index=str, columns={'Method': 'Remarks'})
BUtable = BUtable[['#MO', 'B/eV', 'U/eV', 'N', 'Q', 'DblIon', 'Special', 'Remarks']]
# last-minute check of electron count
if BUtable['N'].sum() != (nalp + nbet):
    # We've lost or gained electrons!
    print('*** Error: there are {:d} electrons but should be {:d}! ***'.format(BUtable['N'].sum(), nalp+nbet))
print()
# print tab-delimited and with energies rounded to .01 eV 
s = BUtable.to_csv(sep='\t', float_format='%.2f', index=False)
print(s)
fbun = '{:s}.bun'.format(mol)
fh = open(fbun, 'w')
fh.write('# Results from automated beb_g09build.py and beb_g09parse.py\n')
fh.write("# Molecule is '{:s}' ({:s} {:s})\n".format(mol, spinname(mult), stoich))
fh.write('# Double-ionization threshold = {:.2f} eV from {:s} (to {:s} dication)\n'.format(VIE2, VIE2_method, spinname(mult_dicat)))
fh.write(s)
fh.close()
print('BUN data file {:s} written and closed.'.format(fbun))
