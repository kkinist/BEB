#!/usr/bin/python3
# Read G09 output files and create BUN data file for BEB calculation
# Expected files are those created by beb_g09build.py
# Karl Irikura, NIST 
#
import sys
import os
import re
from beb_subs import *
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
    #   or 'beta' or 'both') being ionized, return: 
    #   (1) the spin multiplicity of the resulting ion,
    #   (2) the number of alpha and (3) beta electrons in the ion
    nalp_ion = nalp_neut
    nbet_ion = nbet_neut
    if orb_spin == 'alpha':
        nalp_ion -= 1
    elif orb_spin == 'beta' or orb_spin == 'both':
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
def merge_degenerate(df_orbs, thresh):
    # Orbitals are considered degenerate if their binding and kinetic energies are
    #    within thresh[0] energy units or thresh[1] percent
    # An exact match is required if thresh[0] == 0.
    # Return the modified DataFrame
    df_spin = df_orbs.groupby(['Spin'])
    if len(df_spin) > 1:
        # first combine orbitals of same spin
        dflist = []  # list of [combined orbitals]
        for spin, spinorb in df_spin:
            scombo = merge_degenerate(spinorb, thresh)
            dflist.append(scombo.copy())
        df_orbs = pd.concat(dflist)
    if thresh[0] == 0:
        # exact matching requested
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
    else:
        # Approximate match (of either thresh[0] or thresh[1]) is acceptable.
        # sort by energies (using MO label as tie-breaker) and reset index
        df_orbs = df_orbs.sort_values(by=['Energy', 'KE', 'MO'], ascending=[False, False, True]).reset_index(drop=True)
        diffs = df_orbs[['Energy', 'KE']].diff().abs() < thresh[0]  # absolute difference
        df_orbs['Match'] = diffs['Energy'] & diffs['KE']
        for col in ['Energy', 'KE']:
            diffs[col] = (100 * df_orbs[col].diff().abs() / df_orbs[col]) < thresh[1]  # percentage difference
        df_orbs['Match'] = df_orbs['Match'] | (diffs['Energy'] & diffs['KE']) 
        # loop in reverse, combining and deleting matching rows
        for i in range(len(df_orbs.index)-1, 0, -1):
            if df_orbs.loc[i, 'Match']:
                if False:
                    # [You can change the preceding test to 'True' to print a note about unusual conditions] 
                    if df_orbs.loc[i-1]['Spin'] != df_orbs.loc[i]['Spin']:
                        print('\tcombining orbitals of different spin: {:s} and {:s}'.format(df_orbs.loc[i-1]['MO'], df_orbs.loc[i]['MO']))
                    if abs(df_orbs.loc[i-1]['Orbital'] - df_orbs.loc[i]['Orbital']) > 1:
                        print('\tcombining non-consecutive orbitals: {:s} and {:s}'.format(df_orbs.loc[i-1]['MO'], df_orbs.loc[i]['MO']))
                # combine this orbital into the previous one:
                #   use the N-weighted average energies
                #   add the N-values together
                #   combine text rows reasonably
                #   delete this row
                Nsum = df_orbs.loc[i-1, 'N'] + df_orbs.loc[i, 'N']
                for col in ['Energy', 'KE']:
                    # N-weighted mean energy
                    Esum = df_orbs.loc[i-1, col] * df_orbs.loc[i-1, 'N'] + df_orbs.loc[i, col] * df_orbs.loc[i, 'N']
                    df_orbs.loc[i-1, col] = Esum / Nsum
                df_orbs.loc[i-1, 'N'] = Nsum
                for col in ['Special', 'Method']:
                    if df_orbs.loc[i-1, col] != df_orbs.loc[i, col]:
                        # these text fields are different; concatenate
                        df_orbs.loc[i-1, col] = '{:s} and {:s}'.format(df_orbs.loc[i-1, col], df_orbs.loc[i, col])
                if df_orbs.loc[i-1, 'Spin'] != df_orbs.loc[i, 'Spin']:
                    # designate mixed-spin orbital as 'both'
                    df_orbs.loc[i-1, 'Spin'] = 'both'
                #print('\tcombining orbitals {:s} and {:s}'.format(df_orbs.loc[i-1, 'MO'], df_orbs.loc[i, 'MO']))
                # delete the current row
                df_orbs.drop(df_orbs.index[i], inplace=True)
        # remove the 'Match' column
        del df_orbs['Match']
    # Different spins may have been combined; rebuild MO labels
    df_orbs = combine_MOspin(df_orbs, 'Orbital', 'Spin', 'MO')
    return df_orbs
##
def assign_n(contrib, cumul, degen, nmin):
    # assign the principal quantum number in Mulliken contributions of AOs to MOs
    if contrib > 1:
        # avoid n-skipping by disallowing contributions > 1
        cumul -= contrib - 1
    a = cumul / degen
    i = int(a)
    if a == i:
        n = nmin + i - 1
    else:
        n = nmin + i
    return int(n)
##
def special_n(fgau, Nthresh):
    # Read AO contributions to MOs, using 'read_AOpop_in_MOs()',
    #   return a dict of 'special' strings, e.g., '3p' or '3s' with
    #   MO labels as index
    # 'Nthresh' is the minimum contribution to trigger 'special' treatment
    # Take ECPs into account
    # Allow for multiple contributions toward the high-n threshold
    heavy_core_Z = 19     # Z >= this value: core orbitals may be labeled 'heavy_core'
    heavy_core_n = [1,2]  # n in this list:  orbital may be labeled 'heavy_core'
    df = read_AOpop_in_MOs(fgau)
    df_occ = df[df['Occ'] == 'occ'].copy()
    # Make an 'MO' label that includes any spin information
    df_occ = combine_MOspin(df_occ, 'Orbital', 'Spin', 'MO')
    # for each atom, sum the contributions for each L-value to determine n-value
    csum = df_occ.groupby(by=['Spin', 'Atom#', 'L']).cumsum()['Contrib']
    df_occ = df_occ.assign(Cumul = csum)
    grp = df_occ.groupby(by=['Spin', 'Atom#', 'L'])
    df_occ['n'] = 0
    # Find how many core electrons are replaced by pseudopotentials (on each atom)
    ppe = read_ECP_electrons(fgau)
    for name, group in grp:
        iatno = name[1] # center number in molecule (starting from 1)
        Lval = name[2]  # angular momentum value ('s', 'p', 'd', etc.)
        elem = group['Elem'].iloc[0]
        nstart = starting_n(Lval, ppe[iatno])
        degen = L_degeneracy(Lval)
        if ppe[iatno] == 0:
        #if nstart == starting_n(Lval, 0) and ppe[iatno] == 0:
            # routine 'assign_n' determines the values of n
            nvals = group.apply(lambda row: assign_n(row['Contrib'], row['Cumul'], degen, nstart), axis=1)
            df_occ['n'].update(nvals)
        else:
            # this atom has an ECP, so this AO is not 'special'
            pass
    # For each MO, see if it is dominated by one n.  If so, report that n.
    gmo = df_occ.groupby(by=['MO'])
    ndomin = {}  # dominant value of principal quantum number (n) for each MO
    ldomin = {}  # dominant angular momentum within ndomin[mo]
    for mo, group in gmo:
        nsum = {}  # sums for values of n
        for nval in set(group['n']):
            nsum[nval] = group[group['n'] == nval]['Contrib'].sum()
            if nsum[nval] > Nthresh:
                # dominant value of n for this MO
                ndomin[mo] = nval
        if not mo in ndomin:
            # no dominant value of n; combine "high" values
            hi = 0
            nwt = 0
            for n in nsum:
                if n >= 3:
                    # these are 'high' values of n
                    hi += nsum[n]
                    nwt += n * nsum[n]
            if hi > Nthresh:
                # threshold exceeded for mixed high-n 
                # choose the value of n nearest the weighted mean
                ndomin[mo] = int(round(nwt / hi))
            else:
                ndomin[mo] = ldomin[mo] = 0  # no dominant value of n
    # find the dominant value of L
    gmon = df_occ.groupby(by=['MO','n'])
    for (mo, n), group in gmon:
        if ndomin[mo] == n:
            i = group['Contrib'].idxmax()
            ldomin[mo] = group.loc[i, 'L']
    domlbl = {}
    for mo, group in gmo:
        if ndomin[mo] >= 3:
            # mark for (u+1)/n modification
            domlbl[mo] = str(ndomin[mo]) + ldomin[mo]
        else:
            domlbl[mo] = ''
            # check for 'heavy_core' situation
            hvcore = group.apply(lambda x: (elz(x['Elem']) >= heavy_core_Z) and (x['n'] in heavy_core_n), axis=1)
            if hvcore.any():
                # is the total core contribution enough to be 'dominant'?
                corefrac = group[hvcore].Contrib.sum()
                if corefrac > Nthresh:
                    domlbl[mo] = 'heavy_core'
    return domlbl
##
##
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
        elif suff == 'bupp' and not heavies:
            # no heavy atoms, so there will never be a BUPP file
            continue
        else:
            print('\n*** Optional file not found: {:s} ***'.format(fname))
            continue
    print('\n>>> Reading file {:s} <<<'.format(fname))
    nok = count_g09_success(fgau)
    if nok < 1:
        print('*** Calculation did not complete successfully')
        if suff == 'bu':
            sys.exit('\tThis file is required!')
        # Gaussian calculation had a problem, but this is not the essential,
        #   all-electron BU calculation
        continue
    if suff == 'opt':
        # check that geometry optimization converged
        if opt_success(fgau):
            print('Geometry optimization converged: ', end='')
        else:
            sys.exit('\tGeometry optimization FAILED\n')
        # check that there are no imaginary frequencies
        nimag = get_nimag(fgau)
        if nimag > 0:
            print('NOT a minimum! ({:d} imaginary vibrational frequencies)'.format(nimag))
            sys.exit()
        elif nimag < 0:
            print('No vibrational frequencies found')
        else:
            print('local minimum with nimag = 0')
    if suff == 'bu':
        # this is the only file that is essential to get a BEB result
        # extract stoichiometry, charge, and spin multiplicity
        if subminimal_basis(fgau):
            # less than minimal basis set (probably missing basis-set file)
            # Abort!
            sys.exit('Less than minimal basis set!\n'
                     'Check for missing basis set (*.gbs) files.')
        df = read_g09_stoichiometry(fgau)
        stoich = df.iloc[0].Stoich
        elems = df.iloc[0].Elements
        df = read_g09_charge_mult(fgau)
        charge = df.iloc[0].Charge
        mult = df.iloc[0].Mult
        nelec = read_g09_electrons(fgau).iloc[0]
        nalp = nelec.Alpha
        nbet = nelec.Beta
        nelec = nelec.Total
        # extract HF orbital energies and kinetic energies for occupied orbitals
        # orbital energies are negative binding energies
        buae = read_g09_oeke(fgau)[['Orbital', 'Spin', 'Energy', 'KE']]
        buae = combine_MOspin(buae, 'Orbital', 'Spin', 'MO')
        print('HF orbital and kinetic energies (hartree):')
        # In case there are no ECP calculations, cope with n>2 possibility
        zvals = elz(list(elems.keys()))
        heavies = max(zvals) > 10
        if heavies:
            # Find MOs that are dominated by high-n AOs
            specialThresh = 0.75  # AO contribution threshold for AO "special" status (1/n in BEB)
            domlbl = special_n(fgau, specialThresh)
            buae['Special'] = 'deal'
            for mo in domlbl:
                if domlbl[mo] == '':
                    buae.loc[buae['MO'] == mo, 'Special'] = 'none'
                else:
                    buae.loc[buae['MO'] == mo, 'Special'] = domlbl[mo]
            # if there are "special" MOs, print the AO contribution threshold
            m = buae[buae['Special'] != 'none'].shape[0]    # number of 'specials'
            if m > 0:
                print('"Special" AO contribution threshold is', specialThresh)
            print(buae[['MO', 'Energy', 'KE', 'Special']].to_string(index=False))
        else:
            # just create the 'Special' column
            buae['Special'] = 'none'
            print(buae[['MO', 'Energy', 'KE']].to_string(index=False))
        #uhf = ('beta' in buae.Spin.tolist() )   # flag to indicate separate beta orbitals
        uhf = (buae.Spin == 'beta').any()  # flag to indicate separate beta orbitals
        # extract crude molecular ionization threhold
        VIE, iespin, iemeth = weakest_binding(buae)
        mult_ion, nalp_ion, nbet_ion = ion_multiplicity(nalp, nbet, iespin)
        print('Crude VIE = {:.2f} eV to {:s} cation.'.format(VIE, spinname(mult_ion)))
        mulpops_bu = read_AOpop_in_MOs(fgau)   # will be used to match up MOs from different calculations
        mulpops_bu = combine_MOspin(mulpops_bu, 'Orbital', 'Spin', 'MO')  # add MO spin-labels
        mulpops_bu = mulpops_bu[mulpops_bu.Occ == 'occ']  # restrict attention to occupied orbitals
        VIE_method = 'Koopmans'
        VIE2 = 40.0  # guess for vertical double ionization energy (eV)
        VIE2_method = 'guessing'
        mult_dicat = 0  # unknown
    if suff == 'bupp':
        # get valence B/U from HF/ECP calculation
        bupp = read_g09_oeke(fgau)[['Orbital', 'Spin', 'Energy', 'KE']]
        bupp = combine_MOspin(bupp, 'Orbital', 'Spin', 'MO')
        if uhf:
            nelec_pp = len(bupp.index)
        else:
            nelec_pp = len(bupp.index) * 2
        ppcore = nelec - nelec_pp 
        print('ECP (pseudopotential) replaced {:d} electrons.'.format(ppcore))
        bupp['Method'] = 'ECP'
        # Find any MOs that are dominated by high-n AOs
        domlbl = special_n(fgau, specialThresh)
        bupp['Special'] = 'none'
        for mo in domlbl:
            if domlbl[mo]:
                bupp.loc[bupp['MO'] == mo, 'Special'] = domlbl[mo]
        print('HF/ECP orbital and kinetic energies (hartree):')
        # adjust orbital numbers to account for core and re-build MO labels
        bupp['Orbital'] += ppcore // 2
        bupp = combine_MOspin(bupp, 'Orbital', 'Spin', 'MO')
        print(bupp[['MO', 'Energy', 'KE', 'Special']].to_string(index=False))
        # note any light-atom core orbitals remaining
        lightcore = 0
        for el in elems:
            if elz(el) <= 10:
                lightcore += n_core({el: elems[el]})
        if lightcore > 0:
            lightcore //= 2  # convert number of electrons to number of orbitals
            print('({:d} light-atom core orbitals remain)'.format(lightcore))
        # extract VIE but don't replace the all-electron value from the previous file
        VIE_pp, ppspin, ppmeth = weakest_binding(bupp)
        mult_pp, nalp_pp, nbet_pp = ion_multiplicity(nalp-ppcore, nbet-ppcore, ppspin)
        mulpops_bupp = read_AOpop_in_MOs(fgau)
        mulpops_bupp.Orbital += ppcore // 2  # adjust orbital numbers
        mulpops_bupp = combine_MOspin(mulpops_bupp, 'Orbital', 'Spin', 'MO')  # add MO spin-labels
        mulpops_bupp = mulpops_bupp[mulpops_bupp.Occ == 'occ']  # restrict attention to occupied orbitals
        # code for UHF and RHF cases together (delegate separation to 'orbitalPopMatch()'
        momap_bupp = orbitalPopMatch(mulpops_bu, mulpops_bupp)
        #print(momap_bupp)
        if len(momap_bupp) > 0:
            print('Orbitals re-labeled, based mostly upon Mulliken populations:')
            bupp = relabelOrbitals(bupp, momap_bupp)
            print(bupp[['MO', 'Energy', 'KE', 'Special']].to_string(index=False))
        else:
            print('Orbital ordering is consistent with {:s}_bu.out'.format(mol))
        print('Crude VIE = {:.2f} eV to {:s} cation.'.format(VIE_pp, spinname(mult_pp)))
        ecp = True
    if suff == 'ept1':
        # get electron counts for neutral reference
        eneut = read_g09_electrons(fgau).iloc[0]
        nalpneut = eneut.Alpha
        nbetneut = eneut.Beta
        mult_neutref = nalpneut - nbetneut + 1  # the spin mult used for the neutral reference in the EPT1 calc
        # extract correlated orbital energies (negative of binding energies)
        minPS = 0.75
        ept = read_best_ept(fgau, minPS=minPS)
        print('(minimum acceptable pole strength = {:.2f})'.format(minPS))
        # degeneracies will cause a funny orbital ordering--sort before printing
        ept = ept.sort_values(by=['Spin', 'Orbital'])
        # also put columns in the same order as the HF results just displayed
        ept = ept[['Orbital', 'Spin', 'Energy', 'PS', 'Method']]
        ept = combine_MOspin(ept, 'Orbital', 'Spin', 'MO')
        print(ept[['MO', 'Energy', 'Method', 'PS']].to_string(index=False))
        # infer cation ground state
        VIE_ept, iespin, VIE_method = weakest_binding(ept)
        mult_ion, nalp_ion, nbet_ion = ion_multiplicity(nalp, nbet, iespin)
        print('VIE = {:.2f} eV from the {:s} neutral to the {:s} cation from {:s} theory.'.format(VIE_ept,
            spinname(mult_neutref), spinname(mult_ion), VIE_method))
        if abs(VIE_ept - VIE) > 1.5:
            # discrepancy between VIE values; 
            #   maybe EPT calculation has wrong number of orbitals
            print('*** Warning: Large difference of {:.2f} eV between Koopmans and EPT VIE values ***'.format(VIE_ept - VIE))
        VIE = VIE_ept
        mulpops_ept1 = read_AOpop_in_MOs(fgau)
        mulpops_ept1 = mulpops_ept1[mulpops_ept1.Occ == 'occ']  # restrict attention to occupied orbitals
        # adjust orbital numbering if necessary, to match other calculations
        ndispl = mulpops_bu.Orbital.max() - mulpops_ept1.Orbital.max()
        if ndispl > 0:
            # there must be an ECP in the ept1 calculation
            mulpops_ept1.Orbital += ndispl
        mulpops_ept1 = combine_MOspin(mulpops_ept1, 'Orbital', 'Spin', 'MO')
        momap_ept1 = orbitalPopMatch(mulpops_bu, mulpops_ept1)
        minmo = ept.Orbital.min()
        if len(momap_ept1) > 0:
            # possibly relabel some orbitals
            ept = relabelOrbitals(ept, momap_ept1)
            print('Orbitals re-labeled, based mostly upon Mulliken populations:')
            print(ept[['MO', 'Energy', 'Method', 'PS']].to_string(index=False))
        if False:
            keylist = list(momap_ept1.keys())
            keylist = []
            for key in keylist:
                if key < minmo:
                    del momap_ept1[key]
            if len(momap_ept1) > 0:
                print('Orbitals re-ordered, based mostly upon unsigned Mulliken populations:')
                print(momap_ept1)
                # loop once through the rows, changing MO numbers
                for idx in ept.index:
                    imo = ept.loc[idx, 'Orbital'] 
                    if imo in momap_ept1.keys():
                        # change this MO number
                        ept.loc[idx, 'Orbital'] = momap_ept1[imo]
                print(ept[['MO', 'Energy', 'Method', 'PS']].to_string(index=False))
    if suff == 'ept2':
        # get electron counts for cation reference
        ecation = read_g09_electrons(fgau).iloc[0]
        nalpion = ecation.Alpha
        nbetion = ecation.Beta
        mult_ionref = nalpion - nbetion + 1  # the spin multiplicity used for the +1 reference state
        # find ground state of (vertical) dication and (VIE2 - VIE) energy
        ept2 = read_best_ept(fgau, minPS=minPS)
        print('(minimum acceptable pole strength = {:.2f})'.format(minPS))
        print(ept2.to_string(index=False))
        # extract molecular ionization energy
        IEcat, ie2spin, ie2meth = weakest_binding(ept2)
        mult_dicat, nalpdbl, nbetdbl = ion_multiplicity(nalpion, nbetion, ie2spin)
        print('Second IE = {:.2f} eV from the {:s} cation to the {:s} dication from {:s} theory.'.format(IEcat,
            spinname(mult_ionref), spinname(mult_dicat), ie2meth))
        # infer double-ionization threshold, VIE2
        # was this calculation referenced to the ground state of the ion?
        if mult_ionref == mult_ion:
            # yes, the reference is the same multiplicity found in the EPT1 step
            VIE2_ept = VIE + IEcat
            VIE2 = VIE2_ept
            VIE2_method = '({:s} + {:s})'.format(VIE_method, ie2meth)
        else:
            # no, it was referenced to an excited state
            if abs(mult_ionref - mult_neutref) > 1:
                # this reference muliplicity could not be found in the EPT1 step
                print('--- I am confused about the EPT spin states to combine for VIE2 ---')
                print('Try re-running with reference cation spin multiplicity = {:d}'.format(mult_ion))
            else:
                if mult_ionref == mult_neutref - 1:
                    # reference was the low-spin result from EPT1; get that energy
                    targspin = 'alpha'
                if mult_ionref == mult_neutref + 1:
                    # get the high-spin EPT1 energy
                    targspin = 'beta'
                VIE_exc, iespin, exc_method = weakest_binding(ept.loc[ept['Spin'] == targspin])
                VIE2_ept = VIE_exc + IEcat
                VIE2 = VIE2_ept
                VIE2_method = '({:s} + {:s})'.format(exc_method, ie2meth)
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
worry_diff = 0.04   # fractional difference warning threshold for VIE discrepancies
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
        vd = VIE - VIE_ept
        if abs(vd/VIE_ept) > worry_diff:
            # maybe the SCF is an excited state in the CCSD(T) step
            print('*** Warning: Large difference of {:.2f} eV between CCSD(T) and {:s} VIE values ***'.format(vd, VIE_method))
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
        try:
            vd = VIE2 - VIE2_ept
            if abs(vd/VIE2_ept) > worry_diff:
                # maybe the SCF is an excited state in the CCSD(T) step
                print('*** Warning: Large difference of {:.2f} eV between CCSD(T) and {:s} VIE2 values ***'.format(vd, VIE2_method))
        except:
            print('No EPT value for VIE2 is available for comparison.')
        VIE2_method = 'CCSD(T)'
#
# done parsing Gaussian output files
# prepare BUtable for output
#
# combine the different DataFrames; start with the all-electron HF results
# Use the MO labels (e.g., '15', '15a', '15b') for matching up orbitals 
BUtable = buae.copy()
BUtable['Method'] = 'Koopmans'
# change the index to the 'MO' label
BUtable.set_index(['MO'], inplace=True)
if ecp:
    # change the index and use the PP kinetic energies
    print('Taking kinetic energies from ECP calculation.')
    bupp.set_index(['MO'], inplace=True)
    BUtable.update(bupp[['KE', 'Method', 'Special']])
try:
    # use available EPT binding energies
    print('Taking binding energies from EPT calculation.')
    ept.set_index(['MO'], inplace=True)
    BUtable.update(ept['Energy'])  # replace energies
    # If Method == ECP, don't erase that info
    ecprows = BUtable['Method'] == 'ECP'
    BUtable.update(ept['Method'])
    BUtable.loc[ecprows, 'Method'] = 'B ' + ept['Method'] + '; U ECP'
    # any 'Method' listed as NaN (from 'ept') should be 'B Koopmans; U ECP'
    nans = BUtable['Method'].isnull()
    BUtable.loc[nans, 'Method'] = 'B Koopmans; U ECP'
except:
    # no EPT calculations are available
    pass
    
if False:
    BUtable = buae.copy()
    BUtable['Method'] = 'Koopmans'
    # If this is UHF-based, replace orbital numbers with labels like '1a', '1b'
    if uhf:
        # re-label orbitals 
        BUtable = combine_MOspin(BUtable, 'Orbital', 'Spin', 'Label')
    else:
        # RHF case; simple numeric label
        BUtable['Label'] = BUtable['Orbital']
    # make 'Label' the index so merging the DataFrames is easier
    BUtable.set_index(['Label'], inplace=True)
    # Use kinetic energies from ECP calculation, if available
    if ecp:
        if uhf:
            bupp = combine_MOspin(bupp, 'Orbital', 'Spin', 'Label')
        else:
            bupp['Label'] = bupp['Orbital']
        bupp.set_index(['Label'], inplace=True)
        # replace the kinetic energies
        BUtable.update(bupp[['KE', 'Method', 'Special']])
    # If EPT calculations were done, use those values where available
    try:
        if uhf:
            ept = combine_MOspin(ept, 'Orbital', 'Spin', 'Label')
        else:
            ept['Label'] = ept['Orbital']
        ept.set_index(['Label'], inplace=True)
        # Replace with EPT binding energies as available
        BUtable.update(ept['Energy'])
        # If Method == ECP, don't erase that info
        ecprows = BUtable['Method'] == 'ECP'
        BUtable.update(ept['Method'])
        BUtable.loc[ecprows, 'Method'] = 'B ' + ept['Method'] + '; U ECP'
        # any 'Method' listed as NaN (from 'ept') should be 'B Koopmans; U ECP'
        nans = BUtable['Method'].isnull()
        BUtable.loc[nans, 'Method'] = 'B Koopmans; U ECP'
    except:
        # no EPT calculations are available
        pass
# put 'MO' on the same footing as the other columns to cope with degeneracies
BUtable.reset_index(inplace=True)
# add columns: 'N', 'Q', 'DblIon'
BUtable['N'] = 1 if uhf else 2
BUtable['Q'] = 1
BUtable['DblIon'] = 'No'
# make binding energies positive and convert to eV
BUtable['Energy'] = BUtable['Energy'].apply(hartree_eV, args=('to_eV', -1))
BUtable['KE'] = BUtable['KE'].apply(hartree_eV)
# mark 'DblIon' column where B > VIE2
BUtable.ix[BUtable['Energy'] > VIE2, 'DblIon'] = 'Yes'
# combine any degenerate orbitals before installing CCSD(T) binding energy
if True:
    # [Change the test on the preceding line to 'False' if you don't want degenerate
    #   orbitals combined.]
    # Consolidate degenerate orbitals
    # [You can change the thresholds on the following line]
    degen_thresh = (0.05, 0.10)  # (eV, percent) orbital energy differences considered negligible
    if degen_thresh[0] == 0:
        print('Combining any exactly degenerate orbitals')
    else:
        print('Combining any nearly degenerate orbitals: threshold = ({:.2f} eV or {:.2f}%)'.format(*degen_thresh))
    BUtable = merge_degenerate(BUtable, degen_thresh)
# Include any CCSD(T) threshold values
if VIE_method == 'CCSD(T)':
    print('Taking threshold energy from {:s} calculation.'.format(VIE_method))
    # replace binding energy for HOMO
    if uhf and (mult_ion > mult):
        # apply VIE to the highest beta orbital
        homo_spin = 'beta'
    else:
        # apply VIE to the highest alpha orbital
        if uhf:
            homo_spin = 'alpha'
        else:
            homo_spin = 'both'
    #homo = BUtable.loc[BUtable['Spin'] == homo_spin]['Orbital'].idxmax()
    homo = BUtable.loc[BUtable['Spin'] == homo_spin]['Energy'].idxmin()
    BUtable.loc[homo, 'Energy'] = VIE
    if 'U ECP' in BUtable.loc[homo, 'Method']:
        # Don't erase ECP information for HOMO
        BUtable.loc[homo, 'Method'] = 'B ' + VIE_method + '; U ECP'
    else:
        BUtable.loc[homo, 'Method'] = VIE_method
# sort by orbital number and then by spin
BUtable = BUtable.sort_values(by=['Orbital', 'Spin'])
# rename and rearrange columns as needed for input to 'beb_tbl.pl' (BEB program)
BUtable = BUtable.rename(index=str, columns={'MO': '#MO', 'KE': 'U/eV', 'Energy': 'B/eV'})
BUtable = BUtable.rename(index=str, columns={'Method': 'Remarks'})
BUtable = BUtable[['#MO', 'B/eV', 'U/eV', 'N', 'Q', 'DblIon', 'Special', 'Remarks']]
# last-minute check of electron count
if BUtable['N'].sum() != (nalp + nbet):
    # We've lost or gained electrons!
    print('*** Error: there are {:d} electrons but should be {:d}! ***'.format(BUtable['N'].sum(), nalp+nbet))
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
