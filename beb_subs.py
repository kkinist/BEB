# Routines for use in running BEB calculations.
# Karl Irikura 
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
