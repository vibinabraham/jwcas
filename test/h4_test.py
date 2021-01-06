import numpy as np
from sys import argv
from math import factorial
import copy as cp
from jwcas import *

np.set_printoptions(suppress=True, precision=6, linewidth=800)


def test_1():
    molecule = '''
    H      0.00       0.00       0.00
    H      0.00       0.00       2.00
    H      0.00       2.00       0.00
    H      0.00       2.00       2.00
    '''

    charge = 0
    spin  = 0
    basis_set = 'sto-3g'
    orb_basis = 'scf' #scf, lowdin, boys, PM etc
    nelec = 4
    n_orb = 4

    #Integrals from pyscf
    pmol = PyscfHelper()
    pmol.init(molecule,charge,spin,basis_set,orb_basis,cas=False,cas_nel=nelec)

    C = pmol.C
    h = pmol.h
    g = pmol.g
    ecore = pmol.ecore


    # Part 1: JWCAS using JW strings

    # transform the integrals to  spin integrals
    # alpha first beta second ordering
    h_spin, g_spin = spatial_2_spin_ab(n_orb,h, g)

    # Generate the JW strings list
    jw = jordon_wigner_string(n_orb,h_spin,g_spin)

    # Printing
    for i in range(len(jw)):
        print("%12.8f %6d   "%(jw[i].coeff, jw[i].dn),jw[i].string)

    # generate the Hamiltonian using the JW string
    HH = form_jw_string_ham(n_orb,jw)

    # find eigenvalues. since the Full H is diagonalized the eigenvalue can be from any fock space
    ee = np.linalg.eigvalsh(HH)




    # Part 2: Form the Matrix directly using the transformation
    Ham = jordon_wigner_form_mat(n_orb,h_spin,g_spin) #forming the 4^N hamiltonian
    e2 = np.linalg.eigvalsh(Ham)

    assert(np.allclose(e2,ee.T))
    assert(np.allclose(HH,HH.T))
    print("Number of JW strings: ",len(jw))

    from pyscf import fci
    cisolver = fci.direct_spin1.FCI(pmol.mol)
    efci, ci = cisolver.kernel(h, g, n_orb, nelec, ecore=ecore)

    print(efci)
    print(ee[0]+ecore)
    print(e2[0]+ecore)

    assert(abs(efci-ee[0]-ecore)<1e-9)
    assert(abs(efci-e2[0]-ecore)<1e-9)

if __name__== "__main__":
    test_1() 
