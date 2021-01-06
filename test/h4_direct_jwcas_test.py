import numpy as np
from sys import argv
from math import factorial
import copy as cp
from jwcas import *

np.set_printoptions(suppress=True, precision=6, linewidth=800)

def test1():
    molecule = '''
    H 0 0 0
    H 0 0 1
    H 0 2 0
    H 0 2 1
    '''

    charge = 0
    spin  = 0
    basis_set = 'sto-3g'
    orb_basis = 'scf' #scf, lowdin, boys, PM etc
    n_a = 2; n_b = 2
    nelec = n_a + n_b
    n_orb = 4

    #Integrals from pyscf
    pmol = PyscfHelper()
    pmol.init(molecule,charge,spin,basis_set,orb_basis,cas=False,cas_nel=nelec)

    C = pmol.C
    h = pmol.h
    g = pmol.g
    ecore = pmol.ecore
    mol = pmol.mol


    #alpha first beta second ordering
    h_spin, g_spin = spatial_2_spin_ab(n_orb,h, g)

    Ham = jordon_wigner_form_mat(n_orb,h_spin,g_spin) #forming the 4^N hamiltonian
    assert(np.allclose(Ham, Ham.T,rtol=1e-05, atol=1e-05))


    N = form_N(n_orb) #forming the 4^N number operator
    Sz = form_Sz_ab(n_orb) #forming the 4^N Ms operator
    S2 = form_S2_ab(n_orb) #forming the 4^N <S2> operator
    assert(np.allclose(S2, S2.T,rtol=1e-05, atol=1e-05))

    ###Cutting down the matrix based on number of electrons
    Hsmall  = occ_block(n_orb,Ham,N,n_a+n_b)
    Szsmall = occ_block(n_orb,Sz,N,n_a+n_b)
    S2small = occ_block(n_orb,S2,N,n_a+n_b)

    ###Cutting down the matrix based on Ms Symmetry
    Htiny = msblock(n_orb,Hsmall,Szsmall,n_b,n_a)
    S2tiny = msblock(n_orb,S2small,Szsmall,n_b,n_a)

    #running PYSCF FCI for comparison
    from pyscf import fci
    cisolver = fci.direct_spin1.FCI(mol)
    efci, ci = cisolver.kernel(h, g, n_orb, nelec, ecore=ecore)

    #JWCAS_CI eigenvalues and <S2>
    e_val,e_vec = np.linalg.eigh(Htiny)
    S2coupled = e_vec.T @ S2tiny @ e_vec

    print("Ground State  JW FCI     :%16.10f"   %(e_val[0]+ecore)) 
    print("Ground State PYSCF FCI   :%16.10f"   %(efci) )
    print("Error                    :%16.10f"   %(efci - ecore - e_val[0])) 
    assert(abs(efci - ecore - e_val[0]) <1e-10)



    print("\nEigenvalues of Ms = 0 Block:")
    print("  ---------------------------------------")
    print("      Energies                   <S2>")
    print("  ---------------------------------------")
    for i in range(0,e_val.shape[0]):
        print("%16.10f        %16.6f   " %(e_val[i],abs(S2coupled[i,i])))

if __name__== "__main__":
    test1() 
