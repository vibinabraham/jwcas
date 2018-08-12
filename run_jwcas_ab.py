from pyscf import gto, scf, mcscf, fci, ao2mo, lo, molden
import numpy as np
from sys import argv
from math import factorial
import copy as cp
from jwcasci import *

np.set_printoptions(suppress=True, precision=6, linewidth=800)

# set memory requirements
numpy_memory = 2

#Localise the orbitals after SCF??
local = False
local = True

#PYSCF inputs
mol = gto.Mole()
mol.atom = '''
H 0 0 0
H 0 0 1
H 0 1 1
H 1 1 1
'''

mol.max_memory = 1000 # MB
mol.charge = -0
mol.spin = 0
mol.basis = 'sto-3g'

mol.build()

#SCF 
mf = scf.ROHF(mol).run()
C = mf.mo_coeff #MO coeffs
n_b , n_a = mol.nelec #beta and alpha electrons
nel = n_a + n_b
n_orb = mol.nao_nr() #total orbitals



##READING INTEGRALS FROM PYSCF
E_nu = gto.Mole.energy_nuc(mol)
T = mol.intor('int1e_kin_sph')
V = mol.intor('int1e_nuc_sph') 
hcore = T + V
S = mol.intor('int1e_ovlp_sph')
g = mol.intor('int2e_sph')


##AO 2 MO Transformation: local or scf
h = C.T @ hcore @ C

if local == False:
    g = np.einsum("pqrs,pl->lqrs",g,C)
    g = np.einsum("lqrs,qm->lmrs",g,C)
    g = np.einsum("lmrs,rn->lmns",g,C)
    g = np.einsum("lmns,so->lmno",g,C)

    h = C.T @ hcore @ C

if local == True:
    Cl = lo.orth.orth_ao(mol)
    hl = Cl.T @ hcore @ Cl
    gl = np.einsum("pqrs,pl->lqrs",g ,Cl)
    gl = np.einsum("lqrs,qm->lmrs",gl,Cl)
    gl = np.einsum("lmrs,rn->lmns",gl,Cl)
    gl = np.einsum("lmns,so->lmno",gl,Cl)

    g = gl
    h = hl


print("\nNuclear Repulsion        :%16.10f " %E_nu)
print("Electronic SCF energy    :%16.10f " %mf.e_tot)
print("Number of Orbitals       :%16d\n" %(n_orb))


#alpha first beta second ordering
h_spin = spatial_2_spin_oei_ab(n_orb,h)
g_spin = spatial_2_spin_eri_ab(n_orb,g)

#####Operators for JW
ap = np.array([[0, 0], [1, 0]])     #creation operator
am = np.array([[0, 1], [0, 0]])     #annihilation operator
no = np.array([[0, 0], [0, 1]])     #number operator
ho = np.array([[1, 0], [0, 0]])     #hole operator
I2 = np.array([[1, 0], [0, 1]])     #identity operator
Iz = np.array([[1, 0], [0, -1]])    #pauli z operator

####START of JW-FCI
Ham = np.zeros((4**n_orb,4**n_orb))

#uncomment to debug each term
#mat1 =  jordon_wigner_1e(n_orb,h_spin,ap,am,no,ho,I2,Iz)
#mat2 =  jordon_wigner_2e(n_orb,g_spin,ap,am,no,ho,I2,Iz)
#mat3 =  jordon_wigner_3e(n_orb,g_spin,ap,am,no,ho,I2,Iz)
#mat4 =  jordon_wigner_4e(n_orb,g_spin,ap,am,no,ho,I2,Iz)
#
#Ham +=  mat1
#Ham +=  mat2
#Ham +=  mat3
#Ham +=  mat4

Ham = jordon_wigner_form_mat(n_orb,h_spin,g_spin,ap,am,no,ho,I2,Iz) #forming the 4^N hamiltonian

assert(np.allclose(Ham, Ham.T,rtol=1e-05, atol=1e-05))


N = form_N(n_orb,ap,am,no,ho,I2,Iz) #forming the 4^N number operator

Sz = form_Sz_ab(n_orb,ap,am,no,ho,I2,Iz) #forming the 4^N Ms operator

S2 = form_S2_ab(n_orb,ap,am,no,ho,I2,Iz) #forming the 4^N <S2> operator
assert(np.allclose(S2, S2.T,rtol=1e-05, atol=1e-05))

###Cutting down the matrix based on number of electrons
Hsmall  = occ_block(n_orb,Ham,N,n_a+n_b)
Szsmall = occ_block(n_orb,Sz,N,n_a+n_b)
S2small = occ_block(n_orb,S2,N,n_a+n_b)


###Cutting down the matrix based on Ms Symmetry
Htiny = msblock(n_orb,Hsmall,Szsmall,n_b,n_a)
S2tiny = msblock(n_orb,S2small,Szsmall,n_b,n_a)


#running PYSCF FCI for comparison
cisolver = fci.direct_spin1.FCI(mol)
efci, ci = cisolver.kernel(h, g, h.shape[1], mol.nelec, ecore=mol.energy_nuc())

#JWCAS_CI eigenvalues and <S2>
e_val,e_vec = np.linalg.eigh(Htiny)
S2coupled = e_vec.T @ S2tiny @ e_vec

print("Ground State  JW FCI     :%16.10f"   %(e_val[0])) 
print("Ground State PYSCF FCI   :%16.10f"   %(efci - E_nu) )
print("Error from Reality       :%16.10f"   %(efci - E_nu - e_val[0])) 
assert(abs(efci - E_nu - e_val[0]) <1e-10)



print("\nEigenvalues of Ms = 0 Block:")
print("  ---------------------------------------")
print("      Energies                   <S2>")
print("  ---------------------------------------")
for i in range(0,e_val.shape[0]):
    print("%16.10f        %16.6f   " %(e_val[i],abs(S2coupled[i,i])))
