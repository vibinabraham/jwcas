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
#local = True

#PYSCF inputs
mol = gto.Mole()
mol.atom = '''
H 0 0 0
H 0 0 1
H 0 2 0
H 0 2 2.1
'''


mol.max_memory = 1000 # MB
mol.charge = +0
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

################################    ACTIVE ORBITALS
n_start = 0
n_core = n_start 
n_stop  = n_orb 

print("---------------------------")
print("          CAS-CI           ")
print("---------------------------")

print("Active Space Orbitals:")

def active(h,g,n_a,n_b,nel,n_start,n_stop):
# {{{
    n_active = n_stop - n_start

    new_na = n_a - n_start
    new_nb = n_b - n_start
    

    ha = np.zeros((n_active,n_active))
    ga = np.zeros((n_active,n_active,n_active,n_active))

    for j in range(n_start,n_stop):
        J = j - n_start
        for k in range(n_start,n_stop):
            K = k - n_start
            ha[J,K] = h[j,k]
            for l in range(n_start,n_stop):
                L = l - n_start
                for m in range(n_start,n_stop):
                    M = m - n_start
                    ga[J,K,L,M] = g[j,k,l,m]

    return ha,ga,new_na,new_nb,n_active
# }}}


h_active,g_active,n_a,n_b,n_orb = active(h,g,n_a,n_b,nel,n_start,n_stop)

#alpha first beta second ordering

h_spin, g_spin = spatial_2_spin_ab(n_orb,h_active, g_active)


jw = jordon_wigner_string(n_orb,h_spin,g_spin)


for i in range(len(jw)):
    print("%12.8f %6d   "%(jw[i].coeff, jw[i].dn),jw[i].string)


HH = form_jw_string_ham(n_orb,jw)

ee = np.linalg.eigvalsh(HH)

Ham = jordon_wigner_form_mat(n_orb,h_spin,g_spin) #forming the 4^N hamiltonian
e2 = np.linalg.eigvalsh(Ham)
print(e2[:4])
print(ee[:4])
assert(np.allclose(HH,HH.T))
print(len(jw))
