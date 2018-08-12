from pyscf import gto, scf, mcscf, fci, ao2mo, lo, molden
import numpy as np
from sys import argv
from math import factorial
import copy as cp

np.set_printoptions(suppress=True, precision=6, linewidth=800)

def nCr(n, r):
    #{{{
    if n<r:
        return 0
    else:
        return factorial(n) // factorial(r) // factorial(n-r)
    #}}}

def occ_block(n_orb,H,Nmat,occ):
    #{{{
    """
    Function to give a block of the operator with a given specified occupation number

    Parameters:
        n_orb   : total number of orbitals
        H       : the operator that needs to be cut down
        Nmat    : the number operator used to cut down
        occ     : the specified occupancy

    """
    n_occ = nCr(2*n_orb,occ) #size of the cut down matrix
    N = H.shape[0]           #Full Size
    Mb =  np.empty((N,0))    #Empty row matrix to append
    Hsmall = np.empty((0,n_occ)) #Empty column matrix

    for i in range(0,N):
        if  Nmat[i,i] == occ:
            hh = H[:,i].reshape(N,1)
            Mb = np.bmat([Mb,hh])

    for i in range(0,N):
        if  Nmat[i,i] == occ:
            Kp = Mb[i,:].reshape(1,n_occ)
            Hsmall = np.bmat('Hsmall;Kp')

    return Hsmall
    #}}}

def msblock(n,H,Sz,na,nb):
     #{{{
    """
    Takes in a matrix which has a given occupancy and gives a given sz symmetry 
    """
    ms = abs(na - nb)
    
    n_ms = nCr(n,na) * nCr(n,nb) 
    N = H.shape[0]
    size_ms = 0
    Mb =  np.empty((N,0))
    for i in range(0,N):
        if  Sz[i,i] == ms:
            hh = H[:,i].reshape(N,1)
            Mb = np.bmat([Mb,hh])
            size_ms +=1
    H_ms = np.empty((0,size_ms))
    for i in range(0,N):
        if  Sz[i,i] == ms:
            Kp = Mb[i,:].reshape(1,n_ms)
            H_ms = np.bmat('H_ms;Kp')
    return H_ms
    #}}}

def fci_block(n,H,Num,Sz,na,nb):
     #{{{
    """
    Takes in a matrix which has a given occupancy and gives a given sz symmetry 
    """
    ms = abs(na - nb)
    nel = na+nb
    
    n_ms = nCr(n,na) * nCr(n,nb) 
    N = H.shape[0]
    size_ms = 0
    Mb =  np.empty((N,0))
    for i in range(0,N):
        if  Sz[i,i] == ms and Num[i,i] == nel:
            hh = H[:,i].reshape(N,1)
            Mb = np.bmat([Mb,hh])
            size_ms +=1
    H_ms = np.empty((0,size_ms))
    for i in range(0,N):
        if  Sz[i,i] == ms and Num[i,i] == nel:
            Kp = Mb[i,:].reshape(1,n_ms)
            H_ms = np.bmat('H_ms;Kp')
    return H_ms
    #}}}

def form_N(n,ap,am,no,ho,I2,Iz):
# {{{
    N = np.zeros((4**n,4**n))
    for p in range(0,2*n):
        Ia = np.eye(np.power(2,p))
        Ic = np.eye(np.power(2,2*n-p-1))
        N += np.kron(Ia,np.kron(no,Ic))
    return N
# }}}

def spatial_2_spin_eri(n_orb,g):
# {{{
    ##USING LOOPS
    #n = n_orb
    #g_spin = np.zeros((2*n_orb,2*n_orb,2*n_orb,2*n_orb))

    #for p in range(0,n_orb):
    #    for q in range(0,n_orb):
    #        for r in range(0,n_orb):
    #            for s in range(0,n_orb):
    #                g_spin[p,q,r,s] = g[p,q,r,s] 
    #                g_spin[p+n,q+n,r,s] = g[p,q,r,s] 
    #                g_spin[p,q,r+n,s+n] = g[p,q,r,s] 
    #                g_spin[p+n,q+n,r+n,s+n] = g[p,q,r,s] 
    
    g_spin = np.kron(np.eye(2),g)
    g_spin = np.kron(np.eye(2),g_spin.T)

    g_spin = g_spin.swapaxes(1,2)
    return g_spin
# }}}

def spatial_2_spin_oei(n_orb,h):
# {{{
    
    h_spin = np.kron(np.eye(2),h)

    return h_spin
# }}}

def jordon_wigner_1e(n,h_spin,ap,am,no,ho,I2,Iz):
# {{{
    Ham = np.zeros((4**n,4**n))
    for p in range(0,2*n):
        Ia = np.eye(np.power(2,p))
        Ic = np.eye(np.power(2,2*n-p-1))
        
        Ham += h_spin[p,p]*np.kron(Ia,np.kron(no,Ic))

        for q in range(p+1,2*n):

            Zb = np.eye(1)
            for k in range(p+1,q):
                Zb = np.kron(Zb,Iz)

            Ic = np.eye(np.power(2,2*n-q-1))

            Ham += h_spin[p,q]*np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(am,Ic))))
            Ham += h_spin[p,q]*np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(ap,Ic))))
    return Ham
# }}}

def jordon_wigner_2e(n,g_spin,ap,am,no,ho,I2,Iz):
# {{{
    Ham = np.zeros((4**n,4**n))
    for p in range(0,2*n):

        Za = np.eye(1)
        for k in range(0,p):
            Za = np.kron(Za,Iz)

        Ia = np.eye(np.power(2,p))
        
        for q in range(p+1,2*n):

            Zb = np.eye(1)
            for k in range(p+1,q):
                Zb = np.kron(Zb,Iz)

            Ib = np.eye(np.power(2,q-p-1))

            Ic = np.eye(np.power(2,2*n-q-1))
            
            g_val = g_spin[p,q,p,q] - g_spin[p,q,q,p]
            Ham += g_val * np.kron(Ia,np.kron(no,np.kron(Ib,np.kron(no,Ic))))
    return Ham
# }}}

def jordon_wigner_3e(n,g_spin,ap,am,no,ho,I2,Iz):
# {{{
    Ham = np.zeros((4**n,4**n))
    for p in range(0,2*n):

        Za = np.eye(1)
        for k in range(0,p):
            Za = np.kron(Za,Iz)

        Ia = np.eye(np.power(2,p))
        
        for q in range(p+1,2*n):

            Zb = np.eye(1)
            for k in range(p+1,q):
                Zb = np.kron(Zb,Iz)

            Ib = np.eye(np.power(2,q-p-1))
            
            for r in range(q+1,2*n):

                Zc = np.eye(1)
                for k in range(q+1,r):
                    Zc = np.kron(Zc,Iz)

                Ic = np.eye(np.power(2,r-q-1))
                
                Id = np.eye(np.power(2,2*n-r-1))
                    

                g_val1 = g_spin[p,q,r,q] - g_spin[p,q,q,r] 
                g_val2 = g_spin[p,r,q,r] - g_spin[p,r,r,q] 
                g_val3 = g_spin[q,p,r,p] - g_spin[q,p,p,r] 
                g_val4 = g_spin[r,p,q,p] - g_spin[r,p,p,q] 
                g_val5 = g_spin[q,r,p,r] - g_spin[q,r,r,p] 
                g_val6 = g_spin[r,q,p,q] - g_spin[r,q,q,p] 
                print("g_val1",g_val1)
                print("g_val2",g_val2)
                print("g_val3",g_val3)
                print("g_val4",g_val4)
                print("g_val5",g_val5)
                print("g_val6",g_val6)
                         
                Ham += g_val1 * np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(-no,np.kron(Zc,np.kron(am,Id))))))
                #Ham += g_val1 * np.kron(Ia,np.kron(-ap,np.kron(Zb,np.kron(ho,np.kron(Zc,np.kron(ap,Id))))))

                Ham += g_val2 * np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(am,np.kron(Ic,np.kron(no,Id))))))
                #Ham += g_val2 * np.kron(Ia,np.kron(-am,np.kron(Zb,np.kron(-ap,np.kron(Ic,np.kron(ho,Id))))))

                Ham += g_val3 * np.kron(Ia,np.kron(no,np.kron(Ib,np.kron(ap,np.kron(Zc,np.kron(am,Id))))))
                #Ham += g_val3 * np.kron(Ia,np.kron(ho,np.kron(Ib,np.kron(am,np.kron(Zc,np.kron(ap,Id))))))

                Ham += g_val4 * np.kron(Ia,np.kron(no,np.kron(Ib,np.kron(am,np.kron(Zc,np.kron(ap,Id))))))
                #Ham += g_val4 * np.kron(Ia,np.kron(ho,np.kron(Ib,np.kron(ap,np.kron(Zc,np.kron(am,Id))))))

                Ham += g_val5 * np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(ap,np.kron(Ic,np.kron(no,Id))))))
                #Ham += g_val5 * np.kron(Ia,np.kron(-ap,np.kron(Zb,np.kron(-am,np.kron(Ic,np.kron(ho,Id))))))

                Ham += g_val6 * np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(-no,np.kron(Zc,np.kron(ap,Id))))))
                #Ham += g_val6 * np.kron(Ia,np.kron(-ap,np.kron(Zb,np.kron(-ho,np.kron(Zc,np.kron(am,Id))))))
    return Ham
# }}}

def jordon_wigner_4e(n,g_spin,ap,am,no,ho,I2,Iz):
# {{{
    Ham = np.zeros((4**n,4**n))
    for p in range(0,2*n):

        Za = np.eye(1)
        for k in range(0,p):
            Za = np.kron(Za,Iz)

        Ia = np.eye(np.power(2,p))
        
        for q in range(p+1,2*n):

            Zb = np.eye(1)
            for k in range(p+1,q):
                Zb = np.kron(Zb,Iz)

            Ib = np.eye(np.power(2,q-p-1))
            
            for r in range(q+1,2*n):

                Zc = np.eye(1)
                for k in range(q+1,r):
                    Zc = np.kron(Zc,Iz)

                Ic = np.eye(np.power(2,r-q-1))
                
                for s in range(r+1,2*n):

                    Zd = np.eye(1)
                    for k in range(r+1,s):
                        Zd = np.kron(Zd,Iz)

                    Ie = np.eye(np.power(2,2*n-s-1))
                    

                    # Ia P Zb Q Ic (-r) Zd s Ie
                    #print(p,q,r,s)
                    g_val1 = g_spin[p,q,r,s] - g_spin[p,q,s,r]
                    g_val2 = g_spin[p,r,q,s] - g_spin[p,r,s,q]
                    g_val3 = g_spin[p,s,q,r] - g_spin[p,s,r,q]
                    print("gval1    :",g_spin[p,q,s,r])
                    print("gval2    :",g_spin[p,r,s,q])
                    print("gval3    :",g_spin[p,s,r,q])


                    Ham +=  g_val1 * np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(ap,np.kron(Ic,np.kron(am,np.kron(Zd,np.kron(am,Ie)))))))) 
                              
                    Ham +=  g_val1 * np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(am,np.kron(Ic,np.kron(ap,np.kron(Zd,np.kron(ap,Ie)))))))) 
                              
                    Ham +=  g_val2 * np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(am,np.kron(Ic,np.kron(ap,np.kron(Zd,np.kron(am,Ie)))))))) 
                              
                    Ham +=  g_val2 * np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(ap,np.kron(Ic,np.kron(am,np.kron(Zd,np.kron(ap,Ie)))))))) 
                              
                    Ham +=  g_val3 * np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(am,np.kron(Ic,np.kron(am,np.kron(Zd,np.kron(ap,Ie)))))))) 
                              
                    Ham +=  g_val3 * np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(ap,np.kron(Ic,np.kron(ap,np.kron(Zd,np.kron(am,Ie)))))))) 
    return Ham
# }}}

def jordon_wigner_form_mat(n,h_spin,g_spin,ap,am,no,ho,I2,Iz):
# {{{
    Ham = np.zeros((4**n,4**n))
    for p in range(0,2*n):

        Ia = np.eye(np.power(2,p))
        Ic = np.eye(np.power(2,2*n-p-1))

        Za = np.eye(1)
        for k in range(0,p):
            Za = np.kron(Za,Iz)

        Ham += h_spin[p,p]*np.kron(Ia,np.kron(no,Ic))
        
        for q in range(p+1,2*n):

            Zb = np.eye(1)
            for k in range(p+1,q):
                Zb = np.kron(Zb,Iz)

            Ib = np.eye(np.power(2,q-p-1))

            Ic = np.eye(np.power(2,2*n-q-1))
            

            Ham += h_spin[p,q]*np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(am,Ic))))
            Ham += h_spin[p,q]*np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(ap,Ic))))

            g_val = g_spin[p,q,p,q] - g_spin[p,q,q,p]
            Ham += g_val * np.kron(Ia,np.kron(no,np.kron(Ib,np.kron(no,Ic))))

            for r in range(q+1,2*n):

                Zc = np.eye(1)
                for k in range(q+1,r):
                    Zc = np.kron(Zc,Iz)

                Ic = np.eye(np.power(2,r-q-1))

                Id = np.eye(np.power(2,2*n-r-1))
                    

                g_val1 = g_spin[p,q,r,q] - g_spin[p,q,q,r] 
                g_val2 = g_spin[p,r,q,r] - g_spin[p,r,r,q] 
                g_val3 = g_spin[q,p,r,p] - g_spin[q,p,p,r] 
                g_val4 = g_spin[r,p,q,p] - g_spin[r,p,p,q] 
                g_val5 = g_spin[q,r,p,r] - g_spin[q,r,r,p] 
                g_val6 = g_spin[r,q,p,q] - g_spin[r,q,q,p] 
                         
                Ham += g_val1 * -np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(no,np.kron(Zc,np.kron(am,Id))))))

                Ham += g_val2 * np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(am,np.kron(Ic,np.kron(no,Id))))))

                Ham += g_val3 * np.kron(Ia,np.kron(no,np.kron(Ib,np.kron(ap,np.kron(Zc,np.kron(am,Id))))))

                Ham += g_val4 * np.kron(Ia,np.kron(no,np.kron(Ib,np.kron(am,np.kron(Zc,np.kron(ap,Id))))))

                Ham += g_val5 * np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(ap,np.kron(Ic,np.kron(no,Id))))))

                Ham += g_val6 * -np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(no,np.kron(Zc,np.kron(ap,Id))))))
                
                for s in range(r+1,2*n):

                    Zd = np.eye(1)
                    for k in range(r+1,s):
                        Zd = np.kron(Zd,Iz)

                    Ie = np.eye(np.power(2,2*n-s-1))
                    

                    g_val1 = g_spin[p,q,r,s] - g_spin[p,q,s,r]
                    g_val2 = g_spin[p,r,q,s] - g_spin[p,r,s,q]
                    g_val3 = g_spin[p,s,q,r] - g_spin[p,s,r,q]

                    Ham +=  g_val1 * np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(ap,np.kron(Ic,np.kron(am,np.kron(Zd,np.kron(am,Ie)))))))) 
                              
                    Ham +=  g_val1 * np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(am,np.kron(Ic,np.kron(ap,np.kron(Zd,np.kron(ap,Ie)))))))) 
                              
                    Ham +=  g_val2 * np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(am,np.kron(Ic,np.kron(ap,np.kron(Zd,np.kron(am,Ie)))))))) 
                              
                    Ham +=  g_val2 * np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(ap,np.kron(Ic,np.kron(am,np.kron(Zd,np.kron(ap,Ie)))))))) 
                              
                    Ham +=  g_val3 * np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(am,np.kron(Ic,np.kron(am,np.kron(Zd,np.kron(ap,Ie)))))))) 
                              
                    Ham +=  g_val3 * np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(ap,np.kron(Ic,np.kron(ap,np.kron(Zd,np.kron(am,Ie)))))))) 
    return Ham
# }}}

def form_Sz(n,ap,am,no,ho,I2,Iz):
# {{{
    It = np.eye(2**n)
    Sz = np.zeros((4**n,4**n))
    for i in range(0,n):
        Ia1 = np.eye(np.power(2,i))
        Ia2 = np.eye(np.power(2,n-i-1))
        Na = np.kron(Ia1,np.kron(no,Ia2))
        Nb = np.kron(Ia1,np.kron(no,Ia2))
        Sz += (np.kron(Na,It) - np.kron(It,Nb))
    return Sz
# }}}

def form_S2(n_orb,ap,am,no,ho,I2,Iz):
# {{{
    S2 = np.zeros((4**n_orb,4**n_orb))
    s2 = np.array([[0,0],[0,0.75]])

    Ialpha = np.eye(np.power(2,n_orb))
    Zalpha = np.eye(1)
    for k in range(0,n_orb):
        Zalpha = np.kron(Zalpha,Iz)

    for i in range(0,n_orb):

        Ia = np.eye(np.power(2,i))
        Za = np.eye(1)
        for k in range(0,i):
            Za = np.kron(Za,Iz)
        assert(Za.shape == Ia.shape)

        Ic = np.eye(np.power(2,n_orb-i-1))
        Zc = np.eye(1)
        for k in range(0,n_orb-i-1):
            Zc = np.kron(Zc,Iz)
        assert(Zc.shape == Ic.shape)
        
        S2_a = np.kron(Ia,np.kron(s2,Ic))
        S2_b = np.kron(Ia,np.kron(s2,Ic))

        S2 += abs(np.kron(S2_a,Ialpha) - np.kron(Ialpha,S2_b))


        for j in range(i+1,n_orb):
            
            Ib = np.eye(np.power(2,j-i-1))
            Zb = np.eye(1)
            for k in range(i+1,j):
                Zb = np.kron(Zb,Iz)

            Ic = np.eye(np.power(2,n_orb-j-1))
            Zc = np.eye(1)
            for k in range(j,n_orb-1):
                Zc = np.kron(Zc,Iz)
            
            assert(Zc.shape == Ic.shape)
            assert(Zb.shape == Ib.shape)

            #Case A
            # Ia * ap * Zb * am * Ic 
            # Ia *-am * Zb * ap * Ic
            aiaj_up = np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(am,Ic))))
            aiaj_dn = np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(ap,Ic))))


            aiaj = np.kron(aiaj_up,aiaj_dn)

            S2  +=  (-aiaj)

            #Case B
            # Ia * am * Zb *-ap * Ic
            # Ia * ap * Zb * am * Ic
            aiaj_up = np.kron(Ia,np.kron(am,np.kron(Zb,np.kron(ap,Ic))))
            aiaj_dn = np.kron(Ia,np.kron(ap,np.kron(Zb,np.kron(am,Ic))))

            aiaj = np.kron(aiaj_up,aiaj_dn)

            S2  +=  (-aiaj)

            #Case C
            aiaj_up = np.kron(Ia,np.kron(no,np.kron(Ib,np.kron(I2,Ic))))
            aiaj_dn = np.kron(Ia,np.kron(I2,np.kron(Ib,np.kron(no,Ic))))

            aiaj = np.kron(aiaj_up,aiaj_dn)

            S2  -= 0.5 * (aiaj)

            #Case D
            aiaj_up = np.kron(Ia,np.kron(I2,np.kron(Ib,np.kron(no,Ic))))
            aiaj_dn = np.kron(Ia,np.kron(no,np.kron(Ib,np.kron(I2,Ic))))

            aiaj = np.kron(aiaj_up,aiaj_dn)

            S2  -= 0.5 * (aiaj)

            #Case E
            aiaj_up = np.kron(Ia,np.kron(no,np.kron(Ib,np.kron(no,Ic))))
            aiaj_dn = np.kron(Ia,np.kron(I2,np.kron(Ib,np.kron(I2,Ic))))

            aiaj = np.kron(aiaj_up,aiaj_dn)

            S2  += 0.5 * (aiaj)

            #Case F
            aiaj_up = np.kron(Ia,np.kron(I2,np.kron(Ib,np.kron(I2,Ic))))
            aiaj_dn = np.kron(Ia,np.kron(no,np.kron(Ib,np.kron(no,Ic))))

            aiaj = np.kron(aiaj_up,aiaj_dn)
            
            S2  += 0.5 * (aiaj)

            
    return S2
    # }}}


#for alpha beta alternate ordering: use these
def form_Sz_ab(n,ap,am,no,ho,I2,Iz):
# {{{
    It = np.eye(2**n)
    Sz = np.zeros((4**n,4**n))
    for i in range(0,n):
        Ia1 = np.eye(np.power(2,2*i))
        Ic1 = np.eye(np.power(2,2*n-2*i-1))
        Ia2 = np.eye(np.power(2,2*i+1))
        Ic2 = np.eye(np.power(2,2*n-2*i-2))
        Na = np.kron(Ia1,np.kron(no,Ic1))
        Nb = np.kron(Ia2,np.kron(no,Ic2))
        Sz += Na - Nb
    return Sz
# }}}

def form_S2_ab(n_orb,ap,am,no,ho,I2,Iz):
# {{{
    S2 = np.zeros((4**n_orb,4**n_orb))
    s2 = np.array([[0,0],[0,0.75]])


    for i in range(0,n_orb):


        bfor  = 2*i
        aftr  = 2*n_orb-2*i-2

        Ia = np.eye(np.power(2,bfor))
        Ib = np.eye(np.power(2,aftr))
        a_temp = np.kron(s2,I2)
        b_temp = np.kron(I2,s2)
        S2a = np.kron(Ia,np.kron(a_temp,Ib))
        S2b = np.kron(Ia,np.kron(b_temp,Ib))

        S2 += abs(S2a -S2b)
        

        for j in range(i+1,n_orb):
            
            intr = 2*j-2*i-2 
            aftr = 2*n_orb-2*j-2

            Ib = np.eye(np.power(2,intr))
            Zb = np.eye(1)
            for k in range(2*i+2,2*j):
                Zb = np.kron(Zb,Iz)

            Ic = np.eye(np.power(2,aftr))
            Zc = np.eye(1)
            for k in range(2*j,2*n_orb-2):
                Zc = np.kron(Zc,Iz)
            
            assert(Zc.shape == Ic.shape)
            assert(Zb.shape == Ib.shape)

            
            Sptemp = np.kron(ap,am) 
            Smtemp = np.kron(am,ap) 
            ANtemp = np.kron(no,I2)
            BNtemp = np.kron(I2,no)

            ##CASE A
            aiaj = np.kron(Ia,np.kron(Sptemp,np.kron(Ib,np.kron(Smtemp,Ic))))
            S2  +=  (aiaj)

            ##CASE B
            aiaj = np.kron(Ia,np.kron(Smtemp,np.kron(Ib,np.kron(Sptemp,Ic))))
            S2  +=  (aiaj)

            ##CASE C
            aiaj = np.kron(Ia,np.kron(ANtemp,np.kron(Ib,np.kron(BNtemp,Ic))))
            S2  -= 0.5 * (aiaj)

            ##CASE D
            aiaj = np.kron(Ia,np.kron(BNtemp,np.kron(Ib,np.kron(ANtemp,Ic))))
            S2  -= 0.5 * (aiaj)

            ##CASE E
            aiaj = np.kron(Ia,np.kron(ANtemp,np.kron(Ib,np.kron(ANtemp,Ic))))
            S2  += 0.5 * (aiaj)

            ##CASE F
            aiaj = np.kron(Ia,np.kron(BNtemp,np.kron(Ib,np.kron(BNtemp,Ic))))
            S2  += 0.5 * (aiaj)

    return S2
    # }}}

def spatial_2_spin_eri_ab(n_orb,g):
# {{{
    
    g_spin = np.kron(g,np.eye(2))
    g_spin = np.kron(g_spin.T,np.eye(2))

    g_spin = g_spin.swapaxes(1,2)
    return g_spin
# }}}

def spatial_2_spin_oei_ab(n_orb,h):
# {{{
    
    h_spin = np.kron(h,np.eye(2))

    return h_spin
# }}}

