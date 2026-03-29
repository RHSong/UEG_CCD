import numpy as np
from itertools import product
from pyscf.lib import einsum

def get_w(ki, kj, kk, ka, kb, kc, kconserv, eri, t2):
    '''Wijkabc intermediate as described in Scuseria paper before Pijkabc acts'''
    km = kconserv[ki, ka, kj]
    kf = kconserv[kk, kc, kj]
    nvir, nocc = t2.shape[:2]

    term1 = 0.0
    if kf >= nocc:
        term1 = t2[kc-nocc,kk,kj] * eri[kf, kb]

    term2 = 0.0
    if km < nocc and km >= 0:
        term2 = t2[kb-nocc,km,kk] * eri[kj, km]
    return term1-term2

def get_permuted_w(ki, kj, kk, ka, kb, kc, kconserv, eri, t2):
    '''Pijkabc operating on Wijkabc intermediate as described in Scuseria paper'''
    out = get_w(ki, kj, kk, ka, kb, kc, kconserv, eri, t2)
    out = out + get_w(kj, kk, ki, kb, kc, ka, kconserv, eri, t2)
    out = out + get_w(kk, ki, kj, kc, ka, kb, kconserv, eri, t2)
    out = out + get_w(ki, kk, kj, ka, kc, kb, kconserv, eri, t2)
    out = out + get_w(kk, kj, ki, kc, kb, ka, kconserv, eri, t2)
    out = out + get_w(kj, ki, kk, kb, ka, kc, kconserv, eri, t2)
    return out

def get_rw(ki, kj, kk, ka, kb, kc, kconserv, eri, t2):
    '''R operating on Wijkabc intermediate as described in Scuseria paper'''
    ret = (4. * get_permuted_w(ki,kj,kk,ka,kb,kc,kconserv,eri,t2) +
           1. * get_permuted_w(kj,kk,ki,ka,kb,kc,kconserv,eri,t2) +
           1. * get_permuted_w(kk,ki,kj,ka,kb,kc,kconserv,eri,t2) -
           2. * get_permuted_w(ki,kk,kj,ka,kb,kc,kconserv,eri,t2) -
           2. * get_permuted_w(kk,kj,ki,ka,kb,kc,kconserv,eri,t2) -
           2. * get_permuted_w(kj,ki,kk,ka,kb,kc,kconserv,eri,t2))
    return ret

def ccdt(kpts, kconserv, eri, t2, mo_energy, nocc, nao):
    UEGk = UEGMomentumTracker(kpts)
    energy_t = 0.0

    for a in range(nocc, nao):
        for b in range(nocc, a + 1):
            
            for i, j, k in product(range(nocc), repeat=3):
                
                c = UEGk.get_kconserv3(i, j, k, a, b)
                
                if c is None or c < nocc or c >= nao:
                    continue
    
                # Enforce the final permutation symmetry condition: a >= b >= c
                if not (a >= b and b >= c):
                    continue
    
                # Determine the symmetry prefactor to account for restricted loop bounds
                if a == b and b == c:
                    symm_factor = 1.0
                elif a == b or b == c:
                    symm_factor = 3.0
                else:
                    symm_factor = 6.0
    
                # Calculate the scalar energy denominator
                e_denom = (mo_energy[i] + mo_energy[j] + mo_energy[k] - 
                           mo_energy[a] - mo_energy[b] - mo_energy[c])
                
                # Skip infinitesimally small denominators to avoid division by zero
                if abs(e_denom) < 1e-12:
                    continue
    
                # Calculate the permuted W and R intermediates (both return scalars)
                pwijk = get_permuted_w(i, j, k, a, b, c, kconserv,eri, t2)
                rwijk = get_rw(i, j, k, a, b, c, kconserv,eri, t2)
                
                # Accumulate the scalar energy contribution
                energy_t += symm_factor * (pwijk * np.conj(rwijk)) / e_denom
    
    # Standard CCSD(T) scaling factor, normalized by the number of k-points (or volume)
    energy_t *= (1.0 / 3.0)
    print("E(T)=", energy_t)
    return energy_t

class UEGMomentumTracker:
    def __init__(self, kpts):
        self.kpts = kpts
        self.nkpts = kpts.shape[0]
        
        # Build a fast lookup dictionary: tuple(kx, ky, kz) -> index
        # We round and cast to int to guarantee perfectly stable dictionary keys
        self.kpt_dict = {
            tuple(np.rint(k).astype(int)): idx for idx, k in enumerate(kpts)
        }

    def get_kconserv(self):
        '''Replaces pyscf.get_kconserv()'''
        kconserv = np.full((self.nkpts, self.nkpts, self.nkpts), -1, dtype=int)
        
        # For small to medium grids, this nested loop with dict lookup 
        # is often faster and uses way less memory than 4D tensor broadcasting.
        for k in range(self.nkpts):
            for l in range(self.nkpts):
                for m in range(self.nkpts):
                    target_k = self.kpts[k] - self.kpts[l] + self.kpts[m]
                    target_tuple = tuple(np.rint(target_k).astype(int))
                    
                    # Use .get() so it safely returns -1 if out of bounds
                    kconserv[k, l, m] = self.kpt_dict.get(target_tuple, -1)
                    
        return kconserv

    def get_kconserv3(self, ki, kj, kk, ka, kb):
        '''Replaces kpts_helper.get_kconserv3() for the energy loop.
        Finds kc such that k_i + k_j + k_k = k_a + k_b + k_c
        '''
        target_k = (self.kpts[ki] + self.kpts[kj] + self.kpts[kk] 
                  - self.kpts[ka] - self.kpts[kb])
        
        target_tuple = tuple(np.rint(target_k).astype(int))
        
        # Returns the index of kc, or -1 if kc is outside the basis
        return self.kpt_dict.get(target_tuple, -1)
