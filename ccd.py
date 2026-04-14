import numpy as np
from ueg_ccd import ueg_ccd
from pyscf import lib

def ccd_solver(eri, denom, nocc, nao, t2=None, tol=1e-6, max_cycle=50, level_shift=0.0, 
               DoMosaic=True, DoRing=True, DoLadder=True, DoxRing=True):
    if t2 is None: t2 = np.zeros([nao-nocc, nocc, nocc])
    r2 = ueg_ccd.ccd_res(eri, denom, t2, nocc, DoMosaic, DoRing, DoLadder, DoxRing, nao)
    r2 = 2 * r2 + r2.transpose(0,2,1)
    Ed, Ex = ueg_ccd.ccd_ene(eri, t2, nocc, nao)
    E_corr = Ed + Ex
    norm_r2 = np.max(np.abs(r2))
    precond = 1/3 / (denom + level_shift)
    niter = 0
    mydiis = diis(dtype=r2.dtype)
    print("E_corr(init) =", E_corr)
    print("===UEG CCD===")
    while norm_r2 > tol:
        niter += 1
        t2 = t2 - precond * r2

        if niter >= 5:
            mydiis.update_errvec(t2, r2)
            t2 = mydiis.DoDIIS().reshape(t2.shape)

        r2 = ueg_ccd.ccd_res(eri, denom, t2, nocc, DoMosaic, DoRing, DoLadder, DoxRing, nao)
        r2 = 2 * r2 + r2.transpose(0,2,1)
        Ed, Ex = ueg_ccd.ccd_ene(eri, t2, nocc, nao)
        E_corr = Ed + Ex
        norm_r2 = np.max(np.abs(r2))
        print(f"{niter:5d} {E_corr:12.8f} {norm_r2:12.8f}", flush=True)
        if niter > max_cycle:
            break
    return E_corr, Ed, Ex, t2

def drccd_solver(eri, denom, nocc, nao, t2=None, tol=1e-6, max_cycle=50, level_shift=0.0):
    if t2 is None: t2 = np.zeros([nao-nocc, nocc, nocc])
    r2 = ueg_ccd.drccd_res(eri, denom, t2, nocc, nao)
    Ed, Ex = ueg_ccd.ccd_ene(eri, t2, nocc, nao)
    E_corr = Ed + Ex
    norm_r2 = np.max(np.abs(r2))
    precond = 1 / (denom + level_shift)
    niter = 0
    mydiis = diis(dtype=r2.dtype)
    print("E_corr(init) =", E_corr)
    print("===UEG drCCD===")
    while norm_r2 > tol:
        niter += 1
        t2 = t2 - precond * r2

        if niter >= 5:
            mydiis.update_errvec(t2, r2)
            t2 = mydiis.DoDIIS().reshape(t2.shape)

        r2 = ueg_ccd.drccd_res(eri, denom, t2, nocc, nao)
        Ed, Ex = ueg_ccd.ccd_ene(eri, t2, nocc, nao)
        E_corr = Ed + Ex
        norm_r2 = np.max(np.abs(r2))
        print(f"{niter:5d} {E_corr:12.8f} {norm_r2:12.8f}", flush=True)
        if niter > max_cycle:
            break
    return E_corr, Ed, Ex, t2

def build_denom(qconserv, moe, nocc, nao):
    denom = 1e8 * np.ones([nao-nocc, nocc, nocc])
    for i in range(nocc):
        for a in range(nocc, nao):
            for j in range(nocc):
                b = qconserv[i,a,j]
                if b >= nocc:
                    denom[a-nocc,i,j] = moe[a] + moe[b] - moe[i] - moe[j]
    return denom

class diis:
    def __init__(self, max_space=5, dtype=np.float64):
        self.max_space = max_space
        self.n_vecs = 0
        self.current_pos = 0  # Ring Buffer pointer

        self.B_matrix = np.zeros((max_space, max_space), dtype=dtype)

        self.fdiis = lib.H5TmpFile()

    def update_errvec(self, trial_vec, err_vec):
        idx = self.current_pos

        if f'T_{idx}' in self.fdiis:
            del self.fdiis[f'T_{idx}']
            del self.fdiis[f'E_{idx}']

        self.fdiis.create_dataset(f'T_{idx}', data=trial_vec.ravel())
        self.fdiis.create_dataset(f'E_{idx}', data=err_vec.ravel())

        for j in range(self.n_vecs):
            if j == idx and self.n_vecs == self.max_space:
                continue

            ej = self.fdiis[f'E_{j}'][:]
            dot_val = np.vdot(err_vec.ravel(), ej)
            self.B_matrix[idx, j] = dot_val
            self.B_matrix[j, idx] = np.conj(dot_val)

        self.B_matrix[idx, idx] = np.vdot(err_vec.ravel(), err_vec.ravel())

        if self.n_vecs < self.max_space:
            self.n_vecs += 1
        self.current_pos = (self.current_pos + 1) % self.max_space

    def DoDIIS(self):
        if self.n_vecs < 2:
            last_idx = (self.current_pos - 1) % self.max_space
            return self.fdiis[f'T_{last_idx}'][:]

        active_indices = list(range(self.n_vecs))
        B_active = self.B_matrix[np.ix_(active_indices, active_indices)]

        pulay = np.zeros((self.n_vecs + 1, self.n_vecs + 1), dtype=B_active.dtype)
        pulay[:self.n_vecs, :self.n_vecs] = B_active
        pulay[:self.n_vecs, -1] = -1.0
        pulay[-1, :self.n_vecs] = -1.0

        rhs = np.zeros(self.n_vecs + 1, dtype=B_active.dtype)
        rhs[-1] = -1.0

        c = np.linalg.solve(pulay, rhs)

        extrapolated_vec = None

        for i in range(self.n_vecs):
            scaled_T = c[i] * self.fdiis[f'T_{i}'][:]
            if extrapolated_vec is None:
                extrapolated_vec = scaled_T
            else:
                extrapolated_vec += scaled_T

        return extrapolated_vec

#class diis:
#    def __init__(self, errvec_size, trivec_size, max_space=10, dtype=np.float64):
#        self.max_space = max_space
#        self.trial_vecs = np.zeros((trivec_size, max_space), dtype=dtype)
#        self.err_vecs = np.zeros((errvec_size, max_space), dtype=dtype)
#        self.n_vecs = 0
#
#    def update_errvec(self, trial_vec, err_vec):
#        self.trial_vecs[:, 1:] = self.trial_vecs[:, :-1]
#        self.err_vecs[:, 1:] = self.err_vecs[:, :-1]
#
#        self.trial_vecs[:, 0] = trial_vec.ravel()
#        self.err_vecs[:, 0] = err_vec.ravel()
#
#        if self.n_vecs < self.max_space:
#            self.n_vecs += 1
#
#    def DoDIIS(self):
#        if self.n_vecs < 2:
#            return self.trial_vecs[:, 0]
#
#        E = self.err_vecs[:, :self.n_vecs]
#        T = self.trial_vecs[:, :self.n_vecs]
#
#        B = E.conj().T @ E
#        pulay = np.zeros((self.n_vecs + 1, self.n_vecs + 1), dtype=B.dtype)
#        pulay[:self.n_vecs, :self.n_vecs] = B
#        pulay[:self.n_vecs, -1] = -1.0
#        pulay[-1, :self.n_vecs] = -1.0
#
#        rhs = np.zeros(self.n_vecs + 1, dtype=B.dtype)
#        rhs[-1] = -1.0
#
#        c = np.linalg.solve(pulay, rhs)
#
#        return T @ c[:self.n_vecs]
