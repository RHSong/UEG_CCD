import numpy as np
import time
from ueg import *
from hf import *
from ccd import *
from strucfac import *
from ueg_ccd import qconserv, ueg_ccd, ueg_ccd_t

nelec = 57 * 2
nbas = 500
rs = 1.0

# hf
#tw = np.array([0.62, 0.77, 0.12])
my_ueg = UEG(nelec, nbas, rs, twist=True, verbose=True)
ehf, ex, moe, idx = hf_energy(my_ueg)
print("E(HF) =", ehf, ex)
print("HOMO-LUMO gap =", moe[my_ueg.nocc] - moe[my_ueg.nocc-1])

P = np.zeros((my_ueg.nbas, my_ueg.nbas), dtype=int)
P[np.arange(my_ueg.nbas), idx] = 1
my_ueg.rgvecs = P @ my_ueg.rgvecs
my_ueg.eri = None
my_ueg.qconserv = None
my_ueg.build_eri()

qconserv.init_g_map(my_ueg.rgvecs, my_ueg.nbas)
denom = ueg_ccd.build_denom(moe, my_ueg.nocc, my_ueg.nbas)
ls = 0.1

# mp2
t_start = time.perf_counter()
E_corr, Ed, Ex, t2 = ccd_solver(my_ueg.eri, denom, my_ueg.nocc, my_ueg.nbas, level_shift=ls,
                                DoMosaic=False, DoRing=False, DoLadder=False, DoxRing=False,)
print("E(MP2) =", E_corr / nelec, Ed / nelec, Ex / nelec)
t_end = time.perf_counter()
print(f"mp2 time: {t_end-t_start} sec", flush=True)

q, sq, sqd, sqx = strucfac(t2)
q, avgSq = SphereAvg(q, sq)
print("===MP2 Sq===")
print((np.array([q, avgSq / nelec]).T)[:20])

# rpa
t_start = time.perf_counter()
E_corr, Ed, Ex, t2 = drccd_solver(my_ueg.eri, denom, my_ueg.nocc, my_ueg.nbas, level_shift=ls)
print("E(drpa) =", E_corr / nelec, Ed / nelec, Ex / nelec)
t_end = time.perf_counter()
print(f"drpa time: {t_end-t_start} sec", flush=True)

q, sq, sqd, sqx = strucfac(t2)
q, avgSq = SphereAvg(q, sq)
print("===dRPA Sq===")
print((np.array([q, avgSq / nelec]).T)[:20])

# ccd
t_start = time.perf_counter()
E_corr, Ed, Ex, t2 = ccd_solver(my_ueg.eri, denom, my_ueg.nocc, my_ueg.nbas, level_shift=ls)
print("E(CCD) =", E_corr / nelec, Ed / nelec, Ex / nelec)
t_end = time.perf_counter()
print(f"ccd time: {t_end-t_start} sec", flush=True)

q, sq, sqd, sqx = strucfac(t2)
q, avgSq = SphereAvg(q, sq)
print("===CCD Sq===")
print((np.array([q, avgSq / nelec]).T)[:20])

# (T)
t_start = time.perf_counter()
et = ueg_ccd_t.t_ene(my_ueg.eri, t2, moe, my_ueg.nocc, my_ueg.nbas)
print("E(CCD(T)) =", et / nelec)
t_end = time.perf_counter()
print(f"(T) time: {t_end-t_start} sec", flush=True)

t_start = time.perf_counter()
Sq = ueg_ccd_t.structfac_t3(my_ueg.eri, t2, moe, my_ueg.nocc, my_ueg.nbas)
nvir, nocc = t2.shape[:2]
Lq = np.zeros([nvir, nocc, 3])
for a in range(nvir):
    for i in range(nocc):
        Lq[a,i] = my_ueg.rgvecs[nocc+a] - my_ueg.rgvecs[i]
Lq = Lq.reshape((-1,3))
Sq = Sq.reshape((-1,1)) / 3
q, idx = np.unique(Lq.round(decimals=6), axis=0, return_inverse=True)
sq = np.bincount(idx, weights=Sq[:,0])
q, avgSq = SphereAvg(q, sq)
print("===CCD(T) Sq===")
print((np.array([q, avgSq / nelec]).T)[:20])
t_end = time.perf_counter()
print(f"(T) strucfac time: {t_end-t_start} sec", flush=True)

# (bT)
t_start = time.perf_counter()
Foo, Fvv = ueg_ccd.ccfock(my_ueg.eri, t2, my_ueg.nocc, my_ueg.nbas)
#assert np.all(Foo < 0), "Foo raise moe_occ"
#assert np.all(Fvv > 0), "Fvv lower moe_vir"
moe = moe + np.concatenate((Foo, Fvv))
print("Renormalized gap =", np.min(moe[my_ueg.nocc:]) - np.max(moe[:my_ueg.nocc]))
et = ueg_ccd_t.t_ene(my_ueg.eri, t2, moe, my_ueg.nocc, my_ueg.nbas)
print("E(CCD(bT)) =", et / nelec)
t_end = time.perf_counter()
print(f"(bT) time: {t_end-t_start} sec", flush=True)

t_start = time.perf_counter()
Sq = ueg_ccd_t.structfac_t3(my_ueg.eri, t2, moe, my_ueg.nocc, my_ueg.nbas)
nvir, nocc = t2.shape[:2]
Lq = np.zeros([nvir, nocc, 3])
for a in range(nvir):
    for i in range(nocc):
        Lq[a,i] = my_ueg.rgvecs[nocc+a] - my_ueg.rgvecs[i]
Lq = Lq.reshape((-1,3))
Sq = Sq.reshape((-1,1)) / 3
q, idx = np.unique(Lq.round(decimals=6), axis=0, return_inverse=True)
sq = np.bincount(idx, weights=Sq[:,0])
q, avgSq = SphereAvg(q, sq)
print("===CCD(bT) Sq===")
print((np.array([q, avgSq / nelec]).T)[:20])
t_end = time.perf_counter()
print(f"(bT) strucfac time: {t_end-t_start} sec", flush=True)

## bmp2
#t_start = time.perf_counter()
#Foo, Fvv = ueg_ccd.ccfock(my_ueg.qconserv, my_ueg.eri, t2, my_ueg.nocc, my_ueg.nbas)
#moe = moe + np.concatenate((Foo, Fvv))
#denom = ueg_ccd.build_denom(my_ueg.qconserv, moe, my_ueg.nocc, my_ueg.nbas)
#E_corr, Ed, Ex, t2 = ccd_solver(my_ueg.qconserv, my_ueg.eri, denom, my_ueg.nocc, my_ueg.nbas, level_shift=ls,
#                                DoMosaic=False, DoRing=False, DoLadder=False, DoxRing=False,)
#print("E(BMP2) =", E_corr / nelec, Ed / nelec, Ex / nelec)
#t_end = time.perf_counter()
#print(f"bmp2 time: {t_end-t_start} sec")
#
#q, sq, sqd, sqx = strucfac(my_ueg, t2)
#q, avgSq = SphereAvg(q, sq)
#print("===BMP2 Sq===")
#print((np.array([q, avgSq / nelec]).T)[:20])
