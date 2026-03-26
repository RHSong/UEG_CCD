import numpy as np
from ueg import *
from hf import *
from ccd import *
from strucfac import *
from ueg_ccd import ueg_ccd, ueg_ccd_t

nelec = 213 * 2
nbas = 500
rs = 30.0

# hf
my_ueg = UEG(nelec, nbas, rs, verbose=True)
my_ueg.bdrp = True
ehf, ex, moe, idx = hf_energy(my_ueg)
print("E(HF) =", ehf, ex)

P = np.zeros((my_ueg.nbas, my_ueg.nbas), dtype=int)
P[np.arange(my_ueg.nbas), idx] = 1
my_ueg.rgvecs = P @ my_ueg.rgvecs
my_ueg.build_eri()
my_ueg.build_qconserv()

denom = ueg_ccd.build_denom(my_ueg.qconserv, moe, my_ueg.nocc, my_ueg.nbas)
ls = 0.1

# mp2
E_corr, Ed, Ex, t2 = ccd_solver(my_ueg.qconserv, my_ueg.eri, denom, my_ueg.nocc, my_ueg.nbas, level_shift=ls,
                                DoMosaic=False, DoRing=False, DoLadder=False, DoxRing=False,)
print("E(MP2) =", E_corr / nelec, Ed / nelec, Ex / nelec)

q, sq, sqd, sqx = strucfac(my_ueg, t2)
q, avgSq = SphereAvg(q, sq)
print("===MP2 Sq===")
print((np.array([q, avgSq / nelec]).T)[:20])

# rpa
E_corr, Ed, Ex, t2 = drccd_solver(my_ueg.qconserv, my_ueg.eri, denom, my_ueg.nocc, my_ueg.nbas, level_shift=ls)
print("E(drpa) =", E_corr / nelec, Ed / nelec, Ex / nelec)

q, sq, sqd, sqx = strucfac(my_ueg, t2)
q, avgSq = SphereAvg(q, sq)
print("===dRPA Sq===")
print((np.array([q, avgSq / nelec]).T)[:20])

# ccd
E_corr, Ed, Ex, t2 = ccd_solver(my_ueg.qconserv, my_ueg.eri, denom, my_ueg.nocc, my_ueg.nbas, level_shift=ls)
print("E(CCD) =", E_corr / nelec, Ed / nelec, Ex / nelec)

q, sq, sqd, sqx = strucfac(my_ueg, t2)
q, avgSq = SphereAvg(q, sq)
print("===CCD Sq===")
print((np.array([q, avgSq / nelec]).T)[:20])

#et = ueg_ccd_t.t_ene(my_ueg.rgvecs, my_ueg.qconserv, my_ueg.eri, t2, moe, my_ueg.nocc, my_ueg.nbas)
#print("E(CCD(T)) =", et)

# bmp2
Foo, Fvv = ueg_ccd.ccfock(my_ueg.qconserv, my_ueg.eri, t2, my_ueg.nocc, my_ueg.nbas)
moe = moe + np.concatenate((Foo, Fvv))
denom = ueg_ccd.build_denom(my_ueg.qconserv, moe, my_ueg.nocc, my_ueg.nbas)
E_corr, Ed, Ex, t2 = ccd_solver(my_ueg.qconserv, my_ueg.eri, denom, my_ueg.nocc, my_ueg.nbas, level_shift=ls,
                                DoMosaic=False, DoRing=False, DoLadder=False, DoxRing=False,)
print("E(BMP2) =", E_corr / nelec, Ed / nelec, Ex / nelec)

# structfac
q, sq, sqd, sqx = strucfac(my_ueg, t2)
q, avgSq = SphereAvg(q, sq)
print("===BMP2 Sq===")
print((np.array([q, avgSq / nelec]).T)[:20])
