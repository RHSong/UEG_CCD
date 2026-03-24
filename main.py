import numpy as np
from ueg import *
from hf import *
from ccd import *
from strucfac import *

nelec = 14
nbas = 19
rs = 1.0

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
# mp2
#emp2, emp2_d, emp2_x, t2 = mp2(my_ueg, moe)
#print("E(MP2) =", emp2, emp2_d, emp2_x)
#denom = build_denom(my_ueg.qconserv, moe, my_ueg.nocc, my_ueg.nbas)

denom = ueg_ccd.build_denom(my_ueg.qconserv, moe, my_ueg.nocc, my_ueg.nbas)
ls = 0.1

# mp2
E_corr, Ed, Ex, t2 = ccd_solver(my_ueg.qconserv, my_ueg.eri, denom, my_ueg.nocc, my_ueg.nbas, level_shift=ls,
                                DoMosaic=False, DoRing=False, DoLadder=False, DoxRing=False)
print("E(MP2) =", E_corr/nelec, Ed/nelec, Ex/nelec)

# rpa
E_corr, Ed, Ex, t2 = drccd_solver(my_ueg.qconserv, my_ueg.eri, denom, my_ueg.nocc, my_ueg.nbas, level_shift=ls)
print("E(drpa) =", E_corr/nelec, Ed/nelec, Ex/nelec)

# ccd
E_corr, Ed, Ex, t2 = ccd_solver(my_ueg.qconserv, my_ueg.eri, denom, my_ueg.nocc, my_ueg.nbas, level_shift=ls)
print("E(CCD) =", E_corr/nelec, Ed/nelec, Ex/nelec)

Lq, Sqd, Sqx = strucfac(my_ueg, t2)
print((np.array([Lq,Sqd / nelec]).T)[:20])

# (T)
from ueg_ccd import ueg_ccd_t
et = ueg_ccd_t.t_ene(my_ueg.rgvecs, my_ueg.qconserv, my_ueg.eri, t2, moe, my_ueg.nocc, my_ueg.nbas)
print("E(CCD(T)) =", et)
