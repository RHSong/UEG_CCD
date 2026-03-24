import numpy as np
def strucfac(ueg, t2):
    nvir, nocc = t2.shape[:2]
    Lq = np.zeros([nvir, nocc])
    Sqd = np.zeros([nvir, nocc])
    Sqx = np.zeros([nvir, nocc])
    ene = 0
    for a in range(nvir):
        for i in range(nocc):
            Lq[a,i] = np.linalg.norm(ueg.rgvecs[nocc+a] - ueg.rgvecs[i])
            for j in range(nocc):
                b = ueg.qconserv[i,nocc+a,j]
                if b >= ueg.nocc:
                    ene += 2 * t2[a,i,j] * 4 * np.pi / (ueg.qmin() * Lq[a,i])**2 / ueg._volume
                    Sqd[a,i] += 2 * t2[a,i,j] 
                    Sqx[a,i] -= t2[a,j,i]
    Lq = Lq.reshape(nvir*nocc)
    Sqd = Sqd.reshape(nvir*nocc)
    Sqx = Sqx.reshape(nvir*nocc)
    uni_Lq, idx = np.unique(Lq.round(decimals=6), axis=0, return_inverse=True)
    count = np.bincount(idx)
    Sqd_sum = np.bincount(idx, weights=Sqd)
    Sqx_sum = np.bincount(idx, weights=Sqx)
    #uni_Lq *= ueg.qmin()
    print("test:", ene)
    return uni_Lq, Sqd_sum/count, Sqx_sum/count
