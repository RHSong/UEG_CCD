import numpy as np

def hf_energy(ueg):
    ehf = 0
    ex = 0
    qvec = np.zeros(3)
    mo_occ = np.zeros(ueg.nbas, dtype=int)
    mo_occ[:ueg.nocc] = 1
    mo_occ_old = np.zeros(ueg.nbas, dtype=int)
    while np.prod(mo_occ == mo_occ_old) < 1:
        fock = np.zeros(ueg.nbas)
        for p in range(ueg.nbas):
            qvec[:] = ueg.rgvecs[p]
            if ueg.bdrp:
                qvec += np.array([0.25,0.25,0.25])
            qvec *= ueg.qmin()
            fock[p] += 0.5 * np.dot(qvec, qvec)
            for i in (np.where(mo_occ == 1)[0]):
                fock[p] -= ueg.eri[p,i]
        idx = np.argsort(fock)
        moe_idx = idx[:ueg.nocc]
        fock = fock[idx]
        mo_occ_old = np.zeros(ueg.nbas, dtype=int)
        mo_occ_old[moe_idx] = 1
        mo_occ, mo_occ_old = mo_occ_old, mo_occ
    for i in moe_idx:
        qvec[:] = ueg.rgvecs[i]
        if ueg.bdrp:
            qvec += np.array([0.25,0.25,0.25])
        qvec *= ueg.qmin()
        ehf += np.dot(qvec, qvec)
        for j in moe_idx:
            ehf -= ueg.eri[i,j]
            ex -= ueg.eri[i,j]
    ehf /= 2 * ueg.nocc
    ex /= 2 * ueg.nocc
    return ehf, ex, fock, idx

def mp2(ueg, moe):
    ed = 0
    ex = 0
    qconserv = ueg.qconserv
    nocc = ueg.nocc
    nbas = ueg.nbas
    t2 = np.zeros([nbas-nocc, nocc, nocc])
    for i in range(nocc):
        for a in range(nocc, nbas):
            viajb = ueg.eri[a,i]
            for j in range(nocc):
                b = qconserv[i,a,j]
                if b >= nocc:
                    denom = moe[i] + moe[j] - moe[a] - moe[b]
                    t2[a-nocc,i,j] += viajb / denom
                    vibja = ueg.eri[b,i]
                    ed += 2 * viajb * viajb / denom
                    ex -= viajb * vibja / denom
    e_mp2 = ed + ex
    return e_mp2, ed, ex, t2

def check_eris(ueg, moe):
    import pickle
    from pyscf import lib
    nocc = ueg.nocc
    nbas = ueg.nbas
    nvir = nbas - nocc
    qconserv = ueg.qconserv
    ovov = np.zeros([nocc,nvir,nocc,nvir])
    denom = 1e8 * np.ones([nocc,nvir,nocc,nvir])
    for i in range(nocc):
        for a in range(nocc, nbas):
            viajb = ueg.eri[a,i]
            for j in range(nocc):
                b = qconserv[i,a,j]
                if b >= nocc:
                    denom[i,a-nocc,j,b-nocc] = moe[i] + moe[j] - moe[a] - moe[b]
                    ovov[i,a-nocc,j,b-nocc] = viajb
    dov = moe[:nocc,None] - moe[nocc:]
    denom = lib.direct_sum('ia+jb->iajb', dov, dov)
    pickle.dump([denom, ovov], open( "ref.p", "wb" ))
    print(np.sum(2*ovov * ovov / denom))
    return denom, ovov
