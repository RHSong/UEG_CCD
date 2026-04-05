import numpy as np
from ueg_ccd import ueg_ccd
def strucfac(t2):
    nvir, nocc = t2.shape[:2]
    nao = nvir + nocc
    Sqd ,Sqx,Lq = ueg_ccd.strucfac_t2(t2,nocc,nao)
    Lq = Lq.reshape((nvir*nocc, 3))
    Sqd = Sqd.reshape((nvir*nocc, 1))
    Sqx = Sqx.reshape((nvir*nocc, 1))
    uni_Lq, idx = np.unique(Lq.round(decimals=6), axis=0, return_inverse=True)
    Sqd_sum = np.bincount(idx, weights=Sqd[:,0])
    Sqx_sum = np.bincount(idx, weights=Sqx[:,0])
    L2 = np.linalg.norm(uni_Lq, axis=1)
    L2 = L2 * L2
    return uni_Lq, (Sqd_sum + Sqx_sum), Sqd_sum, Sqx_sum

def SphereAvg(Lq, Sq):
    from scipy.interpolate import RBFInterpolator, LinearNDInterpolator
    interp = LinearNDInterpolator(Lq, Sq.real)
    fs_grid = fibonacci_sphere()
    ngrid = len(fs_grid)
    L2 = np.linalg.norm(Lq, axis=1)
    uni_L2, idx = np.unique(L2.round(decimals=10), axis=0, return_inverse=True)
    r = uni_L2[:len(Lq)//2]
    avgSq = np.zeros(len(r))
    for i, rad in enumerate(r):
        Grid = fs_grid * rad
        sq = 0
        for (x, y, z) in Grid:
            sq += interp(x,y,z)
        avgSq[i] = sq / ngrid
    return r, avgSq

def fibonacci_sphere(samples=36):
    points = []
    phi = np.pi * (np.sqrt(5.) - 1.)

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2
        radius = np.sqrt(1 - y * y)
        theta = phi * i
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        points.append((x, y, z))
    return np.array(points)

def checkene(ueg, Sq, Lq):
    L2 = ueg.qmin() * np.linalg.norm(Lq, axis=1)
    L2 = L2 * L2
    ene = np.sum(Sq * 4 * np.pi / L2 / ueg._volume)
    return ene
