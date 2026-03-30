import numpy as np
def strucfac(ueg, t2):
    nvir, nocc = t2.shape[:2]
    Lq = np.zeros([nvir, nocc, 3])
    Sqd = np.zeros([nvir, nocc])
    Sqx = np.zeros([nvir, nocc])
    ene = 0
    for a in range(nvir):
        for i in range(nocc):
            Lq[a,i] = ueg.rgvecs[nocc+a] - ueg.rgvecs[i]
            for j in range(nocc):
                b = ueg.qconserv[i,nocc+a,j]
                if b >= ueg.nocc:
                    ene += 2 * t2[a,i,j] * 4 * np.pi / (ueg.qmin() * np.linalg.norm(Lq[a,i]))**2 / ueg._volume
                    Sqd[a,i] += 2 * t2[a,i,j] 
                    Sqx[a,i] -= t2[a,j,i]
    Lq = Lq.reshape((nvir*nocc, 3))
    Sqd = Sqd.reshape((nvir*nocc, 1))
    Sqx = Sqx.reshape((nvir*nocc, 1))
    uni_Lq, idx = np.unique(Lq.round(decimals=6), axis=0, return_inverse=True)
    count = np.bincount(idx)
    Sqd_sum = np.bincount(idx, weights=Sqd[:,0])
    Sqx_sum = np.bincount(idx, weights=Sqx[:,0])
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
