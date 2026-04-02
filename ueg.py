'''
Helper functions for the uniform electron gas

Author: Timothy C. Berkelbach
        Kyle Bystrom
'''

import numpy as np
import time
import math
import sys
from pyscf import lib

MAGIC_NUMBERS = [1, 7, 19, 27, 33, 57, 81, 93, 123, 147, 171, 179, 203, 251, 257, 305, 341, 365, 389, 437, 461, 485, 515, 587, 619, 691, 739, 751, 799, 847, 895, 925, 949, 1021, 1045, 1141, 1189, 1213, 1237, 1309, 1357, 1365, 1419, 1503, 1551, 1575, 1647, 1743, 1791, 1839, 1863, 1935, 2007, 2103, 2109, 2205, 2301, 2325, 2373, 2469, 2517, 2553, 2601, 2721, 2777, 2801, 2897, 2945, 2969, 3071, 3119, 3191, 3239, 3287, 3407, 3431, 3575, 3695, 3743, 3791, 3887, 3911, 3959, 4067, 4139, 4169, 4337, 4385, 4457, 4553, 4625, 4697, 4729, 4801, 4945, 5041, 5137, 5185, 5257, 5377, 5449, 5497, 5575, 5695, 5743, 5887, 6031, 6043, 6187, 6235, 6355, 6403, 6451, 6619, 6667, 6763, 6859, 6931, 6979, 7075, 7123, 7153, 7249, 7441, 7497, 7521, 7689, 7809, 7881, 8025, 8121, 8217, 8289, 8385, 8409, 8601, 8709, 8733, 8829, 8925, 9045, 9093, 9171, 9315, 9435, 9459, 9627, 9771, 9795, 9843, 9939, 10059, 10131, 10251, 10395, 10443, 10635, 10779, 10827, 11019, 11067, 11075, 11123, 11363, 11459, 11513, 11633, 11753, 11837, 11981, 12053, 12149, 12197, 12293, 12533, 12557, 12797, 12893, 12965, 13037, 13133, 13205, 13301, 13397, 13517, 13613, 13805, 13949, 13997, 14147, 14243, 14363, 14411, 14531, 14771, 14795, 14939, 15155, 15203, 15275, 15419, 15515, 15659, 15791, 15895, 15967, 16135, 16279, 16375, 16519, 16663, 16831, 16879, 17071, 17077, 17269, 17365, 17461, 17557, 17773, 17845, 17941, 18037, 18277, 18325, 18349, 18613, 18805, 18853, 18949, 19093, 19213, 19309, 19381, 19549, 19597, 19837, 19933, 20005, 20197, 20341, 20377, 20479, 20719, 20815, 20863, 21079, 21247, 21367, 21559, 21631, 21823, 21879, 21975, 22119, 22143, 22335, 22575, 22647, 22743, 22887, 22983, 23031, 23127, 23439, 23583, 23703, 23847, 23871, 24111, 24207, 24303, 24405, 24573, 24837, 24885, 25173, 25269, 25341, 25413, 25533, 25677, 25725, 25821, 26001, 26145, 26193, 26529, 26745, 26865, 26961, 27081, 27201, 27369, 27609, 27633, 27825, 28017, 28113, 28257, 28353, 28425, 28545, 28671, 28887, 28991, 29039, 29279, 29423, 29711, 29855, 30047, 30095, 30215, 30551, 30647, 30839, 31031, 31103, 31343, 31439, 31463, 31559, 31799, 31919, 31967, 32231, 32423, 32531, 32675, 32795, 32987, 33059, 33131, 33371, 33401, 33641, 33833, 33881, 34049, 34265, 34457, 34505, 34697, 34889, 35033, 35273, 35513, 35585, 35729, 35825, 36041, 36137, 36257, 36377, 36449, 36785, 37073, 37121, 37193, 37385, 37529, 37561, 37705, 37993, 38089, 38161, 38401, 38497, 38641, 38911, 39007, 39127, 39223, 39607, 39847, 40099, 40243, 40339, 40483, 40651, 40747, 40843, 41155, 41347, 41395, 41755, 41851, 41923, 42115, 42211, 42379, 42499, 42691, 42931, 43003, 43147, 43387, 43507, 43723, 43819, 43867, 44059, 44299, 44395, 44473, 44713, 45025, 45145, 45385, 45553, 45769, 45817, 45961, 46297, 46585, 46681, 46753, 46897, 47089, 47257, 47401, 47497, 47833, 47937, 48297, 48489, 48501, 48693, 48885, 49029, 49173, 49317, 49509, 49557, 49941, 50061, 50181, 50301, 50541, 50685, 50733, 50883, 51219, 51435, 51483, 51627, 51867, 52035, 52179, 52299, 52515, 52635, 52923, 52971, 53355, 53643, 53715, 53811, 54171, 54339, 54435, 54531, 54795, 54891, 54963, 55179, 55467, 55515, 55707, 55803, 56019, 56115, 56259, 56619, 56667, 57051, 57243, 57363, 57555, 57747, 57777, 57873, 58077, 58269, 58365, 58701, 58893, 59085, 59373, 59589, 59757, 59813, 60005, 60245, 60269, 60557, 60941, 61037, 61205, 61349, 61445, 61565, 61805, 62093, 62213, 62525, 62669, 62741, 62933, 63077, 63317, 63461, 63581, 63989, 64085, 64229, 64373, 64493, 64589, 64973, 65117, 65267, 65699, 65795, 65867, 66299, 66539, 66635, 66875, 67043, 67283, 67451, 67691, 67715, 68051, 68243, 68315, 68507, 68699, 68891, 68999, 69239, 69599, 69791, 69815, 69983, 70319, 70415, 70655, 70751, 71015, 71111, 71327, 71591, 71711, 71999, 72359, 72455, 72599, 72743, 72791, 72935, 73223, 73447, 73525, 73885, 74125, 74269, 74509, 74653, 74773, 74893, 75037, 75421, 75445, 75925, 76117, 76237, 76405, 76693, 76813, 76957, 77053, 77365, 77605, 78013, 78205, 78229, 78517, 78805, 78949, 78997, 79117, 79501, 79597, 79885, 80173, 80269, 80389, 80581, 80725, 80797, 80989, 81217, 81313, 81433, 81793, 82057, 82201, 82519, 82663, 82951, 83119, 83599, 83647, 83887, 84127, 84247, 84439, 84727, 84823, 84967, 85159, 85471, 85687, 85735, 86119, 86407, 86551, 86791, 87079, 87271, 87391, 87655, 87703, 88183, 88327, 88423, 88663, 88951, 88959, 89199, 89583, 89727, 89775, 90087, 90447, 90687, 90879, 91047, 91287, 91383, 91623, 91911, 91965, 92157, 92349, 92469, 92589, 92973, 93165, 93285, 93381, 93885, 93981, 94341, 94533, 94617, 95049, 95193, 95433, 95577, 95769, 96105, 96177, 96561, 96969, 97137, 97233, 97377, 97521, 97569, 97713, 98049, 98289, 98385, 98745, 98985, 99225, 99561, 99705, 99873, 100137, 100377, 100401, 100737, 100929, 101073, 101313, 101505, 101673, 101769, 101943]

MAGIC_NUMBERS = np.asarray(MAGIC_NUMBERS, dtype=np.int32)

def sorter_function(inarr):
    val = inarr[ 0 ]*inarr[ 0 ] + inarr[ 1 ]*inarr[ 1 ] + inarr[ 2 ]*inarr[ 2 ]
    return val

class UEG(object):
    def __init__(self, nelec, nbasis, rs, verbose=False):
        self.nelec = nelec
        # Limited to closed shell 
        if ((self.nelec % 2) != 0):
            sys.exit("Need an even number of electrons!")
        self.nocc = self.nelec // 2
        self.nbas = nbasis
        self.dim = 3
        self.rs = rs
        self._ulim = 0
        self._llim = 0
        self._rgvecs = np.zeros(0)
        self.verbose = verbose

        self._volume = self.nelec * (4.0/3.0) * np.pi * self.rs**3
        self._length = self._volume ** (1./ 3.)
        self.madelung = 2.83729747948149 / self._length

        self.vcut = False
        # use Baldereschi point or custom shift
        self.twist = False

        if isinstance(self.twist, np.ndarray):
            self.shift = self.twist
        elif isinstance(self.twist, bool):
            if self.twist:
                self.shift = np.array([0.25,0.25,0.25])
            else:
                self.shift = np.array([0.0,0.0,0.0])

        self.create_gvecs()
        self.build_eri()
        self.build_qconserv()
        if self.verbose:
            self.print_info()

    def create_gvecs(self):
        converged=False
        icur = int( math.ceil( ( 3. / ( 4. * np.pi ) * self.nbas ) ** ( 1. / 3. ) ) )
        Gvecs = []
        while not converged :
            curbas = 0
            self.gnorm = []
            if self.verbose:
                print("...finding gvecs in sphere ix = [ %4d ]" % icur)
            for ix in range(-icur, icur+1):
                ix2 = ix*ix
                for iy in range(-icur, icur+1):
                    iy2 = iy*iy
                    for iz in range(-icur, icur+1):
                        iz2 = iz*iz
                        Gvecs.append(ix)
                        Gvecs.append(iy)
                        Gvecs.append(iz)
                        self.gnorm.append( ix2 + iy2 + iz2 )
                        curbas = curbas + 1
            self.gnorm = sorted( self.gnorm )
            if icur < np.sqrt( self.gnorm[ self.nbas - 1] ):
                converged = False
                icur = icur + 1
            else:
                converged = True
        upper_limit = self.nbas - 1
        same_value = True
        while same_value:
            upper_limit = upper_limit + 1
            if( self.gnorm[ upper_limit ] - self.gnorm[ upper_limit - 1 ] > 0 ):
                same_value = False
            else:
                same_value = True

        lower_limit = self.nbas
        same_value = True
        while same_value:
            lower_limit = lower_limit - 1
            if( self.gnorm[ lower_limit ] - self.gnorm[ lower_limit - 1 ] > 0 ):
                same_value = False
            else:
                same_value = True

        self.ulim = upper_limit
        self.llim = lower_limit
        outbas = self.ulim
        self.nbas = outbas
        self.rgvecs = np.zeros( ( outbas, 3 ) )
        count = 0
        maxnorm  = self.gnorm[ outbas - 1 ]
        for i in range( 0, int(len(Gvecs)/3) ):
            ix = Gvecs[ 3 * i + 0 ]
            iy = Gvecs[ 3 * i + 1 ]
            iz = Gvecs[ 3 * i + 2 ]
            if (ix*ix + iy*iy + iz*iz <= maxnorm ):
                self.rgvecs[ count ][ 0 ] = ix
                self.rgvecs[ count ][ 1 ] = iy
                self.rgvecs[ count ][ 2 ] = iz
                count = count + 1
        zrgvecs = self.rgvecs
        zrgvecs = sorted( zrgvecs, key=sorter_function )
        self.rgvecs = np.array( zrgvecs )

    def print_k(self):
        twopidl = 2.0 * np.pi / self._length
        for i in range(0,len(self.gnorm)):
            print("K%3d : %14.8f" % (i,twopidl*np.sqrt(self.gnorm[ i ])))

    def qmin(self):
        return 2. * np.pi / self._length

    def build_eri(self):
        t1 = time.perf_counter()
        self.eri = np.zeros([self.nbas, self.nbas])
        for p in range(self.nbas):
            for q in range(self.nbas):
                gvec = self.rgvecs[p] - self.rgvecs[q]
                gvec *= self.qmin()
                gvec = np.dot(gvec, gvec)
                if np.isclose(gvec,0, atol=1e-8):
                    self.eri[p,q] = self.madelung
                else:
                    self.eri[p,q] = 4*np.pi / (self._volume*gvec)
        t2 = time.perf_counter()
        print(f"eri time: {t2-t1} sec, mem usage: {lib.current_memory()[0]/1e3} GB", flush=True)

    def build_qconserv(self):
        t1 = time.perf_counter()
        from scipy.spatial import cKDTree
        self.qconserv = np.empty((self.nbas, self.nbas, self.nbas), dtype=np.int32)
        tree = cKDTree(self.rgvecs)
        G_a = self.rgvecs[None, :, None, :]
        G_j = self.rgvecs[None, None, :, :]

        bytes_per_slice = 64 * (self.nbas ** 2)
        max_bytes = 16 * 1024**3  # 16 GB limit
        blk_size = max(1, int(0.9 * max_bytes / bytes_per_slice))

        for i0, i1 in lib.prange(0, self.nbas, blk_size):
            G_i_chunk = self.rgvecs[i0:i1, None, None, :]

            G_b_chunk = G_i_chunk + G_j - G_a
            G_b_flat = G_b_chunk.reshape(-1, 3)

            _, indices = tree.query(G_b_flat, distance_upper_bound=1e-5)

            indices[indices == self.nbas] = -1

            self.qconserv[i0:i1, :, :] = indices.reshape(i1 - i0, self.nbas, self.nbas)
        t2 = time.perf_counter()
        print(f"qconserv time: {t2-t1} sec, mem usage: {lib.current_memory()[0]/1e3} GB", flush=True)

    def print_info(self):
        print("Uniform Electron Gas Parameters")
        print(" - Dimension of System       = %14d " % self.dim)
        print(" - Madelung constant         = %14.8f " % self.madelung)
        print(" - Volume of Box             = %14.8f " % self._volume)
        print(" - Length of Box             = %14.8f " % self._length)
        print(" - Number of Electrons       = %14d " % self.nelec)
        print(" - Number of Basis Functions = %14d " % self.nbas)
        print(" - rs                        = %14d " % self.rs)
