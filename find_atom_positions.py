import createc
import matplotlib.pyplot as plt
import numpy as np
import pdb
from skimage import morphology, measure
from numpy import empty, sqrt, square
from scipy.linalg import lstsq
from scipy.spatial import distance_matrix
from scipy.stats import pearsonr
from dataclasses import dataclass
from multiprocessing import Process, Queue, Array
import multiprocessing
from sklearn.preprocessing import normalize
from math import cos, sin



# DEFINING CONSTANTS
a = 0.409 # nm, lattice constant of silver

d = np.sqrt(6)/4*a # height of triangles
b = np.sqrt(2)*a/2 # width of triangles in nm

dpath = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/"

#corral with 14 atoms:
c1 = "Ag 2021-07-29 corral built/Createc2_210730.105015.dat"
c2 = "Ag 2021-08-10 2p5 nm radius/2p5 nm radius pm20mV line spectrum/Createc2_210810.090437.dat"

# image_file = createc.DAT_IMG(dpath + )
image_file = createc.DAT_IMG(dpath + "Ag 2021-08-10 2p5 nm radius/2p5 nm radius pm20mV line spectrum/Createc2_210810.090437.dat")

def worker(c, v, i, queue, lattice_sites, theta_offset, bbox):
    """
    we only care about the lattice sites that are within bounding box
    defined by atoms arranged in the circle
    """
    x0, y0, x1, y1 = bbox
    xn = v[0] + c.nm_to_pix(np.sqrt(2)*a/2*np.cos(i*2*np.pi/6+theta_offset))
    yn = v[1] + c.nm_to_pix(np.sqrt(2)*a/2*np.sin(i*2*np.pi/6+theta_offset))
    #if xn > 0 and yn > 0 and xn < c.xPix and yn < c.yPix:
    if xn > x0 and yn > y0 and xn < x1 and yn < y1:
        new_site = [xn, yn]
        # in lattice? (i.e. has the site been visited? )
        il = np.any([np.allclose(new_site, l, atol=1e-3) for l in lattice_sites])
        if not il:
            lattice_sites.append(new_site)
            queue.put(new_site) # explore next

@dataclass
class Vector:
    x: float
    y: float
    norm: float
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.arr = np.array([x,y])
        self.norm = np.linalg.norm(self.arr)
        self.normed = self.arr/self.norm

    def __sub__(self, other):
        return Vector(self.x-other.x, self.y-other.y)

    def __add__(self, other):
        return Vector(self.x+other.x, self.y+other.y)

    def rot(self, th):
        rmatrix = np.array([[np.cos(th), -np.sin(th)],[np.sin(th), np.cos(th)]])
        return Vector(*np.dot(rmatrix,self.arr))

@dataclass
class CircCorralData:
    file: str
    def __init__(self, file):
        self.file = file
        self.image_file = createc.DAT_IMG(file)
        # topography,
        self.im = image_file._crop_img(image_file.img_array_list[0][:][:])
        self.imshape = self.im.shape
        self.xPix = self.imshape[0]
        self.yPix = self.imshape[1]
        self.ang_ppx_x = self.image_file.nom_size.x / self.image_file.xPixel
        self.ang_ppx_y = self.image_file.nom_size.y / self.image_file.yPixel
        # nchannels = image_file.channels

    def nm_to_pix(self, nm):
        scale = self.ang_ppx_x# if direction=="x" else self.ang_ppx_y
        return nm/scale*10

    def pix_to_nm(self, pix):
        scale = self.ang_ppx_x #if direction=="x" else self.ang_ppx_y
        return pix*scale/10

    def subtract_plane(self):
        X1, X2 = np.mgrid[:self.xPix, :self.yPix]
        nxny = self.xPix*self.yPix
        X = np.hstack( (np.reshape(X1, (nxny, 1)), np.reshape(X2, (nxny, 1))))
        X = np.hstack((np.ones((nxny, 1)), X))
        YY = np.reshape(self.im, (nxny,1))
        theta = np.dot(np.dot(np.linalg.pinv(np.dot(X.transpose(), X)), X. transpose()), YY)
        plane = np.reshape(np.dot(X, theta), (self.xPix, self.yPix))
        self.im -= plane
        # return plane

    def get_region_centroids(self, diamond_size=4):
        diamond = morphology.diamond(diamond_size)
        # max = morphology.local_maxima(im, connectivity=50)
        maxima = morphology.h_maxima(self.im, np.std(self.im))
        r = morphology.binary_dilation(maxima, selem=diamond)
        plt.imshow(maxima)
        plt.show()
        xim = morphology.label(r)
        # plt.imshow(xim)
        regions = measure.regionprops(xim)
        regions_areas = [r.area for r in regions]
        regions_area_max = max(regions_areas)

        # all regions might be same size, in which case we're spot on
        allsame = np.all([r==regions_areas[0] for r in regions_areas])

        # if we have the 'background' as a region, remove it
        if not allsame:
            regions = [r for r in regions if (r.area != regions_area_max)]
        self.centroids = [list(reversed(r.centroid)) for r in regions]
        # return centroids

    def remove_central_atom(self):
        # this way does not work for the smallest corrals
        # get the distance matrix
        # distmat = distance_matrix(self.centroids, self.centroids)
        #
        # # nearest neighbor distances for every centroid
        # dists = np.ma.masked_equal(distmat,0).min(axis=1)
        #
        # # centroid w largest nearest neighbor distance is the central atom
        # center_idx = np.argmax(dists)
        #
        # # create a copy since we want to save central atom

        r, center = self.nsphere_fit(self.centroids)
        center_idx = np.argmin([np.linalg.norm(center-o) for o in self.centroids])
        ccopy = self.centroids.copy()

        # remove outlier
        ccopy.pop(center_idx)
        return ccopy

    def get_central_atom(self):
        # get the distance matrix
        # distmat = distance_matrix(self.centroids, self.centroids)
        #
        # # nearest neighbor distances for every centroid
        # dists = np.ma.masked_equal(distmat,0).min(axis=1)
        # print(dists)
        # # centroid w largest nearest neighbor distance is the central atom
        # center_idx = np.argmax(dists)
        r, center = self.nsphere_fit(self.centroids)
        center_idx = np.argmin([np.linalg.norm(center-o) for o in self.centroids])
        # remove outlier

        return self.centroids[center_idx]

    def nsphere_fit(self, x, axis=-1, scaling=False):
        r"""
        Fit an n-sphere to ND data.
        The center and radius of the n-sphere are optimized using the Coope
        method. The sphere is described by
        .. math::
           \left \lVert \vec{x} - \vec{c} \right \rVert_2 = r
        Parameters
        ----------
        x : array-like
            The n-vectors describing the data. Usually this will be a nxm
            array containing m n-dimensional data points.
        axis : int
            The axis that determines the number of dimensions of the
            n-sphere. All other axes are effectively raveled to obtain an
            ``(m, n)`` array.
        scaling : bool
            If `True`, scale and offset the data to a bounding box of -1 to
            +1 during computations for numerical stability. Default is
            `False`.
        Return
        ------
        r : scalar
            The optimal radius of the best-fit n-sphere for `x`.
        c : array
            An array of size `x.shape[axis]` with the optimized center of
            the best-fit n-sphere.
        References
        ----------
        - [Coope]_ "\ :ref:`ref-cfblanls`\ "
        """
        # x = preprocess(x, float=True, axis=axis)
        x = np.array(x)
        n = x.shape[-1]
        x = x.reshape(-1, n)
        m = x.shape[0]

        B = empty((m, n + 1), dtype=x.dtype)
        X = B[:, :-1]
        X[:] = x
        B[:, -1] = 1

        if scaling:
            xmin = X.min()
            xmax = X.max()
            scale = 0.5 * (xmax - xmin)
            offset = 0.5 * (xmax + xmin)
            X -= offset
            X /= scale

        d = square(X).sum(axis=-1)

        y, *_ = lstsq(B, d, overwrite_a=True, overwrite_b=True)

        c = 0.5 * y[:-1]
        r = sqrt(y[-1] + square(c).sum())

        if scaling:
            r *= scale
            c *= scale
            c += offset
        self.r, self.c = r, c
        return r, c

    def circle(self, npoints=100):
        theta = 2*np.pi*np.arange(0,1,1/npoints)
        x = self.c[0] + self.r*np.cos(theta)
        y = self.c[1] + self.r*np.sin(theta)
        return x, y

    def plot_circle_fit(self):
        xc, yc = self.circle(100)
        # plt.figure(figsize=(8,8))
        plt.scatter(*np.array(self.centroids).T)
        plt.scatter(xc, yc, alpha=0.5)
        plt.show()

    def bbox(self):
        """
        return (x0, y0, x1, y1) defined by the min and max x, y coords of atoms making up the corral
        """
        xs, ys = np.array(self.remove_central_atom()).T
        return [min(xs), min(ys), max(xs), max(ys)]

    def make_lattice(self, theta):
        theta = -theta
        origin = self.get_central_atom()
        bbox = self.bbox()

        width = 1.1*(bbox[2] - bbox[0])
        height = 1.1*(bbox[3] - bbox [1])

        natoms_per_row = int(np.ceil(self.pix_to_nm(width)/b/2)*2)
        #make this an even # so things are easier
        nrows = int(np.ceil(self.pix_to_nm(height)/d/2.)*2)
        offset = [width/2, height/2]

        ls = []
        for n in np.arange(-nrows/2,nrows/2+1, 1):
            for m in np.arange(-natoms_per_row/2,natoms_per_row/2+1, 1):
                if n%2==0:
                    ls.append(np.array([n*self.nm_to_pix(d), m*self.nm_to_pix(b)]))
                else:
                    ls.append(np.array([n*self.nm_to_pix(d), (m*self.nm_to_pix(b) + self.nm_to_pix(b/2))]))
        rot = np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
        # plt.triplot(*np.array(ls).T)

        ls = np.dot(ls, rot)
        # ls += np.dot(rot, offset)
        ls += np.array(origin)

        # plt.triplot(*np.array(ls).T)
        return ls

    # def get_lattice_pts_over_corral(self, theta_offset):
    #     theta_offset*=-1
    #     pdb.set_trace()
    #
    #     x1 = self.nm_to_pix(np.sqrt(2)*a/2*np.cos(0*2*np.pi/6-theta_offset))
    #     y1 = self.nm_to_pix(np.sqrt(2)*a/2*np.sin(0*2*np.pi/6-theta_offset))
    #     # x2 = self.nm_to_pix(np.sqrt(2)*a/2*np.cos(1*2*np.pi/6+theta_offset))
    #     # y2 = self.nm_to_pix(np.sqrt(2)*a/2*np.sin(1*2*np.pi/6+theta_offset))
    #
    #     # first lattice vector
    #     a1 = Vector(x1,y1)
    #
    #     # atom farthest left
    #     l = sorted(self.remove_central_atom(), key=lambda x: x[0])[0]
    #     others = self.remove_central_atom()
    #     others.remove(l)
    #
    #     # largest dot product, giving "width" of lattice at this theta_offset
    #     o = [np.dot(a1.normed, (Vector(*o)-Vector(*l)).normed) for o in others]
    #     wmax = Vector(*others[np.argmax(o)])
    #     wmin = Vector(*others[np.argmin(o)])
    #
    #     width = np.dot(a1.normed, wmax.arr-l)
    #     a1perp = a1.rot(np.pi/2)
    #
    #     # largest & smallest dot product, giving 'height' for generating lattice
    #     perps = [np.dot((Vector(*o)-Vector(*l)).arr, a1perp.arr) for o in others]
    #     hmax = Vector(*others[np.argmax(perps)])
    #     hmin = Vector(*others[np.argmin(perps)])
    #     hma = np.dot(a1perp.normed, (hmax-Vector(*l)).arr)
    #     hmi = np.dot(a1perp.normed, (hmin-Vector(*l)).arr)
    #
    #     """
    #     - generate lattice over square defined by width, hmi and hma
    #     - rotate the lattice by the defined angle
    #     """
    #     ls = []
    #     natoms_per_row = int(np.ceil(c.pix_to_nm(width)/b)) + 2
    #     nrows = int(np.ceil(c.pix_to_nm(hma-hmi)/d)) + 3
    #     for n in range(nrows):
    #         for m in range(natoms_per_row):
    #             if n%2==0:
    #                 ls.append([n*self.nm_to_pix(d), m*self.nm_to_pix(b)])
    #             else:
    #                 ls.append([n*self.nm_to_pix(d), (m*self.nm_to_pix(b) + self.nm_to_pix(b/2))])
    #
    #     rot = np.array([[cos(theta_offset), -sin(theta_offset)], [sin(theta_offset), cos(theta_offset)]])
    #     # plt.triplot(*np.array(ls).T)
    #     ret = np.array(np.dot(ls, rot))
    #     ret += l
    #     ret += [(-hmi)*sin(-theta_offset),-(-hmi)*cos(-theta_offset)]
    #     # ret += l-[min([r[0] for r in ret])]
    #     # ret += [l[0], hma+(hma-hmi)/2]#hma+3]
    #     # pdb.set_trace()
    #     plt.triplot(*ret.T)
    #     # plt.show()
    #     return ret
    def correlate_lattice_to_atom_positions(self, angle, spread):
        lat = self.make_lattice(angle)
        g2d = np.sum([gauss2d(lat, cen, 0.4) for cen in self.centroids], axis=0)
        return sum(g2d)

def gauss2d(grid, mu, sigma):
    # pdb.set_trace()
    return 1/(sigma*np.sqrt(2*np.pi))**2*np.exp(-(np.sum((grid-mu)**2,axis=-1)/(2*sigma)**2))

if __name__=="__main__":
    c = CircCorralData(dpath + c1   )
    c.subtract_plane()
    c.get_region_centroids()

    radius, center = c.nsphere_fit(c.remove_central_atom())
    angs = np.pi/3*np.arange(0,1.01,0.001)
    corrs = [c.correlate_lattice_to_atom_positions(ang,2) for ang in angs]
    plt.plot(corrs)
    plt.show()

    plt.figure(figsize=(6,6))
    plt.imshow(c.im)
    plt.triplot(*np.array(c.make_lattice(angs[np.argmax(corrs)])).T)
    c.plot_circle_fit()
    ##TO DO:
    """
    - add another layer of refinement to fitting atom positions:
        - get region around local maximum & fit 2D gaussian to that region
    - look at Gaussians and see where it starts to resemble scan
    - implement shift into fitting

    """
    exit(0)


    # breadth first search to create the lattice grid around an atom position

    # theta offset
    to = 0

    origin = c.centroids[0]
    # lattice_sites.append(origin)
    # queue.append(origin)

    pool = multiprocessing.Pool(processes=6)
    m = multiprocessing.Manager()

    # lattice sites
    ls = m.list()
    ls.append(origin)
    q = m.Queue()
    q.put(origin)

    bbox = c.bbox()
    while not q.empty():
        v = q.get()
        workers = [pool.apply_async(worker, (c, v, i, q, ls, to, bbox)) for i in range(6)]
        [w.get() for w in workers]
    plt.triplot(*np.array(ls).T)
    plt.show()
    exit(0)


"""
We want to only generate the lattice over the sites where we need
to optimize the computation to solve for the best lattice fit over the sites
"""
##################################




a1perp.normed


s(Vector(*w)-Vector(*l))
####################################







# atom farthest from the left atom
r = c.centroids[np.argmax(np.linalg.norm(np.array(c.centroids)-np.array(l), axis=1))]
r = Vector(*r)

n= Vector(*[xn, yn])
(r-n).normed

np.linalg.norm([xn, yn])

    # exit(0)
# breadth first search to create the lattice grid around an atom position
theta_offset = 0
lattice_sites = []
queue = []
origin = c.centroids[0]
lattice_sites.append(origin)
queue.append(origin)
#plt.scatter(*np.array([[v[0] + np.sqrt(2)*a/2*np.cos(i*2*np.pi/6+theta_offset),v[1] + np.sqrt(2)*a/2*np.sin(i*2*np.pi/6+theta_offset)] for i in range(6)]).T)
c.yPix

np.any(np.allclose(lattice_sites[0], new_site))
new_site
while len(queue) !=0:
    v = queue.pop(0)
    for i in range(6): # each site has 6 nearest neighbors
        xn = v[0] + c.nm_to_pix(np.sqrt(2)*a/2*np.cos(i*2*np.pi/6+theta_offset))
        yn = v[1] + c.nm_to_pix(np.sqrt(2)*a/2*np.sin(i*2*np.pi/6+theta_offset))
        if xn > 0 and yn > 0 and xn < c.xPix and yn < c.yPix:
            new_site = [xn, yn]

            # in lattice?
            il = np.any([np.allclose(new_site, l, atol=1e-3) for l in lattice_sites])

            #in queue?
            # if len(queue)==0:
            #     iq = False
            # else:
            # iq = np.any([np.allclose(new_site, l, atol=1e-3) for l in queue])
            if not il:
                lattice_sites.append(new_site)
                queue.append(new_site)


lattice_sites
plt.triplot(*np.array(lattice_sites).T)
[np.any(np.isclose(c, lattice_sites)) for c in centroids]
queue
new_site
origin

natoms_per_row = int(ang_ppx_x * xPixels/(b*10))
nrows = int(ang_ppx_y * yPixels/(d*10))
for n in range(nrows):
    for m in range(natoms_per_row):
        if n%2==0:
            lattice_sites.append([n*d/ang_ppx_y*10, m*b/ang_ppx_x*10])
        else:
            lattice_sites.append([n*d/ang_ppx_y*10, (m*b + b/2)/ang_ppx_x*10])





ang_ppx_x
plt.figure(figsize=(10,10))
lines, markers = plt.triplot(*np.array(lattice_sites[:]).T)
lines.set_alpha(0.3)
lines.set_color("green")
plt.imshow(im)
plt.scatter(*np.array(centroids).T)


plt.imshow()
plt.imshow(morphology.diameter_closing(im,diameter_threshold=16, connectivity=10))
plt.imshow(cv2.dilate(cv2.erode(im, kernel, iterations=6), kernel, iterations=6))

plt.imshow()
help(image_file)

def generate_lattice(image_shape, lattice_vectors) :
    center_pix = np.array(image_shape) // 2
    # Get the lower limit on the cell size.
    dx_cell = max(abs(lattice_vectors[0][0]), abs(lattice_vectors[1][0]))
    dy_cell = max(abs(lattice_vectors[0][1]), abs(lattice_vectors[1][1]))
    # Get an over estimate of how many cells across and up.
    nx = image_shape[0]//dx_cell
    ny = image_shape[1]//dy_cell
    # Generate a square lattice, with too many points.
    # Here I generate a factor of 4 more points than I need, which ensures
    # coverage for highly sheared lattices.  If your lattice is not highly
    # sheared, than you can generate fewer points.
    x_sq = np.arange(-nx, nx, dtype=float)
    y_sq = np.arange(-ny, nx, dtype=float)
    x_sq.shape = x_sq.shape + (1,)
    y_sq.shape = (1,) + y_sq.shape
    # Now shear the whole thing using the lattice vectors
    x_lattice = lattice_vectors[0][0]*x_sq + lattice_vectors[1][0]*y_sq
    y_lattice = lattice_vectors[0][1]*x_sq + lattice_vectors[1][1]*y_sq
    # Trim to fit in box.
    mask = ((x_lattice < image_shape[0]/2.0)
             & (x_lattice > -image_shape[0]/2.0))
    mask = mask & ((y_lattice < image_shape[1]/2.0)
                    & (y_lattice > -image_shape[1]/2.0))
    x_lattice = x_lattice[mask]
    y_lattice = y_lattice[mask]
    # Translate to the centre pix.
    x_lattice += center_pix[0]
    y_lattice += center_pix[1]
    # Make output compatible with original version.
    out = np.empty((len(x_lattice), 2), dtype=float)
    out[:, 0] = y_lattice
    out[:, 1] = x_lattice
    return out


x1 = c.nm_to_pix(np.sqrt(2)*a/2*np.cos(0*2*np.pi/6+theta_offset))
y1 = c.nm_to_pix(np.sqrt(2)*a/2*np.sin(0*2*np.pi/6+theta_offset))

x2 = c.nm_to_pix(np.sqrt(2)*a/2*np.cos(1*2*np.pi/6+theta_offset))
y2 = c.nm_to_pix(np.sqrt(2)*a/2*np.sin(1*2*np.pi/6+theta_offset))

plt.scatter([0, x1, x2],[0, y1, y2])

plt.scatter(*np.array(generate_lattice((120, 120), ((x1,y1),(x1,y2)))).T)
