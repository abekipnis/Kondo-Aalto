import createc
import matplotlib.pyplot as plt
import numpy as np
import pdb
from skimage import morphology, measure
from numpy import empty, sqrt, square, meshgrid, linspace, dot, argmax, argmin, reshape, array
from numpy.linalg import norm, pinv, lstsq
from scipy.spatial import distance_matrix
from scipy.optimize import leastsq, least_squares, minimize
from dataclasses import dataclass
from multiprocessing import Process, Queue, Array
import multiprocessing
from sklearn.preprocessing import normalize
from math import cos, sin
import pandas as pd


# DEFINING CONSTANTS
a = 0.409 # nm, lattice constant of silver

d = np.sqrt(6)/4.*a # height of triangles
b = np.sqrt(2)*a/2. # width of triangles in nm

dpath = "/Users/akipnis/Desktop/Aalto Atomic Scale Physics/Summer 2021 Corrals Exp data/"

#corral with 14 atoms:
c1 = "Ag 2021-07-29 corral built/Createc2_210730.105015.dat"
c2 = "Ag 2021-08-10 2p5 nm radius/2p5 nm radius pm20mV line spectrum/Createc2_210810.090437.dat"
c3 = "Ag 2021-08-13 3p8 nm radius/Createc2_210813.102220.dat"
c4 = "Ag 2021-08-13 2p5 nm radius/2p5nm radius pm 20mV line spectrum/Createc2_210813.161840.dat"
c6 = 'Ag 2021-08-11/3p8nm_radius line spectra pm20mV/Createc2_210811.134245.dat'
c5 = "Ag 2021-08-13 2p5 nm radius/pm 100 mV 2p5 nm radius line spectrum/Createc2_210813.172359.dat"
# image_file = createc.DAT_IMG(dpath + )
# image_file = createc.DAT_IMG(dpath + "Ag 2021-08-10 2p5 nm radius/2p5 nm radius pm20mV line spectrum/Createc2_210810.090437.dat")

@dataclass
class Vector:
    x: float
    y: float
    norm: float
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.arr = array([x,y])
        self.norm = norm(self.arr)
        self.normed = self.arr/self.norm

    def __sub__(self, other):
        return Vector(self.x-other.x, self.y-other.y)

    def __add__(self, other):
        return Vector(self.x+other.x, self.y+other.y)

    def rot(self, th):
        rmatrix = array([[np.cos(th), -np.sin(th)],[np.sin(th), np.cos(th)]])
        return Vector(*dot(rmatrix,self.arr))

@dataclass
class CircCorralData:
    file: str
    def __init__(self, file, label):
        self.file = file
        self.label = label
        self.image_file = createc.DAT_IMG(self.file)
        # topography,
        self.im = self.image_file._crop_img(self.image_file.img_array_list[0][:][:])
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
        X = np.hstack((reshape(X1, (nxny, 1)), reshape(X2, (nxny, 1))))
        X = np.hstack((np.ones((nxny, 1)), X))
        YY = np.reshape(self.im, (nxny,1))
        theta = dot(dot(pinv(dot(X.transpose(), X)), X. transpose()), YY)
        plane = np.reshape(dot(X, theta), (self.xPix, self.yPix))
        self.im -= plane
        # return plane

    def get_region_centroids(self, diamond_size=5, sigmaclip=1.5):
        diamond = morphology.diamond(diamond_size)
        maxima = morphology.h_maxima(self.im, sigmaclip*np.std(self.im))
        r = morphology.binary_dilation(maxima, selem=diamond)
        plt.imshow(maxima)
        plt.title(self.label + "\nLocal maxima")
        xim = morphology.label(r)
        regions = measure.regionprops(xim)
        plt.figure()
        plt.imshow(xim)
        plt.show()
        plt.close()
        regions_areas = [r.area for r in regions]
        regions_area_max = max(regions_areas)

        # all regions might be same size, in which case we're spot on
        allsame = np.all([r==regions_areas[0] for r in regions_areas])

        # if we have the 'background' as a region, remove it
        # if not allsame:
        #     regions = [r for r in regions if (r.area != regions_area_max)]

        c = [list(reversed(r.centroid)) for r in regions]
        def is_far_from_edge(c,p):
            a = 0.5 #nanometers, distance from edge
            x = (p[0]>c.nm_to_pix(a) and p[0]<c.imshape[0]-self.nm_to_pix(a))
            y = (p[1]>c.nm_to_pix(a) and p[1]<c.imshape[1]-self.nm_to_pix(a))
            return (x and y)
        # remove centroids close to the edge
        d = [is_far_from_edge(self,d) for d in c]
        c = [d for d in c if is_far_from_edge(self,d)]
        self.centroids = c
        print("%d centroids" %(len(self.centroids)))

    def remove_central_atom(self, data):
        # Check two ways
        # 1:
        # get the distance matrix
        distmat = distance_matrix(data, data)

        # nearest neighbor distances for every centroid
        dists = np.ma.masked_equal(distmat,0).min(axis=1)

        # centroid w largest nearest neighbor distance is the central atom
        center_idx_1 = np.argmax(dists)

        # # create a copy since we want to save central atom

        # 2: get the atom closest to the center in a circular fit of all atoms
        r, center = self.nsphere_fit(data)
        center_idx_2 = argmin([norm(center-o) for o in data])

        if center_idx_1==center_idx_2:
            ccopy = data.copy()
            # remove outlier
            try: #two different options
                ccopy.pop(center_idx_1)
            except:
                ccopy = np.delete(ccopy,center_idx_1,axis=0)
            return ccopy
        else:
            plt.imshow(self.im)
            plt.scatter(*np.array(data).T)
            plt.show()
            raise Exception("Something went wrong removing central atom")

    def get_central_atom(self, data):
        # get the distance matrix
        distmat = distance_matrix(data, data)

        # nearest neighbor distances for every centroid
        dists = np.ma.masked_equal(distmat,0).min(axis=1)
        # print(dists)
        # centroid w largest nearest neighbor distance is the central atom
        center_idx_1 = np.argmax(dists)

        r, center = self.nsphere_fit(data)
        center_idx_2 = argmin([norm(center-o) for o in data])

        if center_idx_1==center_idx_2:
            # remove outlier
            return data[center_idx_1]
        else:
            raise Exception("Something went wrong getting central atom in corral")

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
        x = array(x)
        n = x.shape[-1]
        x = x.reshape(-1, n)
        m = x.shape[0]

        B = np.empty((m, n + 1), dtype=x.dtype)
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
        y, *_ = lstsq(B, d, rcond=None)#, overwrite_a=True, overwrite_b=True)

        c = 0.5 * y[:-1]
        r = sqrt(y[-1] + square(c).sum())

        if scaling:
            r *= scale
            c *= scale
            c += offset
        self.r, self.c = r, c
        return r, c

    def circle(self, r, c, npoints=100):
        theta = 2*np.pi*np.arange(0,1,1/npoints)
        x = c[0] + r*np.cos(theta)
        y = c[1] + r*np.sin(theta)
        return x, y

    def plot_circle_fit(self, points, radius, center, label):
        xc, yc = self.circle(radius, center, npoints=1000)
        plt.scatter(*array(points).T,label=label)
        plt.scatter(xc, yc, alpha=0.5, s=5, label=label)
        # plt.show()

    def bbox(self):
        """
        return (x0, y0, x1, y1) defined by the min and max x, y coords of atoms making up the corral
        """
        if self.occupied:
            xs, ys = array(self.remove_central_atom(self.gauss_fit_locs.T)).T
        else:
            xs, ys = array(self.gauss_fit_locs)
        return [min(xs), min(ys), max(xs), max(ys)]

    def make_lattice(self, theta, offset=0):
        """
        Given a corrals data object, theta, and lattice vector offset value in pixels,
        return an array of shape (N,2) (where N is the # of lattice sites)
        with the locations of the lattice points for that theta, lattice vector offset.
        The
        """
        theta = -theta
        rot = array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
        offset = dot([0,offset], rot)
        if not self.corral:
            origin = array(self.imshape)/2 #the middle
        elif self.occupied:
            origin = self.get_central_atom(self.gauss_fit_locs.T)+offset
        elif not self.occupied:
            origin = self.c + offset
        bbox = self.bbox()

        mult_factor = 2
        width = mult_factor*(bbox[2] - bbox[0])
        height = mult_factor*(bbox[3] - bbox [1])

        natoms_per_row = round_to_even(self.pix_to_nm(width)/b)
        nrows = round_to_even(self.pix_to_nm(height)/d)
        ls = []
        for n in np.arange(-nrows/2,nrows/2+1, 1):
            for m in np.arange(-natoms_per_row/2,natoms_per_row/2+1, 1):
                if n%2==0:
                    ls.append(array([n*self.nm_to_pix(d), m*self.nm_to_pix(b)]))
                else:
                    ls.append(array([n*self.nm_to_pix(d), (m*self.nm_to_pix(b) + self.nm_to_pix(b/2))]))
        ls = array(ls)
        rot = array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
        ls = dot(ls, rot)
        ls += array(origin)
        return ls

    def correlate_lattice_to_atom_positions(self, angle, offset):
        """
        Given a set of atom positions in self.gauss_fit_locs and
        an angle and lattice offset for the triangular lattice,
        return a value that represents the fit between the simulated lattice points and
        the fit locations of the atom. This value should be minimized.

        In this case, return a vector showing the distances between each atom and its nearest lattice site
        """
        lat = self.make_lattice(angle, offset)
        gloc = self.gauss_fit_locs
        mindists = np.min(distance_matrix(gloc.T, lat),axis=1)
        return np.mean(mindists) #np.sum()

        # p = self.gauss_fit_params.T
        # height = np.mean(p[0])
        # width = np.mean(self.gauss_fit_params.T[-2:])
        # p[0] = height
        # p[-2] = width
        # p[-1] = width
        # l = [l(*array(lat).T) for l in list(map(gaussian, *p))]
        # eval = np.sum(l, axis=0)
        #
        # def show():
        #     # lines, markers = plt.triplot(*array(lat).T, label="lattice")#, s=np.sum(l,axis=0))
        #     plt.scatter(*array(lat).T, s=np.sum(l,axis=0))
        #     # markers.set_color("black")
        #     plt.scatter(*self.gauss_fit_locs, label="gauss fit")
        #     plt.scatter(*array(self.centroids).T, label ="max fit")
        #     plt.title("Atom positions, potential lattice sites and Gaussians evaluated at lattice sites")
        #     # g2d = np.sum([gaussian(1, *cen, 23, 23)(*array(lat).T) for cen in positions], axis=0)
        #     plt.legend()
        #     plt.show()
        # # show()
        # return -np.sum(eval)/height

    def get_im_square(self, x, y, sidelen):
        # return the image around x,y with sides sidelen
        return self.im[int(y)-sidelen//2:int(y)+sidelen//2,int(x)-sidelen//2:int(x)+sidelen//2]

    def fit_lattice(self):
        """
        Given a corral data set, fit the triangular lattice to the
        pre-solved Gaussian-fitted atom positions, return the optimal angle
        and lattice offset and plot the topo map, atom positions & lattice together
        """
        def fix_self(s):
            return lambda args: s.correlate_lattice_to_atom_positions(*args)

        # maximum value of shift/offset of lattice
        m = self.nm_to_pix(b)

        # do fitting to find the best
        init = [0.5*np.pi/3,0.5*self.nm_to_pix(b)]
        bounds= ((0,0),(np.pi/3, self.nm_to_pix(b)))
        # ranges = (slice(0,np.pi/3), slice(0, self.nm_to_pix(b), self.nm_to_pix(b)/20))

        #one way to do it
        # result = basinhopping(fix_self(self), init )
        result = least_squares(fix_self(self), init, bounds= bounds,verbose=2, max_nfev=1600, method="dogbox", ftol=1e-11, xtol=1e-10)

        # def print_fun(x, f, accepted):
        #     print("at minimum, x: %1.2lf %1.2lf, %1.2lf accepted %d" %(*x, f, int(accepted)) )
        #
        # class MyBounds:
        #     def __init__(self, xmax=[np.pi/3, self.nm_to_pix(b)], xmin=[0,0] ):
        #         self.xmax = np.array(xmax)
        #         self.xmin = np.array(xmin)
        #
        #     def __call__(self, **kwargs):
        #         x = kwargs["x_new"]
        #         tmax = bool(np.all(x <= self.xmax))
        #         tmin = bool(np.all(x >= self.xmin))
        #         return tmax and tmin
        # mybounds = MyBounds()
        # minimizer_kwargs = {"method": "BFGS"}
        # from scipy.optimize import basinhopping
        # r = basinhopping(fix_self(self), init, minimizer_kwargs=minimizer_kwargs,
        #     niter=200, callback=print_fun, accept_test=mybounds,stepsize=0.05)
        # print(r.x)

        angle, offset = result.x
        print("init:", init)
        print("angle, offset:", angle, offset)

        from  scipy.stats import sigmaclip
        plt.figure(figsize=(6,6))
        new_im = self.im.copy()
        new_im[self.im>np.mean(self.im)+3*np.std(self.im)] = np.inf
        plt.imshow(new_im)
        plt.triplot(*array(self.make_lattice(angle,offset)).T)
        plt.scatter(*self.gauss_fit_locs)

        bb = self.bbox()
        pad = int(self.nm_to_pix(1))
        plt.xlim(max(0,bb[0]-pad), min(bb[2]+pad,self.imshape[0]))
        plt.ylim(max(0,bb[1]-pad),min(self.imshape[1],bb[3]+pad))
        plt.title(self.label+"\nTopography image, lattice fit, atom sites")
        plt.text(bb[0],bb[1],
            "offset: %1.3lf nm\n"
            "angle: %1.4lf degrees"
                %(self.pix_to_nm(offset),angle*180/np.pi),
                bbox={'facecolor':'w', 'alpha':0.5, 'pad':5})#,
        plt.savefig(self.label.split("/")[-1] + "topography_fit.png")
        plt.show()

    def fit_atom_pos_gauss(self, box_size):
        """
        Given a CircCorralData object with first-guess centroids,
        get square of side length box_size and fit 2D Gaussian to the atom shape
        to get a better guess for atom positions. Return a 'reconstruction'
        of the original image using the Gaussian fit parameters for every atom
        and a list of the atoms and their fit parameters
        """
        full_im = np.zeros(self.im.shape)
        fit_params = []
        for n, cen in enumerate(self.centroids):
            # Fit a Gaussian over the atom topography
            f = self.get_im_square(*cen, box_size)
            params = fitgaussian(f,self)
            fitc = gaussian(*params)
            plt.matshow(f);
            plt.contour(fitc(*np.indices(f.shape)), cmap=plt.cm.copper)

            plt.scatter([params[2]], [params[1]], )
            plt.scatter([box_size/2],[box_size/2],c="red")
            plt.title("Fitting centroid %d" %(n))
            plt.show()
            plt.close()
            # plt.show()

            # add the original 'box' back to the center
            params[1] += cen[1] - box_size/2
            params[2] += cen[0] - box_size/2
            fit_params.append(params)

            # add the Gaussian fit to the 'reconstruction' image
            full_im = full_im + gaussian(*params)(*np.indices(self.im.shape))

        # we only really care about the locations
        fp = array([array(fit_params).T[2],array(fit_params).T[1]])

        #dist (Å) 'max height' guess is off from Gaussian fit
        d = np.mean(norm(self.pix_to_nm(array(self.centroids))-self.pix_to_nm(fp.T),axis=1)*10)
        print("Max height guess different from Gaussian fit on average by: %1.2lf Å" %(d))
        self.gauss_fit_params = np.array(fit_params)
        self.gauss_fit_locs = fp
        return full_im, fp

    def compare_fits(self):
        plt.figure(figsize=(9,6))
        plt.imshow(self.im)#, extent=[0, self.pix_to_nm(self.xPix),0 ,self.pix_to_nm(self.yPix)])
        self.plot_circle_fit(self.centroids, self.r_n, self.c_n, "naive")
        self.plot_circle_fit(self.gauss_fit_locs.T, self.r_g, self.c_g, "Gaussian")
        plt.legend()
        plt.text(10,30,
            "r (naive): %1.2lf nm\n"
            "c (naive): x: %1.2lf nm, y: %1.2lf nm\n"
            "r (gauss): %1.2lf nm\n"
            "c (gauss): x: %1.2lf nm, y: %1.2lf nm"
                %(self.pix_to_nm(self.r_n), *self.pix_to_nm(self.c_n),
                self.pix_to_nm(self.r_g), *self.pix_to_nm(self.c_g)),
                bbox={'facecolor':'w', 'alpha':0.5, 'pad':5})#,
        plt.title(self.label + "\nFits to circle, naive & gaussian fit positions")
        plt.savefig(self.label.split("/")[-1] +"_circle_fits.png")
def round_to_even(n):
    # return n rounded up to the nearest even integer
    return int(np.ceil(n/2.)*2)

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments

    first guess for fit parameters
    """
    print(data.shape)
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    width_x, width_y = [data.shape[0]/2,data.shape[1]/2]
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data, c):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)

    s = data.shape
    #height, center_x, center_y, width_x, width_y
    bounds = [[np.mean(data)-3*np.std(data), 0, 0, 0, 0],[np.inf, s[0],s[1],s[0],s[1]]]
    params = [data.max(),data.shape[0]/2,data.shape[1]/2, c.nm_to_pix(.5), c.nm_to_pix(.5)]
    try:
        p = least_squares(errorfunction, params, bounds=bounds).x
    except ValueError as e:
        print(e)
        print("Something went wrong, fix it!")
        pdb.set_trace()
    return p

if __name__=="__main__":
    print("lattice spacing of 111 surface: %1.2lf Å" %(b*10))


    # Clean excel spreadsheet to read .dat file locations and relevant params
    inventory = "/Users/akipnis/Dropbox/papers-in-progress/Small Kondo corrals/Small Kondo corral inventory.xlsx"
    inv = pd.read_excel(inventory, header=2)
    scans = inv['Scan file name & path']
    inv = inv[scans.notnull()]
    scans = scans[scans.notnull()]
    scans = [s.replace('\\',"/") for s in scans]
    inv.drop(labels="Scan file name & path", axis="columns", inplace=True)
    inv['Scan file name & path'] = scans
    inv = inv.reset_index()

    for n in inv.index.values:
        s = inv.iloc[n]
        p = s["Scan file name & path"]
        c = CircCorralData(dpath + p, p)
        c.occupied = s["occupied"]=="t"
        c.corral = True
        c.subtract_plane()
        c.get_region_centroids(diamond_size=5, sigmaclip=2)

        # the box size to fit atom positions
        box_size_nm = 1.5
        box_size_pix = int(c.nm_to_pix(box_size_nm))
        while True:
            try:
                full_im, fit_params = c.fit_atom_pos_gauss(box_size=box_size_pix)
                break
            except ValueError as e:
                print(e)
                print("trying with smaller box size")
                box_size_nm-=0.1
                box_size_pix = int(c.nm_to_pix(box_size_nm))

        # if the corral is occupied, remove central atoms
        if c.occupied:
            atoms_n = c.remove_central_atom(array(c.centroids))
            atoms_g = c.remove_central_atom(c.gauss_fit_locs.T)
        else: # corral is unoccupied, don't try to remove central atom
            atoms_n = array(c.centroids)
            atoms_g = c.gauss_fit_locs.T

        # naive fit from maximum points
        c.r_n, c.c_n = c.nsphere_fit(atoms_n)
        pdb.set_trace()

        # better fit from gaussian fits to atoms
        c.r_g, c.c_g = c.nsphere_fit(atoms_g)

        c.compare_fits()
        plt.savefig(c.file+"_circle_fits.png")
        c.fit_lattice()

    ##TO DO:
    """
    - save circle sizes to excel sheet.
    - calculate mean of radii
    - put file title on top of plots
    - save figures to files
    - put the circle fit data (radius, center) to the plot
    - put the code and data to Triton
    - push like hell !
    """
    exit(0)
