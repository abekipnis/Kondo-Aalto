import datetime
import os
from find_atom_positions import CircCorralData
from read_vertfile import Spec
import numpy as np
import scipy
import scattering_model as sm
import pdb
import matplotlib.pyplot as plt
from scattering_model import ScatteringModel
class ScatteringSimulator():
	def __init__(self, linespec_dir, path):
		self.now_string = str(datetime.datetime.today())
		self.datfile = os.path.split(path)[-1].strip('.dat')

		self.linespec_dir = linespec_dir
		self.fname_head = "%s_%s" %(self.datfile, self.now_string)
		self.fname_head = self.fname_head.replace(" ","_")

		self.c = CircCorralData(path, path.split("/")[-1])
		self.c.subtract_plane()
		self.c.get_region_centroids(diamond_size=5, sigmaclip=2)

		# the box size to fit atom positions
		box_size_nm = 1
		box_size_pix = int(self.c.nm_to_pix(box_size_nm))
		self.c.fit_atom_pos_gauss(box_size=box_size_pix)
		self.c.corral = True
		self.c.occupied = True

		centroids = self.c.centroids
		gauss_fit_locs = self.c.gauss_fit_locs.T
		try:
			self.atoms_n, self.center_atom_loc = c.remove_central_atom(array(centroids))
			self.atoms_g, self.center_atom_loc = c.remove_central_atom(gauss_fit_locs)
		except:
			self.atoms_n = centroids
			self.atoms_g = gauss_fit_locs

		self.spec_files = os.listdir(linespec_dir)
		self.spec_files = sorted([f for f in self.spec_files if f[-4:] =="VERT"])
		print("Found %d VERT files in line spectrum directory" %(len(self.spec_files)))
		self.specs = [Spec(linespec_dir+f) for f in self.spec_files]
		self.xlocs = [s.XPos_nm for s in self.specs]
		self.ylocs = [s.YPos_nm for s in self.specs]


		# bias values
		# biases from Createc are in V, scattering model takes mV
		self.biases = self.specs[0].bias_mv
		self.bias_cut = self.biases>-67
		self.biases = self.biases[self.bias_cut]
		self.biases/=1000.

		# self.biases = self.biases[0::4] # divide the amount of data by 4

		self.x_nm = np.round(self.c.image_file.size[0]/10.)
		self.atoms_g_nm = self.c.pix_to_nm(self.atoms_g)

		self.xlocs = np.array(self.xlocs)-self.c.image_file.offset[0]/10.+self.x_nm/2.
		self.ylocs = np.array(self.ylocs)-self.c.image_file.offset[1]/10.

		# line spectrum points
		self.lsp = np.array([self.xlocs, self.ylocs]).T
		self.spec_dist = np.linalg.norm(self.lsp[-1]-self.lsp[0])

		self.scattering_model = ScatteringModel()

	def simulate_and_save_line_spectrum(self, n_spectra=5, n_bias=5):
		print("simulating line spectrum at %d points %d biases" %(len(self.lsp), len(self.biases)))

		ebias = int(len(self.biases)/n_bias)
		espectra = int(len(self.lsp)/n_spectra)

		biases = self.biases[::ebias]
		lsp = self.lsp[::espectra]

		self.ls = sm.line_spectrum_at_points(lsp, self.atoms_g_nm, biases)
		self.ls = np.array(self.ls)
		self.d = self.ls[~np.isnan(self.ls).any(axis=1)]

		plt.imshow(self.d, extent=[0,self.spec_dist, self.biases[0], self.biases[-1]], aspect='auto')
		plt.imshow(self.ls)
		plt.savefig("%s_line_spectrum.pdf" %(self.fname_head))
		np.save("%s_line_spectrum.npy" %(self.fname_head), self.ls)
		# plt.show()

	def fit_scattering_phase_shift(self, n_bias=5, n_spectra=5):

		ebias = int(len(self.biases)/n_bias)
		espectra = int(len(self.lsp)/n_spectra)

		biases = self.biases[::ebias]
		lsp = self.lsp[::espectra]
		spectra = [s.dIdV[self.bias_cut][::ebias] for s in self.specs]

		# function to calculate spectrum at points
		self.spec = lambda pts, d0: sm.line_spectrum_at_points(pts, self.atoms_g_nm, biases, d0)

		def resid(d0, spectra):
			ls = self.spec(lsp, d0)
			np.save("d0=%1.2lf" %(d0), ls)
			spectra_r = spectra[::espectra]
			line = lambda x, m, b: m*x+b

			# subtract linear fit from spectra
			for s in spectra_r:
				(m,b), _ = scipy.optimize.curve_fit(line, biases, s)
				s -= (biases*m+b)

			# normalize arrays by subtracting minimum and dividing by new maximum
			# https://se.mathworks.com/matlabcentral/answers/196798-how-to-normalize-values-in-a-matrix-to-be-between-0-and-1
			spectra_r -= np.min(spectra_r)
			spectra_r /= np.max(spectra_r)
			spectra_r = spectra_r.T

			ls -= np.min(ls)
			ls /= np.max(ls)
			return np.linalg.norm(ls[:,:,0]-spectra_r)

		# function to calculate spectrum and fit residual
		self.r = lambda d0: resid(d0, spectra)

		# optimizing to fit the scattering phase shift
		# minimize(function, init_vals, options)
		print("FITTING SCATTERING PHASE SHIFT from %d points %d biases" %(len(lsp), len(biases)))
		self.ret = scipy.optimize.minimize(self.r, np.pi/4, options={"disp":True})

		# create the spectra of the final parameter and show it
		self.ls = spec(self.lsp[::espectra], ret.x)
		self.d = ls[~np.isnan(ls).any(axis=1)]

		spec_dist = np.linalg.norm(lsp[-1]-lsp[0])
		plt.imshow(d, extent=[0,spec_dist, biases[0], biases[-1]], aspect='auto')
		plt.imshow(ls)
		plt.savefig("%s_line_spectrum.pdf" %(fname_head))
		np.save("%s_line_spectrum.npy" %(fname_head), ls)
