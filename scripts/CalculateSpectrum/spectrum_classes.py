import numpy as np
from numpy import linalg as la

class Atom:
	"""Atom class for GROMACS normal mode calculations"""

	def __init__(self, squared_mass, charge_mass_ratio):
		"""Initialize an Atom object"""
		self.__squared_mass      = squared_mass
		self.__charge_mass_ratio = charge_mass_ratio

	@property
	def squared_mass(self):
		"""Return the squared atomic mass of an Atom"""
		return self.__squared_mass

	@property
	def charge_mass_ratio(self):
		"""Return the ratio consisting of the atomic charge divided by the square root of the atomic mass of an Atom"""
		return self.__charge_mass_ratio

class NormalMode:
	"""NormalMode class for GROMACS normal mode calculations"""

	def __init__(self, eigenfrequency, eigenvector):
		"""Initialize a NormalMode object"""
		self.__eigenfrequency                       = eigenfrequency
		self.__eigenvector                          = eigenvector
		self.__eigenvector_converted_and_normalized = False
		self.__intensity                            = None

	@property
	def eigenfrequency(self):
		"""Return the eigenfrequency of the NormalMode"""
		return self.__eigenfrequency

	@property
	def intensity(self):
		"""Return the intensity of the NormalMode"""
		return self.__intensity

	def convert_from_cartesian_and_normalize(self, atoms):
		"""Convert from Cartesian back to mass-weighted coordinates and normalize the eigenvectors"""
		if not self.__eigenvector_converted_and_normalized:
			for i, atom in enumerate(atoms):
				self.__eigenvector[i,:] = self.__eigenvector[i,:]*atom.squared_mass
			self.__eigenvector = self.__eigenvector/la.norm(self.__eigenvector)
			self.__eigenvector_converted_and_normalized = True

	def calculate_intensity(self, atoms):
		"""Calculate intensity of a NormalMode"""
		if self.__intensity:
			raise Exception("the intensity is already calculated for this NormalMode")
		else:
			if not self.__eigenvector_converted_and_normalized:
				self.convert_from_cartesian_and_normalize(atoms)
			# get charge mass ratios
			charge_mass_ratios = np.empty([0,0])
			for atom in atoms:
				charge_mass_ratios = np.append(charge_mass_ratios, atom.charge_mass_ratio)
			# iterate over dimensions
			intensity = 0
			for k in range(3):
				intensity += (charge_mass_ratios.dot(self.__eigenvector[:,k]))**2
			self.__intensity = intensity

class Molecule:
	"""Molecule class for GROMACS normal mode calculations"""

	def __init__(self, linear, atoms, normal_modes):
		"""Initialize a Molecule object"""
		self.__atoms   = atoms
		self.__normal_modes = None
		if linear:
			self.__normal_modes = normal_modes[5:]
		else:
			self.__normal_modes = normal_modes[6:]

	@property
	def atoms(self):
		"""Return a list of Atom objects in Molecule"""
		return self.__atoms

	@property
	def normal_modes(self):
		"""Return a list of NormalMode objects in Molecule"""
		return self.__normal_modes
