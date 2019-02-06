""" Average Cost Functions for Horton to determine Charges for Molecular Dynamics 
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import numpy as np
import h5py

# We want to average over cost functions.
# These cost functions contain the three objects of the cost function: A, B, C
# A is a quadratic matrix (97x97), B a vector of the same dimension, and C is a constant.
# In the end, we are interested in the best-fit charges Q which are the solution to
# Q = A^-1 B
# We have multiple snapshots of a molecule, with corresponding As and Bs.
# This script averages over A and B, and outputs a HDF5 file containing these averaged objects.
# Later, these files will be used to obtain corresponding charges.

def read_h5(timesteps):
	"""
	Import cost functions.
	"""
	# Define the names and places of the HDF5 files that HORTON outputs.
	files = ['system300.cost.lnrhoref.-9.h5', 'system600.cost.lnrhoref.-9.h5']

	# The 97x97 matrices corresponding to each timestep
	A_matrices = dict()
	# The corresponding B vectors
	B_vectors = dict()


	# Extract the values for each timestep
	for timestep in timesteps:
		# build the filename
		filename = 'system' + timestep + '.cost.lnrhoref.-9.h5'
		# load the objects (read-only) from HDF5 file
		f = h5py.File(filename, 'r')
		# Extract the A matrix
		A = np.array(f['cost']['A'])
		A_matrices[timestep] = A
		# Extract the B vector
		B = np.array(f['cost']['B'])
		B_vectors[timestep] = B

	return(A_matrices, B_vectors)

def average(A_matrices, B_vectors, timesteps):
	"""
	Average over matrices
	"""

	# Initialize empty
	time = list(B_vectors.keys())[0]
	A = A_matrices[time] * 0
	B = B_vectors[time] * 0

	# Average by adding all objects and dividing by their number
	for timestep in timesteps:
		A += A_matrices[timestep]
		B += B_vectors[timestep]

	# Divide
	number_snapshots = len(timesteps)
	A /= number_snapshots
	B /= number_snapshots

	return(A, B)


if __name__ == '__main__':

	timesteps = ('300', '600')
	A_matrices, B_vectors = read_h5(timesteps)
	A, B = average(A_matrices, B_vectors, timesteps)
	print(A[:10, :10], B[:10])
