""" Average Cost Functions for Horton to determine Charges for Molecular Dynamics.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import numpy as np
import h5py
import shutil

# We want to average over cost functions.
# These cost functions contain the three objects of the cost function: A, B, C
# A is a quadratic matrix (97x97), B a vector of the same dimension, and C is a constant.
# In the end, we are interested in the best-fit charges Q which are the solution to
# Q = A^-1 B
# We have multiple snapshots of a molecule, with corresponding As and Bs.
# This script averages over A and B, and outputs a HDF5 file containing these averaged objects.
# Later, these files will be used to obtain corresponding charges.

def read_h5(work_dir, timesteps):
	"""
	Import cost functions.

	Arguments:
	work_dir: the top-level directory the data is located in.
	timesteps: timestep labels of the snapshot subfolders in the work_dir.

	Returns:
	A_matrices: a dictionary of matrices, indexed by their timestep.
	B_vectors: a dictionary of vectors, indexed by their timestep.
	"""

	# Define the names and places of the HDF5 files that HORTON outputs.
	files = ['system300.cost.lnrhoref.-9.h5', 'system600.cost.lnrhoref.-9.h5']

	# The 97x97 matrices corresponding to each timestep
	A_matrices = dict()
	# The corresponding B vectors
	B_vectors = dict()

	# keep one HDF5 file as a template for writing into later
	template = False

	# Extract the values for each timestep
	print('loading data:')
	for timestep in timesteps:
		# build the path to the raw data
		sys_path = work_dir + '/' + timestep + 'ps'
		# go down into the subdir
		sys_path += '/smampppp_' + timestep + 'ps_new_constrains'
		# build the filename
		filename = 'system' + timestep + '.cost.lnrhoref.-9.h5'

		# special case:	1000ps are saved as 1ns
		if timestep == '1000':
			sys_path = work_dir + '/' + timestep + 'ps'
			# go down into the subdir
			sys_path += '/smampppp_' + timestep + 'ps_new_constrains'
			# build the filename
			filename = 'system' + timestep + '.cost.lnrhoref.-9.h5'

		path = sys_path + '/' + filename
		print(path, '... Done')
		
		# load the objects (read-only) from HDF5 file
		f = h5py.File(path, 'r')
		# Extract the A matrix
		A = np.array(f['cost']['A'])
		A_matrices[timestep] = A
		# Extract the B vector
		B = np.array(f['cost']['B'])
		B_vectors[timestep] = B

		if not template:
			shutil.copyfile(path, './average_cost.h5')
			template = True

	return A_matrices, B_vectors

def average(A_matrices, B_vectors, timesteps):
	""" 
	Average over matrices.

	Arguments:
	A_matrices: a dictionary of NxN matrices, indexed with their timestep.
	B_vectors: a dictionary of vectors with len N, indexed with their timestep.
	timesteps: a list or tuple of timesteps.

	Returns:
	A: the average of all A_matrices.
	B: the average of all B_matrices.
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

	return A, B

def export(A, B, template_path='./average_cost.h5'):
	"""
	Export&save numpy-matrices to HDF5 objects.
	
	Arguments:
	A: Averaged A-matrix.
	B: Averaged B-vector.

	Keyword Arguments:
	template_path: the path to the template file A and B are written into.
	"""
	# Open the template file
	f = h5py.File(template_path, 'r+')
	# Load the template A matrix
	A_old = f['cost/A']
	# Assign the averaged A
	A_old[...] = A
	# Do the same for the B-vectors
	B_old = f['cost/B']
	B_old[...] = B
	# Save changes
	f.close() 

	# Make sure that the changes were written into the template
	f = h5py.File(template_path, 'r')
	assert np.allclose(f['cost/A'][()], A)
	assert np.allclose(f['cost/B'][()], B)

	print('\nData has been written to {}:'.format(template_path))

def main():
	"""
	Run the script.
	"""
	# Uncomment to define an absolute working dir path
	# WORK_DIR = '/work/ws/nemo/fr_jh1130-smamp_shared-0/for_lukas'
	# The WORK_DIR is the top-level directory, all timesteps are subfolders here
	WORK_DIR = '.'
	TIMESTEPS = [str(time) for time in range(100, 1100, 100)] 
	A_matrices, B_vectors = read_h5(work_dir=WORK_DIR, timesteps=TIMESTEPS)

	print('Calculating averages')
	A, B = average(A_matrices, B_vectors, timesteps=TIMESTEPS)
	
	export(A, B)
	print('Done.')

if __name__ == '__main__':
	main()
