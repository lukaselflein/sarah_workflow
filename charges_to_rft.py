""" Transfer Charges from CSV table to .rft file
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import pandas as pd
import shutil

csv_filename = './average_cost_function_check/fitted_points_charges.csv'
charges = pd.read_csv(csv_filename)

print(charges)

def main():
	"""
	Run the script.
	"""

if __name__ == '__main__':
	main()
