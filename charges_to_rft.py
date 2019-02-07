""" Transfer Charges from CSV table to .rft file
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import pandas as pd
import argparse
import shutil


def import_charges(csv_path):
	"""
	Read and filter the best-fit charges.

	Args:
	csv_path: path to the charge table
	
	Return:
	charges: pandas DataFrame with atom names, residue names, charges
	"""
	charges = pd.read_csv(csv_path)

	# select atom names, the residue and q column
	charges = charges[['atom', 'residue', 'q']]
	print('Charge table sucessfully imported.')

	return charges

def parse_rft(rtp_path, charges=None):
	"""
	Substitute the fitted charges for the original charges in the RTF file.

	Args:
	rtp_path: A string containing the original GROMACS topology file path
	charges: A pandas DataFrame containing the best fit charges

	Returns:
	fitted_rtp_text: The original topoly file, but with updated charges
	"""

	with open(rtp_path, 'r') as rtp_file:
		print('Successfully loaded topolgy file {}'.format(rtp_path))
		rtp_text = rtp_file.readlines()

		# save all possible atom names
		atom_names = charges.atom.unique()
		# and residuum names
		residuum_names = charges.residue.unique()
		
		# We will append all modified lines to sub_text
		sub_text = ''
		
		# Keep track of how often we substituted for debugging
		nr_substitutions = 0

		current_residuum = None		
		for line in rtp_text:
			# Atom names are only unique inside one residuum
			# Thus, specify which res we are currently in
			for residuum in residuum_names:
				if residuum in line:
					current_residuum = residuum
					break
			# Now, we can look up the atom name in the charge table.
			# First, select the lines with exactly one atom name
			for atom_name in atom_names:
				# Select lines with at least one atom name
				if atom_name in line[0:7]:
					second_entry = line[8:18].replace('+', '')
					second_entry = second_entry.replace('-', '').strip()
					# Select lines with no atom name in second column
					if not second_entry in atom_names:
						# Now we can substitute the charge with the fitted one
						line = substitute(line, charges, current_residuum)
						nr_substitutions += 1
						break

			# We substituted the charge, if applicable.
			# Now we can append the line to the text file
			sub_text += line

	# Make sure that all substitutions were executed
	assert len(charges.index) == nr_substitutions
	print('{} of {} charges substituted.'.format(nr_substitutions, len(charges.index)))

	return sub_text

def substitute(line, charges, current_residuum):
	"""
	Look up the charge for an atom, and write into a string.

	Args:
	line: String containing an atom name, residuum name, and original charge.
	charges: DataFrame containing atom name, residuum names, new charges.

	Returns:
	line: the original line, with the charge updated from the `charges` table.
	"""
	current_atom = line[0:7].strip()

	# look up the current atom name in the residuum we are currently in
	mask = (charges.atom == current_atom) & (charges.residue == current_residuum)
	# Exract the charge q
	new_charge = charges[mask].q.values[0]
	charge_string = '{: 1.6f}'.format(new_charge) + '  '
	modified_line = line[:24] + charge_string + line[34:]

	return modified_line

def export_rft(text, out_path):
	"""
	Write the modified text to a file.
	"""
	with open(out_path, 'w') as outfile:
		outfile.write(text)
	print('Modified topoly file written to {}'.format(out_path))

def cmd_parser():
	parser = argparse.ArgumentParser(prog='',
					 description='Transfer charges from a csv file to an .rtp file')
	parser.add_argument('-rtp',
        help='The location of the GROMACS topology file (.rtp)',
	default='./n7nh2.rtp', metavar='./n7nh2.rtp')

	parser.add_argument('-csv', metavar='./charges.csv',
        help='The location of the atomname-charge table, in the .csv format.', 
	default='./average_cost_function_check/fitted_points_charges.csv')
	
	parser.add_argument('-out', metavar='./modified.rtp',
        help='The location where the modified rtp file should be saved.',
	default='modified.rtp')

	args = parser.parse_args()

	return args.rtp, args.csv, args.out
	
	
def main():
	"""
	Run the script.
	"""
	rtp_path, csv_path, out_path = cmd_parser()
	charges = import_charges(csv_path)
	sub_text = parse_rft(rtp_path=rtp_path, charges=charges) 
	export_rft(text=sub_text, out_path=out_path)
	print('Done.')

if __name__ == '__main__':
	main()
