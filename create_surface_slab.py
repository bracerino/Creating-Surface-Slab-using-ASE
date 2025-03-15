import os
import sys
import numpy as np
from ase import io
from ase import Atoms
from ase.build import sort
from ase.constraints import FixAtoms
from ase.build import surface, cut, add_vacuum
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

#------------------- USER'S SETTINGS ---------------------
input_file = "TiN.cif" #Set the name of the input structure that will be used for the slab creation
surface_output_file_name = "Out" #Set output file name for the CIF output format
surface_output_file_name_POSCAR = "POSCAR_235" #Set output file name for the POSCAR output format

indices = (2, 3, 5)  # Set Miller index of the surface to be created
layers = 20  # Set the number of layers for the surface
vacuum = 10.0  # Set the vacuum length (in z-direction, Angstrom units) 
supercell = (1, 1, 1) # Set the length of the supercell for the created surface slab
thickness_to_fix = 4 # Thickness from the bottom atomic layer for the atoms to be set as fixed in the POSCAR output file
#---------------------------------------------------------



atoms_input = io.read(input_file)
atoms = surface(atoms_input, indices, layers, vacuum) #vacuum=None(Ang), periodic = False
atoms = atoms.repeat(supercell)

print(f"\n--------------------------------------------------------------------\n")
print(f"Succesfully created surface with Miller indices of {indices}, {layers} layers, {vacuum} Å vacuum, and supercell sizes of {supercell}\n from the input file structure {input_file}.\n\n" )



# Get all atomic positions
positions = atoms.get_positions()
symbols = atoms.get_chemical_symbols()
# Find the minimum z-coordinate (bottom-most atom)
min_z = np.min(positions[:, 2])
arg_min_z = np.argmin(positions[:, 2])
bottom_element = symbols[arg_min_z]

max_z = np.max(positions[:, 2])
arg_max_z = np.argmax(positions[:, 2])
top_element = symbols[arg_max_z]
print(f"The bottom-most element is {bottom_element} with the z-coordinate of: {min_z}")
print(f"The top-most element is {top_element} with the z-coordinate of: {max_z}")



io.write(filename=surface_output_file_name+'.cif', images=atoms, format="cif")

cif_r = io.read(surface_output_file_name+'.cif')
unique_symbols = sorted(set(cif_r.get_chemical_symbols()), key=lambda x: cif_r.get_chemical_symbols().index(x))
print(unique_symbols)
print(f"Elements in structure {unique_symbols}")
print(f"\n--------------------------------------------------------------------\n\n")
# I am using now Pymatgen to standardize the created 'cif' file for the conversion to POSCAR
#Otherwise with only using ASE, the lines in POSCAR defining the number of atoms would look like:
# Ti  N   Ti  N   Ti  N   Ti  N   Ti  N  
#   4   4   4   4   4   4   4   4   4   4
#Instead of:
# Ti  N  
#  20  20
str = Structure.from_file(surface_output_file_name+'.cif')
str = str.get_reduced_structure()
ciff = AseAtomsAdaptor.get_atoms(str)

fixed_indices = [atom.index for atom in ciff if atom.position[2] < min_z + thickness_to_fix]
constraint = FixAtoms(indices=fixed_indices)
ciff.set_constraint(constraint)

io.write(filename=surface_output_file_name_POSCAR, images=ciff, format="vasp")
io.write(surface_output_file_name_POSCAR+".png", atoms, format="png")
#Visualize the created surface structure
from ase.visualize import view
view(atoms)
