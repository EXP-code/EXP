import pandas as pd
import sys
import json
'''
This script generates a makeIon.config file with which to run makeIon.
Supply the elements you wish to include as command line arguments.
For example to generate a makeIon.config file which includes hydrogen, helium
and oxygen the command would be:

python genMakeIonconfig.py H He O

The script uses relative abundance data from Allen, C.W., 1973, Astrophysical
Quantities (see the begining of the chapter on atoms) to calculate the typical
expected mass fractions of the elements the user choses to include. These mass
fractions are outputted to the makeIon.config which can be manually edited should
the user wish to do so.

The default particle type this file is used to generate is trace particles. This
can be changed by changing the particle_type variable in this script
'''

# Check the user has specified at least one element to include
if len(sys.argv) < 2:
      sys.exit('You must specify at least one element to include. See genMakeIonconfig.py docstring.')

# Set the particle type to generate
particle_type = 'Trace'

# Create a dataframe of data related to the relative abundances of different species
species_data = [[1, 'H', 1.0079, 12.00],
                [2, 'He', 4.0026, 10.93],
                [3, 'Li', 6.941, 1.60],
                [4, 'Be', 9.0122, 2.00],
                [5, 'B', 10.811, 3.50],
                [6, 'C', 12.0107, 8.52],
                [7, 'N', 14.0067, 7.96],
                [8, 'O', 15.9994, 8.82],
                [9, 'F', 18.9984, 4.60],
                [10, 'Ne', 20.1797, 7.92],
                [11, 'Na', 22.9897, 6.25],
                [12, 'Mg', 24.305, 7.42],
                [13, 'Al', 26.9815, 6.39],
                [14, 'Si', 28.0855, 7.52],
                [15, 'P', 30.9738,	5.52],
                [16, 'S', 32.065, 7.2],
                [17, 'Cl', 35.453, 5.6],
                [18, 'Ar', 39.948, 6.8],
                [19, 'K', 39.0983, 4.95],
                [20, 'Ca', 40.078, 6.3],
                [21, 'Sc', 44.9559,	3.22],
                [22, 'Ti', 47.867, 5.13],
                [23, 'V', 50.9415, 4.4],
                [24, 'Cr', 51.9961,	5.85],
                [25, 'Mn', 54.938, 5.4],
                [26, 'Fe', 55.845, 7.6],
                [27, 'Ni', 58.6934,	6.3],
                [28, 'Co', 58.9332,	5.1],
                [29, 'Cu', 63.546, 4.5],
                [30, 'Zn', 65.39, 4.2],
                [31, 'Ga', 69.723, 4.20],
                [32, 'Ge', 72.64, 4.80],
                [33, 'As', 74.9216,	4.20],
                [34, 'Se', 78.96, 5.10],
                [35, 'Br', 79.904, 4.50],
                [36, 'Kr', 83.8, 5.10],
                [37, 'Rb', 85.4678,	4.30],
                [38, 'Sr', 87.62, 4.79],
                [39, 'Y',  88.9059, 3.80],
                [40, 'Zr', 91.224, 4.50],
                [41, 'Nb', 92.9064,	4.00],
                [42, 'Mo', 95.94, 3.90],
                [44, 'Ru', 101.07, 3.60],
                [45, 'Rh', 102.9055, 3.20],
                [46, 'Pd', 106.42, 3.48],
                [47, 'Ag', 107.8682, 2.83],
                [48, 'Cd', 112.411,	3.80],
                [49, 'In', 114.818,	3.50],
                [50, 'Sn', 118.71, 3.60],
                [51, 'Sb', 121.76, 3.10],
                [52, 'Te', 127.6, 4.10],
                [53, 'I',  126.9045, 3.50],
                [54, 'Xe', 131.293,	4.10],
                [55, 'Cs', 132.9055, 3.20],
                [56, 'Ba', 137.327,	4.10],
                [57, 'La', 138.9055, 3.70],
                [58, 'Ce', 140.116,	3.95],
                [59, 'Pr', 140.9077, 3.55],
                [60, 'Nd', 144.24, 3.94],
                [62, 'Sm', 150.36, 3.63],
                [63, 'Eu', 151.964,	2.93],
                [64, 'Gd', 157.25, 3.28],
                [65, 'Tb', 158.9253, 2.50],
                [66, 'Dy', 162.5, 3.29],
                [67, 'Ho', 164.9303, 2.70],
                [68, 'Er', 167.259,	3.04],
                [69, 'Tm', 168.9342, 2.50],
                [70, 'Yb', 173.04, 3.40],
                [71, 'Lu', 174.967,	2.80],
                [72, 'Hf', 178.49, 3.00],
                [73, 'Ta', 180.9479, 2.60],
                [74, 'W',  183.84, 3.30],
                [75, 'Re', 186.207,	2.30],
                [76, 'Os', 190.23, 3.20],
                [77, 'Ir', 192.217,	3.10],
                [78, 'Pt', 195.078,	4.20],
                [79, 'Au', 196.9665, 2.89],
                [80, 'Hg', 200.59, 3.20],
                [81, 'Tl', 204.3833, 2.50],
                [82, 'Pb', 207.2, 4.10],
                [83, 'Bi', 208.9804, 3.00],
                [90, 'Th', 232.0381, 3.10],
                [92, 'U',  238.0289, 2.40]]

# Make the species data into a dataframe and apply headings to it
df = pd.DataFrame(species_data, columns = ['Atomic_number', 'Element', 'Atomic_mass', 'Log_abundence'])

# Use only the elements specified by the user
df = df.loc[df['Element'].isin(sys.argv)]

# Convert the log abundances into abundances
df['Abundence'] = df.apply(lambda row: 10**row.Log_abundence, axis=1)

# Get the total mass of each species in untis of atomic masses
df['Species_mass'] = df.apply(lambda row: row.Abundence * row.Atomic_mass, axis=1)

# Get the total mass of all the species
total_mass = df['Species_mass'].sum()

# Get the fraction of the total mass contributed by each species
df['Mass_fraction'] = df.apply(lambda row: row.Species_mass / total_mass, axis=1)

# Create a dictonary to hold the data to output to makeIon.config
output = {}
output['type'] =  particle_type
elements = {}

# Add each element to the elements dictionary
for element in sys.argv[1:]:

    row = df.loc[df['Element'].isin([element])]
    number = str(row['Atomic_number'].values[0])
    mass_frac = row['Mass_fraction'].values[0]
    elements[number] = {"mfrac": mass_frac}

# Place th elements data in the output dictionary
output['elements'] = elements

# Write the output to makeIon.config
file = open("makeIon.config","w")
with open('makeIon.config', 'w') as f:
    json.dump(output, f)
