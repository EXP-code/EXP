import numpy as np
import sys


'''
This script takes a file outputted by the testRatio command called testRatio_(((element_number))).ratio

In the first column of that file is temperature.
In each sucessive column is the ratio between Ion generated recombination coefficients
to the Chianti supplied one for each ionisation level. So if running this code for
hydrogen there will be 1 additional coulmn, for oxygen there is 8.

This script calculates the correction factor needed to make the ratio 1 for each
temperature and ionisation state. It saves them to files and as a  C++ readable header.
'''

def read_data(filemane):

    # Load the file
    data = np.loadtxt(filemane)
    data = data.T

    T = data[0]
    ratios = data[1:]

    return T, ratios

# Create a list to hold the correction factors for every Z, T and ionisation level
all_correction_factors = []

Z_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16]
max_Z = Z_list[-1]
for Z in range(1, max_Z + 1):

    if Z in Z_list:


        # Read in the data
        filename = 'testRatio_' + str(Z) +'.ratio'
        T, ratios = read_data(filename)

        # Get the output filename
        out_file = 'correction_factors_' + str(Z) + '.data'

        # Get and save the correction_factors as a function of temperature, and append them to the masterlist
        correction_factors  = 1/ratios.T
        T = np.reshape(T, (len(T),1))
        T_correction_factors = np.hstack((T, correction_factors))
        header = 'The first column is temperature. Each sucessive column is (for increasing ionisation level) \n the correction factor to the recombination coeficient for Z = ' + str(Z) + ' to make our code agree with chianti.'
        np.savetxt(out_file, T_correction_factors, header=header)

        # Want all the elements to have the same size arrays, so add extra 1s to pad out the rows that wouldn't exits
        padding = np.ones((200, max_Z - Z))
        correction_factors = np.hstack((correction_factors, padding))
        all_correction_factors.append(correction_factors)
    else:
        correction_factors = [1.]
        all_correction_factors.append(correction_factors)


# Save the correction factors, and if it hasn't been done yet temperatures,
# as a C++ header. Faffy formatting stuff writes it in a form C++ will recognise
correct_cpp_header = open('correction_factors.h', 'a')
T_as_string = '{'
for temp in T:
    T_as_string += str(temp[0]) +', '
T_as_string += '};\n\n'
correct_cpp_header.write('double correct_temps[' + str(len(T)) +'] = ' + T_as_string)
#ifndef INCLUDED_A_H \n #define INCLUDED_A_H \n extern
correct_as_string = '{'
for Z in range(1, max_Z + 1):
    Z_index = Z - 1 #Because python indexing starts from 0
    data_Z_string = '{'
    if Z in Z_list:
        for temp in range(len(T)):
            data_T_string = '{'
            for ionisation_level in range(max_Z):
                data_T_string += str(all_correction_factors[Z_index][temp][ionisation_level]) + ', '
            data_T_string += '}, '
            data_Z_string += data_T_string
    data_Z_string += '}, '
    correct_as_string += data_Z_string
correct_as_string += '};\n\n'
correct_cpp_header.write('double correct_facts[' + str(len(all_correction_factors)) + '][' + str(len(T)) +'][' + str(Z) + '] = ' + correct_as_string)
#correct_cpp_header.write('\n #endif') extern
