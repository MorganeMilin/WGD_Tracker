"""
Script that retrieves Ka, Ks and ratio information for each codeml output file (NG)
Author: Morgane MILIN
Last Update: October 2023
Input: All output of codeml analysis
"""

##################
# Modules import #
##################

import os
import sys
import Ks_functions as fct


##################
# Ks Compilation #
##################

outfile = open(sys.argv[1], 'w')
for param in sys.argv[2:]:
    data_path, data_file = fct.split_path_and_file_name(param)
    #print('data_path =', data_path, '; data_file =', data_file)

    ######
    # Ka #
    ######
    inputfile_Ka = open(data_path + data_file + '2NG.dN')
    dico_data = {}
    i = 0
    for line in inputfile_Ka:
        i += 1
        dico_data[i] = line
    SP1 = dico_data[2].split()[0].replace('\n', '')
    SP2 = dico_data[3].split()[0].replace('\n', '')
    Ka = dico_data[3].split()[1].replace('\n', '')
    inputfile_Ka.close()

    ######
    # Ks #
    ######
    inputfile_Ks = open(data_path + data_file + '2NG.dS')
    dico_data = {}
    i = 0
    for line in inputfile_Ks:
        i += 1
        dico_data[i] = line
    if SP1 == dico_data[2].split()[0].replace('\n', '') and SP2 == dico_data[3].split()[0].replace('\n', ''):
        Ks = dico_data[3].split()[1].replace('\n', '')
    else:
        raise ValueError('Error !!!')
    inputfile_Ks.close()

    #########
    # ratio #
    #########
    inputfile_ratio = open(data_path + data_file + '2NG.t')
    dico_data = {}
    i = 0
    for line in inputfile_ratio:
        i += 1
        dico_data[i] = line
    if SP1 == dico_data[2].split()[0].replace('\n', '') and SP2 == dico_data[3].split()[0].replace('\n', ''):
        ratio = dico_data[3].split()[1].replace('\n', '')
    else:
        raise ValueError('Error !!!')
    inputfile_ratio.close()
    outfile.write('\t'.join([SP1, SP2, Ka, Ks, ratio]) + '\n')
outfile.close()
