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
    alt_file = next((i for i in ''.join([data_path, data_file]).split('/') if i.startswith('pairwise')), None)
    SP1, SP2, Ka, Ks, ratio, statut = None, None, None, None, None, None
    #print('data_path =', data_path, '; data_file =', data_file)
    #print('alt_file =', alt_file)

    ######
    # Ka #
    ######
    if os.path.isfile(data_path + data_file + '2NG.dN'):
        with open(data_path + data_file + '2NG.dN') as inputfile_Ka:
            dico_data = {}
            i = 0
            for line in inputfile_Ka:
                i += 1
                dico_data[i] = line
            if "2" in dico_data[1]:
                SP1 = dico_data[2].split()[0].replace('\n', '')
                SP2 = dico_data[3].split()[0].replace('\n', '')
                Ka = dico_data[3].split()[1].replace('\n', '')
                statut = True
                pass
            else:
                print('2NG.dN only 1 seq')
                print(data_path + data_file + alt_file + '.fa')
                if os.path.isfile(data_path + data_file + alt_file + '.fa'):
                    seq_name = []
                    for name, seq in fct.buff_fas_reader(data_path + data_file + alt_file + '.fa'):
                        seq_name.append(name)
                    if len(seq_name) != 2:
                        raise ValueError(f"Problem regarding the number of sequences analysed: {seq_name}")
                    else:
                        SP1, SP2, Ka, Ks, ratio, statut = seq_name[0], seq_name[1], 'Na', 'Na', 'Na', False
    else:
        print('2NG.dN False')
        print(data_path + data_file + alt_file + '.fa')
        if os.path.isfile(data_path + data_file + alt_file + '.fa'):
            seq_name = []
            for name, seq in fct.buff_fas_reader(data_path + data_file + alt_file + '.fa'):
                seq_name.append(name)
            if len(seq_name) != 2:
                raise ValueError(f"Problem regarding the number of sequences analysed: {seq_name}")
            else:
                SP1, SP2, Ka, Ks, ratio, statut = seq_name[0], seq_name[1], 'Na', 'Na', 'Na', False

    if statut:
        ######
        # Ks #
        ######
        if os.path.isfile(data_path + data_file + '2NG.dS'):
            with open(data_path + data_file + '2NG.dS') as inputfile_Ks:
                dico_data = {}
                i = 0
                for line in inputfile_Ks:
                    i += 1
                    dico_data[i] = line
                if SP1 == dico_data[2].split()[0].replace('\n', '') and SP2 == dico_data[3].split()[0].replace('\n', ''):
                    Ks = dico_data[3].split()[1].replace('\n', '')
                else:
                    raise ValueError('Error !!!')
                pass

        #########
        # ratio #
        #########
        if os.path.isfile(data_path + data_file + '2NG.t'):
            with open(data_path + data_file + '2NG.t') as inputfile_ratio:
                dico_data = {}
                i = 0
                for line in inputfile_ratio:
                    i += 1
                    dico_data[i] = line
                if SP1 == dico_data[2].split()[0].replace('\n', '') and SP2 == dico_data[3].split()[0].replace('\n', ''):
                    ratio = dico_data[3].split()[1].replace('\n', '')
                else:
                    raise ValueError('Error !!!')
                pass
    #print(SP1, SP2, Ka, Ks, ratio)
    outfile.write('\t'.join([SP1, SP2, Ka, Ks, ratio]) + '\n')
outfile.close()
