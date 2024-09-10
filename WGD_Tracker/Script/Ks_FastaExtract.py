"""
Ks - Step 1 - Fasta extraction
Author: Morgane MILIN
Last Update: July 2024
"""

##################
# Modules import #
##################

import os
import sys
import Ks_functions as fct


############################
# Get Files and Parameters #
############################

dict_param = {param.split('=')[0]: param.split('=')[1] for param in sys.argv[1:]}
data_path, data_file = fct.split_path_and_file_name(dict_param['data_file'])
folder_name, sp1, sp2, intragenomic = data_path.split('RBBH')[0], data_file.split('_vs_')[0].split('Flt_')[1], data_file.split('_vs_')[1].split('.')[0], eval(dict_param['intragenomic'])
fasta1 = folder_name + sp1 + '.cds' if os.path.exists(folder_name + sp1 + '.cds') else None
fasta2 = folder_name + sp2 + '.cds' if os.path.exists(folder_name + sp2 + '.cds') else None
print(f'\ndata file: {data_path + data_file}\nsp1: {sp1}, sp2: {sp2}, intragenomic: {intragenomic}\nfasta1: {fasta1}\nfasta2: {fasta2}')


#######################
# Dictionary Creation #
#######################

tot_line, dict_cds = 0, {}
inputfile_prep = open(data_path + data_file)
for line in inputfile_prep:
    line = line.replace('\n', '').split('\t')
    tot_line += 1
    if line[0] not in dict_cds:
        dict_cds[line[0]] = True
    if line[1] not in dict_cds:
        dict_cds[line[1]] = True
inputfile_prep.close()
print(f'Total number of extracted fasta file: {tot_line}')

dict_cds = fct.dico_creation(infile=fasta1, dico_data=dict_cds)
if not intragenomic:
    dict_cds = fct.dico_creation(infile=fasta2, dico_data=dict_cds)
print('len(dict_cds) =', len(dict_cds), '\n\n')


####################
# Fasta extraction #
####################

line_count = 0
inputfile = open(data_path + data_file)
for line in inputfile:
    line = line.replace('\n', '').split('\t')
    line_count += 1
    check = len(str(tot_line)) - len(str(line_count))
    nb_file = '0' * check + str(line_count)

    outfile = open(folder_name + '/Ks/fasta/pairwise_' + nb_file + '.fa', 'w')

    if line[0] in dict_cds:
        outfile.write('>' + line[0] + '\n')
        outfile.write(dict_cds[line[0]] + '\n')
    else:
        raise ValueError(f'{line[0]} should be in dict_cds !!!')

    if line[1] in dict_cds:
        outfile.write('>' + line[1] + '\n')
        outfile.write(dict_cds[line[1]] + '\n')
    else:
        raise ValueError(f'{line[1]} should be in dict_cds !!!')

    outfile.close()
inputfile.close()
