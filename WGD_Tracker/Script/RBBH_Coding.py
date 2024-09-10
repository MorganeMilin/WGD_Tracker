"""
RBBH Pipeline - Keep sequences located in a coding region
Author: Morgane MILIN
Last Update: July 2024
Input: (1) Filtered BLAST output text file (2) gff file
"""

##################
# Modules import #
##################

import os
import sys
import RBBH_function as fct


############################
# Get Files and Parameters #
############################

dict_param = {param.split('=')[0]: param.split('=')[1] for param in sys.argv[1:]}
data_path, data_file = fct.split_path_and_file_name(dict_param['data_file'])
folder_name, sp1, sp2 = data_path.split('RBBH')[0], data_file.split('_vs_')[0].split('Flt_')[1], data_file.split('_vs_')[1].split('.')[0]
gff1 = folder_name + sp1 + '.gff' if os.path.exists(folder_name + sp1 + '.gff') else None
gff2 = folder_name + sp2 + '.gff' if os.path.exists(folder_name + sp2 + '.gff') else None
intragenomic, coding_type = eval(dict_param['intragenomic']), dict_param['coding_type']
infos1 = eval(dict_param['SP1_infos']) if dict_param['SP1_infos'] != '' else None
infos2 = eval(dict_param['SP2_infos']) if dict_param['SP2_infos'] != '' else None
if infos1 and infos2:
    sp1_limit = float(infos1[2]) if sp1.startswith(infos1[0]) else float(infos2[2]) if sp1.startswith(infos2[0]) else None
    sp1_ft_gff = infos1[1] if sp1.startswith(infos1[0]) else infos2[1] if sp1.startswith(infos2[0]) else None
    sp2_limit = float(infos1[2]) if sp2.startswith(infos1[0]) else float(infos2[2]) if sp2.startswith(infos2[0]) else None
    sp2_ft_gff = infos1[1] if sp2.startswith(infos1[0]) else infos2[1] if sp2.startswith(infos2[0]) else None
elif infos1 and not infos2:
    sp1_limit = float(infos1[2]) if sp1.startswith(infos1[0]) else None
    sp1_ft_gff = infos1[1] if sp1.startswith(infos1[0]) else None
    sp2_limit = float(infos1[2]) if sp2.startswith(infos1[0]) else None
    sp2_ft_gff = infos1[1] if sp2.startswith(infos1[0]) else None
elif infos2 and not infos1:
    sp1_limit = float(infos2[2]) if sp1.startswith(infos2[0]) else None
    sp1_ft_gff = infos2[1] if sp1.startswith(infos2[0]) else None
    sp2_limit = float(infos2[2]) if sp2.startswith(infos2[0]) else None
    sp2_ft_gff = infos2[1] if sp2.startswith(infos2[0]) else None
else:
    raise ValueError('Error !!!')
print(f"\nCoding region filtration type: {coding_type}")
print(f'data file: {data_path + data_file}\ngff1: {gff1}\ngff2: {gff2}\nsp1: {sp1}, sp2: {sp2}')
print(f'sp1_limit: {sp1_limit}, sp1_ft_gff: {sp1_ft_gff}, sp2_limit: {sp2_limit}, sp2_ft_gff: {sp2_ft_gff}')


#######################
# Dictionary Creation #
#######################

dict_coding = {}
dict_coding = fct.dico_coding_generator(dico_coding=dict_coding, gff=gff1, ft_type=sp1_ft_gff, species=sp1)
if coding_type == 'double' and not intragenomic:
    dict_coding = fct.dico_coding_generator(dico_coding=dict_coding, gff=gff2, ft_type=sp2_ft_type, species=sp2)


##########################
# Keeping Coding Regions #
##########################

input_file = open(data_path + data_file)
outfile_coding = open(folder_name + 'RBBH/cds_' + data_file, 'w')
outfile_non_coding = open(folder_name + 'RBBH/non_coding_' + data_file, 'w')
nb_line, kept_line, noCDS, status, non_coding = 0, 0, 0, None, None
for line in input_file:
    line = line.replace('\n', '').split()
    nb_line += 1

    ##########
    # Double #
    ##########
    if coding_type == 'double':

        if sp1 in dict_coding:
            if line[0] in dict_coding[sp1]:
                status, non_coding = fct.coding_check(species=sp1, col=0, col_start=6, col_end=7, limit=sp1_limit)
                if not status:
                    continue

        ################
        # Intragenomic #
        ################
        if intragenomic:
            if line[1] in dict_coding[sp1]:
                status, non_coding = fct.coding_check(species=sp1, col=1, col_start=8, col_end=9, limit=sp1_limit)

        ################
        # Intergenomic #
        ################
        else:
            if sp2 in dict_coding:
                if line[1] in dict_coding[sp2]:
                    status, non_coding = fct.coding_check(species=sp2, col=1, col_start=8, col_end=9, limit=sp2_limit)

    #########################
    # Intergenomic - simple #
    #########################
    else:
        if sp1_limit:
            sp_simple, sp_col, start, end, sp_limit = sp1, 0, 6, 7, sp1_limit
        elif sp2_limit:
            sp_simple, sp_col, start, end, sp_limit = sp2, 1, 8, 9, sp2_limit
        else:
            raise ValueError('PROBLEM !!!')

        if sp_simple in dict_coding:
            if line[sp_col] in dict_coding[sp_simple]:
                status, non_coding = fct.coding_check(species=sp_simple, col=sp_col, col_start=start, col_end=end, limit=sp_limit)

    ###########
    # Outfile #
    ###########
    if status:
        outfile_coding.write('\t'.join(line) + '\n')
        kept_line += 1
    if non_coding:
        outfile_non_coding.write('\t'.join(line) + '\n')
        noCDS += 1

input_file.close()
outfile_coding.close()
outfile_non_coding.close()
print(f'For a total of {nb_line} lines: {nb_line - kept_line} were deleted. However {kept_line} have been kept.')
print(f'A total of {noCDS} lines are not overlaping a coding region.')
