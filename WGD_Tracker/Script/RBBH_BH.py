"""
RBBH Pipeline - Best Blast Hits
Author: Morgane MILIN
Last Update: July 2024
Input: Filtered BLAST output text file

This script aims to keep a predefined number of best alignments for each query sequences
and/or for each predefined kb window.
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
folder_name, sp1, sp2, intragenomic = data_path.split('RBBH')[0], data_file.split('_vs_')[0].split('Flt_')[1], data_file.split('_vs_')[1].split('.')[0], eval(dict_param['intragenomic'])
limit = dict_param['limit'] if dict_param['limit'] != '' else 1
interval = dict_param['interval'] if dict_param['interval'] != '' else 50000
q_col, q_start, q_end, s_col, s_start, s_end = 0, 6, 7, 1, 8, 9
print(f'\ndata file: {data_path + data_file}\nq_col = {q_col}, q_start = {q_start}, q_end = {q_end}\ns_col = {s_col}, s_start = {s_start}, s_end = {s_end}')


############################################
# Column of interest & Dictionary creation #
############################################

if sp1.endswith('_genomic'):
    dico_query = fct.dict_genomic(datafile=data_path + data_file, col=q_col, start=q_start, end=q_end)
else:
    dico_query = fct.dict_CDS(datafile=data_path + data_file, col=q_col, start=q_start, end=q_end)
print('len(dico_query) =', len(dico_query))


############################
# Search for the Best Hits #
############################

outfile_BH = open(folder_name + 'RBBH/BH_' + data_file, 'w')
nb_line, kept_line, res_split_bh_4 = 0, 0, []
for key in sorted(dico_query.keys()):
    if sp1.endswith('_genomic'):
        mini, maxi = 1, int(interval)
        tempo_query_file = fct.tmp_data_creation(input_data=data_path + data_file, col=q_col, key_ref=key)
        while mini < int(dico_query[key]):
            tempo_res_sorted = fct.read_through_predefined_windows(input_data=tempo_query_file, pos_min=mini, pos_max=maxi, start=q_start, end=q_end)
            if tempo_res_sorted:
                result_bh = fct.best_hits_search(res=tempo_res_sorted, query_col=q_col, subject_col=s_col, start=s_start, end=s_end, limit=limit)
                for BH_4 in result_bh:
                    outfile_BH.write('\t'.join(BH_4) + '\n')
                    kept_line += 1
            mini += int(interval)
            maxi += int(interval)
    else:
        tempo_query_file = fct.tmp_data_creation(input_data=data_path + data_file, col=q_col, key_ref=key)
        if tempo_query_file:
            result_bh = fct.best_hits_search(res=tempo_query_file, query_col=q_col, subject_col=s_col, start=s_start, end=s_end, limit=limit)
            for BH_4 in result_bh:
                outfile_BH.write('\t'.join(BH_4) + '\n')
                kept_line += 1
outfile_BH.close()

print(f'A total of {kept_line} lines have been kept.')
