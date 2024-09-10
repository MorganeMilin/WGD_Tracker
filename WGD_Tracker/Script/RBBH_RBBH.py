"""
RBBH Pipeline - Reciprocal Best Blast Hits
Author: Morgane MILIN
Last Update: July 2024
Input: BH output text files
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
folder_name, file_name, sp1, sp2 = dict_param['folder_name'], dict_param['file_name'], dict_param['SP1'], dict_param['SP2']
file1 = folder_name + '/RBBH/BH_' + file_name + sp1 + '_vs_' + sp2 + '.txt' if os.path.exists(folder_name + '/RBBH/BH_' + file_name + sp1 + '_vs_' + sp2 + '.txt') else None
file2 = folder_name + '/RBBH/BH_' + file_name + sp2 + '_vs_' + sp1 + '.txt' if os.path.exists(folder_name + '/RBBH/BH_' + file_name + sp2 + '_vs_' + sp1 + '.txt') else None
BH1_path, BH1_file = fct.split_path_and_file_name(file1)
BH2_path, BH2_file = fct.split_path_and_file_name(file2)
print('\nBH1_path =', BH1_path, '; BH1_file =', BH1_file)
print('BH2_path =', BH2_path, '; BH2_file =', BH2_file)
bh1_q_col, bh1_q_start, bh1_q_end, bh1_s_col, bh1_s_start, bh1_s_end = 0, 6, 7, 1, 8, 9
bh2_q_col, bh2_q_start, bh2_q_end, bh2_s_col, bh2_s_start, bh2_s_end = 0, 6, 7, 1, 8, 9
bh_col = 12
print('bh1_q_start =', bh1_q_start, ';bh1_q_end =', bh1_q_end, ';bh1_s_start =', bh1_s_start, ';bh1_s_end =', bh1_s_end)
print('bh2_q_start =', bh2_q_start, ';bh2_q_end =', bh2_q_end, ';bh2_s_start =', bh2_s_start, ';bh2_s_end =', bh2_s_end)
print('bh_col =', bh_col)


###############################################
# Creation of a dictionary with BH2 file data #
###############################################

dico_BH2_data = {}
input_file_BH2_tempo = open(BH2_path + BH2_file)
for line in input_file_BH2_tempo:
    line = line.replace('\n', '').split()
    if line[bh2_q_col] not in dico_BH2_data:
        dico_BH2_data[line[bh2_q_col]] = {}
    if line[bh2_s_col] not in dico_BH2_data[line[bh2_q_col]]:
        dico_BH2_data[line[bh2_q_col]][line[bh2_s_col]] = []
    dico_BH2_data[line[bh2_q_col]][line[bh2_s_col]].append(line)
input_file_BH2_tempo.close()
print('len(dico_BH2_data) =', len(dico_BH2_data))


################################################
# Search for Reciprocal Best Blast Hits - RBBH #
################################################

outfile_rbbh = open(folder_name + '/RBBH/RB' + BH1_file, 'w')

##########################
# Opening the first file #
##########################

kept_line, gene_list_done = 0, []
input_file_BH1 = open(BH1_path + BH1_file)
for l_b1 in input_file_BH1:
    l_b1 = l_b1.replace('\n', '').split()

    b1_q, b1_s = [int(l_b1[bh1_q_start]), int(l_b1[bh1_q_end])], [int(l_b1[bh1_s_start]), int(l_b1[bh1_s_end])]

    ###################################
    # Comparison with the second file #
    ###################################

    if l_b1[bh1_s_col] in dico_BH2_data:
        if l_b1[bh1_q_col] in dico_BH2_data[l_b1[bh1_s_col]]:
            for l_b2 in dico_BH2_data[l_b1[bh1_s_col]][l_b1[bh1_q_col]]:
                b2_q, b2_s = [int(l_b2[bh2_q_start]), int(l_b2[bh2_q_end])], [int(l_b2[bh2_s_start]), int(l_b2[bh2_s_end])]

                #############################
                # si genomique vs genomique #
                #############################
                if sp1.endswith('_genomic') and sp2.endswith('_genomic'):

                    if min(b1_q) in range(min(b2_s), max(b2_s) + 1) or min(b2_s) in range(min(b1_q), max(b1_q) + 1):
                        if min(b1_s) in range(min(b2_q), max(b2_q) + 1) or min(b2_q) in range(min(b1_s), max(b1_s) + 1):

                            ####################################################################################
                            # We extract the overlap info of the alignment between both BH files               #
                            # and we keep the orientation of the alignment for the query_id and the subject_id #
                            ####################################################################################

                            q_start, q_end, s_start, s_end = None, None, None, None

                            if b1_q[0] < b1_q[1]:
                                q_start, q_end = max(min(b1_q), min(b2_s)), min(max(b1_q), max(b2_s))
                            elif b1_q[0] > b1_q[1]:
                                q_start, q_end = min(max(b1_q), max(b2_s)), max(min(b1_q), min(b2_s))

                            if b1_s[0] < b1_s[1]:
                                s_start, s_end = max(min(b1_s), min(b2_q)), min(max(b1_s), max(b2_q))
                            elif b1_s[0] > b1_s[1]:
                                s_start, s_end = min(max(b1_s), max(b2_q)), max(min(b1_s), min(b2_q))

                            q_len = (min(max(b1_q), max(b2_s))) - (max(min(b1_q), min(b2_s)))
                            s_len = (min(max(b1_s), max(b2_q))) - (max(min(b1_s), min(b2_q)))

                            overlap_list = [str(q_len), str(q_start), str(q_end), str(s_len), str(s_start), str(s_end)]

                            ###############
                            # File output #
                            ###############

                            res = l_b1
                            res[bh_col] = max(l_b1[bh_col], l_b2[bh_col])

                            outfile_rbbh.write('\t'.join(res) + '\t' + '\t'.join(overlap_list) + '\n')
                            kept_line += 1

                #################
                # si cds vs cds #
                #################
                elif not sp1.endswith('_genomic') and not sp2.endswith('_genomic'):
                    if [l_b1[bh1_q_col], l_b1[bh1_s_col]] in gene_list_done:
                        continue
                    else:
                        gene_list_done.append([l_b1[bh1_q_col], l_b1[bh1_s_col]])

                        ###############
                        # File output #
                        ###############

                        res = l_b1
                        res[bh_col] = max(l_b1[bh_col], l_b2[bh_col])
                        outfile_rbbh.write('\t'.join(res) + '\n')
                        kept_line += 1

                # si cds vs genomique
                else:
                    pass
                    # A revoir !!!

input_file_BH1.close()
outfile_rbbh.close()

print(f'A total of {kept_line} lines have been kept.')
