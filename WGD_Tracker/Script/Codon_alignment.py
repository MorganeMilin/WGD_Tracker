"""
Codon alignment Script
Author: Morgane MILIN
Last Update: October 2023
Input: locus_1_NT_file.fasta locus_1_AA_file.fasta
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
NT_path, NT_file = fct.split_path_and_file_name(dict_param['NT_file'])
AA_path, AA_file = fct.split_path_and_file_name(dict_param['AA_file'])


################
# Seq Recovery #
################

list_rm = []
NT_name, NT_seq, AA_name, AA_seq = '', '', '', ''
for name1, seq1 in fct.buff_fas_reader(NT_path + NT_file):
    NT_name, NT_seq = name1, seq1
    for name2, seq2 in fct.buff_fas_reader(AA_path + AA_file):
        AA_name, AA_seq = name2, seq2

        ##################################################
        # List creation of codon that need to be removed #
        ##################################################
        if NT_name == AA_name:
            nt_start, nt_end = 0, 3
            for i in range(len(AA_seq)):
                check_list = [AA_seq[i].count('-'), AA_seq[i].count('!'),
                              NT_seq[nt_start:nt_end].count('-'), NT_seq[nt_start:nt_end].count('!'),
                              NT_seq[nt_start:nt_end].count('TAA'), NT_seq[nt_start:nt_end].count('TAG'), NT_seq[nt_start:nt_end].count('TGA')]
                if check_list == [0, 0, 0, 0, 0, 0, 0]:
                    pass
                else:
                    if i not in list_rm:
                        list_rm.append(i)
                nt_start += 3
                nt_end += 3
print('len(list_rm) =', len(list_rm))


################
# Seq Recovery #
################

outfile = open(NT_path + 'Codon_alignment_NT.fasta', 'w')
NT_name, NT_seq, AA_name, AA_seq = '', '', '', ''
for NT_name, NT_seq in fct.buff_fas_reader(NT_path + NT_file):
    for AA_name, AA_seq in fct.buff_fas_reader(AA_path + AA_file):

        ###################
        # Codon alignment #
        ###################
        if NT_name == AA_name:
            if len(NT_seq)/3 != len(AA_seq):
                raise ValueError("len(NT_seq) / 3 != len(AA_seq) while it should be equal !!!")
            NT_seq_codon, AA_seq_codon, nt_start, nt_end = '', '', 0, 3
            for i in range(len(AA_seq)):
                if i in list_rm:
                    pass
                else:
                    NT_seq_codon += NT_seq[nt_start:nt_end]
                    AA_seq_codon += AA_seq[i]
                nt_start += 3
                nt_end += 3
            if '-' in NT_seq_codon or '!' in NT_seq_codon or NT_seq_codon == 'TAA' or NT_seq_codon == 'TAG' or NT_seq_codon == 'TGA':
                raise ValueError("There are still characters not allowed in the codon sequences such as '-', '!' or codon stop")
            outfile.write('>' + NT_name + '\n')
            outfile.write(NT_seq_codon + '\n')
outfile.close()
