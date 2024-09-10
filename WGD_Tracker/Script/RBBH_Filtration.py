"""
RBBH Pipeline - Step 1 - param Filtration
Author: Morgane MILIN
Last Update: July 2024
Input: blast output text file
This script aims to filter blast output based on sequences identity, alignment length, sequences aligned coverage and, if intragenomic analysis by removing alignment against itself.
"""

##################
# Modules import #
##################

import os
import sys
from statistics import mean, median
import RBBH_function as fct


############################
# Get Files and Parameters #
############################

dict_param = {param.split('=')[0]: param.split('=')[1] for param in sys.argv[1:]}
folder_name, blast_file = fct.split_path_and_file_name(dict_param['blast_file'])
sp1, sp2, intragenomic = blast_file.split('_vs_')[0], blast_file.split('_vs_')[1].split('.')[0], eval(dict_param['intragenomic'])
fasta1 = folder_name + sp1 + '.cds' if os.path.exists(folder_name + sp1 + '.cds') else folder_name + sp1 + '.fasta' if os.path.exists(folder_name + sp1 + '.fasta') else None
fasta2 = folder_name + sp2 + '.cds' if os.path.exists(folder_name + sp2 + '.cds') else folder_name + sp2 + '.fasta' if os.path.exists(folder_name + sp2 + '.fasta') else None
identity = int(dict_param['identity']) if dict_param['identity'] != '' else 70
len_align = int(dict_param['len_align']) if dict_param['len_align'] != '' else 60
len_ratio = eval(dict_param['len_ratio']) if dict_param['len_ratio'] != '' else None
corr_intra = eval(dict_param['corr_intra']) if dict_param['corr_intra'] != '' else None
print(f"\nIdentity = {identity}\nLength alignment = {len_align}")
if len_ratio: print(f"Ratio alignment = {len_ratio}")
if intragenomic and corr_intra: print(f"Intragenomic = {intragenomic} & Correction of sequences name = {corr_intra}")
elif intragenomic: print(f"Intragenomic = {intragenomic}")
print(f'blast file: {folder_name + blast_file}\nfasta1: {fasta1}\nfasta2: {fasta2}\nsp1: {sp1}, sp2: {sp2}')


#######################
# Dictionary Creation # {name_data: data_length, etc}
#######################

dict_data = {sp1: {}} if intragenomic else {sp1: {}, sp2: {}}
if len_ratio:
    data_list = [(sp1, fasta1), (sp2, fasta2)] if not intragenomic else [(sp1, fasta1)]
    for dataset in data_list:
        for name, seq in fct.buff_fas_reader(dataset[1]):
            if name not in dict_data[dataset[0]]:
                dict_data[dataset[0]][name] = len(seq)
            else:
                raise ValueError('Having the data name twice is not possible, suggesting a potential issue with the input file or the settings.')
    print('len(dict_data_sp1) =', len(dict_data[sp1]), '; len(dict_data_sp2) =', len(dict_data[sp2]))


####################
# Param Filtration #
####################

inputfile = open(folder_name + blast_file)
outfile = open(folder_name + 'RBBH/pFlt_' + blast_file.replace('.blast', '.txt'), 'w')
nb_line, intra_rm, default_rm, ratio_rm, kept_line, data_kept, data_rm = 0, 0, 0, 0, 0, [], []
for line in inputfile:
    line = line.replace('\n', '').split()
    nb_line += 1

    ################
    # Intragenomic # alignments on itself
    ################
    if intragenomic:
        if corr_intra:
            if corr_intra[0].join(line[0].split(corr_intra[0])[:corr_intra[1]]) == corr_intra[0].join(line[1].split(corr_intra[0])[:corr_intra[1]]):
                intra_rm += 1
                continue
        else:
            if sp1.endswith('_genomic'):
                if line[0] == line[1]:
                    sp1_start, sp1_end = min(int(line[6]), int(line[7])), max(int(line[6]), int(line[7]))
                    sp2_start, sp2_end = min(int(line[8]), int(line[9])), max(int(line[8]), int(line[9]))
                    if sp1_start in range(sp2_start, sp2_end) or sp2_start in range(sp1_start, sp1_end):
                        intra_rm += 1
                        continue
            else:
                if line[0] == line[1]:
                    intra_rm += 1
                    continue

    ######################
    # Default parameters # identity, length
    ######################
    if float(line[2]) < identity or int(line[3]) < len_align:
        default_rm += 1
        continue
    else:
        ###################
        # align == simple #
        ###################
        if len_ratio:
            if len_ratio[1] == 'simple':

                ######################
                # Colomn of interest #
                ######################
                col = 0 if not sp1.endswith('_genomic') else 1 if not sp2.endswith('_genomic') else None
                start = 6 if not sp1.endswith('_genomic') else 8 if not sp2.endswith('_genomic') else None
                end = 7 if not sp1.endswith('_genomic') else 9 if not sp2.endswith('_genomic') else None
                species = sp1 if not sp1.endswith('_genomic') else sp2 if not sp2.endswith('_genomic') else None

                if line[col] in dict_data[species]:
                    blast_min, blast_max = min(int(line[start]), int(line[end])), max(int(line[start]), int(line[end]))
                    blast_len = blast_max - blast_min
                    if (float(blast_len) / float(dict_data[species][line[col]])) * 100 >= float(len_ratio[0]):
                        data_kept.append((float(blast_len) / float(dict_data[species][line[col]])) * 100)
                        kept_line += 1
                        outfile.write('\t'.join(line) + '\n')
                    else:
                        data_rm.append((float(blast_len) / float(dict_data[species][line[col]])) * 100)
                        ratio_rm += 1
                    continue
                else:
                    raise ValueError(f'Unknown data...\nPlease check if the given fasta file is correct.'
                                     f'\nspecies = {species}, col = {col}, start = {start}, end = {end}')

            ###################
            # align == double #
            ###################
            elif len_ratio == 'double':
                blast1_min, blast1_max = min(int(line[6]), int(line[7])), max(int(line[6]), int(line[7]))
                blast2_min, blast2_max = min(int(line[8]), int(line[9])), max(int(line[8]), int(line[9]))
                blast_len = min((blast1_max - blast1_min), (blast2_max - blast2_min))
                if line[0] in dict_data[sp1] and line[1] in dict_data[sp2]:
                    data_len = min(int(dict_data[sp1][line[0]]), int(dict_data[sp2][line[1]]))
                    if (float(blast_len) / float(data_len)) * 100 >= float(len_ratio[0]):
                        data_kept.append((float(blast_len) / float(data_len)) * 100)
                        kept_line += 1
                        outfile.write('\t'.join(line) + '\n')
                    else:
                        data_rm.append((float(blast_len) / float(data_len)) * 100)
                        ratio_rm += 1
                    continue
                else:
                    raise ValueError('Unknown data...\nPlease check: if both given fasta files are correct.')

        ###########
        # Genomic #
        ###########
        else:
            kept_line += 1
            outfile.write('\t'.join(line) + '\n')
            continue

outfile.close()
inputfile.close()

print(f'For a total of {nb_line} lines: {intra_rm + default_rm + ratio_rm} were deleted '
      f'({intra_rm} alignments on itself, {default_rm} default settings and {ratio_rm} alignment seq coverage.')
if len_ratio: print(f'seq coverage rm: min = {round(min(data_rm), 2)} ; mean = {round(mean(data_rm), 2)}; median = {round(median(data_rm), 2)}; max= {round(max(data_rm), 2)}). ')
print(f'However {kept_line} have been kept.')
if len_ratio: print(f'(min = {round(min(data_kept), 2)} ; mean = {round(mean(data_kept), 2)}; median = {round(median(data_kept), 2)}; max= {round(max(data_kept), 2)}).')
