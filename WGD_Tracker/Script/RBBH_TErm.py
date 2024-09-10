"""
RBBH Pipeline - Step 2 - Removing repeated regions
Author: Morgane MILIN
Last Update: July 2024
Input: (1) Filtered BLAST output text file (2) hardmasked fasta file + (if cds) gff file or TEs list txt file
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
masked1 = folder_name + sp1 + '.masked' if os.path.exists(folder_name + sp1 + '.masked') else None
masked2 = folder_name + sp2 + '.masked' if os.path.exists(folder_name + sp2 + '.masked') else None
gff1 = folder_name + sp1 + '.gff' if os.path.exists(folder_name + sp1 + '.gff') else None
gff2 = folder_name + sp2 + '.gff' if os.path.exists(folder_name + sp2 + '.gff') else None
intragenomic, TErm_type = eval(dict_param['intragenomic']), dict_param['TErm_type']
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
    raise ValueError('Error: If you need to remove repeated sequences from both genomes, you must provide information '
                     'for "SP1_TErm" and "SP2_TErm" in the configuration file.')
print(f"\nRepeat region filtration type: {TErm_type}")
print(f'data file: {data_path + data_file}\nmasked1: {masked1}\nmasked2: {masked2}\ngff1: {gff1}\ngff2: {gff2}\nsp1: {sp1}, sp2: {sp2}')
print(f'sp1_limit: {sp1_limit}, sp1_ft_gff: {sp1_ft_gff}, sp2_limit: {sp2_limit}, sp2_ft_gff: {sp2_ft_gff}')


#######################
# Dictionary Creation #
#######################

# dict_fasta = {chrom: seq, etc} or dict_TErm = {cds: TE%, etc}
dict_TErm, dict_fasta = {}, {}
sp, sp_col, start, end, sp_limit, sp_ft, sp_masked, sp_gff = None, None, None, None, None, None, None, None
##########
# double #
##########
if TErm_type == 'double':

    if sp1.endswith('_genomic'):
        dict_fasta = fct.dico_genomic(dico_fa=dict_fasta, datafile=masked1)
    else:
        dict_TErm = fct.dico_CDS(fa_masked=masked1, gff=gff1, dico_TEs=dict_TErm, limit=sp1_limit, target=sp1_ft_gff)

    if not intragenomic:
        if sp2.endswith('_genomic'):
            dict_fasta = fct.dico_genomic(dico_fa=dict_fasta, datafile=masked2)
        else:
            dict_TErm = fct.dico_CDS(fa_masked=masked2, gff=gff2, dico_TEs=dict_TErm, limit=sp2_limit, target=sp2_ft_gff)
##########
# simple #
##########
else:
    if sp1_limit is not None:
        sp, sp_col, start, end, sp_limit, sp_ft, sp_masked, sp_gff = sp1, 0, 6, 7, sp1_limit, sp1_ft_gff, masked1, gff1
    elif sp2_limit is not None:
        sp, sp_col, start, end, sp_limit, sp_ft, sp_masked, sp_gff = sp2, 1, 8, 9, sp2_limit, sp2_ft_gff, masked2, gff2
    else:
        raise ValueError('PROBLEM !!!')
    print(f'sp: {sp}, sp_col: {sp_col}, start: {start}, end: {end}, sp_limit: {sp_limit}, sp_ft: {sp_ft}, '
          f'sp_masked: {sp_masked}, sp_gff: {sp_gff}')

    if sp.endswith('_genomic'):
        dict_fasta = fct.dico_genomic(dico_fa=dict_fasta, datafile=sp_masked)
    else:
        dict_TErm = fct.dico_CDS(fa_masked=sp_masked, gff=sp_gff, dico_TEs=dict_TErm, limit=sp_limit, target=sp_ft)
print(f'\nlen(dict_fasta) = {len(dict_fasta)} \nlen(dict_TErm) = {len(dict_TErm)}')


#############################
# Removing Repeated Regions #
#############################

input_file = open(data_path + data_file)
outfile = open(folder_name + 'RBBH/TErm_' + data_file, 'w')
nb_line, TErm_count, sp1_rm, sp2_rm, kept_line = 0, 0, 0, 0, 0
for line in input_file:
    line = line.replace('\n', '').split()
    nb_line += 1

    ##########
    # DOUBLE #
    ##########
    if TErm_type == 'double':

        #############
        # Species 1 #
        #############
        if sp1.endswith('_genomic'):
            if line[0] in dict_fasta:
                TEs_res1 = fct.TE_recover(start=int(line[6]), end=int(line[7]), seq=dict_fasta[line[0]])
                if TEs_res1 > sp1_limit:
                    TErm_count += 1
                    continue
            else:
                raise ValueError(f"ERROR: {line[0]} isn't in dict_fasta, meaning its absence in the masked fasta file. "
                                 f"Please verify both fasta and config files for potential issues.")
        else:
            if line[0] in dict_TErm:
                TErm_count += 1
                sp1_rm += 1
                continue

        #############
        # Species 2 #
        #############
        if sp2.endswith('_genomic'):
            if line[1] in dict_fasta:
                TEs_res2 = fct.TE_recover(start=int(line[8]), end=int(line[9]), seq=dict_fasta[line[1]])
                if TEs_res2 > sp2_limit:
                    TErm_count += 1
                    continue
            else:
                raise ValueError(f"ERROR: {line[1]} isn't in dict_fasta, meaning its absence in the masked fasta file. "
                                 f"Please verify both fasta and config files for potential issues.")
        else:
            if line[1] in dict_TErm:
                TErm_count += 1
                sp2_rm += 1
                continue

    ##########
    # SIMPLE #
    ##########
    else:
        if sp.endswith('_genomic'):
            if line[sp_col] in dict_fasta:
                TEs_res1 = fct.TE_recover(start=int(line[start]), end=int(line[end]), seq=dict_fasta[line[sp_col]])
                if TEs_res1 > sp_limit:
                    TErm_count += 1
                    continue
            else:
                raise ValueError(f"ERROR: {line[sp_col]} isn't in dict_fasta, meaning its absence in the masked fasta "
                                 f"file. Please verify both fasta and config files for potential issues.")
        else:
            if line[sp_col] in dict_TErm:
                TErm_count += 1
                continue
    outfile.write('\t'.join(line) + '\n')
    kept_line += 1

input_file.close()
outfile.close()

print(f'For a total of {nb_line} lines: {TErm_count} were deleted (with sp1 rm = {sp1_rm} and sp2 rm = {sp2_rm}).\n'
      f'However {kept_line} have been kept.')
