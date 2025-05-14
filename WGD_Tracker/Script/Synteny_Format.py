"""
Step 1 - Search for Syntenic blocks - Formatting dataset
Removing pairwise genes present several time on the same chromosome (tandem duplication or other repeats).
Author: Morgane MILIN
Last Update: July 2024
Output: Filtered dataset file
"""

import sys
import os
import glob
from statistics import mean, median
import Synteny_functions as fct

# py script dir_path=$data_dir SP1=$SP1 SP2=$SP2

############################
# Get Files and Parameters #
############################

dict_param = {param.split('=')[0]: param.split('=')[1] for param in sys.argv[1:]}
wkdir_path, sp1, sp2, intragenomic = dict_param['dir_path'], dict_param['SP1'], dict_param['SP2'], eval(dict_param['Intragenomic'])
RBH_file = glob.glob(wkdir_path + "/RBBH/" + "RBBH_*")[0]
Ks_file = glob.glob(wkdir_path + "/Ks/Res_compil_NG_Ks_total.txt")[0]
gff1_file, gff1_motif = glob.glob(wkdir_path + "/" + sp1 + "*.gff")[0], dict_param['SP1_motif']
gff2_file, gff2_motif = glob.glob(wkdir_path + "/" + sp2 + "*.gff")[0], dict_param['SP2_motif']
Ks_min = float(dict_param['Ks_min']) if 'Ks_min' in dict_param else 0.01
Ks_max = float(dict_param['Ks_max']) if 'Ks_max' in dict_param else 3
repeat_limit = int(dict_param['limit']) if 'limit' in dict_param else 1
print(f"PATH: {wkdir_path}\nRBH_file: {RBH_file}\nKs_file: {Ks_file}\ngff1_file: {gff1_file}\ngff2_file: {gff2_file}")
print()


#######
# gff #
#######

dico_sp = {sp1: {}, sp2: {}}
for gff_file, motif, species in [(gff1_file, gff1_motif, sp1), (gff2_file, gff2_motif, sp2)]:
    print('gff_file =', gff_file)
    inputfile = open(gff_file)
    for line in inputfile:
        if not line.startswith('#'):
            line = line.replace('\n', '').split('\t')
            if line[2] == motif:
                chrom, gene, start, end = line[0], line[8].split('=')[1], int(line[3]), int(line[4])
                if chrom not in dico_sp[species]:
                    dico_sp[species][chrom] = {}
                if (start, end) not in dico_sp[species][chrom]:
                    dico_sp[species][chrom][(start, end)] = [gene]
                else:
                    if type(dico_sp[species][chrom][(start, end)]) is list:
                        dico_sp[species][chrom][(start, end)].append(gene)
                    else:
                        raise ValueError("dictionary format error")
    inputfile.close()
    if intragenomic:
        break
print(f'Number of sequences for {sp1}: {len(dico_sp[sp1])}\nNumber of sequences for {sp2}: {len(dico_sp[sp2])}')

dict_gff = {sp1: {}, sp2: {}}
for species in [sp1, sp2]:
    for chrom in sorted(dico_sp[species]):
        n = 0
        for position in sorted(dico_sp[species][chrom]):
            if type(dico_sp[species][chrom][position]) is list:
                n += 1
                for gene in dico_sp[species][chrom][position]:
                    if gene not in dict_gff[species]:
                        dict_gff[species][gene] = [chrom, n]
            else:
                raise ValueError(f"dictionary format error: {type(dico_sp[species][chrom][position])}")
print(f'Number of coding sequences for {sp1}: {len(dict_gff[sp1])}\nNumber of coding sequences for {sp2}: {len(dict_gff[sp2])}')


#############################
# RBH and Ks dataset Fusion #
#############################

dico_RBH, verif_count, dict_out, data_file = {}, 0, {}, wkdir_path + '/Synteny/dataset_' + os.path.basename(RBH_file).replace(".txt", "_v1.txt")
inputfile = open(RBH_file)
for line in inputfile:
    line = line.replace('\n', '').split('\t')
    sp1_gene, sp2_gene = line[0], line[1]
    if sp1_gene is None or sp2_gene is None:
        raise ValueError(f"That's not possible !!!\nsp1_gene = {sp1_gene}, sp2_gene = {sp2_gene}\nline = {line}")
    if (sp1_gene, sp2_gene) not in dico_RBH:
        dico_RBH[(sp1_gene, sp2_gene)] = float(line[11])
    else:
        if float(line[11]) > dico_RBH[(sp1_gene, sp2_gene)]:
            dico_RBH[(sp1_gene, sp2_gene)] = float(line[11])
inputfile.close()

outfile = open(data_file, 'w')
inputfile = open(Ks_file)
for line in inputfile:
    line = line.replace('\n', '').split('\t')
    sp1_gene, sp2_gene, cds_name = line[0], line[1], line[1]
    bit_score = dico_RBH[(sp1_gene, sp2_gene)] if (sp1_gene, sp2_gene) in dico_RBH else dico_RBH[(sp2_gene, sp1_gene)] if (sp2_gene, sp1_gene) in dico_RBH else None
    if bit_score is None:
        raise ValueError(f"Error !!!\nsp1_gene = {sp1_gene}, sp2_gene = {sp2_gene}\nline = {line}")
    dict_out[(sp1_gene, sp2_gene)] = True
    dict_out[(sp2_gene, sp1_gene)] = True
    outfile.write('\t'.join([sp1_gene, sp2_gene, line[2], line[3], line[4]]) + '\t' + str(bit_score) + '\n' +
                  '\t'.join([sp2_gene, sp1_gene, line[2], line[3], line[4]]) + '\t' + str(bit_score) + '\n')
inputfile.close()
outfile.close()
print(f'len(dico_RBH): {len(dico_RBH)}')


########
# Data #
########

dict_data1, dict_data2 = {}, {}
inputfile = open(data_file)
for line in inputfile:
    line = line.replace('\n', '').split('\t')
    if float(line[3]) < Ks_min or float(line[3]) >= Ks_max:
        continue
    gene1, gene2 = line[0], line[1]
    if gene1 in dict_gff[sp1] and gene2 in dict_gff[sp2]:
        sp1_chrom, sp1_nb = dict_gff[sp1][gene1][0], dict_gff[sp1][gene1][1]
        sp2_chrom, sp2_nb = dict_gff[sp2][gene2][0], dict_gff[sp2][gene2][1]
        fct.prep_data(chrom=(sp1_chrom, sp2_chrom), nb=[sp1_nb, sp2_nb], infos=[gene1, gene2, float(line[3]), float(line[5])], dico=dict_data1)
        fct.prep_data(chrom=(sp1_chrom, sp2_chrom), nb=[sp2_nb, sp1_nb], infos=[gene1, gene2, float(line[3]), float(line[5])], dico=dict_data2)
    elif gene2 in dict_gff[sp1] and gene1 in dict_gff[sp2]:
        sp1_chrom, sp1_nb = dict_gff[sp1][gene2][0], dict_gff[sp1][gene2][1]
        sp2_chrom, sp2_nb = dict_gff[sp2][gene1][0], dict_gff[sp2][gene1][1]
        fct.prep_data(chrom=(sp1_chrom, sp2_chrom), nb=[sp1_nb, sp2_nb], infos=[gene2, gene1, float(line[3]), float(line[5])], dico=dict_data1)
        fct.prep_data(chrom=(sp1_chrom, sp2_chrom), nb=[sp2_nb, sp1_nb], infos=[gene2, gene1, float(line[3]), float(line[5])], dico=dict_data2)
    else:
        raise ValueError(f"That's not possible !!! check gff for the following gene: {gene1, gene2}")
inputfile.close()
print(f'len(dict_data1): {len(dict_data1)}, len(dict_data2): {len(dict_data2)}')


##############################
# Retrieve gene pairs + Stat #
##############################

sp1_count, sp2_count, dict_pairwise, dict_stat = {}, {}, {}, {'sp1_gene': [], 'sp2_gene': [], 'stat_sp1': [], 'stat_sp2': []}
for dataset, count_list, stat, statut in [(dict_data1, sp2_count, 'stat_sp2', True), (dict_data2, sp1_count, 'stat_sp1', False)]:
    for chrom in sorted(dataset):
        for species in sorted(dataset[chrom]):
            dict_stat[stat].append(len(dataset[chrom][species]))
            if len(dataset[chrom][species]) != 1:
                if len(dataset[chrom][species]) not in count_list:
                    count_list[len(dataset[chrom][species])] = 1
                else:
                    count_list[len(dataset[chrom][species])] += 1

            if len(dataset[chrom][species]) != repeat_limit:
                for species2 in sorted(dataset[chrom][species]):
                    if (dataset[chrom][species][species2][0], dataset[chrom][species][species2][1]) not in dict_pairwise:
                        dict_pairwise[(dataset[chrom][species][species2][0], dataset[chrom][species][species2][1])] = True

            if statut:
                for species2 in sorted(dataset[chrom][species]):
                    if dataset[chrom][species][species2][0] not in dict_stat['sp1_gene']:
                        dict_stat['sp1_gene'].append(dataset[chrom][species][species2][0])
                    if dataset[chrom][species][species2][1] not in dict_stat['sp2_gene']:
                        dict_stat['sp2_gene'].append(dataset[chrom][species][species2][1])

print(f"{sp1} gene nb: {len(dict_stat['sp1_gene'])}")
print(f"Stat {sp1}: min = {min(dict_stat['stat_sp1'])}, mean = {round(mean(dict_stat['stat_sp1']), 3)}, median = {round(median(dict_stat['stat_sp1']), 3)}, max = {max(dict_stat['stat_sp1'])}")
print(f'\nNumber of hits for {sp1} genes')
for i in sorted(sp1_count):
    print(i, sp1_count[i])
if not intragenomic:
    print(f"\n{sp2} gene nb: {len(dict_stat['sp2_gene'])}")
    print(f"Stat {sp2}: min = {min(dict_stat['stat_sp2'])}, mean = {mean(dict_stat['stat_sp2'])}, median = {median(dict_stat['stat_sp2'])}, max = {max(dict_stat['stat_sp2'])}")
    print(f'\nNumber of hits for {sp2} genes')
    for i in sorted(sp2_count):
        print(i, sp2_count[i])
print('\nlen(dict_pairwise) =', len(dict_pairwise))


###############################################
# Data Filtering Step - Remove repeated genes #
###############################################

dict_gene, line_count, Ks_rm, gff_rm, repeats_rm, keep_count = {'sp1': {}, 'sp2': {}}, 0, 0, 0, 0, 0
outfile = open(data_file.replace('_v1.txt', '_v2.txt'), 'w')
inputfile = open(data_file)
for line in inputfile:
    line = line.replace('\n', '').split('\t')
    line_count += 1

    if float(line[3]) < Ks_min or float(line[3]) >= Ks_max:
        Ks_rm += 1
        continue

    if line[0] in dict_gff[sp1] and line[1] in dict_gff[sp2]:
        if (line[0], line[1]) in dict_pairwise or (line[1], line[0]) in dict_pairwise:
            repeats_rm += 1
            continue
        else:
            outfile.write('\t'.join(line) + '\n')
            if line[0] not in dict_gene['sp1']:
                dict_gene['sp1'][line[0]] = True
            if line[1] not in dict_gene['sp2']:
                dict_gene['sp2'][line[1]] = True
            keep_count += 1
    else:
        gff_rm += 1

inputfile.close()
outfile.close()

print(f'On a total of {line_count} homologous gene pairs:\n'
      f'\t{Ks_rm} were removed due to Ks value,\n'
      f'\t{gff_rm} were removed because at least of gene was not in the initial gff input,\n'
      f'\t{repeats_rm} were removed because one of both gene at least was matching with several genes from the same chromosome of the other species.\n'
      f'\tAs a results, {keep_count} homologous gene pairs have been kept ({round(keep_count / line_count * 100, 3)}% of the initial dataset)\n')
print("len(dict_gene['sp1']) =", len(dict_gene['sp1']))
print("len(dict_gene['sp2']) =", len(dict_gene['sp2']))

