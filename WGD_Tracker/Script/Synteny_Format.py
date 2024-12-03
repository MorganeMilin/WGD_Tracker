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
import pprint

# py script dir_path=$data_dir SP1=$SP1 SP2=$SP2

############################
# Get Files and Parameters #
############################

dict_param = {param.split('=')[0]: param.split('=')[1] for param in sys.argv[1:]}
wkdir_path, sp1, sp2 = dict_param['dir_path'], dict_param['SP1'], dict_param['SP2']
intragenomic, corr_name = eval(dict_param['Intragenomic']), eval(dict_param['corr_SB'])
RBH_file = glob.glob(wkdir_path + "/RBBH/" + "RBBH_*")[0]
Ks_file = glob.glob(wkdir_path + "/Ks/Res_compil_NG_Ks_total.txt")[0]
gff1_motif, gff2_motif = dict_param['SP1_motif'], dict_param['SP2_motif']
Ks_min = float(dict_param['Ks_min']) if 'Ks_min' in dict_param else 0.01
Ks_max = float(dict_param['Ks_max']) if 'Ks_max' in dict_param else 3
print(f"PATH: {wkdir_path}\nRBH_file: {RBH_file}\nKs_file: {Ks_file}\n")


#######
# gff #
#######

dict_sp, dict_gff = fct.retrieve_gff_infos(wkdir_path=wkdir_path, SP=[sp1, sp2], Motif=[gff1_motif, gff2_motif], corr_name=corr_name, intragenomic=intragenomic)


#############################
# RBH and Ks dataset Fusion #
#############################

dico_RBH, verif_count, dict_out, data_file = {}, 0, {}, wkdir_path + '/Synteny/dataset_' + os.path.basename(RBH_file).replace(".txt", "_v1.txt")
inputfile = open(RBH_file)
for line in inputfile:
    line = line.replace('\n', '').split('\t')

    if intragenomic and corr_name:
        sp1_gene, sp2_gene = corr_name[1].join(line[0].split(corr_name[1])[:corr_name[2]]), corr_name[1].join(line[1].split(corr_name[1])[:corr_name[2]])
    elif corr_name and corr_name[0] == sp1:
        sp1_gene, sp2_gene = corr_name[1].join(line[0].split(corr_name[1])[:corr_name[2]]), corr_name[4].join(line[1].split(corr_name[4])[:corr_name[5]])
    elif corr_name and corr_name[0] == sp2:
        sp1_gene, sp2_gene = corr_name[4].join(line[0].split(corr_name[4])[:corr_name[5]]), corr_name[1].join(line[1].split(corr_name[1])[:corr_name[2]])
    else:
        sp1_gene, sp2_gene = line[0], line[1]

    if sp1_gene is None or sp2_gene is None:
        raise ValueError(f"Gene name not found!\nsp1_gene = {sp1_gene}, sp2_gene = {sp2_gene}\nline = {line}")

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

    if intragenomic and corr_name:
        sp1_gene, sp2_gene = corr_name[1].join(line[0].split(corr_name[1])[:corr_name[2]]), corr_name[1].join(line[1].split(corr_name[1])[:corr_name[2]])
    elif corr_name and corr_name[0] == sp1:
        sp1_gene, sp2_gene = corr_name[1].join(line[0].split(corr_name[1])[:corr_name[2]]), corr_name[4].join(line[1].split(corr_name[4])[:corr_name[5]])
    elif corr_name and corr_name[0] == sp2:
        sp1_gene, sp2_gene = corr_name[4].join(line[0].split(corr_name[4])[:corr_name[5]]), corr_name[1].join(line[1].split(corr_name[1])[:corr_name[2]])
    else:
        sp1_gene, sp2_gene = line[0], line[1]

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


##############################
# Retrieve gene pairs + Stat #
##############################

dict_data, sp1_count, sp2_count, dict_pairwise, bitscore_dict = {sp1: {}, sp2: {}}, {}, {}, {}, {}
dict_stat = {'sp1_gene': [], 'sp2_gene': [], 'stat_sp1': [], 'stat_sp2': []}
inputfile = open(data_file)
for line in inputfile:
    line = line.replace('\n', '').split('\t')
    if float(line[3]) < Ks_min or float(line[3]) >= Ks_max:
        continue
    sp1_gene, sp2_gene = line[0], line[1]
    if sp1_gene in dict_gff[sp1] and sp2_gene in dict_gff[sp2]:
        sp1_chrom, sp1_nb, sp2_chrom, sp2_nb = dict_gff[sp1][sp1_gene][0], dict_gff[sp1][sp1_gene][1], dict_gff[sp2][sp2_gene][0], dict_gff[sp2][sp2_gene][1]
        if (sp1_chrom, sp2_chrom) not in bitscore_dict:
            bitscore_dict[(sp1_chrom, sp2_chrom)] = {}
        if (sp1_nb, sp2_nb) not in bitscore_dict[(sp1_chrom, sp2_chrom)]:
            bitscore_dict[(sp1_chrom, sp2_chrom)][(sp1_nb, sp2_nb)] = [line[0], line[1], float(line[3]), float(line[5])]
        else:
            if float(line[5]) > bitscore_dict[(sp1_chrom, sp2_chrom)][(sp1_nb, sp2_nb)][3]:
                bitscore_dict[(sp1_chrom, sp2_chrom)][(sp1_nb, sp2_nb)] = [line[0], line[1], float(line[3]), float(line[5])]

        fct.prep_data(chrom=(sp1_chrom, sp2_chrom), nb=[sp1_nb, sp2_nb], infos=[sp1_gene, sp2_gene, float(line[3]), float(line[5])], dico=dict_data[sp1])
        fct.prep_data(chrom=(sp1_chrom, sp2_chrom), nb=[sp2_nb, sp1_nb], infos=[sp1_gene, sp2_gene, float(line[3]), float(line[5])], dico=dict_data[sp2])

    else:
        raise ValueError(f"Error! Please, check gff for the following gene: {sp1_gene, sp2_gene}")
inputfile.close()
print(f'len(dict_data[sp1]): {len(dict_data[sp1])}, len(dict_data[sp2]): {len(dict_data[sp2])}')

for dataset, count_list, stat, statut in [(dict_data[sp1], sp2_count, 'stat_sp2', True), (dict_data[sp2], sp1_count, 'stat_sp1', False)]:
    for chrom in sorted(dataset):
        for species in sorted(dataset[chrom]):
            dict_stat[stat].append(len(dataset[chrom][species]))
            if len(dataset[chrom][species]) != 1:
                if len(dataset[chrom][species]) not in count_list:
                    count_list[len(dataset[chrom][species])] = 1
                else:
                    count_list[len(dataset[chrom][species])] += 1

            if len(dataset[chrom][species]) != 1:
                for species2 in sorted(dataset[chrom][species]):
                    if (dataset[chrom][species][species2][0], dataset[chrom][species][species2][1]) not in dict_pairwise:
                        dict_pairwise[(dataset[chrom][species][species2][0], dataset[chrom][species][species2][1])] = True

            if statut:
                for species2 in sorted(dataset[chrom][species]):
                    if dataset[chrom][species][species2][0] not in dict_stat['sp1_gene']:
                        dict_stat['sp1_gene'].append(dataset[chrom][species][species2][0])
                    if dataset[chrom][species][species2][1] not in dict_stat['sp2_gene']:
                        dict_stat['sp2_gene'].append(dataset[chrom][species][species2][1])
print(f"{sp1} gene nb: {len(dict_stat['sp1_gene'])}\n"
      f"Stat {sp1}: min = {min(dict_stat['stat_sp1'])}, mean = {round(mean(dict_stat['stat_sp1']), 3)}, median = {round(median(dict_stat['stat_sp1']), 3)}, max = {max(dict_stat['stat_sp1'])}\n"
      f"Number of hits for {sp1} genes")
for i in sorted(sp1_count):
    print(i, sp1_count[i])
if not intragenomic:
    print(f"\n{sp2} gene nb: {len(dict_stat['sp2_gene'])}\n"
          f"Stat {sp2}: min = {min(dict_stat['stat_sp2'])}, mean = {mean(dict_stat['stat_sp2'])}, median = {median(dict_stat['stat_sp2'])}, max = {max(dict_stat['stat_sp2'])}\n"
          f"Number of hits for {sp2} genes")
    for i in sorted(sp2_count):
        print(i, sp2_count[i])
print(f"\nlen(dict_pairwise) = {len(dict_pairwise)}")


###############################################
# Data Filtering Step - Remove repeated genes #
###############################################

dict_gene, line_count, Ks_rm, gff_rm, repeats_rm, keep_count = {'sp1': {}, 'sp2': {}}, 0, 0, 0, 0, 0
outfile2 = open(data_file.replace('_v1.txt', '_v2.txt'), 'w')
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
            outfile2.write('\t'.join(line) + '\n')
            if line[0] not in dict_gene['sp1']:
                dict_gene['sp1'][line[0]] = True
            if line[1] not in dict_gene['sp2']:
                dict_gene['sp2'][line[1]] = True
            keep_count += 1
    else:
        gff_rm += 1

inputfile.close()
outfile2.close()

# dict_data = {(sp1_chrom, sp2_chrom): {(sp1_rank, sp2_rank): [sp1_gene, sp2_gene, Ks]}}
#out_check, db_rm, check_count = {}, 0, 0
dict_data, keep_line, db_rm = {}, 0, 0
inputfile = open(data_file.replace('_v1.txt', '_v2.txt'))
for line in inputfile:
    line = line.replace('\n', '').split('\t')
    if 0.01 <= float(line[3]) < 3:
        if line[0] in dict_gff[sp1] and line[1] in dict_gff[sp2]:
            chrom_sp1, sp1_nb = dict_gff[sp1][line[0]][0], dict_gff[sp1][line[0]][1]
            chrom_sp2, sp2_nb = dict_gff[sp2][line[1]][0], dict_gff[sp2][line[1]][1]
            if (chrom_sp1, chrom_sp2) not in dict_data:
                dict_data[(chrom_sp1, chrom_sp2)] = {}
            if (sp1_nb, sp2_nb) not in dict_data[(chrom_sp1, chrom_sp2)]:
                dict_data[(chrom_sp1, chrom_sp2)][(sp1_nb, sp2_nb)] = [float(line[5]), '\t'.join(line)]
                keep_line += 1
            else:
                db_rm += 1
                if float(line[5]) > dict_data[(chrom_sp1, chrom_sp2)][(sp1_nb, sp2_nb)][0]:
                    dict_data[(chrom_sp1, chrom_sp2)][(sp1_nb, sp2_nb)] = [float(line[5]), '\t'.join(line)]
        else:
            raise ValueError('GFF PROBLEM !!!')
    else:
        raise ValueError('KS PROBLEM !!!')
inputfile.close()


outfile3 = open(wkdir_path + '/Synteny/dataset_' + os.path.basename(RBH_file).replace(".txt", "_v3.txt"), 'w')
for chrom in dict_data:
    for nb in dict_data[chrom]:
        res_line = dict_data[chrom][nb][1].split('\t')
        print(f'chrom: {chrom}\nnb: {nb}\nres_line: {res_line}')
        outfile3.write('\t'.join(res_line) + '\n')
outfile3.close()


if keep_line != keep_count - db_rm:
    raise ValueError(f'{keep_line} != {keep_count} - {db_rm}')
else:
    print(f'On a total of {line_count} homologous gene pairs, {Ks_rm + gff_rm + repeats_rm + db_rm} were removed:\n'
          f'\t{Ks_rm} were removed due to Ks value,\n'
          f'\t{gff_rm} were removed because one of both gene, at least, was not in the initial gff input,\n'
          f'\t{repeats_rm} were removed because one of both gene, at least, was matching with several genes from the same chromosome of the other species,\n'
          f'\t{db_rm} were removed because the hit was present several time within the data file.\n'
          f'\tAs a results, {keep_line} homologous gene pairs have been kept ({round(keep_line / line_count * 100, 2)}% of the initial dataset)\n')
