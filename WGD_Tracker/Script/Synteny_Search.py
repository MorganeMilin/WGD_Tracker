"""
Step 1 - Search for Syntenic blocks
Searching for Syntenic Blocks
Author: Morgane MILIN
Last Update: August 2024
"""

##################
# Modules import #
##################

import sys
import glob
import Synteny_functions as fct


############################
# Get Files and Parameters #
############################

dict_param = {param.split('=')[0]: param.split('=')[1] for param in sys.argv[1:]}
wkdir_path, sp1, sp2 = dict_param['dir_path'], dict_param['SP1'], dict_param['SP2']
intragenomic, corr_name = eval(dict_param['Intragenomic']), eval(dict_param['corr_SB'])
data_file = glob.glob(wkdir_path + "/Synteny/dataset_*_v3.txt")[0]
gff1_motif, gff2_motif, gap_limit, gene_limit = dict_param['SP1_motif'], dict_param['SP2_motif'], int(dict_param['gap_limit']), int(dict_param['gene_limit'])
print(f"\ndata_file: {data_file}")


#######
# gff #
#######

dict_sp, dict_gff = fct.retrieve_gff_infos(wkdir_path=wkdir_path, SP=[sp1, sp2], Motif=[gff1_motif, gff2_motif], corr_name=corr_name, intragenomic=intragenomic)
print(f'\nNumber of Chromosomes:\nSpecies 1: {len(dict_sp[sp1])}\nSpecies 2: {len(dict_sp[sp2])}')
print(f'Number of genes:\nSpecies 1: {len(dict_gff[sp1])}\nSpecies 2: {len(dict_gff[sp2])}')


###########
# Dataset #
###########

# dict_data = {(sp1_chrom, sp2_chrom): {(sp1_rank, sp2_rank): [sp1_gene, sp2_gene, Ks]}}
dict_data, db_rm = {}, 0
inputfile = open(data_file)
for line in inputfile:
    line = line.replace('\n', '').split('\t')
    if 0.01 <= float(line[3]) < 3:
        if line[0] in dict_gff[sp1] and line[1] in dict_gff[sp2]:
            chrom1, nb1, chrom2, nb2 = dict_gff[sp1][line[0]][0], dict_gff[sp1][line[0]][1], dict_gff[sp2][line[1]][0], dict_gff[sp2][line[1]][1]
            if (chrom1, chrom2) not in dict_data:
                dict_data[(chrom1, chrom2)] = {}
            if (nb1, nb2) not in dict_data[(chrom1, chrom2)]:
                dict_data[(chrom1, chrom2)][(nb1, nb2)] = [line[0], line[1], float(line[3]), float(line[5])]
            else:
                raise ValueError("Error! That's not expected")
        else:
            raise ValueError("Error! That's not expected")
    else:
        raise ValueError("Error! That's not expected")
inputfile.close()


print('\n#####################################################')
print('# Syntenic Blocks Search - Step 1 - Formatting data #')
print('#####################################################\n')

header = '#Sp1Chrom_Sp2Chrom_SyntenicBlock\tSp1Chrom\tSp2Chrom\tSyntenicBlock\torientation\tNb_genes_in_block\tSp1_gene_rank_list\tSp1_gene_name_list\tSp2_gene_rank_list\tSp2_gene_name_list\tKs\n'
outfile = open(wkdir_path + '/Synteny/Syntenic_blocks_STEP_1.txt', 'w')
outfile.write(header)
for chrom in sorted(dict_data):
    m = 0
    for nb in sorted(dict_data[chrom]):
        m += 1
        outfile.write('\t'.join(['-'.join([str(chrom[0]), str(chrom[1]), str(m)]), str(chrom[0]), str(chrom[1]), str(m), 'None', str(1), str(nb[0]), str(dict_data[chrom][nb][0]), str(nb[1]), str(dict_data[chrom][nb][1]), str(dict_data[chrom][nb][2])]) + '\n')
outfile.close()
fct.dotplot_extraction(dataset=wkdir_path + '/Synteny/Syntenic_blocks_STEP_1.txt', outname=wkdir_path + '/Synteny/Res_compil_NG_Syntenic_blocks_STEP_1.txt', datatype='str')


print('\n##################################################################')
print('# Syntenic Blocks Search - Step 2 - Detection of syntenic blocks #')
print('##################################################################\n')

inputfile = open(wkdir_path + '/Synteny/Syntenic_blocks_STEP_1.txt')
inputfile.readline()
chrom, block_nb, vec_sens, vec_sp1, vec_sp2, vec_sp2_flt_sens = None, None, None, None, None, {}
for line in inputfile:
    line = line.replace('\n', '').split('\t')
    if chrom is None:
        chrom, block_nb, vec_sp1, vec_sp2, vec_sens = (line[1], line[2]), int(line[3]), [], [], []
        vec_sp1.append(int(line[6]))
        vec_sp2.append(int(line[8]))
        vec_sens.append(line[4])
    else:
        if chrom == (line[1], line[2]):
            vec_sp1.append(int(line[6]))
            vec_sp2.append(int(line[8]))
            vec_sens.append(line[4])
        else:
            ##############################
            # Check vector for Species 1 #
            ##############################
            if vec_sp1 != sorted(vec_sp1):
                raise ValueError(f'ERROR in vec_sp1 for {chrom}')
            if len(vec_sp1) != len(vec_sp2):
                raise ValueError(f'ERROR !!! len(vec_sp1) != len(vec_sp2)')
            ################################
            # Detection of syntenic blocks #
            ################################
            if vec_sp2:
                vec_sp2_flt = fct.rm_isolated_outlier_data(vec_sp2=vec_sp2, gap_limit=gap_limit)     # remove isolated outlier data
                vec_sp2_flt_sens = fct.generate_oriented_SB(chrom=chrom, vec_sp2_flt=vec_sp2_flt, vec_sp2_flt_sens=vec_sp2_flt_sens)    # Generate oriented syntenic blocks
            ################################################
            # Reset initial vector for new chromosome pair #
            ################################################
            chrom, block_nb, vec_sp1, vec_sp2, vec_sens = (line[1], line[2]), int(line[3]), [], [], []
            vec_sp1.append(int(line[6]))
            vec_sp2.append(int(line[8]))
            vec_sens.append(line[4])
##############################
# Check vector for Species 1 #
##############################
if vec_sp1 != sorted(vec_sp1):
    raise ValueError(f'ERROR in vec_sp1 for {chrom}')
if len(vec_sp1) != len(vec_sp2):
    raise ValueError(f'ERROR !!! len(vec_sp1) != len(vec_sp2)')
################################
# Detection of syntenic blocks #
################################
if vec_sp2:
    vec_sp2_flt = fct.rm_isolated_outlier_data(vec_sp2=vec_sp2, gap_limit=gap_limit)     # remove isolated outlier data
    vec_sp2_flt_sens = fct.generate_oriented_SB(chrom=chrom, vec_sp2_flt=vec_sp2_flt, vec_sp2_flt_sens=vec_sp2_flt_sens)    # Generate oriented syntenic blocks
inputfile.close()
print(f'Number of pair of chromosomes: {len(vec_sp2_flt_sens)}')


####################################
# Transformation into a dictionary #
####################################

inputfile = open(wkdir_path + '/Synteny/Syntenic_blocks_STEP_1.txt')
inputfile.readline()
dict_blocks = {}
for line in inputfile:
    line = line.replace('\n', '').split('\t')
    chrom = (line[1], line[2])
    if chrom in vec_sp2_flt_sens:
        for m in vec_sp2_flt_sens[chrom]:
            if int(line[8]) in vec_sp2_flt_sens[chrom][m][1]:
                orientation = vec_sp2_flt_sens[chrom][m][0]
                if '-'.join([chrom[0], chrom[1], str(m)]) not in dict_blocks:
                    dict_blocks['-'.join([chrom[0], chrom[1], str(m)])] = [1, [line[6]], [line[7]], [line[8]], [line[9]], [line[10]], orientation]
                else:
                    # Number of gene
                    dict_blocks['-'.join([chrom[0], chrom[1], str(m)])][0] += 1
                    # Species 1 genes
                    dict_blocks['-'.join([chrom[0], chrom[1], str(m)])][1].append(line[6])
                    dict_blocks['-'.join([chrom[0], chrom[1], str(m)])][2].append(line[7])
                    # Species 2 genes
                    dict_blocks['-'.join([chrom[0], chrom[1], str(m)])][3].append(line[8])
                    dict_blocks['-'.join([chrom[0], chrom[1], str(m)])][4].append(line[9])
                    # Ks
                    dict_blocks['-'.join([chrom[0], chrom[1], str(m)])][5].append(line[10])
                    # check orientation
                    if orientation != dict_blocks['-'.join([chrom[0], chrom[1], str(m)])][6]:
                        raise ValueError("That's not possible !!!")
print(f'Number of syntenic blocks: {len(dict_blocks)}')
inputfile.close()


###########
# Outfile #
###########

outfile = open(wkdir_path + '/Synteny/Syntenic_blocks_STEP_2.txt', 'w')
outfile.write(header)
for i in dict_blocks:
    scaff, chrom, m = i.split('-')[0], i.split('-')[1], i.split('-')[2]
    sp1_nb, sp1_genes = str(dict_blocks[i][1]), str(dict_blocks[i][2])
    sp2_nb, sp2_genes = str(dict_blocks[i][3]), str(dict_blocks[i][4])
    nb_gene, Ks_values, orientation = str(dict_blocks[i][0]), str(dict_blocks[i][5]), str(dict_blocks[i][6])
    outfile.write('\t'.join([i, scaff, chrom, m, orientation, nb_gene, sp1_nb, sp1_genes, sp2_nb, sp2_genes, Ks_values]) + '\n')
outfile.close()
fct.dotplot_extraction(dataset=wkdir_path + '/Synteny/Syntenic_blocks_STEP_2.txt', outname=wkdir_path + '/Synteny/Res_compil_NG_Syntenic_blocks_STEP_2.txt', datatype='list')


print("\n####################")
print("# Gap verification #")
print("####################\n")

inputfile = open(wkdir_path + '/Synteny/Syntenic_blocks_STEP_2.txt')
inputfile.readline()
dict_synteny, block, sens, vec_sp1, vec_rank_sp1, vec_sp2, vec_rank_sp2, vec_Ks = {}, None, None, None, None, None, None, None
for line in inputfile:
    line = line.replace('\n', '').split('\t')

    ##############################################
    # On vérifie le nombre de gène dans le block #
    ##############################################
    if int(line[5]) >= gene_limit:
        block, sens, vec_sp1, vec_sp2 = line[0], line[4], eval(line[7]), eval(line[9])
        vec_rank_sp1, vec_rank_sp2, vec_Ks = [], [], []
        for x in eval(line[6]): vec_rank_sp1.append(int(x))
        for x in eval(line[8]): vec_rank_sp2.append(int(x))
        for x in eval(line[10]): vec_Ks.append(float(x))

        if block not in dict_synteny:
            dict_synteny[block] = {}
        else:
            raise ValueError(f"That's impossible !!! each line is a different syntenic block !!! {block}")

        # Check vector for Species 1
        if vec_rank_sp1 != sorted(vec_rank_sp1):
            raise ValueError(f'ERROR in vec_rank_sp1 for {block}\n{vec_rank_sp1}\n!=\n{sorted(vec_rank_sp1)}')

        # On break en fonction des gaps
        fct.SB_gap_only(vec_sp1=vec_rank_sp1, vec_sp2=vec_rank_sp2, sp1_gene=vec_sp1, sp2_gene=vec_sp2, Ks_vec=vec_Ks, orientation=sens, dico=dict_synteny, key=block, gap_limit=gap_limit)

inputfile.close()


###########
# Outfile #
###########

ref_verif, gene_verif = {}, []
outfile = open(wkdir_path + '/Synteny/Syntenic_blocks_STEP_3.txt', 'w')
outfile.write(header)
for i in dict_synteny:
    for j in dict_synteny[i]:
        if len(dict_synteny[i][j][3]) >= gene_limit:
            scaff, chrom = i.split('-')[0], i.split('-')[1]
            ref_verif[(scaff, chrom)] = 1 if (scaff, chrom) not in ref_verif else ref_verif[(scaff, chrom)] + 1
            m = str(ref_verif[(scaff, chrom)])
            sp1_nb, sp1_genes, sp2_nb, sp2_genes = str(dict_synteny[i][j][1]), str(dict_synteny[i][j][3]), str(dict_synteny[i][j][2]), str(dict_synteny[i][j][4])
            nb_gene, Ks_values, orientation = str(len(dict_synteny[i][j][3])), str(dict_synteny[i][j][5]), str(dict_synteny[i][j][0])
            outfile.write('\t'.join(['-'.join([scaff, chrom, m]), scaff, chrom, m, orientation, nb_gene, sp1_nb, sp1_genes, sp2_nb, sp2_genes, Ks_values]) + '\n')
        else:
            gene_verif.append(len(dict_synteny[i][j][3]))

outfile.close()
print(f"nb rm = {len(gene_verif)}; min = {min(gene_verif)}; max = {max(gene_verif)}")
fct.dotplot_extraction(dataset=wkdir_path + '/Synteny/Syntenic_blocks_STEP_3.txt', outname=wkdir_path + '/Synteny/Res_compil_NG_Syntenic_blocks_STEP_3.txt', datatype='list')
