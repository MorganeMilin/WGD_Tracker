"""
Dotplot - Step 1 - Data formatting
Author: Morgane MILIN
Last Update: July 2024
"""

##################
# Modules import #
##################

import sys
import os
import glob
import Dotplot_functions as fct

# py Dotplot_formatting_v2.py dir_path="./../data" SP1="Smaritima_cds" SP2="Osativa_cds" name_corr="\"['Spartina_maritima_scaffold_', False]"\" Intragenomic="False" control_file="dotplot_control_file_Smaritima_cds_vs_Osativa_cds.txt" dataset="Ks" gff1_motif="mRNA" gff2_motif="mRNA"

############################
# Get Files and Parameters #
############################

dict_param = {param.split('=')[0]: param.split('=')[1] for param in sys.argv[1:]}
wkdir_path, sp1, sp2 = dict_param['dir_path'], dict_param['SP1'], dict_param['SP2']
name_corr = eval(eval(dict_param['name_corr'])) if 'name_corr' in dict_param and dict_param['name_corr'].replace('"', "").replace("'", "") else False
intragenomic = eval(dict_param['Intragenomic']) if 'Intragenomic' in dict_param else True if sp1 == sp2 else False
control_file, dataset = dict_param['control_file'], dict_param['dataset']
gff1_motif = dict_param['gff1_motif'] if 'gff1_motif' in dict_param else None
gff2_motif = dict_param['gff2_motif'] if 'gff2_motif' in dict_param else None
axis_format = eval(dict_param['axis_format']) if 'axis_format' in dict_param else True
data_format = eval(dict_param['data_format']) if 'data_format' in dict_param else True
data_file = glob.glob(wkdir_path + "/RBBH/RBBH_*")[0] if dataset == "RBBH" else glob.glob(wkdir_path + "/Ks/Res_compil_NG_Ks_total.txt")[0] if dataset == "Ks" else glob.glob(wkdir_path + "/Synteny/Res_compil_*_STEP_3.txt")[0] if dataset == "Synteny" else None
print(f'\ncontrol_file: {control_file}\nname_corr: {name_corr}, {type(name_corr)}\ndataset: {dataset}\ndata_file: {data_file}')
if not data_file:
    raise ValueError("You forgot to mention the dataset to use for the Dotplot representation")


############################
# Dotplot format parameter #
############################

dict_infos, dict_gff, statut = {}, {}, None
input_control = open(control_file)
outfile_formatting = open(wkdir_path + '/Dotplot/dotplot_formatting_file.txt', 'w') if axis_format else None
for line in input_control:
    line = line.replace('\n', '').split('\t')
    if statut:
        break
    for axis in ['x axis', 'y axis']:
        if line[0] == axis:
            species, axis_title, axis_order = line[1], line[2], line[3].split(',')
            correction = name_corr[0] if name_corr and sp1 == line[1] else name_corr[1] if name_corr and sp2 == line[1] else None
            print(f"\naxis: {axis}\nspecies: {species}\ncorrection: {correction}")
            fasta_file = glob.glob(wkdir_path + "/" + species + ".fasta")[0]
            dict_fasta, dict_fasta_v2, axis_lines = fct.formatting_dotplot(fasta=fasta_file, list_chr=axis_order, name_corr=correction)
            print(f"fasta_file: {fasta_file}\nlen(dict_fasta) = {len(dict_fasta)}")
            axis_ticks = fct.ticks_generator(chr_list=axis_order, dico_data=dict_fasta)
            if not species.endswith('_genomic'):
                gff_file = glob.glob(wkdir_path + "/" + species + ".gff")[0]
                gff_motif = gff1_motif if sp1 == species else gff2_motif if sp2 == species else None
                dict_gff = fct.gff_infos(datafile=gff_file, motif=gff_motif, name_corr=correction)
                print(f"gff_file: {gff_file}\nlen(dict_gff) = {len(dict_gff)}")
            if axis not in dict_infos:
                dict_infos[axis] = [species, correction, axis_order, dict_fasta_v2, dict_gff]
                if intragenomic:
                    dict_infos['y axis'] = [species, correction, axis_order, dict_fasta_v2, dict_gff]
            if axis_format:
                outfile_formatting.write(f'{axis}\ttitle\t{axis_title}\n')
                outfile_formatting.write(f'{axis}\tticks\t{axis_ticks}\n')
                outfile_formatting.write(f'{axis}\tlabel\t{axis_order}\n')
                outfile_formatting.write(f'{axis}\tlines\t{axis_lines}\n')
                if intragenomic:
                    outfile_formatting.write(f'y axis\ttitle\t{axis_title}\n')
                    outfile_formatting.write(f'y axis\tticks\t{axis_ticks}\n')
                    outfile_formatting.write(f'y axis\tlabel\t{axis_order}\n')
                    outfile_formatting.write(f'y axis\tlines\t{axis_lines}\n')
                    statut = True
                    break
input_control.close()
if axis_format:
    outfile_formatting.close()

x_sp, x_corr, x_order, dict_fasta_x, dict_gff_x = dict_infos['x axis'][0], dict_infos['x axis'][1], dict_infos['x axis'][2], dict_infos['x axis'][3], dict_infos['x axis'][4]
y_sp, y_corr, y_order, dict_fasta_y, dict_gff_y = dict_infos['y axis'][0], dict_infos['y axis'][1], dict_infos['y axis'][2], dict_infos['y axis'][3], dict_infos['y axis'][4]
print(f"\nx_sp: {x_sp}; x_correction: {x_corr}; x_order: {x_order}\ny_sp: {y_sp}; y_correction: {y_corr}; y_order: {y_order}\n")


#############################################
# Data formatting for dotplot visualisation #
#############################################

statut, x_col, x_start, x_end, y_col, y_start, y_end = False, None, None, None, None, None, None
line_count, Ks_rm, gff_rm, chrom_rm = 0, 0, 0, 0
outfile_dataset = open(wkdir_path + '/Dotplot/dotplot_dataset_file.txt', 'w') if data_format else None
inputfile = open(data_file)
for line in inputfile:
    line = line.replace('\n', '').split('\t')
    line_count += 1

    #################
    # Ks filtration #
    #################
    if dataset != "RBBH":
        if float(line[3]) < 0.01 or float(line[3]) >= 3:
            Ks_rm += 1
            continue

    #######################################
    # define which species for which axis #
    #######################################
    if not statut:
        if intragenomic:
            x_col, x_start, x_end, y_col, y_start, y_end = 0, 6, 7, 1, 8, 9
        else:
            if line[0] in dict_gff_x and line[1] in dict_gff_y:
                x_col, x_start, x_end, y_col, y_start, y_end = 0, 6, 7, 1, 8, 9
            elif line[1] in dict_gff_x and line[0] in dict_gff_y:
                x_col, x_start, x_end, y_col, y_start, y_end = 1, 8, 9, 0, 6, 7
            else:
                gff_rm += 1
                continue
        statut = True

    #########################################
    # define the data type - genomic or cds #
    #########################################
    res = [None, None]
    for pos, species, correction, data_col, dict_gff in [(0, x_sp, x_corr, x_col, dict_gff_x), (1, y_sp, y_corr, y_col, dict_gff_y)]:
        if species.endswith('_genomic'):
            res[pos] = line[data_col].replace(correction, '') if correction else line[data_col]
        else:
            if line[data_col] in dict_gff:
                res[pos] = dict_gff[line[data_col]][0]
            else:
                gff_rm += 1
                break
    x_name, y_name = res[0], res[1]
    if x_name not in x_order or y_name not in y_order:
        chrom_rm += 1
        continue
    else:
        x_start = line[x_start] if x_sp.endswith('_genomic') else dict_gff_x[line[x_col]][1] if line[x_col] in dict_gff_x else None
        x_end = line[x_end] if x_sp.endswith('_genomic') else dict_gff_x[line[x_col]][2] if line[x_col] in dict_gff_x else None
        y_start = line[y_start] if y_sp.endswith('_genomic') else dict_gff_y[line[y_col]][1] if line[y_col] in dict_gff_y else None
        y_end = line[y_end] if y_sp.endswith('_genomic') else dict_gff_y[line[y_col]][2] if line[y_col] in dict_gff_y else None

    ############################
    # Correction des positions #
    ############################
    if x_name in dict_fasta_x and y_name in dict_fasta_y:
        x_corr_start = int(x_start) + int(dict_fasta_x[x_name])
        x_corr_end = int(x_end) + int(dict_fasta_x[x_name])
        y_corr_start = int(y_start) + int(dict_fasta_y[y_name])
        y_corr_end = int(y_end) + int(dict_fasta_y[y_name])
    else:
        chrom_rm += 1
        print(x_name, y_name, 'issue...')
        continue

    #########
    # Color #
    #########
    if dataset == "RBBH":   # color dots with identity values
        data_color = str(line[2])
    else:                   # color dots with Ks values
        data_color = str(line[3])

    ###########
    # Outfile #
    ###########
    outfile_dataset.write(x_name + '\t' + str(x_corr_start) + '\t' + y_name + '\t' + str(y_corr_start) + '\t' + data_color + '\n')

inputfile.close()
if data_format:
    outfile_dataset.close()
print(f"On a total of {line_count} lines, {Ks_rm} were removed due to Ks values, "
      f"{gff_rm} were removed because genes were not in the gff files and "
      f"{chrom_rm} were removed because the chromosome was not in the fasta file or on the chrom order list")