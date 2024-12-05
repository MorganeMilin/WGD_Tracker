"""
Karyotype representation
Author: Morgane MILIN
Last Update: December 2024
"""

import sys
import os
import glob
from svg_turtle import SvgTurtle
import Karyotype_functions as fct
from random import randint


print("\n############################")
print("# Get Files and Parameters #")
print("############################\n")

dict_param = {param.split('=')[0]: param.split('=')[1] for param in sys.argv[1:]}
wkdir_path, sp1, sp2, intragenomic = dict_param['dir_path'], dict_param['SP1'], dict_param['SP2'], eval(dict_param['Intragenomic'])
gff1_motif, gff2_motif = dict_param['SP1_motif'], dict_param['SP2_motif']
data_file, control_file = glob.glob(wkdir_path + "/Synteny/Syntenic_blocks_STEP_3.txt")[0], dict_param['control_file']
chrom_corr_name = eval(eval(dict_param['name_corr'])) if dict_param['name_corr'] != '' else None
#chrom_corr_name=['Spartina_maritima_scaffold_', 'Chr']
gene_corr_name = eval(dict_param['corr_SB']) if dict_param['corr_SB'] != '' else None
corr_size = eval(dict_param['corr_size']) if dict_param['corr_size'] != '' else 500
outname = dict_param['outname']
print(f'outname: {outname}')
print(f'chrom_corr_name: {chrom_corr_name}, {type(chrom_corr_name)}')


#########
# Fasta #
#########

dict_chrom, tmp_chrom_order = {}, {sp1: [], sp2: []} # Chromosome length
for species, correction in [[sp1, chrom_corr_name[0]], [sp2, chrom_corr_name[1]]]:
    print('correction:', correction, type(correction))
    if species not in dict_chrom:
        dict_chrom[species] = {}
    fasta_file = wkdir_path + '/' + species + '.fasta' if os.path.exists(wkdir_path + '/' + species + '.fasta') else None
    print('fasta_file:', fasta_file)
    for name, seq in fct.buff_fas_reader(fasta_file):
        name2 = name.replace(correction, '') if correction else name
        dict_chrom[species][name2] = len(seq)
        if name2 not in tmp_chrom_order[species]:
            tmp_chrom_order[species].append(name2)
    print('len(dict_chrom[species]) =', len(dict_chrom[species]))


#######
# gff #
#######

dict_gff = fct.retrieve_gff_infos(wkdir_path=wkdir_path, SP=[sp1, sp2], Motif=[gff1_motif, gff2_motif], chrom_corr=chrom_corr_name, gene_corr=gene_corr_name, intragenomic=intragenomic)


################
# Control file #
################

chrom_order, dict_color, target_color = {sp1: None, sp2: None}, {}, None
inputfile = open(control_file)
for line in inputfile:
    if not line.startswith('#'):
        line = line.replace('\n', '').split('\t')

        # Chromosome order
        if line[1] == 'order':
            chrom_order[line[0]] = eval(line[2])

        # Color
        elif line[0] == 'target_color':
            target_color = line[1]
        elif line[0] == 'color':
            tmp_color = eval(line[1])
            for key in tmp_color:
                dict_color[key] = fct.rgb_to_hex(int(tmp_color[key].split(',')[0]), int(tmp_color[key].split(',')[1]), int(tmp_color[key].split(',')[2]))

inputfile.close()


####################
# Chromosome order #
####################

if not chrom_order[sp1]:
    chrom_order[sp1] = tmp_chrom_order[sp1]
if not chrom_order[sp2]:
    chrom_order[sp2] = tmp_chrom_order[sp2]


#########
# Color #
#########

if not dict_color:
    chrom_list = []
    inputfile = open(data_file)
    inputfile.readline()
    for line in inputfile:
        line = line.replace('\n', '').split('\t')
        chrom = line[1].replace(chrom_corr_name[0], '') if target_color == sp1 and chrom_corr_name else line[1] if target_color == sp1 else line[2].replace(chrom_corr_name[1], '') if target_color == sp2 and chrom_corr_name else line[2] if target_color == sp2 else None
        if chrom not in chrom_list:
            chrom_list.append(chrom)
    inputfile.close()

    list_hex = []
    for chrom in chrom_list:
        if chrom not in dict_color:
            color_temp = fct.rgb_to_hex(randint(0, 255), randint(0, 255), randint(0, 255))
            while color_temp in list_hex:
                color_temp = fct.rgb_to_hex(randint(0, 255), randint(0, 255), randint(0, 255))
            list_hex.append(color_temp)
            dict_color[chrom] = color_temp
print('dict_color:', dict_color)

##########
# Blocks #
##########

dict_blocks = {sp1: {}, sp2: {}}
inputfile = open(data_file)
inputfile.readline()
for line in inputfile:
    line = line.replace('\n', '').split('\t')
    #print(f'line: {line}')

    color = line[1].replace(chrom_corr_name[0], '') if target_color == sp1 and chrom_corr_name else line[1] if target_color == sp1 else line[2].replace(chrom_corr_name[1], '') if target_color == sp2 and chrom_corr_name else line[2] if target_color == sp2 else None
    chrom1 = line[1].replace(chrom_corr_name[0], '') if chrom_corr_name else line[1]
    chrom2 = line[2].replace(chrom_corr_name[1], '') if chrom_corr_name else line[2]

    for chrom, gene_list, species in [[chrom1, eval(line[7]), sp1], [chrom2, eval(line[9]), sp2]]:
        gene1, gene2 = gene_list[0], gene_list[-1]
        if gene1 not in dict_gff[species] and gene2 not in dict_gff[species]:
            raise ValueError("That's not expected! GFF Error!")

        # Retrieve SB position
        start = min([dict_gff[species][gene1][1], dict_gff[species][gene1][2], dict_gff[species][gene2][1], dict_gff[species][gene2][2]])
        end = max([dict_gff[species][gene1][1], dict_gff[species][gene1][2], dict_gff[species][gene2][1], dict_gff[species][gene2][2]])

        if chrom not in dict_blocks[species]:
            dict_blocks[species][chrom] = [[start, end, color]]
        else:
            dict_blocks[species][chrom].append([start, end, color])

inputfile.close()
print("len(dict_blocks[sp1]):", len(dict_blocks[sp1]), "\nlen(dict_blocks[sp2]):", len(dict_blocks[sp2]))


#############
# Karyotype #
#############

def draw_caryotype(t, corr_size):
    print('\nStart drawing!!')
    largeur, deplacement, rotation = 400, 1000, 90
    deplacement2, position, pos_range = 700, -20000, 0
    t.penup()

    ###########
    # Species #
    ###########
    for species in [sp1, sp2]:
        print(f"species: {species}")
        for chrom in chrom_order[species]:
            if type(chrom) is list:
                pos_range2 = pos_range
                for sub_chrom in chrom:
                    if sub_chrom in dict_chrom[species]:
                        print(f'chrom: {sub_chrom}')
                        longueur = int(dict_chrom[species][sub_chrom]) / corr_size

                        # Syntenic Blocks
                        if species not in dict_blocks:
                            raise KeyError(f'{species} is not defined in dict_block')
                        if sub_chrom not in dict_blocks[species]:
                            continue
                        fct.draw_SB(t=t, species=species, chrom=sub_chrom, dico_block=dict_blocks, dico_color=dict_color, largeur=largeur, rotation=rotation, pos=position, pos_range=pos_range2, corr_size=corr_size)
                        # Le contour du chromosome
                        t.goto(position, pos_range)
                        fct.rectangle(t=t, longueur=longueur, largeur=largeur, rotation=rotation, position=[position, pos_range2], color='black')
                        # Nom des chromosomes
                        fct.draw_chrom_name(t=t, pos=position, pos_range=pos_range2 - 500, color='black', chrom_name=sub_chrom)
                        pos_range2 += (longueur + deplacement2)
                    t.goto(position, pos_range2)
            else:
                if chrom in dict_chrom[species]:
                    print(f'chrom: {chrom}')
                    longueur = int(dict_chrom[species][chrom]) / corr_size

                    # Syntenic Blocks
                    if species not in dict_blocks:
                        raise KeyError(f'{species} is not defined in dict_block')
                    if chrom not in dict_blocks[species]:
                        continue
                    fct.draw_SB(t=t, species=species, chrom=chrom, dico_block=dict_blocks, dico_color=dict_color, largeur=largeur, rotation=rotation, pos=position, pos_range=pos_range, corr_size=corr_size)
                    # Le contour du chromosome
                    t.goto(position, pos_range)
                    fct.rectangle(t=t, longueur=longueur, largeur=largeur, rotation=rotation, position=[position, pos_range], color='black')
                    # Nom des chromosomes
                    fct.draw_chrom_name(t=t, pos=position, pos_range=pos_range - 500, color='black', chrom_name=chrom)
            position += deplacement
            t.goto(position, pos_range)
        t.penup()
        position = -20000
        pos_range = -10500


###################
# SVG output file #
###################

def write_file(draw_funct, filename, width, height):
    t = SvgTurtle(width, height)
    draw_funct(t, corr_size)
    t.save_as(filename)


def main():
    write_file(draw_caryotype, wkdir_path + '/Karyotype/' + outname + '.svg', 50000, 50000)
    print('Done.')


if __name__ == '__main__':
    main()
