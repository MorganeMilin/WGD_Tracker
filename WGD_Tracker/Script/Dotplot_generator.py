"""
Dotplot - Step 2 - Dotplot generator
Author: Morgane MILIN
Last Update: July 2024
"""

##################
# Modules import #
##################

import sys
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import Dotplot_functions as fct


############################
# Get Files and Parameters #
############################

dict_param = {param.split('=')[0]: param.split('=')[1] for param in sys.argv[1:]}
wkdir_path, control_file = dict_param['dir_path'], dict_param['control_file']
data_file = glob.glob(wkdir_path + "/Dotplot/dotplot_dataset_file.txt")[0]
formatting_file = glob.glob(wkdir_path + "/Dotplot/dotplot_formatting_file.txt")[0]
markersize = int(dict_param['markersize']) if 'markersize' in dict_param else 1
ticks_size = dict_param['ticks_size'] if 'ticks_size' in dict_param else 10
label_size = dict_param['label_size'] if 'label_size' in dict_param else 15
x_lines = eval(dict_param['x_lines']) if 'x_lines' in dict_param else True
y_lines = eval(dict_param['y_lines']) if 'y_lines' in dict_param else True
outname = dict_param['outname'] if 'outname' in dict_param else 'Dotplot'
color_bar_label = dict_param['color_bar_label'] if 'color_bar_label' in dict_param else 'Ks'


####################
# Formatting infos #
####################

dict_axis = {}
input_formatting = open(formatting_file)
for line in input_formatting:
    line = line.replace('\n', '').split('\t')
    for axis in ['x axis', 'y axis']:
        if axis == line[0]:
            if axis not in dict_axis:
                dict_axis[axis] = {}
            for infos in ['title', 'ticks', 'label', 'lines']:
                if infos == line[1]:
                    if infos not in dict_axis[axis]:
                        dict_axis[axis][infos] = line[2]
                    else:
                        raise ValueError('Not expected !!!')
input_formatting.close()


#################################
# Retrieve data and Color infos #
#################################

fig_size, dict_color, list_else1, list_else2, list_color = None, {}, [], [], []     # dict_color = {limit: [color, [list1], [list2]], etc}
input_control = open(control_file)
for line in input_control:
    if line.startswith('x axis') or line.startswith('y axis') or line.startswith('#'):
        continue
    elif line.startswith('figsize='):
        fig_size = line.split('=')[1].split('\t')[0]
    else:
        line = line.replace('\n', '').split('\t')
        if line != ['']:
            if line[1] not in dict_color:
                dict_color[line[1]] = [line[0], [], []]
input_control.close()

inputfile = open(data_file)
for line in inputfile:
    line = line.replace('\n', '').split('\t')
    if dict_color:
        for limit in dict_color:
            if float(eval(limit)[0]) <= float(line[4]) < float(eval(limit)[1]):
                dict_color[limit][1].append(int(line[1]))
                dict_color[limit][2].append(int(line[3]))
    else:
        list_color.append(float(line[4]))
    list_else1.append(int(line[1]))
    list_else2.append(int(line[3]))
inputfile.close()


###########
# Dotplot #
###########

cmaps_perso = fct.colormaps_perso()
x = figsize=(int(eval(fig_size)[0]), int(eval(fig_size)[1])) if fig_size else (20, 15) if not dict_color else (15, 15)
plt.figure(figsize=x)
if dict_color:
    plt.plot(list_else1, list_else2, '.', color='lightgrey', markersize=markersize)
    for limit in sorted(dict_color):
        plt.plot(dict_color[limit][1], dict_color[limit][2], '.', color=dict_color[limit][0], markersize=markersize)
else:
    plt.scatter(x=list_else1, y=list_else2, c=list_color, marker='.', cmap=cmaps_perso, s=markersize)
    plt.colorbar(label=color_bar_label, orientation="vertical")
plt.xticks(eval(dict_axis['x axis']['ticks']), eval(dict_axis['x axis']['label']), rotation=45, fontsize=ticks_size)
plt.yticks(eval(dict_axis['y axis']['ticks']), eval(dict_axis['y axis']['label']), fontsize=ticks_size)
plt.xlabel(dict_axis['x axis']['title'], fontsize=label_size, style='italic')
plt.ylabel(dict_axis['y axis']['title'], fontsize=label_size, style='italic')
plt.axvline(0, color='black', linewidth=1)
plt.axhline(0, color='black', linewidth=1)
if x_lines:
    for i in eval(dict_axis['x axis']['lines']):
        plt.axvline(float(i), color='darkgrey', linewidth=1)
if y_lines:
    for i in eval(dict_axis['y axis']['lines']):
        plt.axhline(float(i), color='darkgrey', linewidth=1)
plt.margins(x=0.001, y=0.001)
plt.savefig(f"{wkdir_path}/Dotplot/{outname}.png")
plt.show()
