def split_path_and_file_name(inputfile):
    donnee_path, donnee_file = '', ''
    if '/' in inputfile:
        for i in range(len(inputfile[::-1])):
            if inputfile[::-1][i] == '/':
                donnee_path, donnee_file = inputfile[::-1][i:][::-1], inputfile[::-1][:i][::-1]
                break
    elif inputfile == 'None':
        donnee_path, donnee_file = None, None
    else:
        donnee_path, donnee_file = './', inputfile
    return donnee_path, donnee_file


def buff_fas_reader(file):
    current_seq, current_name = [], None
    file_iter = open(file)
    for ligne in file_iter:
        if not ligne.startswith('>'):
            current_seq.append(ligne.strip())
        else:
            if len(current_seq) != 0:
                yield current_name, ''.join(current_seq)
            # We retrieve everything that is after the '>', then we split with '\t' and keep the info in 1st position
            current_name = ligne[1:].strip().split()[0]
            current_seq = []
    file_iter.close()
    yield current_name, ''.join(current_seq)


def formatting_dotplot(fasta, list_chr, name_corr):
    total, dico_len, dico_chr, dico_chr2, list_lines = 0, {}, {}, {}, []
    for name, seq in buff_fas_reader(fasta):
        nom = name.replace(name_corr, '') if name_corr else name
        dico_len[nom] = len(seq)
    for chrom in list_chr:
        if chrom in dico_len:
            if chrom not in dico_chr:
                dico_chr2[chrom] = total
                total += dico_len[chrom]
                dico_chr[chrom] = total
        else:
            raise ValueError(f"Please check if your fasta file is correct.\n{chrom} is not present in it.")
    for key in dico_chr:
        list_lines.append(dico_chr[key])
    return dico_chr, dico_chr2, list_lines


def ticks_generator(chr_list, dico_data):
    ticks_list = []
    for i in range(0, len(chr_list)-1):
        if i == 0:
            start = 0
            end = dico_data[str(chr_list[i])]
            ticks_list.append((end - start) / 2 + start)
        start = dico_data[str(chr_list[i])]
        end = dico_data[str(chr_list[i+1])]
        ticks_list.append((end - start) / 2 + start)
    return ticks_list


def gff_infos(datafile, motif, name_corr):
    # dico_gff = {name: [chr, start, end]}
    gene, gene_start, gene_end, dico_gff = None, None, None, {}
    input_gff = open(datafile)
    for ligne in input_gff:
        if not ligne.startswith('#'):
            ligne = ligne.replace('\n', '').split('\t')
            if ligne[2] == motif:
                name, start, end = ligne[8].split('=')[1], ligne[3], ligne[4]
                chrom = ligne[0].replace(name_corr, '') if name_corr else ligne[0]
                if name not in dico_gff:
                    dico_gff[name] = [chrom, start, end]
                else:
                    raise ValueError('Error !!!')
    input_gff.close()
    return dico_gff


def colormaps_perso():
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap
    colors = {"pink": ["#ffc8ff", "#c80000"], "red": ["#c80000", "#7d007d"],
              "purple": ["#7d007d", "#6400ff"], "dark_blue": ["#6400ff", "#0000e1"],
              "light_blue": ["#0000e1", "#007dff"], "blue_green": ["#007dff", "#00c8c8"],
              "light_green": ["#00c8c8", "#00af00"], "dark_green": ["#00af00", "#006400"],
              "light_orange": ["#006400", "#e17d00"], "dark_orange": ["#e17d00", "#ff6400"],
              "light_brown": ["#ff6400", "#a07d3e"], "brown": ["#a07d3e", "#6d4600"],
              "darkbrown": ["#6d4600", "#441900"], "darkgrey": ["#441900", "#191919"],
              "black": ["#191919", "#000000"]}
    full_colormap_colors = []
    for key in colors:
        # On génère un dégradé linéaire entre les couleurs d'une gamme
        cmap = LinearSegmentedColormap.from_list("", colors[key])
        # Ajouter les dégradés de chaque game à notre colormap
        full_colormap_colors.extend(cmap(np.linspace(0, 1, 1000 // len(colors))))
    return LinearSegmentedColormap.from_list("custom_cmap", full_colormap_colors)
