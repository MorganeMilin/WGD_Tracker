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


def retrieve_gff_infos(wkdir_path, SP, Motif, chrom_corr, gene_corr, intragenomic):
    import glob

    dico_gff = {SP[0]: {}, SP[1]: {}}
    for species, motif, chr_corr in [(SP[0], Motif[0], chrom_corr[0]), (SP[1], Motif[1], chrom_corr[1])]:

        gff_file = glob.glob(wkdir_path + "/" + species + "*.gff")[0]
        print('gff_file =', gff_file)

        # Complete dico_gff
        inputfile = open(gff_file)
        for line in inputfile:
            if not line.startswith('#'):
                line = line.replace('\n', '').split('\t')
                if line[2] == motif:
                    gene = gene_corr[1].join(line[8].split('=')[1].split(gene_corr[1])[:gene_corr[2]]) if gene_corr and gene_corr[0] == species else gene_corr[4].join(line[8].split('=')[1].split(gene_corr[4])[:gene_corr[5]]) if gene_corr and gene_corr[3] == species else line[8].split('=')[1]
                    chrom = line[0].replace(chr_corr, '') if chrom_corr else line[0]
                    start, end = int(line[3]), int(line[4])
                    if gene not in dico_gff[species]:
                        dico_gff[species][gene] = [chrom, start, end]
        inputfile.close()
        print(f'len(dico_gff) for {species} = {len(dico_gff[species])}')

        if intragenomic:
            break

    return dico_gff


def rgb_to_hex(r, g, b):
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)


def rectangle(t, longueur, largeur, rotation, position, color):
    """Fonction pour tracer un rectangle depuis le coin bas gauche"""
    t.penup()
    t.goto(position[0], position[1])
    t.pendown()
    t.pensize(10)
    t.pencolor(color)
    t.fd(largeur)
    t.rt(rotation)
    t.bk(longueur)
    t.rt(rotation)
    t.fd(largeur)
    t.rt(rotation)
    t.bk(longueur)
    t.rt(rotation)
    t.penup()


def rectangle_fill(t, longueur, largeur, rotation, position, color):
    """Fonction pour tracer un rectangle depuis le coin bas gauche"""
    t.penup()
    t.goto(position[0], position[1])
    t.pendown()
    t.pensize(1)
    t.fillcolor(color)
    t.begin_fill()
    t.pencolor(color)
    t.fd(largeur)
    t.rt(rotation)
    t.bk(longueur)
    t.rt(rotation)
    t.fd(largeur)
    t.rt(rotation)
    t.bk(longueur)
    t.rt(rotation)
    t.end_fill()
    t.penup()


def draw_SB(t, species, chrom, dico_block, dico_color, largeur, rotation, pos, pos_range, corr_size):
    for block in dico_block[species][chrom]:
        if block[2] in dico_color:
            SB_color = dico_color[block[2]]
        else:
            continue
        SB_longueur = (int(block[1]) - int(block[0])) / corr_size
        SB_position = int(block[0]) / corr_size
        rectangle_fill(t=t, longueur=SB_longueur, largeur=largeur, rotation=rotation, position=[pos, pos_range + SB_position], color=SB_color)


def draw_chrom_name(t, pos, pos_range, color, chrom_name):
    t.goto(pos, pos_range)
    t.pendown()
    t.pencolor(color)
    t.write(chrom_name, move=False, align='left', font=('Arial', 200, 'normal'))
    t.penup()
