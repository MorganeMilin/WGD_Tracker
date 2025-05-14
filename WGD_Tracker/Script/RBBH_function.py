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


def dico_coding_generator(gff, dico_coding, ft_type, species):
    # dico_coding = {Species: {Chr: [[start, end], [start, end], etc], etc}, etc}
    if species not in dico_coding:
        dico_coding[species] = {}

    inputfile_gff = open(gff)
    for line in inputfile_gff:
        if not line.startswith('#'):
            line = line.replace('\n', '').split('\t')

            if line[2] == ft_type:

                ###############
                # Chrom infos #
                ###############

                chrom = line[0]
                if chrom not in dico_coding[species]:
                    dico_coding[species][chrom] = []

                ########################
                # start & end position #
                ########################

                pos = [line[3], line[4]]
                if pos not in dico_coding[species][chrom]:
                    dico_coding[species][chrom].append(pos)
                else:
                    continue

    inputfile_gff.close()
    return dico_coding


def coding_check(species, col, col_start, col_end, limit, dico):
    overlap, non_coding = None, False
    for pos in dico[species][col]:
        # Si blast alignment est inclu completement dans un coding region
        if int(pos[0]) <= int(col_start) < int(col_end) <= int(pos[1]):
            overlap = True
            break
        # Si coding region completement inclu dans blast
        elif int(col_start) <= int(pos[0]) < int(pos[1]) <= int(col_end):
            overlap = True
            break
        # Si blast alignment et coding region sont chevauchant:
        elif int(pos[0]) <= int(col_start) < int(pos[1]) or int(pos[0]) < int(col_end) <= int(pos[1]):
            if int(pos[1]) - int(col_start) >= limit and int(col_end) - int(pos[0]) >= limit:
                overlap = True
                break
            else:
                overlap = False
        else:
            overlap = False
            non_coding = True

    return overlap, non_coding


def dico_genomic(dico_fa, datafile):
    print('datafile =', datafile)
    for name, seq in buff_fas_reader(datafile):
        print(name)
        dico_fa[name] = seq
    print("len(dico_fa) =", len(dico_fa))
    return dico_fa


def dico_CDS(fa_masked, gff, dico_TEs, limit, target):
    dico_len, name = {}, None
    inputfile_gff = open(gff)
    for line in inputfile_gff:
        if not line.startswith('#'):
            line = line.replace('\n', '').split('\t')

            chrom = line[0]
            if chrom not in dico_len:
                dico_len[chrom] = {}

            if line[2] == target:
                name = line[8].split('=')[1]
                if name not in dico_len[chrom]:
                    dico_len[chrom][name] = [[line[3], line[4]]]
                else:
                    dico_len[line[0]][name].append([line[3], line[4]])
    inputfile_gff.close()

    for name, seq in buff_fas_reader(fa_masked):
        if name in dico_len:
            for gene in dico_len[name]:
                total_len, N_count = 0, 0
                for position in dico_len[name][gene]:
                    start, end = int(position[0]), int(position[1])
                    interest_sequence = str(seq[start - 1:end].upper())
                    total_len += (end - start + 1)
                    N_count += int(interest_sequence.count('N'))
                res_TEs = round((float(N_count) / float(total_len)) * 100, 3)
                if res_TEs > limit:
                    dico_TEs[gene] = res_TEs

    return dico_TEs


def TE_recover(start, end, seq):
    interest_sequence = str(seq[start - 1:end].upper())
    total_len = (end - start + 1)
    N_count = int(interest_sequence.count('N'))
    res_TEs = round((float(N_count) / float(total_len)) * 100, 3)
    return res_TEs


def dict_genomic(datafile, col, start, end):
    # dico = {chrom: length_max, etc}
    query_dict = {}
    input_compil = open(datafile)
    for ligne in input_compil:
        ligne = ligne.replace('\n', '').split()
        if ligne[col] not in query_dict:
            query_dict[ligne[col]] = max(int(ligne[start]), int(ligne[end]))
        elif ligne[col] in query_dict:
            if max(int(ligne[start]), int(ligne[end])) > query_dict[ligne[col]]:
                query_dict[ligne[col]] = max(int(ligne[start]), int(ligne[end]))
    input_compil.close()
    return query_dict


def dict_CDS(datafile, col, start, end):
    # dico = {CDS: [[start, end], [start, end], etc], etc}
    query_dict = {}
    input_compil = open(datafile)
    for ligne in input_compil:
        ligne = ligne.replace('\n', '').split()
        if ligne[col] not in query_dict:
            query_dict[ligne[col]] = []
        if [ligne[start], ligne[end]] not in query_dict[ligne[col]]:
            query_dict[ligne[col]].append([ligne[start], ligne[end]])
    input_compil.close()
    return query_dict


def tmp_data_creation(input_data, col, key_ref):
    tmp_data = []
    inputfile_compil = open(input_data)
    for ligne in inputfile_compil:
        ligne = ligne.split()
        if ligne[col] == key_ref:
            tmp_data.append(ligne)
    inputfile_compil.close()
    return tmp_data


def read_through_predefined_windows(input_data, pos_min, pos_max, start, end):
    tmp_res = []
    for ligne in input_data:
        # Si blast alignment completement dans l'intervalle
        if pos_min <= int(ligne[start]) < int(ligne[end]) <= pos_max:
            tmp_res.append(ligne)
        # If overlapping, then expect at least 50% within this interval (will limit duplicates)
        elif pos_min <= int(ligne[start]) < pos_max < int(ligne[end]) <= pos_max:
            # If at least 50% is included in the interval
            if (pos_max - int(ligne[start]) + 1) / (int(ligne[end]) - int(ligne[start]) + 1) * 100 >= 50:
                tmp_res.append(ligne)
        elif int(ligne[start]) < pos_min < int(ligne[end]) <= pos_max:
            # If at least 50% is included in the interval
            if (int(ligne[end]) - pos_min + 1) / (int(ligne[end]) - int(ligne[start]) + 1) * 100 >= 50:
                tmp_res.append(ligne)
    tmp_res_sorted = sorted(tmp_res, key=lambda bit_score: float(bit_score[11]), reverse=True)
    return tmp_res_sorted


def best_hits_search(res, query_col, subject_col, start, end, limit):
    tmp_q, tmp_s, res_out_bh, nb = [], {}, [], 0
    for ligne in res:
        tempo_ligne = ligne

        if ligne[query_col] != tmp_q:
            tmp_q = ligne[query_col]
            tmp_s[ligne[subject_col]] = []
            tmp_s[ligne[subject_col]].append([min(int(ligne[start]), int(ligne[end])), max(int(ligne[start]), int(ligne[end]))])
            nb = 1
            tempo_ligne.extend([str(nb)])
            res_out_bh.append(tempo_ligne)
            tempo_ligne = tempo_ligne[:-1] 
            continue

        elif ligne[query_col] == tmp_q and ligne[subject_col] not in tmp_s and nb < int(limit):
            tmp_s[ligne[subject_col]] = []
            tmp_s[ligne[subject_col]].append([min(int(ligne[start]), int(ligne[end])), max(int(ligne[start]), int(ligne[end]))])
            nb += 1
            tempo_ligne.extend([str(nb)])
            res_out_bh.append(tempo_ligne)
            tempo_ligne = tempo_ligne[:-1] 
            continue

        elif ligne[query_col] == tmp_q and ligne[subject_col] in tmp_s and nb < int(limit):
            for tmp_key in range(len(tmp_s[ligne[subject_col]])):
                pos_mini = min(tmp_s[ligne[subject_col]][tmp_key][0], tmp_s[ligne[subject_col]][tmp_key][1])
                pos_maxi = max(tmp_s[ligne[subject_col]][tmp_key][0], tmp_s[ligne[subject_col]][tmp_key][1])
                if pos_mini <= int(ligne[start]) < pos_maxi or pos_mini < int(ligne[end]) <= pos_maxi:
                    continue
                else:
                    tmp_s[ligne[subject_col]].append([min(int(ligne[start]), int(ligne[end])), max(int(ligne[start]), int(ligne[end]))])
                    nb += 1
                    tempo_ligne.extend([str(nb)])
                    res_out_bh.append(tempo_ligne)
                    tempo_ligne = tempo_ligne[:-1]
                    break
            continue

        elif nb == int(limit):
            break
        continue
    return res_out_bh
