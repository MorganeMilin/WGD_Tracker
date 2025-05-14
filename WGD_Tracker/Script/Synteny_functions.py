from collections import Counter
import copy


def prep_data(chrom, nb, infos, dico):
    if chrom not in dico:
        dico[chrom] = {}
    if nb[0] not in dico[chrom]:
        dico[chrom][nb[0]] = {}
    if nb[1] not in dico[chrom][nb[0]]:
        dico[chrom][nb[0]][nb[1]] = infos
    else:
        if infos[3] > dico[chrom][nb[0]][nb[1]][3]:
            dico[chrom][nb[0]][nb[1]] = infos


def dotplot_extraction(dataset, outname, datatype):
    outfile = open(outname, 'w')
    infile = open(dataset)
    infile.readline()
    n,  stat = 0, {}
    for line in infile:
        line = line.replace('\n', '').split('\t')
        n += 1

        if int(line[5]) not in stat:
            stat[int(line[5])] = 1
        else:
            stat[int(line[5])] += 1

        if datatype == 'list':
            sp1_gene_list, sp2_gene_list, ks_list = eval(line[7]), eval(line[9]), eval(line[10])
        else:
            sp1_gene_list, sp2_gene_list, ks_list = [line[7]], [line[9]], [line[10]]

        for i in range(len(sp1_gene_list)):
            sp1, sp2, ks_value = sp1_gene_list[i], sp2_gene_list[i], ks_list[i]
            outfile.write('\t'.join([sp1, sp2, str(n), str(ks_value)]) + '\n')

    infile.close()
    outfile.close()
    print(f'{outname}\nNumber of Syntenic Blocks Detected: {n}.')


def rm_isolated_outlier_data(vec_sp2, gap_limit):
    vec_sp2_flt, start = [], False
    for i in range(len(vec_sp2)):
        ###############################################
        # 1st data and last data == no synteny at all #
        ###############################################
        if not start and i == len(vec_sp2)-1:
            continue
        ############
        # 1st data #
        ############
        if not start:
            mini, maxi = min(vec_sp2[i], vec_sp2[i+1]), max(vec_sp2[i], vec_sp2[i+1])
            if maxi - mini > gap_limit:
                continue
            vec_sp2_flt.append(vec_sp2[i])
            start = True
        #############
        # Last data #
        #############
        elif i == len(vec_sp2)-1:
            mini, maxi = min(vec_sp2[i], vec_sp2[i-1]), max(vec_sp2[i], vec_sp2[i-1])
            if maxi - mini > gap_limit:
                continue
            vec_sp2_flt.append(vec_sp2[i])
        ################
        # Dataset body #
        ################
        else:
            dist_inf, dist_sup, dist_ext = abs(vec_sp2[i] - vec_sp2[i-1]), abs(vec_sp2[i] - vec_sp2[i+1]), abs(vec_sp2[i-1] - vec_sp2[i+1])
            if dist_ext < dist_inf or dist_ext < dist_sup:
                continue
            if not vec_sp2[i-1] > vec_sp2[i] > vec_sp2[i+1] and not vec_sp2[i-1] < vec_sp2[i] < vec_sp2[i+1]:
                if dist_inf > gap_limit and dist_sup > gap_limit:
                    continue
            vec_sp2_flt.append(vec_sp2[i])
    return vec_sp2_flt


def generate_oriented_SB(chrom, vec_sp2_flt, vec_sp2_flt_sens):
    m = 0
    for i in range(len(vec_sp2_flt)):
        ############
        # 1st data #
        ############
        if i == 0:
            m += 1
            if chrom not in vec_sp2_flt_sens:
                vec_sp2_flt_sens[chrom] = {}
            if m not in vec_sp2_flt_sens[chrom]:
                vec_sp2_flt_sens[chrom][m] = [None, [vec_sp2_flt[i]]]
        #############
        # Last data #
        #############
        elif i == len(vec_sp2_flt) - 1:
            ref = vec_sp2_flt_sens[chrom][m][0]
            sens_inf = 'normal' if vec_sp2_flt[i] > vec_sp2_flt[i-1] else 'reverse' if vec_sp2_flt[i] < vec_sp2_flt[i-1] else None
            i_min = abs(vec_sp2_flt[i] - vec_sp2_flt[i - 1])
            if ref is None:  # Fusion, Break, initialisation
                vec_sp2_flt_sens[chrom][m][0] = sens_inf
                vec_sp2_flt_sens[chrom][m][1].append(vec_sp2_flt[i])
            elif ref == sens_inf:  # Fusion, Break, initialisation
                vec_sp2_flt_sens[chrom][m][1].append(vec_sp2_flt[i])
            else:   # Break, initialisation
                m += 1
                if (chrom, m) not in vec_sp2_flt_sens:
                    vec_sp2_flt_sens[chrom][m] = [None, [vec_sp2_flt[i]]]
        ################
        # Dataset body #
        ################
        else:
            ref = vec_sp2_flt_sens[chrom][m][0]
            sens_inf = 'normal' if vec_sp2_flt[i] > vec_sp2_flt[i-1] else 'reverse' if vec_sp2_flt[i] < vec_sp2_flt[i-1] else None
            sens_sup = 'normal' if vec_sp2_flt[i] < vec_sp2_flt[i+1] else 'reverse' if vec_sp2_flt[i] > vec_sp2_flt[i+1] else None
            if sens_inf == sens_sup:
                if ref is None:
                    vec_sp2_flt_sens[chrom][m][0] = sens_inf
                    vec_sp2_flt_sens[chrom][m][1].append(vec_sp2_flt[i])
                elif ref == sens_inf:
                    vec_sp2_flt_sens[chrom][m][1].append(vec_sp2_flt[i])
                else:
                    print('WTF !!!')
            else:
                i_min = abs(vec_sp2_flt[i] - vec_sp2_flt[i - 1])
                i_max = abs(vec_sp2_flt[i] - vec_sp2_flt[i + 1])
                statut = None
                if i_min < i_max: statut = 1
                elif i_min > i_max: statut = 2
                elif i_min == i_max and (ref is None or ref == sens_inf): statut = 1
                else: statut = 2

                if statut == 1:     # Fusion, Break, initialisation
                    if ref is None:
                        vec_sp2_flt_sens[chrom][m][0] = sens_inf
                    vec_sp2_flt_sens[chrom][m][1].append(vec_sp2_flt[i])
                    m += 1
                    if (chrom, m) not in vec_sp2_flt_sens:
                        vec_sp2_flt_sens[chrom][m] = [None, []]
                elif statut == 2:   # Break, initialisation
                    if vec_sp2_flt_sens[chrom][m][1]:
                        m += 1
                        if (chrom, m) not in vec_sp2_flt_sens:
                            vec_sp2_flt_sens[chrom][m] = [sens_sup, [vec_sp2_flt[i]]]
                    else:
                        vec_sp2_flt_sens[chrom][m][0] = sens_sup
                        vec_sp2_flt_sens[chrom][m][1].append(vec_sp2_flt[i])
                else:
                    raise ValueError("That's not possible !!! ERROR !!!")
    return vec_sp2_flt_sens


def SB_gap_only(vec_sp1, vec_sp2, sp1_gene, sp2_gene, Ks_vec, orientation, dico, key, gap_limit):
    m, statut_sp1, statut_sp2 = 0, False, False
    for i in range(len(vec_sp1)):

        ############
        # 1st data #
        ############
        if i == 0:
            m += 1
            if m not in dico[key]:
                dico[key][m] = [orientation, [vec_sp1[i]], [vec_sp2[i]], [sp1_gene[i]], [sp2_gene[i]], [Ks_vec[i]]]  # [orientation, [rank SP1], [rank SP2], [gene SP1], [gene SP2], [Ks]]

        ##############################
        # Dataset body and Last data #
        ##############################
        else:
            ref_sp1, ref_sp2 = dico[key][m][1][-1], dico[key][m][2][-1]
            statut_sp1 = True if ref_sp1 + gap_limit >= vec_sp1[i] else False

            if orientation == 'normal':
                statut_sp2 = True if ref_sp2 + gap_limit >= vec_sp2[i] else False
            else:
                statut_sp2 = True if vec_sp2[i] + gap_limit >= ref_sp2 else False

            if statut_sp1 and statut_sp2:
                dico[key][m][1].append(vec_sp1[i])
                dico[key][m][2].append(vec_sp2[i])
                dico[key][m][3].append(sp1_gene[i])
                dico[key][m][4].append(sp2_gene[i])
                dico[key][m][5].append(Ks_vec[i])
            else:
                m += 1
                dico[key][m] = [orientation, [vec_sp1[i]], [vec_sp2[i]], [sp1_gene[i]], [sp2_gene[i]], [Ks_vec[i]]]

