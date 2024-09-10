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


def dico_creation(infile, dico_data):
    for name, seq in buff_fas_reader(infile):
        if name in dico_data:
            dico_data[name] = seq.upper()
    return dico_data


