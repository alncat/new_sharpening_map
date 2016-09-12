import pandas as pd

def read_data(filename):
    with open(filename, 'r') as f:
        counter = 0
        new_sharp = []
        sharp = []
        pdbs = []
        for line in f:
            if counter % 6 == 0:
                words = line.split()
                pdb_code = words[3].split('_')[1]
                pdbs.append(pdb_code.split('.')[0])
            if counter % 6 == 2:
                words = line.split(':')
                new_sharp.append(float(words[1]))
            if counter % 6 == 5:
                words = line.split(':')
                sharp.append(float(words[1]))
            counter += 1
    new_sharp = pd.Series(new_sharp)
    sharp = pd.Series(sharp)
    return new_sharp, sharp, pdbs
