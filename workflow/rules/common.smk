import pandas as pd

tsv = config['tsv']
df = pd.read_csv(tsv, sep='\t', header=None, names=['ID','R1','R2'])

SAMPLE_DICT = {}
for index, row in df.iterrows():
    ID = str(row['ID'])
    R1 = str(row['R1'])
    R2 = str(row['R2'])
    SAMPLE_DICT[ID] = {'R1': R1, 'R2': R2}
