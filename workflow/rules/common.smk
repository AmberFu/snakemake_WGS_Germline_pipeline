import os
import glob
import pandas as pd
from pathlib import Path

###################
### Raw BAMs:
###################

def file_extension(s: str) -> str:
    if s == ".bam":
        return "BAM"
    elif s == ".cram":
        return "CRAM"
    else:
        return "OTHERS"

### Get input BAM/CRAM full path as a list: (independent list)
BAMCRAM_FILESPATH = glob.glob(config["bam_cram_folder"]+"/UDN*am")
BAMCRAM_FILENAME = list(map(os.path.basename, BAMCRAM_FILESPATH))
BAMCRAM_STEMNAME = []
SAMPLE_ID = []
BAM_CRAM = []
for f in BAMCRAM_FILENAME:
    filePathObject =  Path(f)
    f_stem = filePathObject.stem
    f_extension = filePathObject.suffix
    BAMCRAM_STEMNAME.append(f_stem)
    SAMPLE_ID.append(f_stem.split("-")[0])
    BAM_CRAM.append(file_extension(f_extension))

### Tabular configuration
# https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#tabular-configuration
info_dict = {'BAMCRAM_FILENAME': BAMCRAM_FILENAME, \
             'BAMCRAM_STEMNAME': BAMCRAM_STEMNAME, \
             'SAMPLE_ID': SAMPLE_ID, \
             'BAM_CRAM': BAM_CRAM}  
df_BAM = pd.DataFrame(info_dict).set_index(["BAMCRAM_FILENAME"], drop=False).sort_index()

###################
### Read-in CSV:
###################
sample_fq_csv = config["sample_fq_csv"]
df_Sample_Family = pd.read_csv(sample_fq_csv, header=None, 
                               names=['udnid','family_id','fq_prefix'], 
                               usecols=[0,3,4]).set_index(
    ["udnid"], drop=False).sort_index()
## Add Prefix:
df_Sample_Family['prefix'] = df_Sample_Family.apply(lambda x: "WGS_" + x['family_id'] + "_" + x['udnid'], axis=1)


########################
### Merged: 124 samples
########################
SAMPLES = df_Sample_Family.merge(df_BAM, how='inner', left_on='fq_prefix', right_on='BAMCRAM_STEMNAME')

### Prepare a Dict: {"prifix": {"r1": "xxxx", "R2": "oooo"}}
prefix_dict = {}
for index, row in SAMPLES.iterrows():
    prefix_dict[row.BAMCRAM_STEMNAME] = row.prefix

