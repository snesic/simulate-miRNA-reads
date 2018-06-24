# This script evaluates how miRNA aligners (STAR, microrazers, miraligner, razers and quagmir) perform on simulated reads.
# Simulated reads are initially miRNAs from the miRbase database. They are further edited by cutting bases or adding them from their corresponding hairpin on both 3p and 5p ends, adding random sequences (random gain) on both ends, inserting snips, etc.. The information on simulated reads has to be contained in read names in the following manner:
#hairpinName_mirnaName_sequenceStartInHairpin:sequenceEndInHairpin_gainOrLoss5p:gainOrLoss3p_mut:posAndBase_add:posAndBase
# Example: hsa-let-7a-1_hsa-let-7a-5p_5:29_2:0_mut:null_add:null

import pandas as pd
import re
import math
import numpy as np

#------Read fastq file and transform it to pandas DataFrame

def fastq_to_df(simFile):

    df = pd.read_csv(simFile, header=None)
    readsDF = pd.DataFrame({'NAME'     : df[0][df.index % 4 == 0 ].tolist(),
                            'SEQUENCE' : df[0][df.index % 4 == 1 ].tolist(),
                            'QUALS'    : df[0][df.index % 4 == 3 ].tolist()})

    readsDF.NAME = readsDF.NAME.str.replace('@','')
    return readsDF.drop_duplicates()


#------Extract info from read names into separate DataFrame columns

def extract_info_from_read_name(df):

    d = df.copy()
    d['NAME_MIR'] = d.NAME.apply(lambda x: x.split("_")[1])
    d['NAME_HP'] = d.NAME.apply(lambda x: x.split("_")[0])
    d['ADD'] = d.NAME.apply(lambda x: x.split("_")[5].replace('add:',''))
    d['MUT'] = d.NAME.apply(lambda x: x.split("_")[4].replace('mut:',''))
    d['LEN_READ'] = d.SEQUENCE.str.len()

    pos = d.NAME.apply(lambda x: x.split("_")[2].split(":")[0]) # pos of miRNA in Hairpin (without gain or loss)
    shift = d.NAME.apply(lambda x: x.split("_")[3].split(":")[0]) # substract shift
    d['POS'] = pos.astype('int64') - shift.astype('int64')

    add=d.NAME.apply(lambda x: x.split("_")[5])
    a = add[add.str.startswith('add:loss-5p')].str.replace('add:loss-5p-','').str.len() # number of 5'lost bases
    b = add[add.str.startswith('add:gain-5p')].str.replace('add:gain-5p-','').str.len() * (-1) # number of 5'gained bases
    d.POS = d.POS.add(a.add(b,fill_value=0), fill_value=0) # add a and subtract b

    return d

#------Evaluate aligners from a DataFrame prepaired in each aligner function below

def evaluate(df, col, tool): # col - which column to compare: hairpin or mirna name

    d = df.copy()
    d['TOOL'] = tool
    d['ALIGNED']  = 'yes'
    d['GOOD_ALIGNED']  = 'no'

    d.loc[d.MIRNA.isnull(), 'ALIGNED'] = 'NA'
    d.loc[d.MIRNA==d[col], 'GOOD_ALIGNED'] = 'yes'
    d.loc[d.START!=d.POS,'GOOD_ALIGNED'] = 'no'

    #Multi-aligned reads: Take the first one sort it and mark it as multi-yes or multi-no
    d = d.sort_values(['NAME', 'GOOD_ALIGNED'],ascending=[1,0])
    d.loc[d.NAME.duplicated(keep=False), 'ALIGNED'] = 'multi'
    d.loc[(~d.NAME.duplicated()).multiply(d.ALIGNED=='multi'), 'ALIGNED'] = 'multi-' + d.loc[(~d.NAME.duplicated()).multiply(d.ALIGNED=='multi'),'GOOD_ALIGNED']


    #columns = ['NAME', 'NAME_MIR', 'GOOD_ALIGNED', 'MIRNA', 'ADD', 'MUT', 'LEN_READ', 'QUALS', 'ALIGNED', 'TOOL']
    columns = ['NAME', 'ALIGNED', 'TOOL']

    return d.loc[~d.NAME.duplicated(),columns]


#---------Functions to prepare DataFrame for each aligner and evaluate it

def miraligner(simFile, resFile, tool='miraligner'):

    reads = fastq_to_df(simFile)

    results = pd.read_csv(resFile, sep='\t')
    results = results[['name', 'mir', 'start']]
    results.columns=['NAME', 'MIRNA', 'START']
    results.START=results.START.astype('int')-1

    allDF = pd.merge(reads, results, how='left', on=['NAME'])
    allDF = extract_info_from_read_name(allDF)

    return evaluate(allDF, 'NAME_MIR', tool)



def quagmir(simFile, resFile, tool='quagmir'):

    reads = fastq_to_df(simFile)

    results = pd.read_csv(resFile, sep='\t').reset_index()
    results = results[['MIRNA','SEQUENCE', 'READS', 'LEN_READ']]
    results = results.drop_duplicates()

    allDF = pd.merge(reads, results, how='left', on=['SEQUENCE'])
    allDF = extract_info_from_read_name(allDF)
    allDF['START'] = allDF.POS - 1

    return evaluate(allDF, 'NAME_MIR', tool)



def microrazers(simFile, resfile, tool='microrazers'):

    reads = fastq_to_df(simFile)

    results = pd.read_csv(resfile, header=None, sep='\t')
    results = results[[0,4,5]]
    results.columns = ['NAME', 'MIRNA', 'START']
    results = results.drop_duplicates()

    allDF = pd.merge(reads, results, how='left', on=['NAME'])
    allDF = extract_info_from_read_name(allDF)

    return evaluate(allDF,'NAME_HP', tool)


def sam(simFile, resFile, tool): # star and razers

    reads = fastq_to_df(simFile)

    results = pd.read_csv(resFile, sep='\t', comment='@', header=None)
    results = results[[0,2,4]]
    results.columns=['NAME', 'MIRNA', 'START']
    results = results.drop_duplicates()

    allDF = pd.merge(reads, results, how='left', on=['NAME'])
    allDF = extract_info_from_read_name(allDF)

    return evaluate(allDF, 'NAME_HP', tool)


#---------Read files and score algorithms----------------

files=pd.read_csv('mirna_files1.csv')

df=pd.DataFrame(columns = ['NAME', 'ALIGNED', 'TOOL'])


quag_dir = 'quagmir/'
micr_dir = 'microrazers/'
star_dir = 'star/'
mira_dir = 'miraligner/'
razr_dir = 'razers3/'

for col, row in files.iterrows():
    quag = quagmir(row.simReads, quag_dir+row.quagmir, 'quagmir;' + row.simReads.split('.')[0])
    mcr = microrazers(row.simReads, micr_dir+row.microrazers, 'microrazers;' + row.simReads.split('.')[0])
    star = sam(row.simReads, star_dir+row.star, 'star;'+ row.simReads.split('.')[0])
    mir = miraligner(row.simReads, mira_dir+row.miraligner, 'miraligner;'+ row.simReads.split('.')[0])
    razr = sam(row.simReads, razr_dir+row.razers3, 'razers3;'+ row.simReads.split('.')[0])
    df = pd.concat([df,quag,mcr,star,mir], ignore_index=True)


#------Write results in file-----------

df.to_csv('test2.tsv', index=False, sep='\t', na_rep='NA')
