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
    d['ADD'] = d.NAME.str.find('add:null')
    d['MUT'] = d.NAME.str.find('mut:null')
    d['LEN_READ'] = d.SEQUENCE.str.len()
    
    pos = d.NAME.apply(lambda x: x.split("_")[2].split(":")[0]).astype('int64') # pos of miRNA in Hairpin (without gain or loss)
    shift = d.NAME.apply(lambda x: x.split("_")[3].split(":")[0]).astype('int64') # substract shift
    d['POS'] = pos - shift
    
    add=d.NAME.apply(lambda x: x.split("_")[5])
    a = add[add.str.startswith('add:loss-5p')].str.replace('add:loss-5p-','').str.len() # number of 5'lost bases
    b = add[add.str.startswith('add:gain-5p')].str.replace('add:gain-5p-','').str.len() * (-1) # number of 5'gained bases
    d.POS = d.POS.add(a.add(b,fill_value=0), fill_value=0) # add a and subtract b
    
    return d

#------Evaluate aligners from a DataFrame prepaired in each aligner function below

def evaluate(df, col, tool): # col - which column to compare: hairpin or mirna name
    
    d = df.copy()
    d['TOOL'] = tool
    d['ALIGNED']  = 'NA'
    d['GOOD_ALIGNED']  = 'no'
    
    d.loc[d.MIRNA==d[col], 'ALIGNED'] = 'yes'
    d.loc[d.MIRNA==d[col], 'GOOD_ALIGNED'] = 'yes'
    d.loc[d.START!=d.POS, 'GOOD_ALIGNED'] = 'no'
    
    columns = ['NAME', 'NAME_MIR', 'GOOD_ALIGNED', 'MIRNA', 'ADD', 'MUT', 'LEN_READ', 'QUALS', 'ALIGNED', 'TOOL']
    
    return d[columns]


#---------Functions to prepare DataFrame for each aligner and evaluate it

def miraligner(simFile, resFile, tool='mirexpress'):
    
    readsDF = fastq_to_df(simFile)
    
    resultsDF = pd.read_csv(resFile, sep='\t')
    resultsDF = resultsDF[['name', 'mir','start']]
    resultsDF.columns=['NAME', 'MIRNA','START']
    resultsDF.START=resultsDF.START.astype('int') - 1
    
    allDF = pd.merge(readsDF, resultsDF, how='left', on=['NAME'])
    allDF = extract_info_from_read_name(allDF)
    
    return evaluate(allDF, 'NAME_MIR', tool)



def quagmir(simFile, resFile, tool='quagmir'):
    
    readsDF = fastq_to_df(simFile)
    
    resultsDF = pd.DataFrame.from_csv(resFile, sep='\t').reset_index()
    resultsDF = resultsDF[['MIRNA','SEQUENCE', 'READS', 'LEN_READ']]
    resultsDF=resultsDF.drop_duplicates()
    
    allDF = pd.merge(readsDF, resultsDF, how='left', on=['SEQUENCE'])
    allDF = extract_info_from_read_name(allDF)
    allDF['START'] = allDF.POS - 1
    
    return evaluate(allDF, 'NAME_MIR', tool)



def microrazers(simFile, resfile, tool='microrazers'):
    
    readsDF = fastq_to_df(simFile)
    
    resultsDF = pd.read_csv(resfile, header=None, sep='\t')
    resultsDF = resultsDF[[0,4,5]]
    resultsDF.columns = ['NAME', 'MIRNA','START']
    resultsDF = resultsDF.drop_duplicates()

    allDF = pd.merge(readsDF, resultsDF, how='left', on=['NAME'])
    allDF = extract_info_from_read_name(allDF)
    
    return evaluate(allDF,'NAME_HP', tool)


def sam(simFile, resFile, tool): # star and razers
    
    readsDF = fastq_to_df(simFile)
    
    resultsDF=pd.read_csv(resFile, sep='\t', comment='@', header=None)
    resultsDF = resultsDF[[0,2,4]]
    resultsDF.columns=['NAME', 'MIRNA','START']
    resultsDF = resultsDF.drop_duplicates()
    
    allDF = pd.merge(readsDF, resultsDF, how='left', on=['NAME'])
    allDF = extract_info_from_read_name(allDF)
    
    return evaluate(allDF,'NAME_HP', 'microrazers')


#---------Read files and score algorithms----------------

files=pd.read_csv('mirna_files1.csv')

df=pd.DataFrame(columns = ['NAME', 'NAME_MIR', 'GOOD_ALIGNED', 'MIRNA', 'ADD', 'MUT', 'LEN_READ', 'QUALS', 'ALIGNED', 'TOOL'])
df.ADD=df.ADD.astype('int')
df.MUT=df.MUT.astype('int')
df.LEN_READ=df.LEN_READ.astype('int')


quag_dir = 'quagmir/'
micr_dir = 'microrazers/'
star_dir = 'star/'
mira_dir = 'miraligner/'
for col, row in files.iterrows():
    quag = quagmir(row.simReads, quag_dir+row.quagmir, 'quagmir_' + row.simReads.split('.')[0])
    mcr = microrazers(row.simReads, micr_dir+row.microrazers, 'microrazers' + row.simReads.split('.')[0])
    star = sam(row.simReads, star_dir+row.star, 'star'+ row.simReads.split('.')[0])
    mir = miraligner(row.simReads, mira_dir+row.miraligner, 'miraligner'+ row.simReads.split('.')[0])
    df = pd.concat([df,quag,mcr,star,mir], ignore_index=True)


#------Write results in file-----------

df.to_csv('test2.tsv', index=False, sep='\t', na_rep='NA')





