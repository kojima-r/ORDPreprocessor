import pandas as pd
import pickle

import warnings
warnings.filterwarnings("ignore", message='not removing hydrogen atom without neighbors')


df2=pd.read_csv("all_ord_reaction.tsv",sep="\t",header=None)

df3=df2.drop_duplicates(subset=[1,2])

df4=df3.sort_values(1)


fp=open("all_ord_reaction_uniq.tsv","w")
for i,r in df4.iterrows():
    fp.write(r[0]+"\t"+r[1]+"\t"+r[2]+"\n")

import numpy as np
with open("all_ord_reaction_uniq_prep.tsv","w") as fp:
    for i,r in df4.iterrows():
        len_out=len(r[1])
        len_in =len(r[2])
        if len_in < 256:
            if len_out >= 15 and len_out<256:
                fp.write(r[0]+"\t"+r[1]+"\t"+r[2]+"\n")



