
import pandas as pd
import pickle

import warnings
warnings.filterwarnings("ignore", message='not removing hydrogen atom without neighbors')


df = pickle.load(open("All_ord_reaction-data_nitpysrt.pkl", 'rb'))
#df = pickle.load(open("All_ord_reaction-data_nitpysrt_allsmiles_dotsep.pkl", 'rb'))


# In[54]:
print(df.columns)

from rdkit import Chem
def canonicalize(s):
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m)
    
fp=open("all_ord_reaction.tsv","w")
for i,r in df.iterrows():
    skip=False
    p=r["product_smiles"]
    if p =="NoData":
        continue
    try:
        p=canonicalize(p)
    except:
        continue
    rid=r["reaction_id"]
    arr=[]
    try:
        for j in range(len(df.columns)):
            rs="role"+str(j)
            ss="smiles"+str(j)
            if r[rs]=='REACTANT' and r[ss]!='NoData':
                smi=canonicalize(r[ss])
                arr.append(smi)
    except:
        continue
    if len(arr)>0:
        fp.write(rid+"\t"+p+"\t"+".".join(arr))
        fp.write("\n")
    


