
import pandas as pd
import pickle

import numpy as np
from rdkit import Chem

obj = pickle.load(open("All_ord_reaction-data_nitpysrt_allsmiles_dotsep.dict0.pkl", 'rb'))
n=-1
print(obj.keys())
for k in obj.keys():
    print(k,len(obj[k]))
for i,rid in enumerate(obj["reaction_id"]):
    if rid=="ord-0a0551b5bf3049c0a58006b54583054c":
        n=i

print(obj["products"][0])
print(obj["products"][n-1])
print(obj["products"][n])
print(obj["products"][n+1])
#print(df[df['reaction_id']=="ord-0a0551b5bf3049c0a58006b54583054c"]["products0"])

