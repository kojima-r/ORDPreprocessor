#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import pickle

import warnings
warnings.filterwarnings("ignore", message='not removing hydrogen atom without neighbors')


# In[3]:


df = pickle.load(open("All_ord_reaction-data_nitpysrt_allsmiles_dotsep.pkl", 'rb'))

#for r in df.iterrows():
r = df.iterrows()


# In[12]:


rr=next(r)


# In[20]:


prods=[rr[1]["products"+str(i)] for i in range(10)]
prods=[x for x in prods if x != 'NoData']


# In[76]:


import numpy as np
from rdkit import Chem

def canonicalize(mol):
    mol = Chem.MolToSmiles(Chem.MolFromSmiles(mol),True)
    return mol
in_smi=[rr[1]["smiles"+str(i)] for i in range(104)]
in_role=[rr[1]["role"+str(i)] for i in range(104)]
in_smi=[(x2,canonicalize(x1)) for x1,x2 in zip(in_smi,in_role) if x1 != 'NoData']
in_smi

in_smi_=sorted(list(set(in_smi)))
print(in_smi_)


def run(args):
    i,rr = args

    in_smi=[rr["smiles"+str(i)] for i in range(104)]
    in_role=[rr["role"+str(i)] for i in range(104)]
    prods=[rr["products"+str(i)] for i in range(10)]
    flag=False
    skip_flag=False
    reaction_id=rr['reaction_id']
    name=rr['file_name']
    
    
    out_info_list=[]
    if rr['yield']!="NoData":
        out_info_list.append(("YIELD","{:0.0f}".format(rr['yield'])))
    if rr['temp']!="NaN":
        out_info_list.append(("TEMP","{:0.1f}".format(rr['temp'])))
    try:
        prods_=[("PRODUCT",canonicalize(x)) for x in prods if x != 'NoData' and x != '']
    except:
        print("skip",reaction_id, "product:",prods)
        skip_flag=True
    in_smi_=[]
    try:
        for x1,x2 in zip(in_smi,in_role):
            if x1 != 'NoData' and x1 != '':
                in_smi_.append((x2,canonicalize(x1)))
                if (x2=="UNSPECIFIED") or (x2=="NoData"):
                    flag=True
    except:
        print("skip",reaction_id,"input:", in_smi)
        skip_flag=True
    a=None
    b=None
    if not skip_flag:
        in_smi_=sorted(list(set(in_smi_)))
        if flag:
            #out_data_unspecified_role.append((name, reaction_id, in_smi_,prods_, out_info_list))
            a=(name, reaction_id, in_smi_,prods_, out_info_list)
        else:
            #out_data_role.append((name, reaction_id, in_smi_,prods_, out_info_list))
            b=(name, reaction_id, in_smi_,prods_, out_info_list)
    return a,b
    #except:
    #    print("skip",reaction_id)
    #    pass

from multiprocessing import Pool
out_data_role=[]
out_data_unspecified_role=[]

args_v=[(i,rr) for i,rr in  df.iterrows()]
p = Pool(128)
results=p.map(run, args_v)
for a,b in results:
    if a is not None:
        out_data_unspecified_role.append(a)
    if b is not None:
        out_data_role.append(b)

#print(en(out_data_unspecified_role),len(out_data_role))


# In[75]:


in_set=set()
prods_set=set()
info_set=set()
for name, reaction_id, in_smi_,prods_, out_info_list in out_data_role:
        in_set|=set([x1 for x1,x2 in in_smi_]) 
        prods_set|=set([x1 for x1,x2 in prods_]) 
        info_set|=set([x1 for x1,x2 in out_info_list]) 
print(in_set)
print(prods_set)
print(info_set)


# In[74]:


with open("ord_simple_data.tsv","w") as fp:
    for name, reaction_id, in_smi_,prods_, out_info_list in out_data_role:
        in_str=".".join([x1+":"+x2 for x1,x2 in in_smi_]) 
        prods_str=".".join([x1+":"+x2 for x1,x2 in prods_]) 
        info_str=".".join([x1+":"+x2 for x1,x2 in out_info_list]) 
        fp.write("\t".join([reaction_id,in_str,prods_str, info_str]))
        fp.write("\n")


# In[32]:


import re
data=[]
for line in open("ord_simple_data.tsv","r"):
#for line in open("all_ord_reaction_uniq_with_attr_v1.tsv","r"):
    arr=line.strip("\n").split("\t")
    reaction_id,in_str,prods_str, info_str=arr[0],arr[1],arr[2],arr[3]
    v=[]
    if (len(in_str)>0):
        pairs=[v.split(":") for v in in_str.split(".")]
        v.extend(pairs)
    if (len(prods_str)>0):
        pairs=[v.split(":") for v in prods_str.split(".")]
        v.extend(pairs)
    if (len(info_str)>0):
        pairs=[v.split(":") for v in re.split('\.[TY]', info_str)]
        new_pairs=[]
        for pair in pairs:
            if pair[0]=="EMP":
                new_pairs.append(("TEMP",pair[1]))
            elif pair[0]=="IELD":
                new_pairs.append(("YIELD",pair[1]))
            else:
                new_pairs.append(pair)
        v.extend(new_pairs)
    ##to dictionary
    reaction_dict={}
    for pair in v:
        if pair[0] not in reaction_dict:
            reaction_dict[pair[0]]=[]
        reaction_dict[pair[0]].append(pair[1])
    data.append((reaction_id,reaction_dict))


# In[35]:


roles=set()
for reaction_id,inp in data:
    for pair in inp.items():
        roles.add(pair[0])
roles


# In[43]:


r_list=list(roles)
fp=open("ord_simple_data_v2.tsv","w")
for reaction_id,reaction_dict in data:
    reaction_dict_out={}
    line=[]
    for r in r_list:
        if r in reaction_dict:
            reaction_dict_out[r]=".".join(reaction_dict[r])
        else:
            reaction_dict_out[r]="None"
        line.append(r+":"+reaction_dict_out[r])
    fp.write("\t".join(line))
    fp.write("\n")
        

# In[ ]:


