#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pickle

df = pickle.load(open("All_ord_reaction-data_nitpysrt_allsmiles_dotsep.pkl", 'rb'))

import numpy as np
from rdkit import Chem
s_smi=[]
s_role=[]
s_prod=[]
for e in df.columns:
    if "smiles" in e:
        s_smi.append(e)
    if "role" in e:
        s_role.append(e)
    if "products" in e:
        s_prod.append(e)

def canonicalize(mol):
    mol = Chem.MolToSmiles(Chem.MolFromSmiles(mol),True)
    return mol


out_data_role=[]
out_data_unspecified_role=[]

for i,rr in  df.iterrows():
    in_smi=[rr[s] for s in s_smi]
    in_role=[rr[s] for s in s_role]
    prods=[rr[s] for s in s_prod]
    flag=False
    skip_flag=False
    reaction_id=rr['reaction_id']
    name=rr['file_name']
    prods_=[]
    #####
    out_info_list=[]
    if rr['yield']!="NoData":
        out_info_list.append(("YIELD","{:0.0f}".format(rr['yield'])))
    if rr['temp']!="NaN":
        out_info_list.append(("TEMP","{:0.1f}".format(rr['temp'])))
    #####
    try:
        prods_=[("PRODUCT",canonicalize(x)) for x in prods if x != 'NoData' and x != '']
    except:
        print("skip",reaction_id, "product:",prods)
        skip_flag=True
    #####
    if reaction_id=="ord-0a0551b5bf3049c0a58006b54583054c":
        print(prods)
        print(prods_)
    #####
    in_smi_=[]
    try:
        for x1,x2 in zip(in_smi,in_role):
            if x1 != 'NoData' and x1 != '':
                in_smi_.append((x2,canonicalize(x1)))
                #if (x2=="UNSPECIFIED") or (x2=="NoData"):
                #    flag=True
    except:
        print("skip",reaction_id,"input:", in_smi)
        skip_flag=True
    if not skip_flag:
        if flag:
            out_data_unspecified_role.append((name, reaction_id, in_smi_,prods_, out_info_list))
        else:
            out_data_role.append((name, reaction_id, in_smi_,prods_, out_info_list))
   
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

with open("all_ord_reaction_uniq_with_attr_v1.tsv","w") as fp:
    for name, reaction_id, in_smi_,prods_, out_info_list in out_data_role:
        in_str=".".join([x1+":"+x2 for x1,x2 in in_smi_]) 
        prods_str=".".join([x1+":"+x2 for x1,x2 in prods_]) 
        info_str=".".join([x1+":"+x2 for x1,x2 in out_info_list]) 
        fp.write("\t".join([reaction_id,in_str,prods_str, info_str]))
        fp.write("\n")


print("all_ord_reaction_uniq_with_attr_v1.tsv")
