#!/usr/bin/env python
# coding: utf-8

import glob
import json
import pickle
import pandas as pd

import rdkit
from rdkit import rdBase, Chem
from rdkit.Chem import PandasTools, AllChem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import rdMolDescriptors
RDLogger.DisableLog('rdApp.*')#エラーの非表示

from multiprocessing import Pool

def run(filename):
    print('{}...'.format(filename))
    json_open = open(filename)
    json_obj = json.load(json_open)
    return json_obj

#ord_dataの読み込み
files = glob.glob("./ord_data_json/*")
print('Number of Files :', len(files))
p = Pool(64)
jsons=list(p.map(run, files))

print('FINISH')


'''SMILES'''
def pick_inputSMILESandRole(reactions,i):
    smiles = []
    roles = []
    for key in reactions[i]['inputs'].keys():
        if 'components' in reactions[i]['inputs'][key]:
            for n in range(len(reactions[i]['inputs'][key]['components'])):
                check_keys = reactions[i]['inputs'][key]['components'][n].keys()
                role_check = False
                for role_key in  check_keys:
                    big_key = role_key.upper()
                    if 'ROLE' in big_key: #reaction_roleもしくはroleをヒット
                        roles.append(reactions[i]['inputs'][key]['components'][n][role_key])
                        role_check = True
                if role_check is False:
                    roles.append('NoData')
                smile = 'NoData'
                for m in range(len(reactions[i]['inputs'][key]['components'][n]['identifiers'])):
                    check = reactions[i]['inputs'][key]['components'][n]['identifiers'][m]['type']
                    if check == 'SMILES':
                        smile = reactions[i]['inputs'][key]['components'][n]['identifiers'][m]['value']
                    else:
                        continue 
                smiles.append(smile)
        else:
            print(reactions[i]['inputs'][key])
    return smiles, roles

def pick_productSMILES(reactions,r_idx):
    smiles = []
    smile = 'NoData'
    if 'products' not in reactions[r_idx]['outcomes'][0]:
        print(reactions[r_idx]['outcomes'][0])
        return [smile]
    for m in range(len(reactions[r_idx]['outcomes'][0]['products'][0]['identifiers'])):
        check = reactions[r_idx]['outcomes'][0]['products'][0]['identifiers'][m]['type']
        if check == 'SMILES':
            smile = reactions[r_idx]['outcomes'][0]['products'][0]['identifiers'][m]['value']
        else:
            continue 
        smiles.append(smile)
    return smiles

def pick_smiles(reactions, r_idx):
    input_list = []
    role_list = []
    product_list = []
    input_smiles_list, input_role_list = pick_inputSMILESandRole(reactions,r_idx)
    for input_smiles, role in zip(input_smiles_list, input_role_list):
        smiles_list = split_smarts(input_smiles)
        input_list.extend(smiles_list)
        for _ in range(len(smiles_list)):
            role_list.append(role)
    product_smiles_list = pick_productSMILES(reactions,r_idx)
    for prod in product_smiles_list:
        prod_list = split_smarts(prod)
        product_list.extend(prod_list)
    return input_list,role_list,product_list

'''SMARTS'''
def split_smarts(smarts):
    smarts_list = []
    dot_num = smarts.count('.')
    if dot_num > 0:
        for i in range(dot_num):
            dot_idx = smarts.find('.')
            smart = smarts[:dot_idx]
            smarts = smarts[dot_idx+1:]
            smarts_list.append(smart)
    smarts_list.append(smarts)
    return smarts_list


def split_reaction(reactions, r_idx):
    smarts = reactions[r_idx]['identifiers'][0]['value']
    input_smarts_list = []
    input_list = []
    if '>>' in smarts:
        darrow_idx = smarts.find('>>')
        input_smarts = smarts[:darrow_idx]
        product_smarts = smarts[darrow_idx+2:]    
        input_smarts_list.append(input_smarts)

    elif '>' in smarts:
        arrow_num = smarts.count('>')
        for i in range(arrow_num):
            arrow_idx = smarts.find('>')
            smart = smarts[:arrow_idx]
            smarts = smarts[arrow_idx+1:]
            input_smarts_list.append(smart)
        product_smarts = smarts
        
    #input
    for smarts in input_smarts_list:
        smarts_list = split_smarts(smarts)
        input_list.extend(smarts_list)
    
    #product
    product_list = split_smarts(product_smarts)
    
    """
    #input
    input_list = input_smarts_list
    #product
    product = product_smarts
    """
    #role
    input_num = len(input_list)
    roles = ['NoData']*input_num
    
    return input_list, roles, product_list

'''その他'''
def pick_Yield(reactions,i):
    Yield = 'NoData' 
    try:
        check = len(reactions[i]['outcomes'][0]['products'][0]['measurements'])
    except:
        check = 0
    if check != 0:
        for n in range(check):
            check2 = reactions[i]['outcomes'][0]['products'][0]['measurements'][n]['type']
            if check2 == 'YIELD':
                if 'percentage' in reactions[i]['outcomes'][0]['products'][0]['measurements'][n].keys():
                    if 'details' in reactions[i]['outcomes'][0]['products'][0]['measurements'][n]:
                        if reactions[i]['outcomes'][0]['products'][0]['measurements'][n]['details'] == 'CALCULATEDPERCENTYIELD':
                            continue
                    Yield = reactions[i]['outcomes'][0]['products'][0]['measurements'][n]['percentage']['value']
                else:
                    continue
            
            else:
                continue
    return Yield

def pick_temprature(reactions,i):
    try:
        temp = reactions[i]['conditions']['temperature']['setpoint']['value']
    except:
        temp = 'NaN'
    return temp

'''まとめ'''
def pick_reaction(reactions,i):
    first_check = len(reactions[i]['inputs'])
    if first_check == 0:  #smarts形式で入っているデータを判定
        input_list, role_list, product_list = split_reaction(reactions,i) #smarts
    else :
        input_list, role_list, product_list = pick_smiles(reactions, i) #smiles
    Yield = pick_Yield(reactions,i)
    temp = pick_temprature(reactions,i)
    reaction_id = None
    if 'reaction_id' in reactions[i]:
        reaction_id = reactions[i]['reaction_id'] 
    else:
        reaction_id = reactions[i]['reactionId'] 
    ###
    if reaction_id=="ord-0a0551b5bf3049c0a58006b54583054c":
        print(product_list)
    ###
    return reaction_id, input_list, role_list, temp, product_list, Yield


name_list = []
reaction_id_list = []
input_list = []
role_list = []
temp_list = []
product_list = []
yield_list = []

for num, json in enumerate(jsons):
    if 'name' not in json:
        file_name=""
    else:
        file_name = json['name']
    print('{}: {}...'.format(num, file_name), end='')
    json_reactions = json['reactions']
    for n in range(len(json_reactions)):
        #name
        name_list.append(file_name)
        
        #reactions
        reaction_id, inputs, roles, temp, products, Yield= pick_reaction(json_reactions,n)
        reaction_id_list.append(reaction_id)
        input_list.append(inputs)
        role_list.append(roles)
        temp_list.append(temp)
        product_list.append(products)
        yield_list.append(Yield)
    print('FINISH')


#input_df = pd.DataFrame({'smiles':input_list,
#                         'role':role_list})
print("... first step")
print('length of input_list:',len(input_list))
print('length of role_list:',len(role_list))
print('length of temp_list:',len(temp_list))
print('length of product_list:',len(product_list))
print('length of yield_list:',len(yield_list))

def sort_products(prod_list):
    prod_list_w=[]
    for el in prod_list:
        try:
            #mol = Chem.MolFromSmarts(el)
            mol = Chem.MolFromSmiles(el)
            Chem.SanitizeMol(mol)
            prod_MW = rdMolDescriptors._CalcMolWt(mol)
            prod_list_w.append((prod_MW,el))
        except:
            print(">>",el)
            pass
    new_prod_list=[smi for _, smi in sorted(prod_list_w,reverse=True)]
    return new_prod_list

#empty = {
#        'smiles':['NoData' for i in range(len(input_df.columns))],
#        'tanimoto':['NoData' for i in range(len(input_df.columns))]}
#empty_df = pd.DataFrame(empty)
new_prod_list=[]
new_input_list=[]
new_role_list=[]
new_tanimoto_list=[]
for n in range(len(input_list)): #reactionごとに計算
    if n % 10000 == 0:
        print(n,'...',end='')
    input_smiles = input_list[n] #１つのreactionに用いるinputのSMILES
    input_roles = role_list[n]   #１つのreactionに用いるinputのrole
    #Tanimoto係数
    mols = []   #１つのreactionに用いるinputのMolオブジェクト
    for smile in input_smiles:
        if smile is None or smile == 'NoData':
            mols.append('NoData')
        else:
            mols.append(Chem.MolFromSmiles(smile))
    new_prods=[]
    prod_fp=None
    if len(product_list[n]) > 0:
        new_prods = sort_products(product_list[n])
        if len(new_prods)>0:
            prod = Chem.MolFromSmiles(new_prods[0])
            if prod is not None:
                prod_fp = AllChem.GetMorganFingerprintAsBitVect(prod, 2, 2048)
        else:
            print("skip:",product_list[n],"=>",new_prods)
    if prod_fp is not None:
        if 'NoData' in mols:
            morgan_fps = ['NoData' if mol == 'NoData' else AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in mols]
            false_num = [i for i, x in enumerate(morgan_fps) if x == 'NoData']
            tanimoto = [-1.0 if l in false_num else DataStructs.TanimotoSimilarity(prod_fp,morgan_fps[l]) for l in range(len(morgan_fps))]
        else:
            morgan_fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in mols]
            tanimoto = DataStructs.BulkTanimotoSimilarity(prod_fp, morgan_fps)
        #reaction_df['tanimoto'] = tanimoto
        if not (len(input_smiles)==len(input_roles) and len(tanimoto)==len(input_smiles)):
            print("Alert1:",n)
        input_tuple_list=sorted([(a,b,c) for a,b,c in zip(tanimoto, input_smiles,input_roles)], reverse=True)
        ###
        new_input_list.append(   [b for a,b,c in input_tuple_list])
        new_role_list.append(    [c for a,b,c in input_tuple_list])
        new_tanimoto_list.append([a for a,b,c in input_tuple_list])
    else:
        if not (len(input_smiles)==len(input_roles)):
            print("Alert2:",n)
        new_input_list.append(input_smiles)
        new_role_list.append(input_roles)
        new_tanimoto_list.append(["NoData" for _ in range(len(input_smiles))])
    new_prod_list.append(new_prods)

    #sort_reaction_df = reaction_df.sort_values('tanimoto', ascending=False).reset_index(drop=True)
    #input_dic[n] = sort_reaction_df

print("... second step")
print('length of input_list:',len(new_input_list))
print('length of role_list:',len(new_role_list))
print('length of tanimoto_list:',len(new_tanimoto_list))
print('length of prod_list:',len(new_prod_list))

reaction_dict={'file_name': name_list,
            'reaction_id': reaction_id_list,
            'smiles': new_input_list,
            'role': new_role_list,
            'tanimoto': new_tanimoto_list,
            'products': new_prod_list,
            'temp': temp_list,
            'yield': yield_list}

with open('All_ord_reaction-data_nitpysrt_allsmiles_dotsep.dict0.pkl', 'wb') as f:
    pickle.dump(reaction_dict, f)


