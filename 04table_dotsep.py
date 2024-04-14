#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import json
import pandas as pd
#rdkit関連のimport
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

# In[2]:
#ord_dataの読み込み
files = glob.glob("./ord_data_json/*")
print('Number of Files :', len(files))
p = Pool(64)
jsons=list(p.map(run, files))

print('FINISH')


# In[6]:


'''SMILES'''
def Pick_inputSMILESandRole(reactions,i):
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

def Pick_productSMILES(reactions,i):
    smiles = []
    smile = 'NoData'
    if 'products' not in reactions[i]['outcomes'][0]:
        print(reactions[i]['outcomes'][0])
        return smile
    for m in range(len(reactions[i]['outcomes'][0]['products'][0]['identifiers'])):
        check = reactions[i]['outcomes'][0]['products'][0]['identifiers'][m]['type']
        if check == 'SMILES':
            smile = reactions[i]['outcomes'][0]['products'][0]['identifiers'][m]['value']
        else:
            continue 
        smiles.append(smile)
    return smiles

def PicK_smiles(reactions, i):
    input_list = []
    role_list = []
    product_list = []
    input_smiles_list, input_role_list = Pick_inputSMILESandRole(reactions,i)
    for input_smiles, role in zip(input_smiles_list, input_role_list):
        smiles_list = Dot_separate(input_smiles)
        input_list.extend(smiles_list)
        for i in range(len(smiles_list)):
            role_list.append(role)
    product_smiles_list = Pick_productSMILES(reactions,i)
    for prod in product_smiles_list:
        prod_list = Dot_separate(prod)
        product_list.extend(prod_list)
    return input_list,role_list,product_list

'''SMARTS'''

def Dot_separate(smarts):
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


def  Reaction_separate(reactions, i):
    smarts = reactions[i]['identifiers'][0]['value']
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
        smarts_list = Dot_separate(smarts)
        input_list.extend(smarts_list)
    
    #product
    product_list = Dot_separate(product_smarts)
    
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
def Pick_Yield(reactions,i):
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

def Pick_temprature(reactions,i):
    try:
        temp = reactions[i]['conditions']['temperature']['setpoint']['value']
    except:
        temp = 'NaN'
    return temp

'''まとめ'''
def Pick_reaction(reactions,i):
    first_check = len(reactions[i]['inputs'])
    if first_check == 0:  #smarts形式で入っているデータを判定
        input_list, role_list, product_list = Reaction_separate(reactions,i) #smarts
    else :
        input_list, role_list, product_list = PicK_smiles(reactions, i) #smiles
    Yield = Pick_Yield(reactions,i)
    temp = Pick_temprature(reactions,i)
    reaction_id = None
    if 'reaction_id' in reactions[i]:
        reaction_id = reactions[i]['reaction_id'] 
    else:
        reaction_id = reactions[i]['reactionId'] 
    
    return reaction_id, input_list, role_list, temp, product_list, Yield


# In[7]:


#j_reactions = jsons[34]['reactions']
#print(jsons[34])
#Pick_reaction(j_reactions,0)


# In[8]:


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
        reaction_id, inputs, roles, temp, products, Yield= Pick_reaction(json_reactions,n)
        reaction_id_list.append(reaction_id)
        input_list.append(inputs)
        role_list.append(roles)
        temp_list.append(temp)
        product_list.append(products)
        yield_list.append(Yield)
    print('FINISH')


# In[9]:


input_df = pd.DataFrame({'smiles':input_list,
                         'role':role_list})
print('length of raction_list:',len(input_df.index))
print('length of temp_list:',len(temp_list))
print('length of product_list:',len(product_list))
print('length of yield_list:',len(yield_list))


# In[10]:


def Product_selection(product_list): #分子量の最も大きいものを選出
    product = product_list[0]
    if len(product_list) > 1:
        max_MW = 0
        for prod in product_list:
            try:
                mol = Chem.MolFromSmarts(prod)
                Chem.SanitizeMol(mol)
                prod_MW = rdMolDescriptors._CalcMolWt(mol)
                if prod_MW > max_MW:
                    product = prod
                    max_MW = prod_MW
            except:
                continue  
    return product


# In[11]:


input_dic = {}
empty = {'smiles':['NoData' for i in range(len(input_df.columns))],'tanimoto':['NoData' for i in range(len(input_df.columns))]}
empty_df = pd.DataFrame(empty)
for n in range(len(input_df.index)): #reactionごとに計算
    if n % 10000 == 0:
        print(n,'...',end='')
    smiles = input_df['smiles'][n] #１つのreactionに用いるinputのSMILES
    input_roles = input_df['role'][n]   #１つのreactionに用いるinputのrole
    reaction_df= pd.DataFrame({'smiles':smiles, 'roles':input_roles}) #DataFrame化
    #Tanimoto係数
    mols = []   #１つのreactionに用いるinputのMolオブジェクト
    for smile in smiles:
        if smile is None or smile == 'NoData':
            mols.append('NoData')
        else:
            mols.append(Chem.MolFromSmiles(smile))
    if len(product_list[n]) == 0:
        input_dic[n] = empty_df
        continue
    max_prod = Product_selection(product_list[n])
    prod = Chem.MolFromSmiles(max_prod)
    if prod is None:
        input_dic[n] = empty_df
        continue
    prod_fp = AllChem.GetMorganFingerprintAsBitVect(prod, 2, 2048)
    if 'NoData' in mols:
        morgan_fps = ['NoData' if mol == 'NoData' else AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in mols]
        false_num = [i for i, x in enumerate(morgan_fps) if x == 'NoData']
        tanimoto = [-1.0 if l in false_num else DataStructs.TanimotoSimilarity(prod_fp,morgan_fps[l]) for l in range(len(morgan_fps))]
    else:
        morgan_fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in mols]
        tanimoto = DataStructs.BulkTanimotoSimilarity(prod_fp, morgan_fps)
    reaction_df['tanimoto'] = tanimoto
    #sort
    sort_reaction_df = reaction_df.sort_values('tanimoto', ascending=False).reset_index(drop=True)
    input_dic[n] = sort_reaction_df


# In[12]:


#input
max_len = 0
max_num = 0
for n, key in enumerate(input_dic.keys()):
    list_len = len(input_dic[key])
    if max_len < list_len:
        max_len = list_len
print('max input number: ', max_len)


# In[13]:


#products
p_max_len = 0
p_max_num = 0
for n, prods in enumerate(product_list):
    list_len = len(prods)
    if p_max_len < list_len:
        p_max_len = list_len
print('max product number: ', p_max_len)


# In[14]:


#input
smiles_list = [ [] for i in range(max_len)]
roles_list = [ [] for i in range(max_len)]
tanimoto_list = [ [] for i in range(max_len)]
for key in input_dic.keys():
    if key % 10000 == 0:
        print(key,'...',end='')
    for i in range(max_len):
        try:
            smiles_list[i].append(input_dic[key]['smiles'][i])
        except:
            smiles_list[i].append('NoData')
        try:
            roles_list[i].append(input_dic[key]['roles'][i])
        except:
            roles_list[i].append('NoData')
        try:
            tanimoto_list[i].append(input_dic[key]['tanimoto'][i])
        except:
            tanimoto_list[i].append('NoData')


# In[15]:


#products
prods_list = [ [] for i in range(p_max_len)]
for n, prods in enumerate(product_list):
    if n % 10000 == 0:
        print(n,'...',end='')
    for i in range(p_max_len):
        try:
            prods_list[i].append(prods[i])
        except:
            prods_list[i].append('NoData')

print(len(smiles_list[0]))
print(len(prods_list[0]))
reaction_dict={'file_name': name_list,
                             'reaction_id': reaction_id_list,
                             'temp': temp_list,
                             'yield': yield_list}


for n in range(len(prods_list)):
    prods_name = 'products'+str(n)
    reaction_dict[prods_name] = prods_list[n]
for n in range(len(smiles_list)):
    smiles_name = 'smiles'+str(n)
    role_name = 'role'+str(n)
    tanimoto_name = 'tanimoto'+str(n)
    reaction_dict[smiles_name] = smiles_list[n]
    reaction_dict[role_name] = roles_list[n]
    reaction_dict[tanimoto_name] = tanimoto_list[n]


# In[19]:

import pickle
with open('All_ord_reaction-data_nitpysrt_allsimles_dotsep.dict.pkl', 'wb') as f:
    pickle.dump(reaction_dict, f)

reactions_df = pd.DataFrame(reaction_dict)
print(reactions_df.shape)
print(reactions_df)
#pickle化
reactions_df.to_pickle('All_ord_reaction-data_nitpysrt_allsmiles_dotsep.pkl')


