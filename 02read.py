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
from multiprocessing import Pool


RDLogger.DisableLog('rdApp.*')#エラーの非表示

'''SMILES'''
def pick_input_role(reactions,i):
    smiles = []
    roles = []
    for key in reactions[i]['inputs'].keys():
        if 'components' not in reactions[i]['inputs'][key]:
            print("skip inputs:",reactions[i]['inputs'][key])
            continue
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
    return smiles, roles

def pick_product(reactions,i):
    smiles = []
    smile = 'NoData'
    for m in range(len(reactions[i]['outcomes'][0]['products'][0]['identifiers'])):
        check = reactions[i]['outcomes'][0]['products'][0]['identifiers'][m]['type']
        if check == 'SMILES':
            smile = reactions[i]['outcomes'][0]['products'][0]['identifiers'][m]['value']
        else:
            continue 
        smiles.append(smile)
    return smiles

def pick_smiles(reactions, i):
    input_list = []
    role_list = []
    product_list = []
    input_smiles_list, input_role_list = pick_input_role(reactions,i)
    for input_smiles, role in zip(input_smiles_list, input_role_list):
        smiles_list = dot_separate(input_smiles)
        input_list.extend(smiles_list)
        for i in range(len(smiles_list)):
            role_list.append(role)
    product_smiles_list = pick_product(reactions,i)
    for prod in product_smiles_list:
        prod_list = dot_separate(prod)
        product_list.extend(prod_list)
    return input_list,role_list,product_list

'''SMARTS'''

def dot_separate(smarts):
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


def reaction_separate(reactions, i):
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
        smarts_list = dot_separate(smarts)
        input_list.extend(smarts_list)
    
    #product
    product_list = dot_separate(product_smarts)
    
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
def pick_yield(reactions,i):
    Yield = 'NoData' 
    try:
        check = len(reactions[i]['outcomes'][0]['products'][0]['measurements'])
    except:
        check = 0
    if check != 0:
        for n in range(check):
            obj=reactions[i]['outcomes'][0]['products'][0]['measurements'][n]
            check2 = obj['type']
            if check2 == 'YIELD':
                if 'percentage' in obj.keys():
                    if "details" in obj and obj['details'] == 'CALCULATEDPERCENTYIELD':
                        continue
                    Yield = obj['percentage']['value']
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
        input_list, role_list, product_list = reaction_separate(reactions,i) #smarts
    else :
        input_list, role_list, product_list = pick_smiles(reactions, i) #smiles
    yield_val = pick_yield(reactions,i)
    temp = pick_temprature(reactions,i)
    reaction_id = reactions[i]['reactionId'] 
    
    return reaction_id, input_list, role_list, temp, product_list, yield_val


def run(filename):
    print('{}...'.format(filename))
    json_open = open(filename)
    json_obj = json.load(json_open)

    if 'name' not in json_obj:
        name=""
    else:
        name = json_obj['name']
    #print('{}: {}...'.format(num, name), end='')
    json_reactions = json_obj['reactions']
    out=[]
    for n in range(len(json_reactions)):
        #reactions
        reaction_id, inputs, roles, temp, products, yield_val= pick_reaction(json_reactions,n)
        obj = { "ID":reaction_id,
                "NAME":name,
                "IN":inputs,
                "ROLE":roles,
                "TEMP":temp,
                "PRODUCT":products,
                "YIELD":yield_val}
        out.append(obj)
    return out
    
#ord_dataの読み込み
files = glob.glob("./ord_data_json/*")

#jsonの読み込み
print('Number of Files :', len(files))


p = Pool(32)
result=p.map(run, files)

with open("ord.jsonl","w") as fp:
    for res in result:
        for obj in res:
            fp.write(json.dumps(obj))
            fp.write("\n")


