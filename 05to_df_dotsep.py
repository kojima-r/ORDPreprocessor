
import pandas as pd
import pickle

import warnings
warnings.filterwarnings("ignore", message='not removing hydrogen atom without neighbors')


with open('All_ord_reaction-data_nitpysrt_allsmiles_dotsep.dict0.pkl', 'rb') as f:
    reaction_dict0=pickle.load(f)

name_list=reaction_dict0['file_name']
reaction_id_list=reaction_dict0['reaction_id']
new_input_list=reaction_dict0['smiles']
new_role_list=reaction_dict0['role']
new_tanimoto_list=reaction_dict0['tanimoto']
new_prod_list=reaction_dict0['products']
temp_list=reaction_dict0['temp']
yield_list=reaction_dict0['yield']


#input
max_len = max([len(el) for el in new_input_list])
print('max input number: ', max_len)

#products
p_max_len = max([len(el) for el in new_prod_list])
print('max product number: ', p_max_len)

def get_nth_from_list(n,target_list):
    out=[]
    for e in target_list:
        if n <len(e):
            out.append(e[n])
        else:
            out.append("NoData")
    return out

reaction_dict={'file_name': name_list,
            'reaction_id': reaction_id_list,
            'temp': temp_list,
            'yield': yield_list}


for n in range(p_max_len):
    prods_name = 'products'+str(n)
    reaction_dict[prods_name] = get_nth_from_list(n,new_prod_list)
for n in range(max_len):
    smiles_name = 'smiles'+str(n)
    role_name = 'role'+str(n)
    tanimoto_name = 'tanimoto'+str(n)
    reaction_dict[smiles_name] = get_nth_from_list(n,new_input_list)
    reaction_dict[role_name] = get_nth_from_list(n,new_role_list)
    reaction_dict[tanimoto_name] = get_nth_from_list(n,new_tanimoto_list)


with open('All_ord_reaction-data_nitpysrt_allsmiles_dotsep.dict.pkl', 'wb') as f:
    pickle.dump(reaction_dict, f)

reactions_df = pd.DataFrame(reaction_dict)
print(reactions_df.shape)
print(reactions_df)
#pickle化
with open('All_ord_reaction-data_nitpysrt_allsmiles_dotsep.pkl', 'wb') as f:
    pickle.dump(reactions_df, f)
