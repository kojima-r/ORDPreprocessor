import pickle
import pandas as pd
with open('All_ord_reaction-data_nitpysrt.dict.pkl', 'rb') as f:
    reaction_dict=pickle.load(f)
reactions_df = pd.DataFrame(reaction_dict)
print(reactions_df.shape)
#pickle化
reactions_df.to_pickle('All_ord_reaction-data_nitpysrt.pkl')

