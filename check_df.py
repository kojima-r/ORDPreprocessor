
import pandas as pd
import pickle

import numpy as np
from rdkit import Chem

df = pickle.load(open("All_ord_reaction-data_nitpysrt_allsmiles_dotsep.pkl", 'rb'))
print(df[df['reaction_id']=="ord-0a0551b5bf3049c0a58006b54583054c"]["products0"])

