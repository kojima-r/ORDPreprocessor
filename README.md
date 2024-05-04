# ORDPreprocessor
#### 00download
- output: `ord-data/`

#### 01build_json
- input: `ord-data/data/`
- output: `ord_data_json/`

#### 02read
- input: `ord_data_json/`
- output: `ord.jsonl`

`ord.jsonl`は以下のキーの値をもつ
  - "ID":reaction_id 反応ID
  - "NAME":name DB名
  - "IN":inputs　入力物質のリスト
  - "ROLE":roles　役割のリスト
  - "TEMP":temp　温度
  - "PRODUCT":products　生成物のリスト
  - "YIELD":yield_val　収率
                
#### 03stat
- input: `ord.jsonl`
統計値の計算

DB名の統計，役割の統計，入力物質の数と生成物の数の統計，温度の統計，収率の統計

#### 04table
- input: `ord_data_json/`
- output: `All_ord_reaction-data_nitpysrt.dict.pkl`
04table_to_df.py
- output: `All_ord_reaction-data_nitpysrt.pkl`

出力されるdictは以下のキーの値をもつ(pandasに入力するため，キーは列名に対応する)

- 'file_name': name_list
- 'reaction_id': reaction_id_list
- 'temp': temp_list
- 'product_smiles': product_list
- 'yield': yield_list
- 'smiles1'
- 'smiles2'

...

- 'role1'
- 'role2'

...

- 'tabuniti1'
- 'tabuniti2'

...


#### 04table_dotsep
- input: `ord_data_json/`
- output: `All_ord_reaction-data_nitpysrt_allsimles_dotsep.dict.pkl`
- output: `All_ord_reaction-data_nitpysrt_allsmiles_dotsep.pkl`

##### 05to_tsv
- input: `All_ord_reaction-data_nitpysrt_allsmiles_dotsep.pkl`
- output: `ord_simple_data.tsv`
- output: `ord_simple_data_v2.tsv`

##### 05to_tsv_with_attr
- input: `All_ord_reaction-data_nitpysrt_allsmiles_dotsep.pkl`
- output: `all_ord_reaction_uniq_with_attr_v1.tsv`

> all_ord_reaction_uniq_with_attr202404_v1.tsv

##### 06to_tsv2
- input: `All_ord_reaction-data_nitpysrt.pkl`
- output: `all_ord_reaction.tsv`

##### 07to_tsv3
-input: `all_ord_reaction.tsv`
-output: `all_ord_reaction_uniq.tsv"`
-output: `all_ord_reaction_uniq_prep.tsv`


