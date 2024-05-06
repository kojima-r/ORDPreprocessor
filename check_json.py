import json

#filename="ord_data_json/ord_dataset-e984b2d4813f44d59867e1771acc9b66.pb.gz.json"
#filename="ord_data_json/ord_dataset-06d4002fc4d34860a0688cba690e12dc.pb.gz.json"
filename="ord_data_json/ord_dataset-171b840ae6e84e45bab43b987d09f5c7.pb.gz.json"
print('{}...'.format(filename))
json_open = open(filename)
json_obj = json.load(json_open)
print(json_obj.keys())
rs=json_obj['reactions']
for i in range(len(rs)):
    #print(rs[i].keys())
    if 'reactionId' in rs[i]:
        rid = rs[i]['reactionId'] 
    elif 'reaction_id':
        rid = rs[i]['reaction_id'] 
    r = rs[i]
    #if rid =="ord-f4de5ee1862b463da5529745e2a17356":
    if rid =="ord-0a0551b5bf3049c0a58006b54583054c":
        print(r)
        print(r["outcomes"][0]["products"][0]["identifiers"][2]["value"])
#ord_dataの読み込み
#files = glob.glob("./ord_data_json/*")

