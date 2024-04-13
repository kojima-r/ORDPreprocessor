import ord_schema.message_helpers
import json
from google.protobuf.json_format import MessageToJson
from ord_schema.proto import dataset_pb2
import glob
import os
from multiprocessing import Pool

out_path="ord_data_json"
os.makedirs(out_path,exist_ok=True)
def run(filename):
    print(filename)
    data = ord_schema.message_helpers.load_message(filename, dataset_pb2.Dataset)
    json_s=MessageToJson(data)
    #o=json.loads(s)
    name=os.path.basename(filename)
    with open(out_path+"/"+name+".json","w") as fp:
        fp.write(json_s)

filename_list=[filename for filename in glob.glob("ord-data/data/**/*.*")]
p = Pool(16)
p.map(run, filename_list)


