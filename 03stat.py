import json
name_count={}
name2_count={}

role_count={}
in_count={}
out_count={}
temp_count={}
yield_count={}
for line in open("ord.jsonl"):
    x=json.loads(line)
    if x["NAME"] not in name_count:
        name_count[x["NAME"]]=0
    name_count[x["NAME"]]+=1
    ###
    name2=x["NAME"]
    if "uspto" in name2:
        name2="uspto"
    if name2 not in name2_count:
        name2_count[name2]=0
    name2_count[name2]+=1
    ###
    for r in x["ROLE"]:
        if r not in role_count:
            role_count[r]=0
        role_count[r]+=1
    ### check: One record contains only one molecule
    for e in x["IN"]:
        if "." in e:
            print("error: multiple molecules")
    ##
    l=len(x["IN"])
    if l not in in_count:
        in_count[l]=0
    in_count[l]+=1
    l=len(x["PRODUCT"])
    if l not in out_count:
        out_count[l]=0
    out_count[l]+=1
    ###
    if x["TEMP"]!="NaN":
        if name2 not in temp_count:
            temp_count[name2]=0
        temp_count[name2]+=1
    if x["YIELD"]!="NoData":
        if name2 not in yield_count:
            yield_count[name2]=0
        yield_count[name2]+=1

print(name2_count)
print(role_count)
print(in_count)
print(out_count)

print(temp_count)
print(yield_count)

for k,v in name2_count.items():
    t=0
    if k in temp_count:
        t=temp_count[k]
    y=0
    if k in yield_count:
        y=yield_count[k]
    print(k,":",t,"/",v,"  ",y,"/",v)



