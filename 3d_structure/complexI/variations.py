from Bio.PDB import *
import numpy as np
import pandas as pd

parser=PDBParser()
structure=parser.get_structure('ComplexI', '5xtd.pdb')
chains=structure[0]

# just to keep in mind
# ppb=PPBuilder()
# pp_list=ppb.build_peptides(chains["A"])
# for pp in pp_list:
#     print(pp)

# pp=pp_list[0]

# function to establish if a residue is in the membrane, intermembrane or matrix space
# 14.7 and -14.7 are the coordinates of the planes that model the membranes
def find_location (coord):
    z=coord[2]
    if z>=-14.7 and z<=14.7:
        return "membrane"
    elif z>14.7:
        return "intermembrane"
    else:
        return "matrix"
    
# dictionary of mitochondrial protein with coordinates of CA atoms
mt_pp_dic={}

chain2protein = {
    "s": "MT-ND1",
    "i": "MT-ND2",
    "j": "MT-ND3",
    "r": "MT-ND4",
    "k": "MT-ND4L",
    "l": "MT-ND5",
    "m": "MT-ND6"
}

for chain_id in chain2protein.keys():
        protein_name = chain2protein[chain_id]
        mt_pp_dic[protein_name] = {}
        for res in chains[chain_id]:
                if res.id[0]==" ":
                        mt_pp_dic[protein_name][res.id[1]] = (res["CA"].coord, find_location(res["CA"].coord))

# read the df containing variations produced with the R script 
var_df=pd.read_csv("../../vcf file/mt_var_01")
# extract only variations related to complexI
cI_var=var_df[var_df["SYMBOL"].isin(chain2protein.values())]

locations=[]
for i, row in cI_var.iterrows():
    locations.append(mt_pp_dic[row["SYMBOL"]][row["Protein_position"]][1])
cI_var["Location"]=locations

# count number of aminoacid per location per protein
# for normalisation

dic_count={"MT-ND1": [0,0,0],
"MT-ND2": [0,0,0],
"MT-ND3":[0,0,0],
"MT-ND4":[0,0,0],
"MT-ND4L":[0,0,0],
"MT-ND5":[0,0,0],
"MT-ND6": [0,0,0]}

for prot_id in chain2protein.values(): 
    for values in mt_pp_dic[prot_id].values():
        if values[1]=="membrane":
            dic_count[prot_id][0]+=1
        elif values[1]=="matrix":
            dic_count[prot_id][1]+=1
        else:
            dic_count[prot_id][2]+=1

# rearrange df and add counts to the cI_var df in order to have 
# everything in one df

count_df=pd.DataFrame(dic_count).transpose()
count_df.columns=["Membrane", "Matrix", "Intermembrane"]
count_df.reset_index(inplace=True)
count_df.rename(columns={'index': 'SYMBOL'}, inplace=True)
# also add the chain for visualization in pymol
protein2chain = {v: k for k, v in chain2protein.items()}
protein2chain = pd.DataFrame(list(protein2chain.items()), columns=['SYMBOL', 'Chain'])
merged_df=pd.merge(cI_var, protein2chain, on='SYMBOL', how='left')
def_df = pd.merge(merged_df, count_df, on='SYMBOL', how='left')

def_df.to_csv("complexI.csv", index=False)
