from  pymol import cmd
import pandas as pd
import numpy as np
mt_chains=["s","i", "j", "r", "k", "l", "m"]

cmd.load("5xtd.pdb")
cmd.spectrum()
cmd.select("nuc_sub", "not chain " + "+".join(mt_chains))
cmd.set("cartoon_transparency", 0.8, "nuc_sub")
# cmd.select("hetatm", "hetatm")
cmd.set("sphere_transparency", 0.8, "hetatm")
cmd.set("stick_transparency", 0.8, "hetatm")
cmd.select("membrane", "chain 6+7")
cmd.hide("sphere", "membrane")
cmd.show("nb_spheres", "membrane")

data=pd.read_csv("complexI.csv")

pos_dict=data.groupby("Chain")["Protein_position"].apply(list).to_dict()

for chain in np.unique(data["Chain"]):
    ress="+".join(str(v) for v in data["Protein_position"][data["Chain"]==chain].to_list())
    cmd.color("Green", "chain " + chain + " and resi " + ress)




