#from pymol import cmd, stored
import re
import sys
import os

suffix = "_TopDes_"

def create_output_states(name, folder, suffix, num_des, res_list = []):
    num_des = int(num_des)
    res_list = [str(x) for x in res_list]
    for x in range(1, num_des+1):
        x = (4-len(str(x)))*"0" + str(x)
        cmd.create("%s"%name, "*%s%s*"%(suffix, x), 1, x)
            
    cmd.delete("*%s*"%suffix)
    cmd.show_as("cartoon", "%s"%name)
    cmd.show("spheres", "het")
    cmd.show("sticks", "resn CRO")
    cmd.select("Imp_Res", "resi %s"%("+".join(res_list)))
    cmd.show("sticks", "Imp_Res")
    
    cmd.save("%s/%s.pse"%(folder, name), "%s"%name, 0, "pse")

for root, dirs, files in os.walk("GFP/DomInsertion/ShortDomain/RosRemodelTop3/Top3Relax/FinalRelax/D229L239_Best", topdown=False):
    print(root)
    name = root.split("/")
    #name[0] is name for create_output fxn
    #print(name)
    
    for value in files:
        if suffix in value:
            #print(value[:-4])
            print(value)
            cmd.load("%s/%s"%(root, value))
            pdb_num = int(value.split("_")[-1][:-4])
            if pdb_num > 1:
                #print(value, value[:-8])
                cmd.align("%s"%value[:-4], "%s_0001"%value[:-9])
    
    create_output_states(name[-1], "/".join(name[:-1]), suffix, 300, [9, 20, 22, 40, 42, 44, 51, 144])
    #cmd.delete("all")
    
    
    
    
