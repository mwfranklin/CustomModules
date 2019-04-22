import os
import subprocess
import numpy as np
import pandas as pd
import glob

def ligands(folder_name, geom_id):
    ligands = np.zeros(4) #N, O, S, other
    ligand_lines = subprocess.check_output(["sed", "-e", "1,/TER/d", "%s/%s.pdb"%(folder_name, geom_id)])
    ligand_lines = ligand_lines.decode("utf-8").strip().split("\n")[1:]
    #print(ligand_lines)
    for line in ligand_lines:
        atom_id = line[12:16].strip()
        if "N" in atom_id:
            ligands[0] += 1
        elif "O" in atom_id:
            ligands[1] += 1
        elif "S" in atom_id:
            ligands[2] += 1
        else:
            ligands[3] += 1            
    #print(ligands)
    return ligands

def angles(rotated_matrix):
    these_angles = []
    #print(rotated_matrix)
    for x in range(0, len(rotated_matrix)):
        for y in range(x+1, len(rotated_matrix)):
            #print("new angles:")
            #print(rotated_matrix[x])
            #print(rotated_matrix[y])
            v1_u = rotated_matrix[x]/np.linalg.norm(rotated_matrix[x])
            v2_u = rotated_matrix[y]/np.linalg.norm(rotated_matrix[y])
            these_angles.append( np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))) )
    return(these_angles)
    
def angle_rmsd(folder_name, geom_id):
    #the [geom_id]_orig.pdb files contain the site that has already been rotated and aligned with the template
    #therefore, by calculating the angles in the same order, we can guarantee that the angles are corresponding
    #we can also grab the actual coords to use for calculating valence and nVESCUM later
    angle_rmsd = 0
    template_lines = subprocess.check_output(["sed", "-e", "/TER/,$d", "%s/%s_orig.pdb"%(folder_name, geom_id)])
    template_lines = template_lines.decode("utf-8").strip().split("\n")[1:] #fist line will have coords of [0,0,0]
    template_coords = []
    for entry in template_lines:
        template_coords.append([ float(entry[30:38].strip()), float(entry[38:46].strip()) , float(entry[46:54].strip())])
    template_coords = np.asarray(template_coords)
    template_angles = angles(template_coords)
    #print(template_angles)
    
    ligand_coords = []
    ligand_ids = []
    ligand_lines = subprocess.check_output(["sed", "-e", "1,/TER/d", "%s/%s_orig.pdb"%(folder_name, geom_id)])
    ligand_lines = ligand_lines.decode("utf-8").strip().split("\n")[1:] #fist line will have coords of [0,0,0]
    for entry in ligand_lines:
        ligand_coords.append([ float(entry[30:38].strip()), float(entry[38:46].strip()) , float(entry[46:54].strip())])
        ligand_ids.append(str(entry[12:16].strip()))
    ligand_coords = np.asarray(ligand_coords)
    ligand_angles = angles(ligand_coords)
    #print(ligand_angles)
    
    return np.sum( np.abs(np.asarray(template_angles) - np.asarray(ligand_angles)) ), np.max( np.abs(np.asarray(template_angles) - np.asarray(ligand_angles)) ), ligand_ids, ligand_coords #calculate the difference between each pair of corresponding angles and add up differences    

def get_orig_charge(metal, filename):
    try:
        this_charge = subprocess.check_output(["grep", "^FORMUL", filename])
        this_charge = this_charge.decode("utf-8").strip().split("\n")
        #print(this_charge)
        for line in this_charge: 
            #print(line)
            if metal in line:
                if line[19].isdigit():
                    charge = line.replace(")", "(").split("(")[1][3:].strip()
                else:
                    charge = line[21:25].strip()
                #print(charge)
                this_charge = int(charge[0])
                if charge[-1] == "-": this_charge * -1
                break
        return(this_charge)
    except subprocess.CalledProcessError:
        return(9)

def bond_valences(ligands, coords, metal, charge, bond_params):
    valence = 0
    distances = np.sqrt((coords*coords).sum(axis = 1))
    #print(distances)
    for x in range(0, len(ligands)):
        this_ligand = ligands[x]
        if "N" in this_ligand:
            this_ligand = "N"
        elif "O" in this_ligand:
            this_ligand = "O"
        elif "S" in this_ligand:
            this_ligand = "S"
        else:
            this_ligand = this_ligand      
        this_dist = distances[x]
        #print(this_ligand, this_dist)
        this_param = bond_params[(bond_params.Metal.str.upper() == metal) & (bond_params.Charge == charge) & (bond_params.CoordAtom == this_ligand)].values[0]
        r0 = this_param[4]
        b = this_param[5]
        #print(bond_params[(bond_params.Metal.str.upper() == metal) & (bond_params.Charge == charge) & (bond_params.CoordAtom == this_ligand)] )
        #print(r0, b)
        if r0 is np.nan: #if one is np.nan then both r0 and b will be undefined
            print("Check this metal:", metal, charge, ligands)
            return 1000
        else:
            #print (np.exp((r0 - this_dist)/b) )
            valence += np.exp((r0 - this_dist)/b)
            #print(valence)
    return(valence)

def vescum(coords, bond_valence):
    return(np.linalg.norm(coords.sum(axis = 0))/bond_valence)
    
def calc_cmm_params(metal_name, geom_name, pdb_id):
    good_metals = ["CU", "CO", "FE", "MN", "MG", "MO", "NI", "ZN"]
    if metal_name[0:2] not in good_metals:
        return(0)
        
    if geom_name != "Irr":
        gRMSD, max_dev, ligand_ids, ligand_coords = angle_rmsd(metal_name, geom_name)
        these_ligands = ligands(metal_name, geom_name)
        #print(ligand_ids, ligand_coords) 
    else:
        aligned_pdbs = glob.glob("%s/*.out"%metal_name)
        aligned_pdbs = [ x.split("/")[-1][:-4] for x in aligned_pdbs ]
        if aligned_pdbs[0] != "findgeo":
            new_geom = aligned_pdbs[0]
        else:
            new_geom = aligned_pdbs[1]
        #print("Fake geom use:", new_geom)
        gRMSD, max_dev, ligand_ids, ligand_coords = angle_rmsd(metal_name, new_geom)
        these_ligands = ligands(metal_name, new_geom)
        #print(ligand_ids, ligand_coords)
        gRMSD = 1000
                
    bond_params = pd.read_csv("/panfs/pfs.local/work/sluskylab/MSEAL/data/MetalParamsFiltered.txt", header = 0)
    this_charge = get_orig_charge(metal_name[0:2], "%s.pdb"%pdb_id)
    #print(this_charge, ligand_ids, ligand_coords)
    valence = bond_valences(ligand_ids, ligand_coords, metal_name[0:2], this_charge, bond_params)
    #print(valence)
    if valence != 1000:
        nVESCUM = vescum(ligand_coords, valence)
    else:
        nVESCUM = 1000
    #print(these_ligands, gRMSD, valence, nVESCUM)    
    cmm_params = these_ligands.tolist()
    cmm_params.extend([gRMSD, max_dev, valence, nVESCUM])
    #print(cmm_params)
    cmm_params = "\t".join(map(str, cmm_params))
    return(cmm_params)
    
#bond_params = pd.read_csv("MetalParamsFiltered.txt", header = 0)    
#gRMSD, ligand_ids, ligand_coords = angle_rmsd("SampleData/1muw_A/MG_453_5942_A", "spv")    
#this_charge = get_orig_charge("MG", "SampleData/1muw_A/1muw_A.pdb")
#valence = bond_valences(ligand_ids, ligand_coords, "MG", this_charge, bond_params)
#if valence != np.nan:
 #   nVESCUM = vescum(ligand_coords, valence)
#print(glob.glob("SampleData/1muw_A/MG_453_5942_A/*.out"))

#valence = bond_valences(ligand_ids, ligand_coords, "MG", bond_params)

#angle_rmsd("SampleData/1r1v_A/ZN_501_1533_A", "lin")
#ligands("SampleData/1muw_A/MG_454_5943_A", "tev")    
#ligands("SampleData/1muw_A/MG_454_6761_A", "trv")    
#ligands("SampleData/6rxn_A/FE_53_671_A", "tri")    