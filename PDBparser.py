import numpy as np
import os
import glob
import decimal
import subprocess
import scipy.spatial
import math
import itertools

aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE", "SEC"]
oneletAA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "M", "C"]
noncanAA = ["MSE", "IAS", "SEC"] #selenomethionene, L-asp acid for crosslinking, selenocysteine
metals = ["CUA", "CU", "FE", "MG", "ZN"]

class Protein:
    def __init__(self, res_list, res_nums):
        self.sequence = ""
        self.NumRes = len(res_list)
        self.num_atoms = 0
        self.Coords = []
        self.Coords_index = []
        self.metals = []
        self.residues = res_list
        self.res_nums = []
        
        for x in res_list:
            self.sequence += x.onelet
            self.num_atoms += len(x.Atoms)
            self.Coords.extend(x.Coords)
            new_resnum = str(x.resnum) + x.chain
            self.res_nums.append(new_resnum)
            if x.type == "metal": self.metals.append(new_resnum)
            for atom in x.Atoms:
                self.Coords_index.append(new_resnum)
        self.Coords = np.asarray(self.Coords)
        self.KDTree = scipy.spatial.KDTree(self.Coords)
        
    def get_neighbor_res(self, atom_coords, r):
        neighbor_atoms = self.KDTree.query_ball_point(atom_coords, r)
        neighbor_atoms = [y for x in neighbor_atoms for y in x]
        neighbor_res = [self.Coords_index[x] for x in neighbor_atoms]
        neighbor_res = sorted(set(neighbor_res))
        return neighbor_res
        
    def get_metal_shells(self, metals, cutoff1, cutoff2):
        print("There are this many metals:", len(metals), metals)
        metal_first_shell = [[] for x in metals]
        for x in range(0, len(metals)):
            first_shell = self.get_neighbor_res(self.residues[self.res_nums.index(metals[x])].Coords, cutoff1)
            #print(metals[x], first_shell)
            metal_first_shell[x] = first_shell #sorted(set(first_shell).difference([metals[x]])) #don't actually want to remove metal for the mononuclear cases
            #print(metals[x], metal_first_shell[x])
        #print(metal_first_shell)
        
        metal_sites = []
        for i, x in enumerate(metal_first_shell):
            #print("\n", i, x)
            for y in range(len(metals)):
                #check if metals are w/in 3A, then check if any res are in common between metals
                if metals[y] in x: 
                    metal_first_shell[y].extend(x)
                    metal_first_shell[i].extend(metal_first_shell[y])
                elif len(set(x).intersection(metal_first_shell[y])) != 0:
                    metal_first_shell[y].extend(x)
                    metal_first_shell[i].extend(metal_first_shell[y])
            metal_first_shell = [sorted(set(x)) for x in metal_first_shell]
            #print(metal_first_shell)
        metal_first_shell.sort()
        #print("finished first shell")
        #print(metal_first_shell)
        metal_first_shell = list(metal_first_shell for metal_first_shell, _ in itertools.groupby(metal_first_shell))
        for value in metal_first_shell:
            #print(value)
            metal_sites.append([sorted(set(metals).intersection(value)) , value])
        
        #print("Metal ligands are: ", metal_sites)
        for site in metal_sites:
            metal = site[0]
            site_res = site[1]
            #print(metal, site_res)
            site_coords = [self.residues[self.res_nums.index(x)].Coords for x in site_res]
            site_coords = [y for x in site_coords for y in x]
            site_coords = np.asarray(site_coords)    
            new_res = self.get_neighbor_res(site_coords, cutoff2)
            site[1] = new_res
            #pymol_res = sorted([int(x[:-1]) for x in new_res])
            #print("+".join(map(str, pymol_res)))
        return metal_sites #,pymol_res
                    
class Residue:
    def __init__(self, resnum, name, chain, atom_list, res_index):
        #print("atomlist is ", len(atom_list))
        self.name = name
        self.resnum = resnum #from PDB
        self.res_index = res_index #total count
        if self.name in metals:
            self.type = "metal"
        else:
            self.type = "protein"
        try:
            self.onelet = oneletAA[aa.index(self.name)]
        except:
            self.onelet = "X"
        self.chain = chain
        self.Atoms = []
        self.Coords = np.zeros((len(atom_list), 3))
        self.num_atoms = len(atom_list)
        self.psi = 500.00
        self.phi = 500.00
        self.chi1 = 500.00   
        for x in range(0, len(atom_list)):
            self.Atoms.append(atom_list[x][0])
            self.Coords[x][0] = float(atom_list[x][1])      
            self.Coords[x][1] = float(atom_list[x][2])
            self.Coords[x][2] = float(atom_list[x][3])

    def set_phi(self, prev_res):
        #phi is prevC-N-C-C
        if prev_res.resnum + 1 == self.resnum and self.name not in metals:
            self.phi = dihedral(prev_res.Coords[2], self.Coords[0], self.Coords[1], self.Coords[2])
        
    def set_psi(self, next_res):
        #psi is N-c-c-nextN
        if self.resnum + 1 == next_res.resnum and self.name not in metals:
            self.psi = dihedral(self.Coords[0], self.Coords[1], self.Coords[2], next_res.Coords[0])
        
def chi1_4th(res_id):
    if (res_id == "VAL" or res_id == "ILE"): return "CG1"
    if res_id == "CYS": return "SG"
    if res_id == "SER": return "OG"
    if res_id == "THR": return "OG1"
    else: return "CG"
        
def dihedral(a, b, c, d):
    p = np.vstack((a, b, c, d))
    #print(p)
    # Calculate vectors between points, b1, b2, and b3 in the diagram
    b = p[:-1] - p[1:]
    #print(b)
    # "Flip" the first vector so that eclipsing vectors have dihedral=0
    b[0] *= -1
    # Use dot product to find the components of b1 and b3 that are not
    # perpendicular to b2. Subtract those components. The resulting vectors
    # lie in parallel planes.
    v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    #print(v)
    b1 = b[1] / np.linalg.norm(b[1])
    #print(b1)
    x = np.dot(v[0], v[1])
    #print(x)
    m = np.cross(v[0], b1)
    y = np.dot(m, v[1])
    return np.degrees(np.arctan2( y, x ))

def set_all_phi_psi(residues, res_nums):
    for value in range(0,len(res_nums)):
        if value == 0:
            residues[value].set_psi(residues[value+1])
        elif value == len(res_nums) -1:
            residues[value].set_phi(residues[value-1])
        else:
            #print("prev ", residues[value-1].Coords)
            #print("next ", residues[value+1].Coords)
            residues[value].set_phi(residues[value-1])
            residues[value].set_psi(residues[value+1])

def calc_RMSD(res_coords_native, res_coords_model):
    if len(res_coords_model) != len(res_coords_native):
        print("Uneven number of atoms")
        return
    rmsd = np.sqrt(np.sum((res_coords_native - res_coords_model)**2)/len(res_coords_native))
    return rmsd
            
class Atom:
    def __init__(self, name, resnum, x, y, z):
        self.name = name
        self.resnum = int(resnum)
        self.coords = [float(x), float(y), float(z)]

def create_res(pdb):
    all_res = []
    res_nums = []
    with open(pdb, "r") as pdb_file:
        prev_res = 0
        res_atoms = []
        resnum = 0
        chain = ""
        name = ""
        res_index = 0
        for line in pdb_file:
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                #print(res_index)
                new_res = int(line[22:26])
                #print(prev_res, resnum, len(res_atoms))
                if prev_res == new_res:
                    res_atoms.append([line[12:16].strip(), float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    #print(res_atoms)
                else:
                    #print("atoms to append: ", len(res_atoms), name)
                    if len(res_atoms) > 0:
                        #print("atoms present")
                        all_res.append(Residue(resnum, name, chain, res_atoms, res_index))
                        #print(resnum, name, chain, res_index)
                        res_nums.append(resnum)
                        res_index += 1
                    res_atoms = []
                    prev_res = int(line[22:26])
                    resnum = int(line[22:26])
                    name = line[17:20].strip()
                    chain = line[21]
                    res_atoms.append([line[12:16].strip(), float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    
            
        if len(res_atoms) != 0: 
            all_res.append(Residue(resnum, name, chain, res_atoms, res_index))
            res_nums.append(resnum)

    return(all_res, res_nums)

def get_CA_coords(pdb, res_list = [], chainID = "A"):
    with open("%s.pdb" %pdb, "r") as coord_file:
        if len(res_list) > 0:
            res_coords = np.zeros((len(res_list), 3))
            for line in coord_file:
                if line[12:16].strip() == "CA" and int(line[22:26]) in res_list and line[21] == chainID:
                    res_index = res_list.index(int(line[22:26]))
                    res_coords[res_index][0] = float(line[30:38])
                    res_coords[res_index][1] = float(line[38:46])
                    res_coords[res_index][2] = float(line[46:54])
        else:
            res_coords = []
            count = 0
            for line in coord_file:
                if line[12:16].strip() == "CA" and line[21] == chainID:
                    res_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    count +=1
                    #res_coords.extend(coord_line)
            res_coords = np.asarray(res_coords)
            #print(res_coords)
            res_coords = np.reshape(res_coords, (count, 3))
    
    return res_coords
    
def get_all_coords(pdb, res_list = [], chainID = "A"):
    res_list = [int(x) for x in res_list]
    with open("%s.pdb" %pdb, "r") as coord_file:
        if len(res_list) > 0:
            res_coords = np.zeros((len(res_list), 3))
            for line in coord_file:
                if int(line[22:26]) in res_list and line[21] == chainID:
                    res_index = res_list.index(int(line[22:26]))
                    res_coords[res_index][0] = float(line[30:38])
                    res_coords[res_index][1] = float(line[38:46])
                    res_coords[res_index][2] = float(line[46:54])
        else:
            res_coords = []
            count = 0
            for line in coord_file:
                if line[21] == chainID:
                    res_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    count +=1
                    #res_coords.extend(coord_line)
            res_coords = np.asarray(res_coords)
            print(res_coords)
            res_coords = np.reshape(res_coords, (count, 3))

    return res_coords

def get_res_centroids(pdb1, res_list = [], chainID = "A"):
    res_list = [int(x) for x in res_list]
    with open("%s.pdb" %pdb1, "r") as coord_file:
        if len(res_list) > 0:
            res_centroids = np.zeros((len(res_list),3))
            res_coords = []
            count = 0
            res_num = res_list[0]
            for line in coord_file:
                if line[0:4] == "ATOM" and line[12:16].strip() == "N" and len(res_coords) != 0:
                    res_coords = np.asarray(res_coords)
                    #print(res_coords)
                    res_coords = np.reshape(res_coords, (count, 3))
                    res_centroids[res_list.index(res_num)] = np.average(res_coords, axis = 0)
                    #print(res_coords)
                    #print(res_centroids[res_list.index(res_num)])
                    count = 0
                    res_coords = []
                
                res_num = int(line[22:26].strip())    
                if res_num in res_list and line[21] == chainID:
                    count += 1
                    res_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        else:
            res_centroids = []
            res_count = 0
            res_coords = []
            count = 0
            res_num = res_list[0]
            for line in coord_file:
                if line[0:4] == "ATOM" and line[12:16].strip() == "N" and len(res_coords) != 0:
                    res_coords = np.asarray(res_coords)
                    #print(res_coords)
                    res_coords = np.reshape(res_coords, (count, 3))
                    res_centroids.append(np.average(res_coords, axis = 0))
                    count = 0
                    res_count += 1
                    res_coords = []
                
                res_num = int(line[22:26].strip())    
                if line[21] == chainID:
                    count += 1
                    res_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    
            res_centroids = np.asarray(res_centroids)
            #print(res_centroids)
            res_centroids = np.reshape(res_centroids, (res_count, 3))

    return res_centroids
    
def remove_H(pdb1):
    with open("%s.pdb" %pdb1, "r") as origPDB:
        pdb_lines = origPDB.readlines()
    
    with open("%s.pdb" %pdb1, "w+") as new_PDB:
        for line in pdb_lines:
            if line[76:78].strip() != "H": new_PDB.write(line)

def remove_het(pdb1):
     with open("%s.pdb" %pdb1, "r") as origPDB:
         pdb_lines = origPDB.readlines()
    
     with open("%s.pdb" %pdb1, "w+") as new_PDB:
         for line in pdb_lines:
             if line[0:6] != "HETATM": new_PDB.write(line)
             
def clean_pdb(pdb1, HET = False, NoH = True):
    #This will remove any header info as well as any HET atoms by default. If Het is set to True, any HETATMs will be included in clean PDB
    with open("%s.pdb" %pdb1, "r") as origPDB:
        pdb_lines = origPDB.readlines()
    print(pdb_lines[0])
    with open("%s_clean.pdb" %pdb1, "w+") as newPDB:
        finished = False
        for line in pdb_lines:
            if line[0:6] == "EXPDTA":
                if "NMR" in line:
                    print("NMR structure")
                    finished = clean_NMR(pdb_lines, pdb1, HET, NoH)
                    #print(finished)
                else:
                    print("Not NMR")
                    finished = clean_xray(pdb_lines, pdb1, HET, NoH)
            elif finished == True:
                continue
            else: continue
            
def clean_xray(pdb_file, pdb1, het, NoH):
    with open("%s_clean.pdb" %pdb1, "w+") as newPDB:
        for line in pdb_file:
            if het == False:
                if line[0:4] == "ATOM":
                    if NoH == True and line[76:78].strip() != "H":
                        newPDB.write(line)
                    else:
                        newPDB.write(line)
            else:
                if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                    if NoH == True and line[76:78].strip() != "H":
                        newPDB.write(line)
                    else:
                        newPDB.write(line)
    return True

def clean_NMR(pdb_file, pdb1, het, NoH):
    #Clean NMR will only keep the first model of 20
    model = True
    with open("%s_clean.pdb" %pdb1, "w+") as newPDB:
        for line in pdb_file:
            #print(line[0:6], model)
            if line[0:6] == "ENDMDL":
                model = False
            elif model == True:
                if het == False and line[0:4] == "ATOM":
                    if NoH == True and line[76:78].strip() != "H":
                        newPDB.write(line)
                    else:
                        newPDB.write(line)
                elif het == True:
                    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                        if NoH == True and line[76:78].strip() != "H":
                            newPDB.write(line)
                        else:
                            newPDB.write(line)
                else: continue
            else:
                continue
    return True
       
def get_res_numbers(pdb1, chainID = "A"):
    res_list = []
    with open("%s.pdb" %pdb1, "r") as pdb_file:
        for line in pdb_file:
            if line[0:4] == "ATOM" and line[12:16].strip() == "CA"and line[21] == chainID:
                res_list.append(line[22:26].strip())
            elif line[0:6] == "HETATM" and line[12:16].strip() == "CA" and line[17:20] in noncanAA and line[21] == chainID:
                res_list.append(line[22:26].strip())
            elif "ENDMDL" in line: break
    return res_list
    
def pdb_to_seq(pdb_file, chainID = "A"):
    seq = ""
    seq_lines = []
    seq_res = False
    for line in pdb_file:
        if line[0:6] == "SEQRES": 
            seq_res = True
            line = line.split()
            #print(line)
            if line[2] == chainID:
                #seq_line = [oneletAA[aa.index(value)] for value in line[4:]]
                for value in line[4:]:    
                    try:
                        seq += oneletAA[aa.index(value)]
                    except ValueError:
                        seq += "X"
                #seq_line = "".join(seq_line)
                #seq += seq_line
            #print(seq)
        else:    
            #print("No seq res")
            if line[0:4] == "ATOM" and (line[12:16].strip() == "CA" or line[12:16].strip()=="CA1") and line[21] == chainID:
                seq_lines.append(line)
            elif line[0:6] == "HETATM" and line[12:16].strip() == "CA" and line[21] == chainID and line[17:20] in noncanAA:
                seq_lines.append(line)
            elif "ENDMDL" in line: break
    #print(seq_lines)
    if seq_res == True: print("SeqRes records found")
    if seq_res != True:
        start_res = int(seq_lines[0][22:26]) - 1
        #start_res = 0
        for row in seq_lines:
            curr_res = int(row[22:26].strip())
            if curr_res != start_res + 1:
                seq += ("-"* (curr_res - start_res - 1))
            try:
                seq += oneletAA[aa.index(row[17:20])]
            except ValueError:
                seq += "X"
            start_res = curr_res
    #print(seq)
        #seq += oneletAA[aa.index(line[17:20])]    
    return seq

def seq_from_struct(pdb_file, chainID = "A"):
    seq = ""
    seq_lines = []
    for line in pdb_file:
        if line[0:4] == "ATOM" and line[12:16].strip() == "CA" and line[21] == chainID:
            seq_lines.append(line)
        elif line[0:6] == "HETATM" and line[12:16].strip() == "CA" and line[21] == chainID and line[17:20] in noncanAA:
            seq_lines.append(line)
        elif "ENDMDL" in line: break
    #print(seq_lines)
        
    start_res = int(seq_lines[0][22:26]) - 1
    for row in seq_lines:
        curr_res = int(row[22:26].strip())
        if curr_res == start_res:
            continue
        elif curr_res != start_res + 1:
            seq += ("-"* (curr_res - start_res - 1))
            #print(seq, start_res, curr_res)
        try:
            seq += oneletAA[aa.index(row[17:20])]
        except ValueError:
            seq += "X"
        
        start_res = curr_res
    #print(seq)  
    return seq

def one_chain_pdb(pdb_file, chainID = "A"):
    one_chain_lines = []
    for line in pdb_file:
        if line[0:4] == "ATOM" and line[21] == chainID:
            one_chain_lines.append(line)
        elif line[0:6] == "HETATM" and line[21] == chainID:
            one_chain_lines.append(line)
        elif line[0:3] == "TER" and line[21] == chainID:
            one_chain_lines.append(line)
    return one_chain_lines

def renumber_pdb_contig(filename, start_value = 1, inplace = True):
    with open(filename, "r") as orig_pdb:
        orig_pdb = orig_pdb.readlines()
    if inplace == True:
        new_filename = filename
    else: 
        new_filename = filename[:-4] + "Renumb.pdb"
    
    with open(new_filename, "w+") as new_pdb:    
        prev_res = -999
        count = start_value - 1
        for line in orig_pdb:
            if "ANISOU" in line: continue
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                this_res = int(line[22:26].strip())
                if this_res == prev_res:
                    new_res_num = (4-len(str(count)))*" " + str(count)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
                    prev_res = this_res
                else:
                    count += 1
                    new_res_num = (4-len(str(count)))*" " + str(count)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
                    prev_res = this_res
            else:
                new_pdb.write(line)
 
def renumber_pdb(pdb_file, filename, start_value):
    with open(filename, "w+") as new_pdb:
        count = 0
        adj_factor = 0
        for line in pdb_file:
            if "ANISOU" in line: continue
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                if count == 0:
                    count +=1
                    old_start = int(line[22:26].strip()) 
                    adj_factor = start_value - old_start
                    #print(adj_factor)
                    new_res = int(line[22:26].strip()) + adj_factor
                    new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
                else:
                    new_res = int(line[22:26].strip()) + adj_factor
                    new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
            else:
                new_pdb.write(line)
                
def renumber_pdb_inplace(pdb_file, start_value):
    new_lines = []
    count = 0
    adj_factor = 0
    for line in pdb_file:
        if "ANISOU" in line: continue
        if line[0:4] == "ATOM":
            if count == 0:
                count +=1
                old_start = int(line[22:26].strip()) 
                adj_factor = start_value - old_start
                #print(adj_factor)
                new_res = int(line[22:26].strip()) + adj_factor
                new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                new_lines.append(line[0:22] + new_res_num + line[26:])
            else:
                new_res = int(line[22:26].strip()) + adj_factor
                new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                new_lines.append(line[0:22] + new_res_num + line[26:])
        else:
            new_lines.append(line)
    return new_lines
                                
def renumber_pdb_Rosetta(pdb_file, filename, start_value):
    #pdb_file is the readlines input
    with open(filename, "w+") as new_pdb:
        count = 0
        adj_factor = 0
        for line in pdb_file:
            if "ANISOU" in line: continue
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                if count == 0:
                    chain = line[21]
                    count +=1
                    old_start = int(line[22:26].strip()) 
                    adj_factor = start_value - old_start
                    #print(adj_factor)
                    new_res = int(line[22:26].strip()) + adj_factor
                    new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
                elif line[21] == chain:
                    new_res = int(line[22:26].strip()) + adj_factor
                    new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
                    end_res = new_res
                elif line[21] != chain:
                    chain = line[21]
                    old_start = int(line[22:26].strip()) 
                    adj_factor = end_res
                    #print(adj_factor)
                    new_res = int(line[22:26].strip()) + adj_factor
                    new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
                    end_res = new_res
            else:
                new_pdb.write(line)

def ID_neighbors(res_list, resnum, target_res, cutoff):
    #this identifies neighbors within cutoff of ALL res in target_res list
    cutoff = int(cutoff)
    all_neighbors = []
    if (isinstance(target_res, list) == False): target_res = [target_res]
    for value in target_res:
        neighbors = []
        id_value = resnum.index(int(value))
        for res in res_list:
            if np.min(scipy.spatial.distance.cdist(res.Coords, res_list[id_value].Coords, "euclidean")) <= cutoff: 
                neighbors.append(res.resnum)
        if len(all_neighbors) == 0: 
            all_neighbors.extend(neighbors)
        else:
            all_neighbors = sorted(set(all_neighbors).intersection(neighbors))
        #print(all_neighbors)
    all_neighbors = list(sorted(set(all_neighbors).difference(target_res)))
    return(all_neighbors)
    
def ID_neighbors_perl(pdb1, target_res, cutoff):
    #this checks for any res with atom w/in cutoff of any target_res in list
    current_dir = os.getcwd()
    pdb_path = os.path.join(current_dir, pdb1)
    #print(type(target_res))
    if (isinstance(target_res, list) == True): target_res = ":".join(map(str, target_res))
    vicinity = subprocess.check_output(["/Users/meghan/mmtsb/perl/vicinity.pl", "-l", "%s" %target_res, "-hard", "%s" %cutoff, "-soft", "%s" %cutoff - 1, "%s.pdb" %pdb_path])
    vicinity = vicinity.decode("utf-8").strip()
    print(vicinity)
    if (len(vicinity) >2): chain = vicinity[0]
    else: chain = ""
    
    #return of vicinity.pl is in form of A32:45=A47=A52:55; "=" is separation between ranges
    vicinity = vicinity.split("=")    
    #print(chain)
    vicinity = [k[1:] for k in vicinity if k != "+1"]
    vicinity_list = []
    for value in vicinity:
        if ":" in value:
            values = value.split(":")
            x = 0
            for x in range(int(values[0]),int(values[1])+1): vicinity_list.append(str(x))
        else: vicinity_list.append(value)
    #print(vicinity_list)
    return(vicinity_list, chain)
            
def align_seq(seq1, seq2, pdb1 = "pdb1", pdb2 = "pdb2"):
    #there is a better version that is more customized in Sequences.py
    with open("%s_%sSeq.txt" %(pdb1, pdb2), "w+") as seq_file:
        seq_file.write(">" + pdb1 + "\n")
        seq_file.write(seq1 + "\n")
        seq_file.write(">" + pdb2 + "\n")
        seq_file.write(seq2 + "\n")
    clustalW = subprocess.check_output(["/Users/meghan/ClustalW2/clustalw2", "-INFILE=%s_%sSeq.txt" %(pdb1, pdb2)])
    
    pdb1_al_seq = []
    pdb2_al_seq = []
    with open("%s_%sSeq.aln" %(pdb1, pdb2), "r") as aligned_seq:
        for line in aligned_seq:
            if pdb1 in line:
                line = line.strip().split()
                #print(line)
                pdb1_al_seq.extend(line[1])
            elif pdb2 in line:
                line = line.strip().split()
                #print(line)
                pdb2_al_seq.extend(line[1])
            else:
                continue
    #print(pdb1_al_seq)
    #print(pdb2_al_seq)
    os.system("rm *.aln")
    os.system("rm *.dnd")
    return pdb1_al_seq, pdb2_al_seq

def get_POR(pdb1, pdb2):
    pdb1_seq, pdb2_seq = align_seq(pdb1, pdb2)
    pdb1_res = get_res_numbers(pdb1)
    pdb2_res = get_res_numbers(pdb2)
    #print(pdb1_res)
    #print(pdb2_res)
    k = 0
    while k <= len(pdb1_seq)-1:
        if pdb1_seq[k] == "-":
            pdb1_res.insert(k, "-")
        if pdb2_seq[k] == "-":
            pdb2_res.insert(k, "-")
        k+= 1
    #print(pdb1_res)
    #print(pdb2_res)
    return pdb1_res, pdb2_res
    
def align_structures_MMTSB(pdb1, pdb1chain, pdb2):
    #this doesn't quite work...
    #mmtsb alligns by res number; therefore, residues must match
    #aligns SECOND pdb to first pdb; second pdb will be renumbered and its chain changed to match
    pdb1_seq, pdb2_seq = align_seq(pdb1, pdb2)
    pdb1_res = get_res_numbers(pdb1)
    pdb2_res = get_res_numbers(pdb2)
    #print(pdb1_res)
    #print(pdb2_res)
    k = 0
    while k <= len(pdb1_seq)-1:
        if pdb1_seq[k] == "-":
            pdb1_res.insert(k, "-")
        if pdb2_seq[k] == "-":
            pdb2_res.insert(k, "-")
        k+= 1
    #print(pdb1_res)
    #print(pdb2_res)
    with open("%s.pdb" %pdb2, "r") as pdb_file:
        old_pdb = pdb_file.readlines()
    with open("%s_renumb.pdb" %pdb2, "w+") as new_pdb:
        for line in old_pdb:
            if line[0:4] == "ATOM":
                res_num = line[22:26].strip()
                new_res = pdb1_res[pdb2_res.index(res_num)]
                new_pdb.write(line[0:21] + pdb1chain + new_res + line[26:])
    os.system("/Users/meghan/mmtsb/perl/lsqfit.pl -sel heavy %s.pdb %s_renumb.pdb 2>&1 1>%s_aligned.pdb" %(pdb1, pdb2, pdb2)) #-resnumonly
    #os.remove("%s_renumb.pdb" %pdb2)
    print("RMS of aligned structures is: ")
    os.system("/Users/meghan/mmtsb/perl/rms.pl %s.pdb %s_aligned.pdb" %(pdb1, pdb2))

def align_structures(pdb1, pdb2):
    #pdb2 will be aligned to pdb2
    os.system("/Users/meghan/TMalign -A %s.pdb -B %s.pdb >TM_%s.txt" %(pdb1, pdb2, pdb1))
    
    with open("TM_%s.txt" %pdb1, "r") as TM_align:
        tm_lines = TM_align.readlines()
    
    x_matrix = tm_lines[20].strip().split()
    x_matrix = x_matrix[1:]
    y_matrix = tm_lines[21].strip().split()
    y_matrix = y_matrix[1:]
    z_matrix = tm_lines[22].strip().split()
    z_matrix = z_matrix[1:]
    rot_matrix = np.zeros((3,3), dtype = float)
    k = 0
    for k in range(0, 3): 
        rot_matrix[0][k] = x_matrix[k+1]
        rot_matrix[1][k] = y_matrix[k+1]
        rot_matrix[2][k] = z_matrix[k+1]
    trans_matrix = [float(x_matrix[0]), float(y_matrix[0]), float(z_matrix[0])]
    #print(rot_matrix, trans_matrix)
    #rot_matrix[0][1] = x_matrix[1]
    
    with open("%s.pdb" %pdb2, "r") as pdb_file:
        old_pdb = pdb_file.readlines()
    
    X_coords = np.zeros(len(old_pdb), dtype = float)
    Y_coords = np.zeros(len(old_pdb), dtype = float)
    Z_coords = np.zeros(len(old_pdb), dtype = float)    
    for i in range(0, len(old_pdb)):
        if old_pdb[i][0:4] == "ATOM":
            X_coords[i] = float(old_pdb[i][30:38])
            Y_coords[i] = float(old_pdb[i][38:46])
            Z_coords[i] = float(old_pdb[i][46:54])
          
    all_coords = np.zeros((len(X_coords),3))
    j= 0
    for j in range(0, len(X_coords)):
        all_coords[j][0] = trans_matrix[0] + rot_matrix[0][0]*X_coords[j] + rot_matrix[0][1]*Y_coords[j] + rot_matrix[0][2]*Z_coords[j]
        all_coords[j][1] = trans_matrix[1] + rot_matrix[1][0]*X_coords[j] + rot_matrix[1][1]*Y_coords[j] + rot_matrix[1][2]*Z_coords[j]
        all_coords[j][2] = trans_matrix[2] + rot_matrix[2][0]*X_coords[j] + rot_matrix[2][1]*Y_coords[j] + rot_matrix[2][2]*Z_coords[j]
    
    #print("%0.3f" %all_coords[0][0]) # decimal.Decimal.from_float(all_coords[0][0]))
    with open("%s_aligned.pdb" %pdb2, "w+") as new_pdb:
        for i in range(0, len(old_pdb)):
            line = old_pdb[i]
            x_coord = "%0.3f" % all_coords[i][0]
            y_coord = "%0.3f" % all_coords[i][1]
            z_coord = "%0.3f" % all_coords[i][2]
            x_coord = x_coord.rjust(8)
            y_coord = y_coord.rjust(8)
            z_coord = z_coord.rjust(8)
            coords = x_coord+y_coord+z_coord
            #print(len(coords))
            new_pdb.write(line[0:30] + coords + line[54:])
    
    