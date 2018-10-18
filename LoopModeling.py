import numpy as np
import os
import glob
import decimal
import subprocess
import pandas

aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
oneletAA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

def model_same_seq(start_at, loop_size, filename, chainID, pdbname, mDis = False):
    start_at = int(start_at)
    done_list = []
    with open(filename, "r") as orig_pdb, open("%s_UndefLoops.pdb" %pdbname, "w+") as insert_pdb:
        for line in orig_pdb:
            if (("ATOM" in line[0:4] and line[21] == chainID) or (line[17:20] == "MSE" and line[21] == chainID) or ("HETATM" in line[0:6] and line[21] == chainID)): 
                if (int(line[22:26]) < start_at or int(line[22:26]) > (start_at + loop_size -1)):
                    insert_pdb.write(line)
                elif int(line[22:26]) == start_at:
                    if line[12:16].strip() == "N" or line[12:16].strip() == "CA":
                        insert_pdb.write(line)
                    else: 
                        continue
                elif int(line[22:26]) == start_at + loop_size -1:
                    if line[12:16].strip() == "C" or line[12:16].strip() == "CA":
                        insert_pdb.write(line)
                    else: 
                        continue
                else:
                    if int(line[22:26]) not in done_list:
                        if mDis == False:
                            insert_pdb.write(line[0:11] + "  CA " + line[16:26] + "       0.000   0.000   0.000  1.00  1.00\n")
                        else:
                            insert_pdb.write(line[0:11] + "  H  " + line[16:26] + "       0.000   0.000   0.000  1.00  1.00\n")
                        done_list.append(int(line[22:26]))
            elif "ENDMDL" in line:
                break
                
def insert_missing_seq(insert_at, seq, filename, pdbname, mDis = False):
    insert_at = int(insert_at)
    inserted = False
    check_lists = False
    res_pairs = {}
    end_res = 0
    with open(filename, "r") as orig_pdb:
        orig_pdb = orig_pdb.readlines()
        
    with open("%s_insert.pdb" %pdbname, "w+") as insert_pdb:
        for line in orig_pdb:
            if "ATOM" in line:
                if int(line[22:26]) < insert_at - 1:
                    insert_pdb.write(line)
                elif int(line[22:26]) == insert_at -1:
                    if line[12:16].strip() == "N" or line[12:16].strip() == "CA":
                        insert_pdb.write(line)
                    else: 
                        continue
                elif int(line[22:26]) >= insert_at + len(seq) + 1:
                    insert_pdb.write(line)
                    
                elif int(line[22:26]) > insert_at:
                    if inserted == False:
                        atom_num = int(line[6:11])+1
                        start_res = insert_at
                        for char in seq:
                            new_atom_num = (5-len(str(atom_num)))*" " + str(atom_num)
                            new_res_num = (4-len(str(start_res)))*" " + str(start_res)
                            if mDis == True:
                                new_line = "ATOM  " + new_atom_num + " H    " + aa[oneletAA.index(char)] + "  " + new_res_num + "       0.000   0.000   0.000  1.00  1.00"
                            else:
                                new_line = "ATOM  " + new_atom_num + " CA   " + aa[oneletAA.index(char)] + "  " + new_res_num + "       0.000   0.000   0.000  1.00  1.00"
                            print(new_line)
                            insert_pdb.write(new_line + "\n")
                            start_res +=1
                            atom_num+=1
                        inserted = True

                    if line[12:16].strip() == "CA" or line[12:16].strip() == "C":
                        insert_pdb.write(line)
                    else:
                        continue
                
            elif "TER" in line:
                continue
            else:
                insert_pdb.write(line)
                    
def insert_seq_renumb(insert_at, seq, filename, pdbname, mDis = False):
    insert_at = int(insert_at)
    inserted = False
    check_lists = False
    res_pairs = {}
    end_res = 0
    with open(filename, "r") as orig_pdb:
        orig_pdb = orig_pdb.readlines()
        
    with open("%s_insert.pdb" %pdbname, "w+") as insert_pdb:
        for line in orig_pdb:
            if "ATOM" in line:
                if check_lists == True:
                    if line[22:26] in res_pairs:
                        insert_pdb.write(line[0:22] + res_pairs[line[22:26]] + line[27:])
                    else:
                        insert_pdb.write(line)
                    
                elif int(line[22:26]) < insert_at - 1:
                    insert_pdb.write(line)
                elif int(line[22:26]) == insert_at -1:
                    if line[12:16].strip() == "N" or line[12:16].strip() == "CA":
                        insert_pdb.write(line)
                    else: 
                        continue
                elif int(line[22:26]) == insert_at:
                    if inserted == False:
                        atom_num = int(line[6:11])+1
                        start_res = insert_at
                        for char in seq:
                            new_atom_num = (5-len(str(atom_num)))*" " + str(atom_num)
                            new_res_num = (4-len(str(start_res)))*" " + str(start_res)
                            if mDis == True:
                                new_line = "ATOM  " + new_atom_num + " H    " + aa[oneletAA.index(char)] + "  " + new_res_num + "       0.000   0.000   0.000  1.00  1.00"
                            else:
                                new_line = "ATOM  " + new_atom_num + " CA   " + aa[oneletAA.index(char)] + "  " + new_res_num + "       0.000   0.000   0.000  1.00  1.00"
                            print(new_line)
                            insert_pdb.write(new_line + "\n")
                            start_res +=1
                            atom_num+=1
                        inserted = True
                
                    if line[12:16].strip() == "CA" or line[12:16].strip() == "C":
                        atom1 = line[6:11]
                        atom_num = int(line[6:11]) + len(seq)
                        new_atom_num = (5-len(str(atom_num)))*" " + str(atom_num)
                    
                        res1 = (line[22:26])
                        res_num = int(line[22:26]) + len(seq)
                        new_res_num = (4-len(str(res_num)))*" " + str(res_num)
                        if res1 not in res_pairs: res_pairs[res1] = new_res_num
                        insert_pdb.write(line[0:6] + new_atom_num + line[11:22] + new_res_num + line[26:])
                        end_res = new_res_num
                    else:
                        continue
                elif check_lists == False: 
                    atom1 = line[6:11]
                    atom_num = int(line[6:11]) + len(seq)
                    new_atom_num = (5-len(str(atom_num)))*" " + str(atom_num)
                    
                    res1 = (line[22:26])
                    res_num = int(line[22:26]) + len(seq)
                    new_res_num = (4-len(str(res_num)))*" " + str(res_num)
                    if res1 not in res_pairs: res_pairs[res1] = new_res_num
                    
                    insert_pdb.write(line[0:6] + new_atom_num + line[11:22] + new_res_num + line[26:])
                
            else:
                if "ENDMDL" in line: 
                    check_lists = True
                    #print(res_pairs)
                insert_pdb.write(line)

def insert_seq_del_seq_renumb(insert_at, seq, filename, pdbname, stop_at = 0, len_del = 0, byLength = True, mDis = False):
    print(insert_at, insert_at+len_del, stop_at)
    insert_at = int(insert_at)
    stop_at = int(stop_at)
    inserted = False
    check_lists = False
    res_pairs = {}
    end_res = 0
    with open(filename, "r") as orig_pdb:
        orig_pdb = orig_pdb.readlines()
        
    with open("%s_insert.pdb" %pdbname, "w+") as insert_pdb:
        for line in orig_pdb:
            if "ATOM" in line:
                if check_lists == True:
                    continue 
                elif int(line[22:26]) < insert_at - 1:
                    insert_pdb.write(line)
                elif int(line[22:26]) == insert_at -1:
                    if line[12:16].strip() == "N" or line[12:16].strip() == "CA":
                        insert_pdb.write(line)
                    else: 
                        continue
                elif int(line[22:26]) == insert_at:
                    if inserted == False:
                        atom_num = int(line[6:11])+1
                        start_res = insert_at
                        for char in seq:
                            new_atom_num = (5-len(str(atom_num)))*" " + str(atom_num)
                            new_res_num = (4-len(str(start_res)))*" " + str(start_res)
                            if mDis == True:
                                new_line = "ATOM  " + new_atom_num + "  H   " + aa[oneletAA.index(char)] + "  " + new_res_num + "       0.000   0.000   0.000  1.00  1.00"
                            else:
                                new_line = "ATOM  " + new_atom_num + "  CA  " + aa[oneletAA.index(char)] + "  " + new_res_num + "       0.000   0.000   0.000  1.00  1.00"
                            print(new_line)
                            insert_pdb.write(new_line + "\n")
                            start_res +=1
                            atom_num+=1
                        inserted = True
                    else:
                        continue
                elif byLength == True and int(line[22:26]) > insert_at and int(line[22:26])< insert_at+len_del: #the length of residues starting at insert point
                    #print(line, "Being replaced")
                    continue
                elif byLength == True and int(line[22:26])== insert_at+len_del:
                    if line[12:16].strip() == "CA" or line[12:16].strip() == "C":
                        atom1 = line[6:11]
                        atom_num = int(line[6:11]) + len(seq)
                        new_atom_num = (5-len(str(atom_num)))*" " + str(atom_num)
                    
                        res1 = (line[22:26])
                        res_num = int(line[22:26]) + len(seq) - len_del
                        new_res_num = (4-len(str(res_num)))*" " + str(res_num)
                        if res1 not in res_pairs: res_pairs[res1] = new_res_num
                        end_res = new_res_num
                        insert_pdb.write(line[0:6] + new_atom_num + line[11:22] + new_res_num + line[26:])
                elif byLength == False and int(line[22:26]) > insert_at and int(line[22:26]) <= stop_at: #the length of residues starting at insert point
                    #print(line, "Being replaced by stopat")
                    if line[12:16].strip() == "CA": len_del += 1
                    continue
                elif byLength == False and int(line[22:26])== stop_at+1:
                    if line[12:16].strip() == "CA": len_del += 1
                    if line[12:16].strip() == "CA" or line[12:16].strip() == "C":
                        atom1 = line[6:11]
                        atom_num = int(line[6:11]) + len(seq)
                        new_atom_num = (5-len(str(atom_num)))*" " + str(atom_num)
                    
                        res1 = (line[22:26])
                        res_num = int(line[22:26]) + len(seq) - len_del
                        new_res_num = (4-len(str(res_num)))*" " + str(res_num)
                        if res1 not in res_pairs: res_pairs[res1] = new_res_num
                        end_res = new_res_num
                        insert_pdb.write(line[0:6] + new_atom_num + line[11:22] + new_res_num + line[26:])
                    
                elif check_lists == False: 
                    res1 = (line[22:26])
                    res_num = int(line[22:26]) + len(seq) - len_del                
                    new_res_num = (4-len(str(res_num)))*" " + str(res_num)
                    if res1 not in res_pairs: res_pairs[res1] = new_res_num
                    
                    insert_pdb.write(line[0:22] + new_res_num + line[26:])
            
            elif check_lists == False:
                if "ENDMDL" in line: 
                    check_lists = True
                    #print(res_pairs)
                insert_pdb.write(line)
    return (insert_at-1), end_res

def process_PETALS(filename):
    scores = subprocess.check_output(["grep", "^LOOP", filename])
    scores = scores.decode("utf-8").strip().split("\n")
    scores = [x.split()[1:] for x in scores]
    scores = pandas.DataFrame.from_records(scores, columns = ("ModelNum", "dfire", "oscar", "dis", "bb"))
    scores.to_csv(filename[:-4] + "_Scores.txt", sep = "\t", index = False)
                    
    with open(filename, "r") as inData:
        raw_petals = inData.readlines()

    #score_index = score_fxns.index(orderby)
    with open(filename[:-4] + ".pdb", "w+") as outData:
        for line in raw_petals:
            if "LOOP" in line:
                line = line.split()
                outData.write("MODEL %s\n"%line[1])
                outData.write(" ".join(line) + "\n")
            elif "END" in line:
                outData.write("ENDMDL\n")
            else:
                outData.write(line.strip() + "  1.00  0.00\n")
      
def continuous_atoms(filename, pdbname):
    with open(filename, "r") as orig_pdb:
        orig_pdb = orig_pdb.readlines()
        
    with open("%s_continuous.pdb" %pdbname, "w+") as insert_pdb:
        count = 0
        for line in orig_pdb:
            if "ATOM" in line:
                count +=1
                atom1 = str(count)
                new_atom_num = (5-len(str(atom1)))*" " + str(atom1)
                
                insert_pdb.write(line[0:6] + new_atom_num + line[11:])

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
                    res_coords.extend(coord_line)
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

def get_res_numbers(pdb1):
    res_list = []
    with open("%s.pdb" %pdb1, "r") as pdb_file:
        for line in pdb_file:
            if line[0:4] == "ATOM" and line[12:16].strip() == "CA":
                res_list.append(line[22:26].strip())
    return res_list
