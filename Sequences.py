import numpy as np
import os
import glob
import decimal
import subprocess
import scipy.spatial
import math
import textwrap 

aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE"]
oneletAA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "M"]
noncanAA = ["MSE", "IAS"]

#Matrices: "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B" (asp), "Z" (glu), "X" (unknown), "*" gap init penalty]
matrix_aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z","X", "*"]
blosum62 = np.array(
[[4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4],
[-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4],
[-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4],
[-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4],
[0,-3,-3,-3,9,3,4,3,3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
[-1,1,0,0,3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4],
[-1,0,0,2,4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
[0,-2,0,-1,3,-2,-2,6,2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4],
[-2,0,1,-1,3,0,0,2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4],
[-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4],
[-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4],
[-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-1,-1,0,-1,-3,-2,-2,0,1,-1,-4],
[-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,-1,-3,-1,-1,-4],
[-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-1,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4],
[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
[1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4],
[0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4],
[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4],
[-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4],
[0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,-1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4],
[-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4],
[-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
[0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4],
[-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1]])

bbTM = np.array([
   [6,-11,-7,-10,-19,-9,-10,1,-12,-2,-2,-9,-5,-5,-8,-1,0,-10,-10,0,-8,-4,-6,-24],
    [-11,10,-7,-11,-13,-4,-10,-12,-5,-12,-11,-1,-12,-11,-14,-6,-8,-14,-11,-12,-9,-11,-8,-24],
    [-7,-7,10,0,-13,-4,-8,-7,-6,-12,-12,-5,-7,-12,-14,0,-2,-13,-10,-11,5,-7,-7,-24],
    [-10,-11,0,11,-13,-5,-2,-10,-10,-13,-11,-9,-12,-14,-14,-4,-8,-14,-14,-12,6,-6,-8,-24],
    [-19,-13,-13,-13,24,-11,-16,-16,-16,-15,-16,-14,-13,-14,-14,-12,-13,-16,-12,-15,-14,-16,-12,-24],
    [-9,-4,-4,-5,-11,10,-1,-11,-2,-9,-10,-4,-7,-13,-12,-4,-5,-13,-11,-9,-4,-6,-6,-24],
    [-10,-10,-8,-2,-16,-1,10,-11,-12,-14,-13,-10,-12,-15,-11,-8,-8,-14,-14,-12,-5,0,-9,-24],
    [1,-12,-7,-10,-16,-11,-11,7,-14,-8,-8,-12,-8,-10,-13,-2,-4,-13,-12,-6,-8,-2,-8,-24],
    [-12,-5,-6,-10,-16,-2,-12,-14,13,-13,-12,-11,-12,-11,-16,-9,-10,-10,-5,-13,-8,-13,-9,-24],
    [-2,-12,-12,-13,-15,-9,-14,-8,-13,7,3,-12,-3,-2,-14,-8,-5,-8,-9,5,-12,-11,-7,-24],
    [-2,-11,-12,-11,-16,-10,-13,-8,-12,3,6,-12,0,1,-13,-9,-6,-6,-6,3,-11,-10,-6,-24],
    [-9,-1,-5,-9,-14,-4,-10,-12,-11,-12,-12,11,-11,-13,-14,-7,-6,-12,-14,-12,-7,-11,-8,-24],
    [-5,-12,-7,-12,-13,-7,-12,-8,-12,-3,0,-11,12,-5,-15,-7,-2,-10,-10,-2,-9,-10,-7,-24],
    [-5,-11,-12,-14,-14,-13,-15,-10,-11,-2,1,-13,-5,8,-11,-10,-9,-3,1,-1,-13,-12,-7,-24],
    [-8,-14,-14,-14,-14,-12,-11,-13,-16,-14,-13,-14,-15,-11,14,-11,-12,-16,-16,-13,-14,-12,-11,-24],
    [-1,-6,0,-4,-12,-4,-8,-2,-9,-8,-9,-7,-7,-10,-11,8,2,-15,-12,-7,-2,-5,-6,-24],
    [0,-8,-2,-8,-13,-5,-8,-4,-10,-5,-6,-6,-2,-9,-12,2,8,-12,-10,-3,-5,-6,-5,-24],
    [-10,-14,-13,-14,-16,-13,-14,-13,-10,-8,-6,-12,-10,-3,-16,-15,-12,11,-1,-6,-13,-13,-10,-24],
    [-10,-11,-10,-14,-12,-11,-14,-12,-5,-9,-6,-14,-10,1,-16,-12,-10,-1,7,-8,-12,-13,-8,-24],
    [0,-12,-11,-12,-15,-9,-12,-6,-13,5,3,-12,-2,-1,-13,-7,-3,-6,-8,5,-11,-9,-6,-24],
    [-8,-9,5,6,-14,-4,-5,-8,-8,-12,-11,-7,-9,-13,-14,-2,-5,-13,-12,-11,11,-4,-7,-24],
    [-4,-11,-7,-6,-16,-6,0,-2,-13,-11,-10,-11,-10,-12,-12,-5,-6,-13,-13,-9,-4,10,-7,-24],
    [-6,-8,-7,-8,-12,-6,-9,-8,-9,-7,-6,-8,-7,-7,-11,-6,-5,-10,-8,-6,-7,-7,1,-24],
    [-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,1]
])
all_matrices_names = ["BLOSUM62", "bbTM"]
all_matrices = [blosum62, bbTM]

def calc_sim_score(seq1, seq2, matrix="BLOSUM62", gap_pen = -1):
    score= 0
    if matrix != "BLOSUM62": matrix = all_matrices[all_matrices_names.index(matrix)]
    else: matrix = blosum62
    #print(matrix)
    
    for pair in zip(seq1, seq2):
        if "-" in pair: score += int(gap_pen)
        else: 
            #print(pair, matrix[matrix_aa.index(pair[0])][matrix_aa.index(pair[1])])
            score += matrix[matrix_aa.index(pair[0])][matrix_aa.index(pair[1])]
    
    return score

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
            if line[0:4] == "ATOM" and (line[12:16].strip() == "CA" or line[12:16].strip()=="CA1")and line[21] == chainID:
                seq_lines.append(line)
            elif line[0:6] == "HETATM" and line[12:16].strip() == "CA" and line[21] == chainID and line[17:20] in noncanAA:
                seq_lines.append(line)
            elif "ENDMDL" in line: break
    #print(seq_lines)
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
    
def align_seq_perl(pdb1, pdb2):
    SeqFile = open("%s_%sSeq.txt" %(pdb1, pdb2), "w+")
    SeqFile.write(">" + pdb1 + "\n")
    SeqFile.close()
    os.system("perl /Users/meghan/CustomModules/pdb2seq1h.pl %s.pdb >>%s_%sSeq.txt" %(pdb1,pdb1, pdb2))
    SeqFile = open("%s_%sSeq.txt" %(pdb1, pdb2), "a")
    SeqFile.write(">" + pdb2 + "\n")
    SeqFile.close()
    os.system("perl /Users/meghan/CustomModules/pdb2seq1h.pl %s.pdb >>%s_%sSeq.txt" %(pdb2, pdb1, pdb2))
    
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
    
def align_seq(pdb1, pdb2, clustalW_path = "/Users/meghan/CompBioPrograms/ClustalW2/clustalw2", DB_path = "", chainID1="A", chainID2="A"):
    pdb1_seq = []
    pdb2_seq = []
    with open("%s%s.pdb"%(DB_path, pdb1), "r") as inFile:
        inFile = inFile.readlines()
    pdb1_seq = pdb_to_seq(inFile, chainID1)
    pdb1_seq = textwrap.wrap(pdb1_seq, 50)
    with open("%s%s.pdb"%(DB_path, pdb2), "r") as inFile:
        inFile = inFile.readlines()
    pdb2_seq = pdb_to_seq(inFile, chainID2)
    pdb2_seq = textwrap.wrap(pdb2_seq, 50)
    
    with open("%s_%s_%s_%s_Seq.txt" %(pdb1, chainID1, pdb2, chainID2), "w+") as SeqFile:
        SeqFile.write(">" + pdb1+ "_" + chainID1 + "\n")
        SeqFile.write("\n".join(pdb1_seq) + "\n")
        SeqFile.write(">" + pdb2+ "_" + chainID2 + "\n")
        SeqFile.write("\n".join(pdb2_seq) + "\n")

    clustalW = subprocess.check_output(["%s"%clustalW_path, "-INFILE=%s_%s_%s_%s_Seq.txt" %(pdb1, chainID1, pdb2, chainID2), "-SCORE=PERCENT", "-STATS=%s_%sScore.txt"%(pdb1, pdb2)])
    
    pdb1_al_seq = []
    pdb2_al_seq = []
    with open("%s_%s_%s_%s_Seq.aln" %(pdb1, chainID1, pdb2, chainID2), "r") as aligned_seq:
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
    with open("%s_%sScore.txt" %(pdb1, pdb2), "r") as score_file:
        for line in score_file:
            if "aln pw-id avg:" in line:
                #print(line)
                score = float(line.strip()[-5:].strip())
                #print(score)
            if "aln len" in line:
                length = int(line.strip()[8:].strip())
    #print(pdb1_al_seq)
    #print(pdb2_al_seq)
    #os.system("rm %s_%s_%s_%s_Seq.aln"%(pdb1, chainID1, pdb2, chainID2))
    os.system("rm *.dnd")
    os.system("rm %s_%sScore.txt"%(pdb1, pdb2))
    os.system("rm %s_%s_%s_%s_Seq.txt"%(pdb1, chainID1, pdb2, chainID2))
    return "".join(pdb1_al_seq), "".join(pdb2_al_seq), score, length
    
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
    
    