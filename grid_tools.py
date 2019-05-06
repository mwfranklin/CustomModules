import os
import sys
import subprocess
import numpy as np
import scipy.spatial 

#added 11/16 MWF from CustomModules files
def rotation_matrix(vector1, origin):
    #vector1 should be a vector between two points; origin should be the vector to which vector1 is being moved
    vector1 = vector1/np.linalg.norm(vector1)
    origin = origin/np.linalg.norm(origin)
    v = np.cross(vector1, origin)
    c = np.dot(vector1, origin)
    s = np.linalg.norm(v)
    I = np.identity(3)
    k = np.matrix([(0, -v[2], v[1]), (v[2], 0, -v[0]), (-v[1], v[0], 0)])#.reshape(3,3)
    #print(k)
    r = I + k + np.matmul(k,k) * ((1 -c)/(s**2))
    return(r)

#added 11/16 MWF from CustomModules files
def translation(new_point, origin):
    return (new_point - origin)

#added 11/16 MWF from CustomModules files
def rotate_trans_coords(a, b, a_orig, b_orig, res_coord_array):
    #a, b are current points; a_orig, b_orig are where they should move to; res_coord_array contains coordinates to be moved
    #print(a_orig - b_orig)
    rot_mat = rotation_matrix(a-b, a_orig-b_orig)
    #print(rot_mat)
    trans_vect = translation(np.asarray(np.dot(rot_mat, a))[0], a_orig )
    #print(trans_vect)
    matrix_coords = np.transpose(res_coord_array)
    matrix_rot = np.matmul(rot_mat, matrix_coords)
    #print(matrix_rot)
    rotated_coords = np.transpose(matrix_rot)
    #print(rotated_coords)
    rotated_coords = rotated_coords - trans_vect
    return(rotated_coords) #returns a matrix of n by 3 of the translated/rotated vectors

#following functions blatantly copied from http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    #print("D:", D)
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2
    C[1,1] = -1
    #print("S:", S)
    #print("C:", C)
    #print( np.linalg.pinv(S) )
    E, V =  np.linalg.eig(np.dot(np.linalg.pinv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a

def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])

def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    if b == 0:
        if a > c:
            return 0
        else:
            return np.pi/2
    else: 
        if a > c:
            return np.arctan(2*b/(a-c))/2
        else:
            return np.pi/2 + np.arctan(2*b/(a-c))/2

def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    #print(up, down1, down2)
    res1=np.sqrt(np.absolute(up/down1)) #added absolute value 1-10-18 MWF
    res2=np.sqrt(np.absolute(up/down2)) #added absolute value 1-10-18 MWF
    #print(res1, res2)
    return np.array([res1, res2])

#added 11/16 MWF

#given an array of points that have been pre-rotated to x,y,z and a pair of z-boundaries in that array
#return features describing that slice of points
#   farthest_point is the max of the distance between all the points in the slice
#   pocket_area is the area of the convex hull; in otherwords, the approx cross-sectional area of the slice
#   pocket_offset is how far the center of the ellipse is from the origin where the points have been rotated to
#   long/short_axis is the ellipse dimensions
def calc_slice_features(pocket_array, z_depth, interval = 1.4):
    upper_slice_boundary = z_depth + interval
    lower_slice_boundary = z_depth - interval
    #the z-coordinate is not considered for any of these calculations; the slice of points is flattened into just x,y space
    pocket_slice = pocket_array[ (pocket_array[:,2] <= upper_slice_boundary) & (pocket_array[:,2] >= lower_slice_boundary)]
    #print(pocket_slice)
    if len(pocket_slice) > 0:
        all_dist = scipy.spatial.distance.cdist(pocket_slice[:,0:2], pocket_slice[:,0:2])
        farthest_point = np.max(all_dist) 
        #print("farthest:", farthest_point)
        if len(pocket_slice) > 2:
            #calculate the area of the convex hull that encloses the points
            hull = scipy.spatial.ConvexHull(pocket_slice[:,0:2], qhull_options="QJ")
            #print(hull)
            #fit an ellipse to the hull points
            pocket_area = hull.area
            #print("pocket area", pocket_area)
            a = fitEllipse(hull.points[hull.vertices][:,0], hull.points[hull.vertices][:,1])
            center = ellipse_center(a)
            pocket_offset = scipy.spatial.distance.cdist([[0,0]], [center])
            phi = ellipse_angle_of_rotation(a)
            axes = ellipse_axis_length(a)
            long_axis = max(axes)
            short_axis = min(axes)
            #print(long_axis, short_axis)
            return [float(farthest_point), float(pocket_area), float(pocket_offset[0][0]), float(long_axis), float(short_axis)]
        else:
            return[farthest_point, 0, 0, 0, 0]
    else:
        return [0,0,0,0,0]

##Takes in a pocket pdb file and parses down to the relevant grid points    
def process_pocket(pocket_file):
    tp = []
    tpb = []
    tpe = []
    tps = []
    stp = []
    sts = []
    stb = []
    ste = []
    ts = []
    whole_pocket = []

    ##Read in the pocket PDB file
    ##Place the points into the respective lists. 
    ##TP should be key target. but, you can often have issues in pits so stp is useful too. 
    with open(pocket_file, "r") as pocketPDB:
        for entry in pocketPDB:
            if (entry[0:6] == "ATOM  "):
                m_id = int(entry[8:11])
                x = float(entry[31:38])
                y = float(entry[39:46]) 
                z = float(entry[47:54])
                resn = entry[13:16].strip() ##???
                t = np.array([x, y, z])

                if (resn == "PR"): #protein
                    pass
                elif (resn == "TP"): #target pocket
                    tp.append(t)
                    whole_pocket.append(t)
                elif (resn == "TPB"): #target pocket buried
                    tpb.append(t)
                    whole_pocket.append(t)
                elif (resn == "TPE"): #target pocket edge
                    tpe.append(t)
                    #whole_pocket.append(t)
                elif (resn == "TPS"): #target pocket surface
                    tps.append(t)
                    whole_pocket.append(t)
                elif (resn == "STP"): #small pocket
                    stp.append(t)
                    whole_pocket.append(t)
                elif (resn == "STS"): #small pocket surface
                    sts.append(t)
                    whole_pocket.append(t)
                elif (resn == "STB"): #small pocket buried
                    stb.append(t)
                    whole_pocket.append(t)
                elif (resn == "TS"): #target surface? 
                    ts.append(t)
                    whole_pocket.append(t)
                elif (resn == "STE"): #small pocket edge
                    ste.append(t)
                    #whole_pocket.append(t)
    return [tp, tpb, tpe, tps, stp, sts, stb, ts, whole_pocket] #each of these is a list of np arrays; these can be converted to an n by 3 array with np.asarray(X)

def get_vectors(pocket_defs, origin, closest_res):
    #pocket_defs should be a n by 3 numpy array
    #origin should be a [x, y, z] array
    #closest res is also in the form of an [x, y, z] array

    #Find average of subset of points selected to be the pocket opening, i.e. TP/STP residues in grid file
    ##Average of all points in pocket_defs 
    z_target = np.average(pocket_defs, axis = 0)
    #print(z_target)
    #Get Z
    z = z_target - origin
    z /= np.linalg.norm(z)
    #Get X
    x = closest_res - origin
    #print(x)
    #get the orthogonal vector from x
    x -= x.dot(z) * z

    ##Now the cross product to get the Y vector
    y = np.cross(x, z)

    #returns unit vectors of x, y, z
    return [x, y, z, z_target - origin]

def vectors_from_pdb(pocket_grid_name, metal_ion_center, m_closest_residue):
    pocket_defs = process_pocket(pocket_grid_name)
    endpoints = pocket_defs[0] + pocket_defs[4]
    vecs = get_vectors(np.asarray(endpoints), metal_ion_center, m_closest_residue)
 
    return vecs

#added 11-20-18 MWF
def calc_pocket_features(pdb_id, site_id, pocket, origin, nearest_res_coords, center_of_mass, directory):
    feature_labels = ["farPt", "PocketArea", "Offset", "LongAx", "ShortAx"]
    labels1 = [x+"Low" for x in feature_labels]
    labels2 = [x+"Mid" for x in feature_labels]
    labels3 = [x+"High" for x in feature_labels]
    all_labels = ["SEPocket", "Depth", "Vol", "LongPath"] + labels1 + labels2 + labels3

    target_pocket = pocket[0] + pocket[4] #TP and STP residue codes in pocket pdb - target pocket
    #print(target_pocket)
    if len(target_pocket) == 0:
        if len(pocket[3] + pocket[5]) > 0: #use TPS/STS - pocket surface
            target_pocket = pocket[3] + pocket[5]
            pocket_type = 2
            print("USING TPS residues")
        #print(len(pocket[1]))
        elif len(pocket[1]) > 0: #if we have no TP or STP residues, use the TPB ones instead
            pocket_type = 1
            target_pocket = farthest_residue(center_of_mass, pocket[1]) #rather than the average of all of the TPB, 
            print("USING TPB residues from C.O.M.", target_pocket)
        else:
            print("NO SURFACE POCKET")
            if len(pocket[7]) > 0 :
                return(all_labels, [0] * 19)
            else:
                return(all_labels, [4] + [0]*18 )
    else:
        pocket_type = 3
        
    vectors = get_vectors(np.asarray(target_pocket), origin, nearest_res_coords)
    #print(vectors)
    depth_in_pocket = np.linalg.norm(vectors[-1])
    volume = subprocess.check_output(["sed", "1q;d", "%s/%s_Vol.txt"%(directory, site_id)] )
    #print(volume)
    volume = volume.decode("utf-8").strip().split("\n")[0].split()[-1]
    volume = float(volume) #volume of pocket from pocket_grid
    longest_path_out = scipy.spatial.distance.cdist([vectors[-1]], [origin], "cityblock" )[0][0] #if you traverse a grid from origin to the target pocket, how far is it?

    #vectors[2] is the unit z-axis as returned from get_vectors, [0,0,0] is the origin for this z-axis; the pocket array must be pre-translated
    rotated_pocket = np.asarray(rotate_trans_coords(np.array([0,0,0]), vectors[2], np.array([0, 0, 0]), np.array([0,0,1]), np.asarray(pocket[-1])-origin))

    #the pocket grid is currently done on 0.5A intervals; 
    #an interval of 1.4 is approximately the distance between 3 stacked grid points on the diagonal; 
    #a slice containing +- an interval should capture a resonable collection of points
    middle_features = calc_slice_features(rotated_pocket, depth_in_pocket/2)
    upper_features = calc_slice_features(rotated_pocket, depth_in_pocket)
    lower_features = calc_slice_features(rotated_pocket, 0)
    
    this_SITE_features = [pocket_type, depth_in_pocket, volume, longest_path_out]
    this_SITE_features.extend(lower_features)
    this_SITE_features.extend(middle_features)
    this_SITE_features.extend(upper_features)
    
    #solvent-exposed pocket, depth, vol, long_path, [farthest_point, pocket_area, pocket_offset[0][0], long_axis, short_axis] for low/mid/upper
    return(all_labels, this_SITE_features)

class edge_table():
    def __init__(self, pocket_list, target_list,  origin, edge_range=1.5):
        ##Centroid, real or contrived, of the center point
        self.origin = origin 
  
        ##This is the list that holds the water/solvents
        self.pocket_list = pocket_list
  
        #And the targets
        self.target_points = target_list

         ##Attach a distance from the centroid to each atom     
        for i in range(len(pocket_list)):
            self.pocket_list[i].d = dist(origin, self.pocket_list[i])
   
        for i in range(len(target_list)):
            self.target_points[i].d = dist(origin, self.target_points[i])

        ##The algorithm used assumes the list is sorted by
        ##distance from the central ion, or centroid average
        self.pocket_list.sort(key = atom.get_d)
        self.target_points.sort(key = atom.get_d)
 
        #The union thereof
        self.pocket_points = self.pocket_list + self.target_points

        ##This governs the radius in A for a viable edge
        self.edge_range = edge_range

        ##Find the potential origin atoms
        ##Which is to say, find things up to the greatest
        ##Distance of edge range from the center
        self.origin_list = []
        c = 0
        g_dist = 0
        while (g_dist <= edge_range):
            self.origin_list.append(self.pocket_points[c])
            c = c + 1
            if (c <= len(self.pocket_points)):
                g_dist=self.pocket_points[c].d
            else: 
                g_dist=999999
 
        self.origin_l = len(self.origin_list)
        self.check_break = False

        ##From the list passed in, create a list of edges
        self.edges = self.create_edges()
  
        ##Combines the key functionalities of this class
        ##path_out is a list of the solvent molecules,
        ##indexed by distance, from the bulk to the ion
        ##Bulk index is the index of that bulk solvent
        ##And distance, is well, the distance b/w them

    def calculate_path(self):
        self.path_out = []
        ##Testing-Uses the origin_list as start points
        for i in range(len(self.origin_list)):
            putative_path_out = self.travel(i)
            print("Testing origin  " + str(i))

        if (((len(putative_path_out)) < len(self.path_out)) and (len(putative_path_out > 0))):
            self.path_out = putative_path_out
 
        ##Get Bulk Index
        if (len(self.path_out) > 0):
            self.bulk_index = self.path_out[0]
            self.distance = dist(self.origin, self.pocket_points[self.bulk_index])
        else:
            self.path_out_id = []
        print(self.path_out)

        ##Get the actual atom numbers for use in pymol
        for i in range(len(self.path_out)):
            self.path_out_id.append(self.path_out[i].m_id)

        ##Grabs the vector from the first bulk encountered to centroid ion
    def vector_to_bulk(self):
        x_c = origin.x - self.pocket_points[bulk_index].x
        y_c = origin.y - self.pocket_points[bulk_index].y
        z_c = origin.z - self.pocket_points[bulk_index].z
        return [x_c, y_c, z_c]

    def vector_to_first(self):
        x_c = origin.x - self.pocket_points[0].x
        y_c = origin.y - self.pocket_points[0].y
        z_c = origin.z - self.pocket_points[0].z
        return [x_c, y_c, z_c]

     ##Lists out which pocket_points have an edge to which other ways, with the following constraints:
     ##     Waters will only path to pocket_points of equal or greater distance. 
     ##     A shortest route will not involve backtracking.
    def display_edges(self):
        for i in range(len(self.edges)):
            print(str(i) + " -> ", end = ' ')
            for j in range(len(self.edges[i])):
                print(self.edges[i][j], end = " ")
        print()

     ##Main subroutine which attempts to find the
     ##Shortest path to the bulk solvent from ion
     ##Recursively attempts to find the shortest
     ##length travel to a bulk solvent
    def travel(self, dindex, current=[]):
        putative = []
        current.append(dindex)

        if (self.pocket_points[dindex] in ["TP", "STP"]):
            print("TP designed sphere reached at " + str(self.pocket_points[dindex].a_id))
            return [dindex]
        ##If a bulk is found, return the index
        else:
            ##If the index evaluated isn't bulk, then iterate over the edges and recursively probe where those edges lead. 
            ##If a bulk is found, good, assess all other outcomes. 
            ##The shortest path is the one that is returned
            for j in range(len(self.edges[dindex])):
                z = self.travel(self.edges[dindex][j], current)
                if (len(z) > 0):
                    if (len(z) < len(putative) or (len(putative) == 0)):
                        putative = z
                        putative.append(dindex)
    
        ##At termination, return the list. No path means null
        return putative

    ##Subroute that populates the edge table used for the DFS
    def create_edges(self):
        edges = []
        break_check = False

        for i in range(len(self.pocket_points)):
            edges.append([])
            for j in range(i+1, len(self.pocket_points)):
                ##water_counter variable
                if (dist(self.pocket_points[i], self.pocket_points[j]) <= self.edge_range):
                    #Check to make sure we don't count the origins as edges
                    if (break_check ==  False):
                        if (self.pocket_points[j] in self.origin_list):
                            continue  

                if (dist(self.pocket_points[i], self.origin) > self.edge_range):
                    break_check =  True
    
                ##Add the edge to the edge list
                edges[i].append(j)

                if (self.pocket_points[i].resn in ["TP", "STP"]):
                    print("Bulk solvent found at " + str(self.pocket_points[i].a_id))
                    #break

        return edges

def closest_residue(origin, residue_list):
    distances = residue_list - origin
    nearest_neighbor = residue_list[np.argmin( np.sum(distances**2, axis=1))] #return the row that corresponds to the minimum distance
    return nearest_neighbor
    
def farthest_residue(origin, residue_list):
    distances = residue_list - origin
    nearest_neighbor = residue_list[np.argmax( np.sum(distances**2, axis=1))] #return the row that corresponds to the minimum distance
    return nearest_neighbor

##Get the four coordinating residues for the pocket mouth
##Usually as defined by the centroid of the TP set
def mouth_define(centroid, res_list, z_vector, dist_cutoff=1.5, search_max = 2.0):
    #centroid is [x, y, z]; res_list is an n x 3 numpy array, z_vector is [x, y, z]
    first_closest = closest_residue(centroid, res_list)
    
    ##Define vectors from mouth to closest
    fc = first_closest - centroid

    #now get the opposite vector
    ov = fc * -1
    
    #opposite residue centroid definition
    counter_centroid = centroid + ov
    
    far_closest = closest_residue(counter_centroid, res_list)
    y = np.cross(fc, z_vector)
    s = 0.1

    #holders for the y near and far
    y_n = far_closest
    y_f = far_closest

    ##Extend a probe along the vector so we dont overshoot!
    while (s <= search_max):
        wx = y_x * s + centroid.x
        wy = y_y * s + centroid.y
        wz = y_z * s + centroid.z
        w = y*s + centroid

        y_closest = closest_residue(w, res_list)
        if ( ( (np.linalg.norm(y_closest - w) < dist_cutoff) or (s == search_max) ) and (y_closest not in [far_closest, first_closest]) ):
            y_n = y_closest
            s = search_max + 0.1
        else:
            s += 0.1
    ov = y_n - centroid
    ov *= -1
    
    #repeat out the other direction
    s = 0.1
    while (s <= search_max):
        wx = y_x * s + centroid.x
        wy = y_y * s + centroid.y
        wz = y_z * s + centroid.z
        w = y*s + centroid

        y_closest = closest_residue(w, res_list)
        if ( ( (np.linalg.norm(y_closest - w) < dist_cutoff) or (s == search_max) ) and (y_closest not in [far_closest, first_closest, y_n]) ):
            y_f = y_closest
            s = search_max + 0.1
        else:
            s += 0.1
    
    return [first_closest, far_closest, y_n, y_f]


