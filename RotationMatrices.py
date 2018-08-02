import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

def rotation_matrix(vector1, origin):
    #vector1 should be a vector between two points; origin should be the vector to which vector1 is being moved
    vector1 = point/np.linalg.norm(vector1)
    origin = origin/np.linalg.norm(origin)
    v = np.cross(vector1, origin)
    c = np.dot(vector1, origin)
    s = np.linalg.norm(v)
    I = np.identity(3)
    k = np.matrix([(0, -v[2], v[1]), (v[2], 0, -v[0]), (-v[1], v[0], 0)])#.reshape(3,3)
    print(k)
    r = I + k + np.matmul(k,k) * ((1 -c)/(s**2))
    return(r)

def translation(new_point, origin):
    return (new_point - origin)

"""a = np.array([107.934,   77.841,  -54.026]) #lower point
b = np.array([87.014,   62.353,  -39.110]) #upper point

a_orig = np.array([0,0,-15])
b_orig = np.array([0,0,15])
#point = a-b
#origin = a_orig-b_orig

rot_mat = rotation_matrix(a-b, a_orig-b_orig)
print(rot_mat)
trans_vect = translation(np.asarray(np.dot(rot_mat, a))[0], a_orig )
print(trans_vect)

new_a = np.asarray(np.dot(rot_mat, a))[0]
new_b = np.asarray(np.dot(rot_mat, b))[0]
new_points = np.array([new_a, new_b])
print(new_points)
new_points = new_points - trans_vect
print(new_points)
old_points = np.array([a, b])
#np.linalg.norm(new_a - new_b)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(old_points[:,0], old_points[:,1], old_points[:,2], color = "red")
ax.plot(new_points[:,0], new_points[:,1], new_points[:,2], color = "blue")
plt.show()"""