#! /usr/bin/python
# Copyright (c) 2010 Robert L. Campbell
#
# A PyMOL script for drawing a CGO plane from the coordinates of three atoms (pk1,pk2,pk3 by default)

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#
# @brief Script to visualize membrane planes in pymol
# @author JKLeman (julia.koehler1982@gmail.com)

from pymol.cgo import *
from pymol import cmd

# Store xyz coordinates
class XYZCoord():
	
  def __init__(self, x=0, y=0, z=0):
    self.x = x
    self.y = y
    self.z = z
		
  def __str__(self):
    return "(%8.3f, %8.3f, %8.3f)" % (self.x, self.y, self.z)

cmd.extend( "XYZCoord", XYZCoord )

################################################################################

# Subtract v2 from v1, represented as XYZcoord classes (see above)
def subtract( v1, v2 ):

  x = v1.x - v2.x
  y = v1.y - v2.y
  z = v1.z - v2.z
				
  vf = XYZCoord( x, y, z)
  return vf

cmd.extend( "subtract", subtract )

################################################################################

# Add v1 and v2
def add( v1, v2 ):

  x = v1.x + v2.x
  y = v1.y + v2.y
  z = v1.z + v2.z
				
  vf = XYZCoord( x, y, z)
  return vf

cmd.extend( "add", add )

################################################################################

# Cross-product
def cross( v1, v2 ):

  x = v1.y * v2.z - v1.z * v2.y
  y = v1.z * v2.x - v1.x * v2.z
  z = v1.x * v2.y - v1.y * v2.x
				
  vf = XYZCoord( x, y, z)
  return vf

cmd.extend( "cross", cross )

################################################################################

# Return length of vector v, represented as XYZcoord class (see above)
def length( v ):

  return math.sqrt( v.x**2 + v.y**2 + v.z**2 )

cmd.extend( "length", length )

################################################################################

# Normalize vector v to length l
def normalize( v, l ):

  v_length = length( v )

  x = v.x * l / v_length
  y = v.y * l / v_length
  z = v.z * l / v_length

  vf = XYZCoord( x, y, z)

  return vf

cmd.extend( "normalize", normalize )

################################################################################

# A helper function for computing the normal to a triangular facet
def compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3):

  nx = (y2-y1)*(z3-z2) - (z2-z1)*(y3-y2)
  ny = (z2-z1)*(x3-x2) - (x2-x1)*(z3-z2)
  nz = (x2-x1)*(y3-y2) - (y2-y1)*(x3-x2)

  return (nx,ny,nz)



def draw_plane_cgo(name,apex1,apex2,apex3,apex4,color,transparency):

  """
DESCRIPTION
    Create a CGO plane from three arbitary coordinates

USAGE
    draw_plane_cgo apex1, apex2, apex3, apex4, color

    where each apex is a 3-element vector and color is a 3-element RGB
    list defining the color of the plane (where each value of R, G
    and B is between 0 and 1 inclusive).

  """

  # Convert args to floating point numbers
  x1,y1,z1 = map(float,apex1)
  x2,y2,z2 = map(float,apex2)
  x3,y3,z3 = map(float,apex3)
  x4,y4,z4 = map(float,apex4)
  if type(color) == type(''):
    color = map(float,color.replace('(','').replace(')','').split(','))
  print("CGO color = ", color)
  # Compute the normal vector for the triangle
  normal1 = compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3)
  normal2 = compute_normal(x1, y1, z1, x3, y3, z3, x4, y4, z4)
  normal3 = compute_normal(x2, y2, z2, x3, y3, z3, x4, y4, z4)

  # Create the CGO objects
  obj = [

    BEGIN, TRIANGLE_STRIP,
	
    COLOR, color[0], color[1], color[2],
	ALPHA, 0.5,
	#COLOR, 0.5, 0.5, 0.5,
    NORMAL, normal1[0], normal1[1], normal1[2],
    VERTEX, x1, y1, z1,
    VERTEX, x2, y2, z2,
    VERTEX, x3, y3, z3,
    VERTEX, x4, y4, z4,

    END
  ]

  # Display them
  cmd.load_cgo(obj,name)

def draw_plane(name,atom1='(pk1)',atom2='(pk2)',atom3='(pk3)',atom4='(pk4)',color=[0,0,0],transparency=0.5):
  """
DESCRIPTION
    Create a CGO plane from four atomic coordinates

USAGE
    draw_plane name, atom1, atom2, atom3, atom4, color

    where each atom is a standard PyMOL selection (defaults to pk1,pk2
    and pk3) and color is a 3-element RGB tuple defining the color
    of the plane (where each value of R, G and B is between 0 and 1
    inclusive).  The color defaults to (1,1,1).

  """

# get coordinates from atom selections
  coor1 = cmd.get_model(atom1).atom[0].coord
  coor2 = cmd.get_model(atom2).atom[0].coord
  coor3 = cmd.get_model(atom3).atom[0].coord
  coor4 = cmd.get_model(atom4).atom[0].coord
  draw_plane_cgo(name,coor1,coor2,coor3,coor4,color,transparency)



def draw_memb(pdb_name):
    cmd.show_as("cartoon", "%s"%pdb_name)
    cmd.select("mem", "resn MEM")
    cmd.select("emb", "resn EMB")
    cmd.hide("everything", "mem")
    cmd.hide("everything", "emb")
    cmd.alter("emb", "vdw=1.5")
    cmd.rebuild()
    """cmd.set_view (\
         0.958278537,   -0.025210002,    0.284718215,\
         0.285831660,    0.085169815,   -0.954486012,\
        -0.000186668,    0.996048152,    0.088822074,\
        -0.000005720,    0.000001855, -301.575897217,\
         3.013955355,   -3.607778311,   -1.902338266,\
       241.575897217,  361.575897217,  -20.000000000 )"""

    myspace= {"with_mem": [], "center_list": [], "normal_list": []}
    cmd.iterate("%s and resn MEM"%pdb_name, 'with_mem.append("MEM")', space=myspace)
    print(myspace["with_mem"])
    if len(myspace["with_mem"]) == 0:
        print("No membrane")
    else:
    	# set membrane plane length
    	width = 100

    	# Read in center
    	#myspace= {}
    	cmd.iterate_state(1, "%s and resn MEM and name CNTR"%pdb_name, "center_list.append((x,y,z))", space=myspace )
        center_list = myspace["center_list"]
    	center = XYZCoord(center_list[0][0], center_list[0][1], center_list[0][2])

    	# Read in normal position
    	#myspace= {"normal_list": []}
    	cmd.iterate_state(1, "%s and resn MEM and name NORM"%pdb_name, "normal_list.append((x, y, z))", space=myspace )
    	normal_list = myspace["normal_list"]
        normalp = XYZCoord(normal_list[0][0], normal_list[0][1], normal_list[0][2])

    	# compute normal vector, leaflet thickness is 15A
    	normal = subtract( normalp, center )
    	normal = normalize( normal, 15 )

    	# get upper and lower center point along normal
    	upper_centerp = add( center, normal )
    	lower_centerp = subtract( center, normal )
    	print(upper_centerp)
    	print(lower_centerp)

    	# get a vector perpendicular (in membrane plane) to normal
    	v1 = XYZCoord()
    	v1.x = normal.z
    	v1.y = normal.z
    	v1.z = -normal.x - normal.y
    	v1n = normalize( v1, width )

    	# get vector perpendicular (in membrane plane) to v1 and normal
    	v2 = XYZCoord()
    	v2 = cross( normal, v1 )
    	v2n = normalize( v2, width )

    	# get 4 points defining upper plane
    	p1 = add( upper_centerp, v1n )
    	p2 = add( upper_centerp, v2n )
    	p3 = subtract( upper_centerp, v2n )
    	p4 = subtract( upper_centerp, v1n )

    	# get 4 points defining the lower plane 
    	q1 = add( lower_centerp, v1n )
    	q2 = add( lower_centerp, v2n )
    	q3 = subtract( lower_centerp, v2n )
    	q4 = subtract( lower_centerp, v1n )
        with open("tmp.pml", "w+") as temp_pseudo:
        	temp_pseudo.write("pseudoatom p1, pos=[" + str(p1.x) + ", " + str(p1.y) + ", " + str(p1.z) + "]\n")
        	temp_pseudo.write("pseudoatom p2, pos=[" + str(p2.x) + ", " + str(p2.y) + ", " + str(p2.z) + "]\n")
        	temp_pseudo.write("pseudoatom p3, pos=[" + str(p3.x) + ", " + str(p3.y) + ", " + str(p3.z) + "]\n")
        	temp_pseudo.write("pseudoatom p4, pos=[" + str(p4.x) + ", " + str(p4.y) + ", " + str(p4.z) + "]\n")

        	temp_pseudo.write("pseudoatom q1, pos=[" + str(q1.x) + ", " + str(q1.y) + ", " + str(q1.z) + "]\n")
        	temp_pseudo.write("pseudoatom q2, pos=[" + str(q2.x) + ", " + str(q2.y) + ", " + str(q2.z) + "]\n")
        	temp_pseudo.write("pseudoatom q3, pos=[" + str(q3.x) + ", " + str(q3.y) + ", " + str(q3.z) + "]\n")
        	temp_pseudo.write("pseudoatom q4, pos=[" + str(q4.x) + ", " + str(q4.y) + ", " + str(q4.z) + "]\n\n")
        cmd.do("run tmp.pml")
    	draw_plane("arbtr_up", "p1", "p2", "p3", "p4")
    	draw_plane("arbtr_lo", "q1", "q2", "q3", "q4")

    	# set transparency
    	cmd.set("cgo_transparency", 0.5, "arbtr_lo")
    	cmd.set("cgo_transparency", 0.5, "arbtr_up")

    	# remove pseudoatoms
    	cmd.delete("p1")
    	cmd.delete("p2")
    	cmd.delete("p3")
    	cmd.delete("p4")
    	cmd.delete("q1")
    	cmd.delete("q2")
    	cmd.delete("q3")
    	cmd.delete("q4")

#cmd.extend("draw_plane", draw_plane)
