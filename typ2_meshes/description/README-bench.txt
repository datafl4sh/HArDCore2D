Website: http://www.i2m.univ-amu.fr/fvca5/benchmark/Meshes/index.html

**********************************************************************
**********************************************************************
****                 Benchmark on finite volume schemes        *******
****                 for anisotropic diffusion problems        *******
**********************************************************************

**********************************************************************
****                           The different meshes           ********
********************************************************************** 
All of these meshes are discretizations of the domain [0,1]^2 
 
      -> mesh1_i is a classical triangular mesh 
         *******
            mesh1_i has been obtained by duplication of mesh1_0

      -> mesh2_i is a uniform rectangular  mesh
         ******* the mesh size of mesh2_i is equal to 2^(-i-1)

      -> mesh3_i is rectangle mesh locally refined in the neighbourhood of (0,0)
         ******* mesh3_i has been obtained by refinement of mesh3_0

      -> mesh4_i is a quadrangle mesh with acute angles
         *******

      -> mesh5_i is a nonconforming rectangular mesh for the vertical fault 
                  problem

      -> mesh6 is  a conforming distorted quadrangle mesh
         *****

      -> mesh7 is  a nonconforming distorted quadrangle mesh
         *****
 

******************************************************************************
******                        The format of the meshes                 *******
******************************************************************************

We provide two different formats for the meshes

 *.typ1  and  *.typ2


Format .typ1 (triangular or quadrangular control volumes) includes
************
     -> the  total number of vertices
     -> the coordinates of the vertices   
     -> the total number of triangles
     -> the number of each of their three vertices ordered clockwise
     -> the total number of quadrangles
     -> the number of each of their four vertices ordered clockwise
     -> the total number of edges on the boundary of the domain
     -> the number of the vertices that define such edge
     -> the total number of all the edges ( intersection of two control
     volumes or of a control volume and the boundary of the domain)
     -> informations on each edge
                + the number of the two vertices of the edge 
                + the number of the two control volumes  defining
		the edge or for the boundary the "numero" of the
		control volume and 0
 

*******************************************************************************

Format .typ2 (polygonal control volumes) gives
************
     -> the  total number of vertices
     -> the coordinates of the vertices
     -> the number of control volumes
     -> the control volumes with
               + the total number of vertices on the boundary of volume
               + the number of each of their vertices ordered
	       counter-clockwise

*******************************************************************************
 
Example (see figure test.eps):
********* 
We consider a mesh named test with 3 control volumes
      [0,0.5]*[0.1], [0.5,1]*[0,0.5],[0.5,1]*[0.5,1]


 -> file  test.typ1
vertices
8
0.0 0.0
0.5 0.0
1.0 0.0
0.5 0.5 
1.0 0.5
0.0 1.0
0.5 1.0
1.0 1.0
triangles
0
quadrangles
3
1 2 7 6
2 3 5 4
4 5 8 7
edges of the boundary
7
1 2
1 6 
2 3
3 5
5 8
6 7
7 8
all edges
1 2 1 0
1 6 1 0
2 3 2 0
2 4 1 2
3 5 2 0
4 5 2 3
4 7 1 3
5 8 3 0
6 7 1 0
7 8 3 0
******************************************************************************     
 -> file  test.typ2
Vertices
8
0.0 0.0
0.5 0.0
1.0 0.0
0.5 0.5 
1.0 0.5
0.0 1.0
0.5 1.0
1.0 1.0
cells
3
5 1 2 4 7 6
4 2 3 5 4
4 4 5 8 7

*******************************************************************************
