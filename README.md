#   Matlab Finite Element routines for the course Computational Geomechanics - Civil 423

This is a set of matlab routines for 2D (plane-strain & axi-symmetric) Finite Element analysis for laplacian, diffusion, elastic & coupled problems.

© All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Lab, 2015-2019

authors:
Brice Lecampion

See the LICENSE.TXT file for more details.

Weak form operators currently coded:
- Elasticity: stiffness,  pre stress field , boundary loads
- Laplacian
- Mass matrix
- poroelastic Coupling operator 

Problem type currently coded up:
- 2D plane-strain 
- 2D axi-symmetry

Available Element type:
- Linear triangular Element (2D) : Tri3 (for flow and elasticity), Tri6 (for elasticity, poroelasticity)


For mesh generation, either ones directly script it in matlab (for simple geometries), or use available matlab libraries e.g. MESH2D - Delaunay-based unstructured mesh-generation 
(https://ch.mathworks.com/matlabcentral/fileexchange/25555-mesh2d-delaunay-based-unstructured-mesh-generation)


Element classes - contained the B-matrices, isoparametric mapping etc.

Assembly routines coded as function
