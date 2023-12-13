# Explanation

## Intro
We begin this with some motivation for the code. In general, it is advantagous for there to be easy and quick way to setup and test cardiovascular flow dynamics. Because of this, we are working on a pipeline from geometry definition to outputed video animation. 

## Geometry 
The first step in the pipeline uses the meshing algorithm gmsh to build and optimize a mesh from a defined geometry (.geo) file. To triangulate the space, Delaunay triangulation is used as a way to ensure that the minimum interior angle witin each element (tet or triangle) is maximized. 

## CFD 
Once the mesh is complete, is can be passed to the flow solver. This uses the open-source finite-element based solver FENICS. More specifically, it uses the fluid specific solver Dolfinx within FENICS. Both velocity and pressure will be solved for from the traditional Navier-Stokes equation. A common method for doing this involves solving for the next timestep for velocity, and then pressure, and then making a correction to the velocity using the new pressure. Two finite element spaces are used, one using linear Lagrangian elements for the pressure and one using quadratic Lagrangian elements for the velocity. Since not all combinations of finite element spaces are stable for N-S, this choise is made deliberately to help ensure solver stability. 
To actually solve the nodal degrees of freedom, a linear system is built from the mesh. Each of the three previously described steps are then solved using various PETSc solvers. The first two use stablized bi-conjugent gradient method with algebraic multigrid preconditioners. The last step uses regular conjugent gradient method with a successive over-relaxation preconditioner. 



