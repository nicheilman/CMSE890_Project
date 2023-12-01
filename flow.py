from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
import pyvista
#import gmsh

from dolfinx.fem import Constant, Function, FunctionSpace, assemble_scalar, dirichletbc, form, locate_dofs_geometrical
from dolfinx.fem.petsc import assemble_matrix, assemble_vector, apply_lifting, create_vector, set_bc
from dolfinx.io import VTXWriter
from dolfinx.mesh import create_unit_cube
from dolfinx.plot import vtk_mesh
from ufl import (FacetNormal, FiniteElement, Identity, TestFunction, TrialFunction, VectorElement,
                 div, dot, ds, dx, inner, lhs, nabla_grad, rhs, sym)

# Spacetime discretization for meshing
nx = 10
mesh = create_unit_cube(MPI.COMM_WORLD, 10, 10, 10)
t = 0
T = 10
num_steps = 500
dt = T / num_steps

# Define function spaces
v_cg2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
s_cg1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, v_cg2)
Q = FunctionSpace(mesh, s_cg1)

# Create the test and trial function
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)


# I/O and Pathing
from pathlib import Path
folder = Path("results")
folder.mkdir(exist_ok=True, parents=True)
vtx_mesh = VTXWriter(mesh.comm, folder / "mesh.bp", mesh, engine="BP4")
vtx_mesh.write(t)


