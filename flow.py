from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
import pyvista
import BC
from fluids import epsilon, sigma
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
print(len(v))
p = TrialFunction(Q)
q = TestFunction(Q)

# Boundary Conditions
wall_dofs = locate_dofs_geometrical(V, BC.walls)
u_noslip = np.array((0,) * mesh.geometry.dim, dtype=PETSc.ScalarType)
bc_noslip = dirichletbc(u_noslip, wall_dofs, V)

inflow_dofs = locate_dofs_geometrical(Q, BC.inflow)
bc_inflow = dirichletbc(PETSc.ScalarType(8), inflow_dofs, Q)

outflow_dofs = locate_dofs_geometrical(Q, BC.outflow)
bc_outflow = dirichletbc(PETSc.ScalarType(0), outflow_dofs, Q)
bcu = [bc_noslip]
bcp = [bc_inflow, bc_outflow]

# Defining Navier-Stokes 
u_n = Function(V)
u_n.name = "u_n"
U = 0.5 * (u_n + u)
n = FacetNormal(mesh)
f = Constant(mesh, PETSc.ScalarType((0, 0, 0)))
print(len(f))
k = Constant(mesh, PETSc.ScalarType(dt))
mu = Constant(mesh, PETSc.ScalarType(1))
rho = Constant(mesh, PETSc.ScalarType(1))

# Define the variational problem for the first step
p_n = Function(Q)
p_n.name = "p_n"
F1 = rho * dot((u - u_n) / k, v) * dx
F1 += rho * dot(dot(u_n, nabla_grad(u_n)), v) * dx
F1 += inner(sigma(U, p_n, mu), epsilon(v)) * dx
F1 += dot(p_n * n, v) * ds - dot(mu * nabla_grad(U) * n, v) * ds
F1 -= dot(f, v) * dx
a1 = form(lhs(F1))
L1 = form(rhs(F1))

# I/O and Pathing
from pathlib import Path
folder = Path("results")
folder.mkdir(exist_ok=True, parents=True)
vtx_mesh = VTXWriter(mesh.comm, folder / "mesh.bp", mesh, engine="BP4")
vtx_mesh.write(t)




