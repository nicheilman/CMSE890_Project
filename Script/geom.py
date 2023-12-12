from pathlib import Path
from mpi4py import MPI
from dolfinx.io import VTXWriter, gmshio
import gmsh

gmsh.initialize()
gmsh.model.add("aorta")
gmsh.open("aorta.geo")

gmsh.model.mesh.generate(3)
#gmsh.model.mesh.refine()
gmsh.write("aorta_mesh_fine.msh")

gmsh_model_rank = 0
mesh_comm = MPI.COMM_WORLD
domain, cell_markers, facet_markers = gmshio.model_to_mesh(gmsh.model, mesh_comm, gmsh_model_rank, gdim=3)
gmsh.finalize()

folder = Path("../results")
folder.mkdir(exist_ok=True, parents=True)
vtx_mesh = VTXWriter(mesh_comm, folder / "test.bp", domain, engine="BP4")
vtx_mesh.write(0.0)
