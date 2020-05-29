from __future__ import print_function
from fenics import *
import matplotlib.pyplot as plt
import Mesh_Create
import uuid
# Setting physical parameters
f        = Constant( (0.0, 0.0) )
rho      = Constant( 0.001230 )
mu       = Constant( 0.000179 )
u_walls  = Constant( (0.0, 0.0) )
p_inlet  = Constant( 8000.0 )
p_outlet = Constant( 0.0 )


txt = input('digite o nome do arquivo que deseja abrir \n')
Saida =f'{uuid.uuid4()}_{txt}.pvd'
txt = f'{txt}.vtu'
mesh = Mesh_Create.create_mesh(txt)
vtkfile = File(Saida)
vtkfile << mesh

V = FunctionSpace(mesh, 'P', 1)
u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)
def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Plot solution and mesh
#plot(u)
#plot(mesh)

# Save solution to file in VTK format
vtkfile = File('PDE/solution.pvd')
vtkfile << u


# Compute error in L2 norm
error_L2 = errornorm(u_D, u, 'L2')

# Compute maximum error at vertices
vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
import numpy as np
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

# Print errors
print('error_L2  =', error_L2)
print('error_max =', error_max)

# Hold plot
plt.show()
