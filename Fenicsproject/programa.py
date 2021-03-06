from  __future__ import print_function
from fenics import *
import Mesh_Create
import uuid
# Setting physical parameters
f        = Constant( (0.0, 0.0) )
rho      = Constant( 0.001230 )
mu       = Constant( 0.000179 )
u_walls  = Constant( (0.0, 0.0) )
p_inlet  = Constant( 8000.0 )
p_outlet = Constant( 0.0 )


txt =input("digite o arquivo que quer abrir \n")
print(txt)
print('Loading mesh')
Saida =f'{uuid.uuid4()}_{txt}.pvd'
my_mesh = Mesh_Create.create_mesh(txt)
vtkfile = File("mesh.pvd")
vtkfile << mesh

# Define function spaces
print('Defining trial and test functions')
P2 = VectorElement('Lagrange', my_mesh.ufl_cell(), 2)
P1 = FiniteElement('Lagrange', my_mesh.ufl_cell(), 1)
TH = P2 * P1  # <- Taylor-Hood mixed element
W  = FunctionSpace(my_mesh, TH)

# Define test and trial functions
(u, p) = TrialFunctions(W)
(v, q) =  TestFunctions(W)

# Define boundary conditions
# --> Lucas, aqui vc deve definir as condições de contorno
# --> eu deixei o código mais ou menos preparado
# --> as paredes devem ter velocidade NULA
print('Defining boundary conditions...')

# --> define entrada de ar: pontos da superfície com x == 0 (corrija)
class loc_inlet(SubDomain):
    def inside(self, x, on_boundary):
        return (on_boundary and near(x[0], 0.0))

# --> define entrada de ar: pontos da superfície com x == 10.0 (corrija)
class loc_outlet(Subdomain):
    def inside(self, x, on_boundary):
        return (on_boundary and near(x[0], 10.0))

# --> define as paredes da laringe: todos os pontos (depois vamos filtrá-los)
# --> não precisa mexer aqui
class loc_walls(Subdomain):
    def inside(self, x, on_boundary):
        return on_boundary

# --> recuperando os pontos da superfície: vamos marcar os pontos da superfície
# --> de acordo com as regras estabelecidas acima
my_surface = MeshFunction('size_t', my_mesh, mesh.topology().dim()-1)
my_surface.set_all(0)

# pontos com 1 => paredes, com 2 => "pulmão", com 3 => "boca"
subdomain_loc = loc_walls()
subdomain_loc.mark(my_surface, 1)
subdomain_loc = loc_inlet()
subdomain_loc.mark(my_surface, 2)
subdomain_loc = loc_outlet()
subdomain_loc.mark(my_surface, 3)

bcu_walls = DirichletBC(W.sub(0), u_walls, my_surface, 1)
bcs = [bcu_walls]

# Define variational problem a(u, V) = L(v)
print('Solving a(u, V) = L(v)')
ds_surface = Measure('ds')[my_surface]
n = FacetNormal(my_mesh)
a = mu*inner(grad(v), grad(u)) - div(v)*p + q*div(u)*dx 
L = rho*inner(v, f)*dx - p_inlet*dot(v, n)*ds_surface(2) - p_outlet*dot(v, n)*ds_surface(3)

# Compute solution using default solver with default parameters
w = Function(W)
my_problem = LinearVariationalProblem(a, L, w, bcs = bcu)
my_solver  = LinearVariationalSolver(my_problem)

my_solver.solve()

# --> se esta linha der erro, ignore ela
list_timings(TimingClear.clear, [TimingType.user])

# Saving solution
print('Saving solution')
(u, p) = w.split()
with XDMFFile(MPI.comm_world, f'{Saida}.res_pressure.xdmf') as my_file:
	my_file.write(p)
with XDMFFile(MPI.comm_world,f'{Saida}res_velocity.xdmf') as my_file:
	my_file.write(u)

# para executar esse código, vamos usar o mpirun (coloquei para usar 2 threads)
# mpirun -np 2 python3 nome_deste_codigo.py
