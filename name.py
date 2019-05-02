import fenics as fe
import dolfin
import numpy as np
from dolfin.fem.norms import errornorm
from dolfin.common.plotting import plot
import matplotlib.pyplot as plt
import sys

EPSILON = 1.0e-14
DEG = 2

mesh = fe.Mesh('step.xml')

# Control pannel
MODEL = False # flag to use SA model. False defaults to a mixing length model.
UNSTEADY = False
STAB = False
b = fe.Expression(('0', '0'), degree=DEG) # forcing
nu = fe.Constant(1e-2)
rho = fe.Constant(1)
RE = 20
dt = 1.0
# Re = 10 / 1e-4 = 1e5

V   = fe.VectorElement("Lagrange", mesh.ufl_cell(), 2)
P   = fe.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
NU   = fe.FiniteElement("Lagrange", mesh.ufl_cell(), 2)
if MODEL: M   = fe.MixedElement([V, P, NU])
else: M   = fe.MixedElement([V, P])
W   = fe.FunctionSpace(mesh, M)
We = fe.Function(W)

W0 = fe.Function(W)
if MODEL: u0, p0, nu_0 = fe.split(We)
else: u0, p0 = fe.split(We)

if MODEL:
   (v, q, nu_test)    = fe.TestFunctions(W)
   (u, p, nu_trial)    = fe.split(W0)
else:
   (v, q, )    = fe.TestFunctions(W)
   (u, p, )    = fe.split(W0)
   nu_trial = fe.Constant(5) # artificial viscosity!!!
   fv1 = fe.Constant(1)



#-------------------------------------------------------
# Defining essential/Dirichlet boundary conditions
# Step 1: Identify all boundary segments forming Gamma_d
#-------------------------------------------------------
#   (-3., 2.5)
#  |
#  |
#  |_______(0, 1.)
#    bc1  |
#      bc2|__________(3., ,0.)
#        (0,0) bc3

# surface before step
def dbc1(x, on_boundary):
    return on_boundary and np.abs(x[1] - 1.) < EPSILON and x[0] < EPSILON

# surface on step side
def dbc2(x, on_boundary):
    return on_boundary and x[0] < EPSILON and x[1] < 1.0

# surface after step
def dbc3(x, on_boundary):
    return on_boundary and x[1] < EPSILON and x[0] > - 1 * EPSILON

# inflow
def dbc_inflow(x, on_boundary):
    return on_boundary and np.abs(x[0] + .2) < EPSILON

# outlet
def dbc_outflow(x, on_boundary):
    return on_boundary and np.abs(x[0] - 10) < EPSILON

# top
def dbc_top(x, on_boundary):
    return on_boundary and np.abs(x[1] - 2.5) < EPSILON

#--------------------------------------------------------
# Defining essential/Dirichlet boundary conditions
# Step 2: Defining what the boundary values will be (u_D)
#--------------------------------------------------------
uD_X0 = fe.Expression(('0', '0'), degree=DEG)
uD_X1 = fe.Expression(('0', '0'), degree=DEG)
uD_Y0 = fe.Expression(('0', '0'), degree=DEG)
#uD_Y1 = fe.Expression(('%s * pow((2.5 - x[1]), 2)' % RE, '0'), degree=DEG)
uD_Y1 = fe.Expression(('%s' % RE, '0'), degree=DEG)
bc_p  = fe.Constant(('0'))
bc_v  = fe.Constant(('0'))
bc_vf  = fe.Constant(5 * nu)

bc_1 = fe.DirichletBC(W.sub(0), uD_X0, dbc1)
bc_2 = fe.DirichletBC(W.sub(0), uD_Y0, dbc2)
bc_3 = fe.DirichletBC(W.sub(0), uD_X1, dbc3)
bc_inflow = fe.DirichletBC(W.sub(0), fe.Expression(('%s * pow((x[1] - 1) / 2.5, 2)' % RE, '0'), degree=DEG), dbc_inflow)
#bc_inflow = fe.DirichletBC(W.sub(0), fe.Constant((RE, '0')), dbc_inflow)
bc_p = fe.DirichletBC(W.sub(1), bc_p, dbc1)
if MODEL:
   bc_v_x1 = fe.DirichletBC(W.sub(2), bc_v, dbc1)
   bc_v_x0 = fe.DirichletBC(W.sub(2), bc_v, dbc2)
   bc_v_y1 = fe.DirichletBC(W.sub(2), bc_v, dbc3)
   bc_v_in = fe.DirichletBC(W.sub(2), bc_vf, dbc_inflow)
   bc_v_out = fe.DirichletBC(W.sub(2), bc_vf, dbc_outflow)
   bc_v_top = fe.DirichletBC(W.sub(2), bc_vf, dbc_top)

def Max(a, b): return (a + b + abs(a-b)) / 2.
def Min(a, b): return (a + b - abs(a-b)) / 2.

if MODEL:
   sigma = fe.Constant(2./3.)
   Cb1 = fe.Constant(0.1355)
   Cb2 = fe.Constant(0.622)
   kappa = fe.Constant(0.41)
   Cw1 = Cb1 / kappa * kappa + (fe.Constant(1) + Cb2) / sigma
   Cw2 = fe.Constant(0.3)
   Cw3 = fe.Constant(2)
   Cv1 = fe.Constant(7.1)
   Ct1 = fe.Constant(1)
   Ct2 = fe.Constant(2)
   Ct3 = fe.Constant(1.1)
   Ct4 = fe.Constant(2)
   
   d1 = fe.Expression('x[0] - 0', degree=1)
   #d2 = fe.Expression('1 - x[0]', degree=1)
   #d3 = fe.Expression('1 - x[1]', degree=1)
   d = fe.Expression('x[1] - 0', degree=1)

   #d = fe.Constant(10)
   #d = Min(d1, d4)

   xi = nu_trial / nu
   fv1 = xi ** 3 / (xi ** 3 + Cv1 * Cv1 * Cv1)
   #fv1 = fe.elem_pow(xi, 3) / (fe.elem_pow(xi, 3) + Cv1 * Cv1 * Cv1)
   fv2 = 1 - xi / (1 + xi * fv1)
   ft2 = Ct3 * fe.exp(-Ct4 * xi * xi)
   Omega = fe.Constant(0.5) * (fe.grad(u) - fe.grad(u).T)
   S = fe.sqrt(EPSILON + fe.Constant(2) * fe.inner(Omega, Omega)) 
   Stilde = S + nu_trial /(kappa * kappa * d * d) * fv2
   ft2 = Ct3 * fe.exp(fe.Constant(-1) * Ct4 * xi * xi)
   r = Min(nu_trial / (Stilde * kappa * kappa * d * d), 10)
   g = r + Cw2 * (r * r * r * r * r * r - r)
   fw = fe.elem_pow(EPSILON + abs(g * ((1 + Cw3 * Cw3 * Cw3 * Cw3 * Cw3 * Cw3) / (g * g * g * g * g * g + Cw3 * Cw3 * Cw3 * Cw3 * Cw3 * Cw3))), fe.Constant(1. / 6.))
   

ns_conv = fe.inner(v, fe.grad(u)*u)*fe.dx
ns_press = p * fe.div(v) * fe.dx
s = fe.grad(u) + fe.grad(u).T
if MODEL: 
   ns_tv = fe.inner((fv1 * nu_trial) * fe.grad(v), fe.grad(u)) * fe.dx
else: 
   #lmx = 0.1
   lmx = fe.Expression(' 0.02 + 0.28 * (x[0] + x[1] - abs(x[1] - x[0])) / 2.0', degree=2)
   sij = 0.5 * (fe.grad(u) + fe.grad(u).T)
   nu_tv = 2 * lmx ** 2 * (EPSILON + 2 * fe.inner(sij, sij)) ** (0.5)
   #ns_tv = 0
   ns_tv = fe.inner((nu_tv) * fe.grad(v), fe.grad(u)) * fe.dx
ns_visc = nu * fe.inner(fe.grad(v), fe.grad(u)) * fe.dx
ns_conti = q * fe.div(u) * fe.dx
ns_forcing = fe.dot(v, b)*fe.dx

NS = ns_conv + ns_press + ns_tv + ns_visc + ns_conti + ns_forcing

N = 5
fe.parameters["form_compiler"]["quadrature_degree"] = N
if MODEL:
   #dx(metadata={'quadrature_degree':N}) # maybe 10. <10?
   tv_adv = fe.inner(nu_test, fe.inner(u, fe.grad(nu_trial))) * fe.dx(metadata={'quadrature_degree':N})
   tv1 = fe.inner(nu_test, Cb1 * (1 - ft2) * 1 * nu_trial) * fe.dx(metadata={'quadrature_degree':N}) # TODO Missing stilde term
   #tv1 = fe.inner(nu_test, Cb1 * (1 - ft2) * Stilde * nu_trial) * fe.dx(metadata={'quadrature_degree':N})
   tv2 = fe.inner((1 / sigma) * fe.grad(nu_test), (nu + nu_trial) * fe.grad(nu_trial)) * fe.dx(metadata={'quadrature_degree':N})
   tv3 = fe.inner(nu_test / sigma * Cb2, fe.dot(fe.grad(nu_trial), fe.grad(nu_trial))) * fe.dx(metadata={'quadrature_degree':N})
   #tv4 = fe.inner(nu_test * (Cw1 * fw - Cb1 / kappa / kappa * ft2), (nu_trial / d) * (nu_trial / d)) * fe.dx
   tv4 = fe.inner(nu_test, (Cw1 - Cb1 / kappa / kappa * ft2) * (nu_trial / d) **2) * fe.dx(metadata={'quadrature_degree':N}) # TODO - missing fw term
   #tv4 = fe.inner(nu_test, (Cw1 * fw - Cb1 / kappa / kappa * ft2) * (nu_trial / d) **2) * fe.dx(metadata={'quadrature_degree':N})
   #tv5 = fe.inner(ft1, DELU) * fe.dx
   
   #TV = fe.inner(fe.grad(nu_trial), fe.grad(nu_test)) * fe.dx
   TV = tv2 - tv3 + tv_adv + tv4 - tv1

   #TV =  tv2 - tv3 + tv4 + tv_adv
   #TV =  tv2 - tv3 + tv4 + tv_adv
   #TV = -1 * tv1 + tv2 - tv3 + tv4 + tv_adv

weakForm  = NS 
if UNSTEADY: weakForm  += (1.0 / dt) * fe.inner(u-u0, v) * fe.dx 
if MODEL: weakForm += TV
if MODEL and UNSTEADY: weakForm += (1.0 / dt) * fe.inner(nu_trial-nu_0, nu_test) * fe.dx


he  = fe.CellDiameter(mesh)
tau = (1.0/3.0)*(he*he)/(4.0 * (nu + nu_trial) * rho)
h  = fe.CellDiameter(mesh)
vnorm    = fe.sqrt(fe.dot(u0, u0))
tau      = ((2.0*1.0/dt)**2 + (2.0*vnorm/h)**2 + (4.0*nu/h**2)**2)**(-0.5)


res = - tau * fe.inner(fe.dot(u, fe.grad(v)), fe.grad(u)*u) * fe.dx
res += - tau * fe.inner(fe.dot(u, fe.grad(v)),  fe.grad(p)) * fe.dx
res += - tau * fe.inner(fe.dot(u, fe.grad(v)),  -1 * fv1 * nu_trial * fe.div(fe.grad(u))) * fe.dx
res += - tau * fe.inner(fe.dot(u, fe.grad(v)), -1 * (nu + nu_trial * fv1) * fe.div(fe.grad(u))) * fe.dx 
res += - tau * fe.inner(fe.dot(u, fe.grad(q)), fe.div(u)) * fe.dx
res += - tau * fe.inner(fe.dot(u, fe.grad(v)), -1 * b) * fe.dx 

if MODEL and False:

   # adv
   res += - tau * fe.inner(fe.dot(u, fe.grad(nu_test)), fe.inner(u, fe.grad(nu_trial))  ) * fe.dx 

   # t1
   res += - tau * fe.inner(fe.dot(u, fe.grad(nu_test)), -1 * Cb1 * (1 - ft2) * Stilde * nu_trial ) * fe.dx 

   # t2
   res += - tau * fe.inner(fe.dot(u, fe.grad(nu_test)), fe.div((1 / sigma) * (nu + nu_trial) * fe.grad(nu_trial))) * fe.dx 

   # t3
   res += - tau * fe.inner(fe.dot(u, fe.grad(nu_test)),   -1 * Cb2 / sigma * fe.inner(fe.grad(nu_trial), fe.grad(nu_trial))) * fe.dx 

   # t4
   res += - tau * fe.inner(fe.dot(u, fe.grad(nu_test)),   
   (Cw1 * fw - Cb1 / kappa / kappa * ft2) * (nu_trial / d) * (nu_trial / d)) * fe.dx 

stab     = -tau*fe.inner(fe.grad(q), fe.grad(p))*fe.dx
if STAB: 
   weakForm += res
#weakForm = weakForm + stab

dW  = fe.TrialFunction(W)
dFdW = fe.derivative(weakForm, W0, dW)

if MODEL: bcSet   = [bc_1, bc_2, bc_3, bc_inflow, bc_p, bc_v_x0, bc_v_x1, bc_v_y1, bc_v_in]
else: bcSet   = [bc_1, bc_2, bc_3, bc_inflow, bc_p]
problem = fe.NonlinearVariationalProblem(weakForm, W0, bcSet, J=dFdW)

solver = fe.NonlinearVariationalSolver(problem)

#usol    = fe.Function(W)
prm = solver.parameters
if not UNSTEADY: prm["newton_solver"]["maximum_iterations"] = 200

t = 0.0
t_end = 6.1
pFile = fe.File('Pressure.pvd')
uFile = fe.File('Velocity.pvd')
nFile = fe.File('nut.pvd')
vFile = fe.File('Vorticity.pvd')
wFile = fe.File('W.pvd')
T = fe.FunctionSpace(mesh, 'CG', 1)
if not UNSTEADY: solver.solve()
else:
 ii = 0
 while t < t_end:
    print("t =",t)
    solver.solve()
    if MODEL: u1,p1,n1 = W0.split()
    else: u1,p1 = W0.split()
    We.assign(W0)
    if True:
      uFile << u1
      pFile << p1
      if MODEL: nFile << n1
    #wFile << W0
      #omega = fe.curl(u)
      #vFile << fe.project(fe.inner(omega, omega), T)
    t += dt
    #if t > 1.5: dt = 0.1
    ii += 1

if MODEL:
   u, p, nu_t = W0.split()
else:
   u, p = W0.split()
#-------------------------------------------------
# Save this solution to a file for post-processing
#-------------------------------------------------
vtkFile = fe.File('u.pvd')
vtkFile << u
vtkFile = fe.File('p.pvd')
vtkFile << p
if MODEL: vtkFile = fe.File('nut.pvd')
if MODEL: vtkFile << nu_t

#T = fe.TensorFunctionSpace(mesh, 'CG', 1)
#vtkFile = fe.File('tau.pvd')
#vtkFile << fe.project(nu * (fe.grad(u) + fe.grad(u).T), T)

if MODEL:  

   T = fe.FunctionSpace(mesh, 'CG', 1)
   #T = fe.TensorFunctionSpace(mesh, 'CG', 1)
   vtkFile = fe.File('fw.pvd')
   vtkFile << fe.project(fw, T)
   vtkFile = fe.File('g.pvd')
   vtkFile << fe.project(g, T)
   vtkFile = fe.File('r.pvd')
   vtkFile << fe.project(r, T)
   vtkFile = fe.File('Stilde.pvd')
   vtkFile << fe.project(Stilde, T)
#vtkFile << fe.project( (nu_trial / (Stilde * kappa ** 2 * d ** 2 )) ** 6, T)
#vtkFile << fe.project( 1. / (fe.elem_pow(kappa, 2) * fe.elem_pow(d, 2) ), T)
#vtkFile << fe.project(g ** (1./6.), T)

#print(errornorm(usol, usol))

#plot(u, mesh=mesh)
#plt.savefig('usol')
#p#lt.clf()
#plot(p, mesh=mesh)
#plt.savefig('Psol')
#plt.clf()
