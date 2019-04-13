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

MODEL = True

V   = fe.VectorElement("Lagrange", mesh.ufl_cell(), 2)
P   = fe.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
NU   = fe.FiniteElement("Lagrange", mesh.ufl_cell(), 2)
if MODEL: M   = fe.MixedElement([V, P, NU])
else: M   = fe.MixedElement([V, P])
W   = fe.FunctionSpace(mesh, M)

W0 = fe.Function(W)

if MODEL:
   (v, q, nu_test)    = fe.TestFunctions(W)
   (u, p, nu_trial)    = fe.split(W0)
else:
   (v, q, )    = fe.TestFunctions(W)
   (u, p, )    = fe.split(W0)
   nu_trial = fe.Constant(2)
   fv1 = fe.Constant(1)


b = fe.Expression(('0', '0'), degree=DEG)
nu = fe.Constant(1)
rho = fe.Constant(1)
RE = 100

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
    return on_boundary and np.abs(x[0] + 3) < EPSILON

# outlet
def dbc_outflow(x, on_boundary):
    return on_boundary and np.abs(x[0] - 3) < EPSILON

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
uD_Y1 = fe.Expression(('%s' % RE, '0'), degree=DEG)
bc_p  = fe.Constant(('0'))
bc_v  = fe.Constant(('0'))
bc_vf  = fe.Constant(3 * nu)

bc_1 = fe.DirichletBC(W.sub(0), uD_X0, dbc1)
bc_2 = fe.DirichletBC(W.sub(0), uD_Y0, dbc2)
bc_3 = fe.DirichletBC(W.sub(0), uD_X1, dbc3)
bc_inflow = fe.DirichletBC(W.sub(0), fe.Constant((RE, '0')), dbc_inflow)
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
   #Cw1 = (1)
   Cw1 = Cb1 / kappa * kappa + (1 + Cb2) / sigma
   Cw2 = fe.Constant(0.3)
   Cw3 = fe.Constant(2)
   Cv1 = fe.Constant(7.1)
   Ct1 = fe.Constant(1)
   Ct2 = fe.Constant(2)
   Ct3 = fe.Constant(1.1)
   Ct4 = fe.Constant(2)
   
   d1 = fe.Expression('x[0] - 0', degree=1)
   d2 = fe.Expression('1 - x[0]', degree=1)
   d3 = fe.Expression('1 - x[1]', degree=1)
   d4 = fe.Expression('x[1] - 0', degree=1)

   #d = fe.Constant(0.1)
   d = Min(Min(d1, d2), Min(d3, d4))

   xi = nu_trial / nu
   #xi = fe.Constant(1.)
   #fv1 = fe.Constant(1.)
   #fv1 = fe.Expression('xi * xi * xi / (xi * xi * xi + Cv1 * Cv1 * Cv1)', degree=DEG, xi=fe.Constant(1), Cv1=Cv1)
   fv1 = xi * xi * xi / (xi * xi * xi + Cv1 * Cv1 * Cv1)
   fv2 = fe.Constant(1.)
   #fv2 = 1 - xi / (1 + xi * fv1)
   ft2 = Ct3 * fe.exp(-Ct4 * xi * xi)
   #Omega = 1.
   Omega = 0.5 * (fe.grad(u) - fe.grad(u).T)
   S = fe.sqrt(2 * fe.inner(Omega, Omega))
   #S = 1.
   Stilde = S + nu_trial /(kappa * kappa * d * d) * fv2
   #ft2 = fe.Constant(1.)
   ft2 = Ct3 * fe.exp(- 1 * Ct4 * xi * xi)
   #r = fe.Constant(1.)
   r = Min(nu_trial / (Stilde * kappa * kappa * d * d), 10)
   g = r + Cw2 * (r * r * r * r * r * r - r)
   fw = g * ((1 + Cw3 * Cw3 * Cw3 * Cw3 * Cw3 * Cw3) / (g * g * g * g * g * g + Cw3 * Cw3 * Cw3 * Cw3 * Cw3 * Cw3)) ** (1. / 6.)
   #fw = fe.Constant(1)
   

ns_conv = fe.inner(v, fe.grad(u)*u)*fe.dx
ns_press = p * fe.div(v) * fe.dx
ns_tv = fe.inner((fv1 * nu_trial) * fe.grad(v), fe.grad(u)) * fe.dx
ns_visc = nu * fe.inner(fe.grad(v), fe.grad(u)) * fe.dx
ns_conti = q * fe.div(u) * fe.dx
ns_forcing = fe.dot(v, b)*fe.dx

NS = ns_conv + ns_press + ns_tv + ns_visc + ns_conti + ns_forcing

if MODEL:
   tv_adv = fe.inner(nu_test, fe.inner(u, fe.grad(nu_trial))) * fe.dx
   tv1 = fe.inner(nu_test, Cb1 * (1 - ft2) * Stilde * nu_trial) * fe.dx 
   tv2 = fe.inner((1 / sigma) * fe.grad(nu_test), (nu + nu_trial) * fe.grad(nu_trial)) * fe.dx 
   tv3 = fe.inner(nu_test / sigma * Cb2, fe.dot(fe.grad(nu_trial), fe.grad(nu_trial))) * fe.dx
   tv4 = fe.inner(nu_test * (Cw1 * fw - Cb1 / kappa ** 2 * ft2), (nu_trial / d) * (nu_trial / d)) * fe.dx
   #tv5 = fe.inner(ft1, DELU) * fe.dx
   
   TV = fe.inner(fe.grad(nu_trial), fe.grad(nu_test)) * fe.dx
   #TV += tv2 + tv3 + tv_adv
   TV = tv2 + tv3 + tv4
   #TV = tv1 + tv2 + tv3 + tv4 + tv_adv

weakForm  = NS

if MODEL: weakForm += TV

STAB = True
he  = fe.CellDiameter(mesh)
tau = (1.0/3.0)*(he*he)/(4.0 * nu * rho)

res = - tau * fe.inner(fe.dot(u, fe.grad(v)), fe.grad(u)*u) * fe.dx
res += - tau * fe.inner(fe.dot(u, fe.grad(v)),  fe.grad(p)) * fe.dx
res += - tau * fe.inner(fe.dot(u, fe.grad(v)),  -1 * nu_trial * fe.div(fe.grad(u))) * fe.dx
res += - tau * fe.inner(fe.dot(u, fe.grad(v)), -1 * nu * fe.div(fe.grad(u))) * fe.dx 
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
if STAB: weakForm = weakForm + stab

#weakForm += res

dW  = fe.TrialFunction(W)
dFdW = fe.derivative(weakForm, W0, dW)

if MODEL: bcSet   = [bc_1, bc_2, bc_3, bc_inflow, bc_p, bc_v_x0, bc_v_x1, bc_v_y1, bc_v_in, bc_v_out, bc_v_top]
else: bcSet   = [bc_1, bc_2, bc_3, bc_inflow, bc_p]
problem = fe.NonlinearVariationalProblem(weakForm, W0, bcSet, J=dFdW)

solver = fe.NonlinearVariationalSolver(problem)

#usol    = fe.Function(W)
prm = solver.parameters
#prm["newton_solver"]["absolute_tolerance"] = 5E-1
solver.solve()

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
if MODEL:  
   vtkFile = fe.File('nut.pvd')
   vtkFile << nu_t
#print(errornorm(usol, usol))

#plot(u, mesh=mesh)
#plt.savefig('usol')
#p#lt.clf()
#plot(p, mesh=mesh)
#plt.savefig('Psol')
#plt.clf()
