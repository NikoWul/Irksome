import pytest
from fenics import *
from math import isclose
from numpy import sqrt, array, allclose
from irksome import LobattoIIIA, GaussLegendre, Dt, TimeStepper
from ufl.algorithms.ad import expand_derivatives


#@pytest.mark.parametrize("butcher_tableau", [GaussLegendre(3), LobattoIIIA(2)])
#def test_1d_heat_dirichletbc(butcher_tableau):

bt = GaussLegendre(2)
A = array([[1 / 4, 1 / 4 - sqrt(3) / 6],
            [1 / 4 + sqrt(3) / 6, 1 / 4]])
b = array([1 / 2, 1 / 2])
    # btilde = array([1 / 2 + sqrt(3) / 6, 1 / 2 - sqrt(3) / 6])
c = array([1 / 2 - sqrt(3) / 6, 1 / 2 + sqrt(3) / (6)])
    # assert allclose(btilde, bt.btilde)
for (X, Y) in zip([A, b, c], [bt.A, bt.b, bt.c]):
    assert allclose(X, Y)

# Boundary values
u_0 = Constant(2.0)
u_1 = Constant(3.0)
T=2.0        
N = 50
x0 = 0.0
x1 = 10.0
nx = ny = 8
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, 'P', 1)
dt = Constant(10.0 / N)
t = 0
#(x,) = SpatialCoordinate(mesh)
    
# Method of manufactured solutions copied from Heat equation demo.
alpha = 3          # parameter alpha
beta = 1.2         # parameter beta
# Note end linear contribution
u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t*t',
                         degree=2, alpha=alpha, beta=beta, t=0)

u_n = interpolate(u_D, V)
u = TrialFunction(V)
v = TestFunction(V)
f =  Expression('beta*t-2-2*alpha',
                         degree=1, alpha=alpha, beta=beta, t=0)
F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
def boundary(x, on_boundary):
            return on_boundary
        
bc = DirichletBC(V, u_D, boundary)
    
luparams = {"mat_type": "aij", "ksp_type": "preonly", "pc_type": "lu"}
#butcher_tableau= ButcherTableaux()

stepper = TimeStepper(
        F, bt, t, dt, u_n, bcs=bc, solver_parameters=luparams
    )
    
t_end = 2.0
while float(t) < t_end:
    if float(t) + float(dt) > t_end:
        dt.assign(t_end - float(t))
        stepper.advance()
        t.assign(float(t) + float(dt))
    # Check solution and boundary values
    print( norm(u - uexact) / norm(uexact) < 10.0 ** -5)
    #assert isclose(u.at(x0), u_0)
    #assert isclose(u.at(x1), u_1)
