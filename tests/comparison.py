# -*- coding: utf-8 -*-
"""
Created on Wed May 26 08:22:18 2021

@author: nikol
"""

from fenics import * 

from irksome import GaussLegendre, RadauIIA, Dt, TimeStepper, LobattoIIIC


  
from ufl.algorithms.ad import expand_derivatives

import numpy as np

#implicit midpoint rule
butcher_tableau =  LobattoIIIC(2)

ns = butcher_tableau.num_stages  


N = 100
x0 = 0.0
x1 = 10.0
y0 = 0.0
y1 = 10.0
nx = ny = 8
msh = UnitSquareMesh(nx, ny)
V = FunctionSpace(msh, "CG", 1)
#print(V)


  
dt = Constant(10.0 / N)
t = Constant(0.0)
 


#x, y = SpatialCoordinate(msh)
#S = Constant(2.0)
#C = Constant(1000.0)
#B = (x-Constant(x0))*(x-Constant(x1))*(y-Constant(y0))*(y-Constant(y1))/C
#R = (x * x + y * y) ** 0.5
#uexact = B * atan(t)*(pi / 2.0 - atan(S * (R - t)))
#rhs = expand_derivatives(diff(uexact, t)) - div(grad(uexact))
alpha = 3          # parameter alpha
beta = 1.2         # parameter beta
u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t*t',
                         degree=2, alpha=alpha, beta=beta, t=0)
        

#initial condition:


u_n = interpolate(u_D, V)

def boundary(x, on_boundary):
            return on_boundary
        
bc = DirichletBC(V, u_D, boundary)

#Now, we will define the semidiscrete variational problem using
#standard UFL notation, augmented by the ``Dt`` operator from Irksome::
u = TrialFunction(V)
v = TestFunction(V)

F = inner(Dt(u), v) * dx + inner(grad(u), grad(v)) * dx
bc = DirichletBC(V, 0, "on_boundary")



luparams = {"mat_type": "aij",
              "ksp_type": "preonly",
              "pc_type": "lu"}

#Most of Irksome's magic happens in the :class:`.TimeStepper`.  It
#transforms our semidiscrete form `F` into a fully discrete form for
#the stage unknowns and sets up a variational problem to solve for the
#stages at each time step.::

stepper = TimeStepper(F, butcher_tableau, t, dt, u_n,
                        solver_parameters=luparams)

#This logic is pretty self-explanatory.  We use the
#:class:`.TimeStepper`'s :meth:`~.TimeStepper.advance` method, which solves the variational
#problem to compute the Runge-Kutta stage values and then updates the solution.::
		
while (float(t) < 1.0):
    if (float(t) + float(dt) > 1.0):
        dt.assign(1.0 - float(t))
    stepper.advance()
    u_e = interpolate(u_D, V)
       
    error = np.abs(u_e.vector().get_local() - u_n.vector().get_local()).max()
    print('t = %.2f: error = %.3g' % (t, error))
    print(float(t))
    t.assign(float(t) + float(dt))

#Finally, we print out the relative :math:`L^2` error::

print()