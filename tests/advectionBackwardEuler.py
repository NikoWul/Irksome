# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 13:31:10 2021

@author: nikol
"""

from fenics import * 

from irksome import GaussLegendre, RadauIIA, Dt, TimeStepper


  
from ufl.algorithms.ad import expand_derivatives

import numpy as np


msh = IntervalMesh(10, 0, 1)

V = FunctionSpace(msh, "CG", 1)


  
T = 2.0            # final time
num_steps = 20     # number of time steps
dt = T / num_steps # time step size

 
c= Constant(2.0)
u_D = Expression('x[0]-c*t', degree=1, c=c,t=0)
        
def boundary(x, on_boundary):
            return on_boundary
        
bc = DirichletBC(V, u_D, boundary)


u_n = interpolate(u_D, V)
u = TrialFunction(V)
k0= Function(V)
v0=TestFunction(V)
u0= u_n + dt*Constant(1.0)*k0     
F = u*v0*dx+u0*v0*dx- dt*c*Dx(u0,0)*v0*k0*dx
a, L = lhs(F), rhs(F)



u = Function(V)
t = 0
arrayY=[]
arrayX=[]
for n in range(num_steps):

    # Update current time
    t += dt
    u_D.t = t

    # Compute solution
    solve(a == L, u, bc)

    u_e = interpolate(u_D, V)
    #error = np.abs(u_e.vector().get_local() - u.vector().get_local()).max()
    #print('t = %.2f: error = %.3g' % (t, error))
    arrayY.append(u.vector().get_local().max())
    arrayX.append(t)

    # Update previous solution
    u_n.assign(u)

with open('advectionBackwardEuler2.txt', 'w') as file:    
    for i in range(num_steps):
        file.write('{},{}\n'.format(arrayX[i],arrayY[i]))
