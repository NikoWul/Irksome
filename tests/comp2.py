# -*- coding: utf-8 -*-
"""
Created on Wed May 26 16:44:16 2021

@author: nikol
"""
from fenics import * 

from irksome import GaussLegendre, RadauIIA, Dt, TimeStepper


  
from ufl.algorithms.ad import expand_derivatives

import numpy as np

N = 100

dt = Constant(10.0 / N)
t = Constant(0.0)
nx = ny = 8
msh = UnitSquareMesh(nx, ny)
V = FunctionSpace(msh, "CG", 1)
#print(V)
E=V.ufl_element()*V.ufl_element()
Vbig= FunctionSpace(V.mesh(),E)
k = Function (Vbig)
#print(k.__repr__())
k0, k1 = split(k)
#print(k0.__repr__())
#print(k1.__repr__())
v0, v1 = TestFunctions(Vbig)
print(v0.__repr__())
print(v1.__repr__())
u = TrialFunction(V)

v = TestFunction(V)

u0 = u + Constant (0.5) * dt * k0 + Constant ( -0.5)  * dt * k1
print(u0)
u1 = u + Constant (0.5) * dt * k0 + Constant (0.5) * dt * k1
print(u1)
F = (inner(k0 , v0) * dx + inner(k1, v1) * dx +inner(grad(u0), grad(v0)) * dx + inner(grad(u1), grad(v1)) * dx)
print(F)