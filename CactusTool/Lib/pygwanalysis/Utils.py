#----------------------------------------------------------
#
# Utility functions for DiscreteFunction
#
# Released under the MIT License.
# (C) Christian Reisswig 2009-2011
#
#----------------------------------------------------------


import sys

from math import *
from numpy import *





#   Integration coefficients for a fourth
#   order scheme (variant Simpson's rule).
#   N+1 = total number of points
#   i = current point
def Simpson_coeff(i, N):

    if (i < 0):
        return


    if (N == 0): return 0.0             # trivial
    if (N == 1): return 0.5             # trapezoid rule
    if (N == 2):                        # Simpson's rule
        if (i == 0): return 1.0/3.0
        if (i == 1): return 4.0/3.0
        if (i == 2): return 1.0/3.0

    if ((N == 3) or (N == 5)):                # trapezoid rule
        if ((i == 0) or (i == N)): return 0.5
        else: return 1.0

    if (N == 4):                            # Simpson's rule
        if (i == 0): return 1.0/3.0
        if (i == 1): return 4.0/3.0
        if (i == 2): return 2.0/3.0
        if (i == 3): return 4.0/3.0
        if (i == 4): return 1.0/3.0

    if (N == 6):                            # Simpson's rule
        if (i == 0): return 1.0/3.0
        if (i == 1): return 4.0/3.0
        if (i == 2): return 2.0/3.0
        if (i == 3): return 4.0/3.0
        if (i == 4): return 2.0/3.0
        if (i == 5): return 4.0/3.0
        if (i == 6): return 1.0/3.0

    # else use variant Simpson's rule
    if ((i == 0) or (i == N  )): return 17.0/48.0
    if ((i == 1) or (i == N-1)): return 59.0/48.0
    if ((i == 2) or (i == N-2)): return 43.0/48.0
    if ((i == 3) or (i == N-3)): return 49.0/48.0
    if ((i >  3) or (i <  N-3)): return 1.0

    return 0.0


def Lagrange_derivative(x, xi, fi, order=2):
    # return first derivative of second-order Lagrange interpolant

    if order==2:
        a0 = (2*x - xi[1] - xi[2]) / ((xi[0]-xi[1]) * (xi[0]-xi[2]))
        a1 = (2*x - xi[0] - xi[2]) / ((xi[1]-xi[0]) * (xi[1]-xi[2]))
        a2 = (2*x - xi[0] - xi[1]) / ((xi[2]-xi[0]) * (xi[2]-xi[1]))
        return a0*fi[0] + a1*fi[1] + a2*fi[2]

    if order==4:
        x0 = xi[0]
        x1 = xi[1]
        x2 = xi[2]
        x3 = xi[3]
        x4 = xi[4]
        a0 = -((-4*x**3 + x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 + 3*x**2*(x1 + x2 + x3 + x4) - 2*x*(x3*x4 + x2*(x3 + x4) + x1*(x2 + x3 + x4)))/((x0 - x1)*(x0 - x2)*(x0 - x3)*(x0 - x4)))
        a1 = (-4*x**3 + x2*x3*x4 + 3*x**2*(x0 + x2 + x3 + x4) + x0*(x3*x4 + x2*(x3 + x4)) - 2*x*(x3*x4 + x2*(x3 + x4) + x0*(x2 + x3 + x4)))/((x0 - x1)*(x1 - x2)*(x1 - x3)*(x1 - x4))
        a2 = -((-4*x**3 + x0*x1*x3 + x0*x1*x4 + x0*x3*x4 + x1*x3*x4 + 3*x**2*(x0 + x1 + x3 + x4) - 2*x*(x3*x4 + x1*(x3 + x4) + x0*(x1 + x3 + x4)))/((x0 - x2)*(x1 - x2)*(x2 - x3)*(x2 - x4)))
        a3 = (-4*x**3 + x1*x2*x4 + 3*x**2*(x0 + x1 + x2 + x4) + x0*(x2*x4 + x1*(x2 + x4)) - 2*x*(x2*x4 + x1*(x2 + x4) + x0*(x1 + x2 + x4)))/((x1 - x3)*(-x0 + x3)*(-x2 + x3)*(x3 - x4))
        a4 = (-4*x**3 + x1*x2*x3 + 3*x**2*(x0 + x1 + x2 + x3) + x0*(x2*x3 + x1*(x2 + x3)) - 2*x*(x2*x3 + x1*(x2 + x3) + x0*(x1 + x2 + x3)))/((x1 - x4)*(-x0 + x4)*(-x2 + x4)*(-x3 + x4))
        return a0*fi[0] + a1*fi[1] + a2*fi[2] + a3*fi[3] + a4*fi[4]
    
    return 0



def find_root(x, f, val, eps):
    # find root of function (f-val) via bisection in intervall x[0],x[1]
    #eps = 1e-4 # accuracy
    a = x[0]
    b = x[1]
    while (2*eps < abs(b-a)):
        midpoint = (a+b)/2
        if (f(a)*f(b) > 0):
            a = midpoint
        else:
            b = midpoint
    return (a+b)/2


