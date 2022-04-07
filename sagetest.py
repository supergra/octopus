#!/usr/bin/env python

'''
This was an experiment with using "sagemath" to do the solving
instead of my own solver octopus

www.sagemath.org
'''

import sage
from sage import var
from math import cos, radians
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

var('rho CL V Sref bank mass L')

eq1 = L == 0.5*rho*V*V*CL*Sref
eq2 = L*cos(bank) == 9.81*mass

# eq3 = Lstall == 0.5*rho*V*V*CL*Sref

'''
So notice, the same *function* can be used multiple times, with
different variables.

vCLstall = var()
vLstall = var()
vvstall = var()
# Note here, we keep rho and Sref the same variables from earlier, but
# specify a different velocity + CL flight condition
e3 = eqn(L, CL=vCLstall, L=vLstall, V=vvstall, rho=vrho, Sref=vsref)

# Solve a single equation for *any* variable, given the others
# Note: the arguments here are the *variable* names created above, *NOT* the
      arguments to the original function
'''

solve([e1], vCL=1.0, vL=122.5, vv=10, vsref=2)

solve([e1], vrho=1.225, vL=122.5, vv=10, vsref=2)

# Solve two equations, given all but two variables specfied:
e2 = eqn(Lbank, mass=vm, Lbank=vL, bankAngleRad=vbank)
x = solve([e1,e2], vrho=1.15, vv=10, vsref=2, vbank=radians(45), vm=10)


solns = sage.solve([eq1,eq2,rho==1.15,V==10,Sref==2,bank==radians(45), mass==10],L,CL, solution_dict=True)
bits = 30
[[s[L].n(bits), s[CL].n(bits)] for s in solns]
# sys.exit()

# However, not all sets of variables are solvable. Here we have simultaneously
# UNDER-specified the lift equation and OVER-specified the bank equation.
# The solver detects this by evaluating the signatures of each equation
# to make sure the system is "full-rank" before trying to solve it.
try:
    x = solve([e1,e2], vrho=1.15, vv=10, vL=2, vbank=radians(45), vm=10)
except Exception as e:
    print("CORRECT BEHAVIOR:")
    print(e)

# However, these are NONLINEAR equations, so even "full-rank" systems do
# not necessarily have a solution with given input data.
# Here, the solvability check passes, but we get no solution, because we
# asked for the bank angle at which we generate LESS lift than in steady
# flight, which of course has no solution.
try:
    x = solve([e1,e2], vrho=1.15, vv=10, vL=2, vsref=20, vm=10)
except Exception as e:
    print("CORRECT BEHAVIOR:")
    print(e)



# Now, let's solve for the V/bank angle relationship at fixed CL,
# and make a plot!
guessDict = {'vv':0.1,'vL':300}
banks = np.arange(0,90,10)
vels = []

sys = System([e1,e2],[vrho,vCL,vbank,vsref,vm])
for bankVal in banks:
    print("Solving for bank",bankVal)
    sys.setConstVals(vrho=1.15, vCL=0.6, vbank=radians(bankVal), vsref=2, vm=10)
    x = sys.solve(guessDict)
    # x = solve([e1,e2], guessDict, vrho=1.15, vCL=0.6, vbank=radians(bankVal), vsref=2, vm=10)
    v = x['vv']
    vels.append(v)

fig = plt.figure()
ax = fig.add_subplot('111')
ax.plot(banks, vels, marker='o')
ax.set_xlabel("Bank Angle (deg)")
ax.set_ylabel("Velocity (m/s)")
# print("Bank ",bankVal, "v=",v)
plt.show()





