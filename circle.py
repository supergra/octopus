#!/usr/bin/env python3
import sys
from math import radians, degrees, sin, cos, atan, sqrt, pi

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from octopus import Variable, Equation, System, solve

# Constants
g = 9.81
rho = 1.15

def L(V,CL,Sref):
    return 0.5*rho*V*V*CL*Sref

def V_omega(omega,r):
    return omega*r

def omegaDef(V,r):
    return V/r

def omegaMatch(V1,r1,V2,r2):
    return V1*r2-V2*r1

def tension(M,th):
    return M*g/cos(radians(th))

# These are the residual equations:
def rw_geom(rm,l,theta):
    return -rm+l*sin(radians(theta))

def L_yforce(M,m,lam):
    return (M+m)*g/cos(radians(lam))


def tension_centrip_top(theta, M, omega, rm):
    # Centripetal balance at top
    return M*(omega**2)*rm/sin(radians(theta))

# Centripetal balance at bottom
def tension_centrip_bot(theta, M, m, L, lam, omega, rw):
    return (L*sin(radians(lam))+m*omega**2*rw)/sin(radians(theta))


# varNames = '''
# Vw
# Vm
# rw
# rm
# CL
# Sref
# M
# m
# theta
# roll
# l
# lift
# t
# omega
# '''

# for vName in varNames.split():
#     # Auto-add create local variables
#     # This is a bit awkward, since it has to be done locally
#     # and error-checked
#     if not vName:
#         continue

#     v = Variable(vName)

#     if vName in globals():
#         raise Exception("Variable name exists in scope: "+vName)
#     try:
#         exec("{} = v".format(vName))
#     except SyntaxError:
#         raise Exception("Invalid variable name: '{}'".format(vName))

# Or explicitly
Vw = Variable(unitlabel="m/s")
Vm = Variable(unitlabel="m/s")
rw = Variable(unitlabel="m")
rm = Variable(unitlabel="m")
CL = Variable()
Sref = Variable(unitlabel="m^2")
M = Variable(unitlabel="kg")
m = Variable(unitlabel="kg")
theta = Variable(unitlabel="deg",cyclicBounds=(-180,180))
roll = Variable(unitlabel="deg",cyclicBounds=(-180,180))
l = Variable(unitlabel="m")
lift = Variable(unitlabel="N")
t = Variable(unitlabel="N")
omega = Variable(unitlabel="rad/s")


# 8 equations
e1 = Equation(L, L=lift, V=Vw, Sref=Sref, CL=CL)

# Omega options:
# e2 = Equation(V_omega, omega=omega, V_omega=Vw, r=rw)
# e3 = Equation(V_omega, omega=omega, V_omega=Vm, r=rm)
# e2 = Equation(omega, omega=omega, V=Vw, r=rw)
# e3 = Equation(omega, omega=omega, V=Vm, r=rm)

e2 = Equation(omegaMatch, funcIsRes=True, V1=Vw,V2=Vm,r1=rw,r2=rm)
e3 = Equation(V_omega, V_omega=Vw, omega=omega, r=rw)

e4 = Equation(tension, tension=t, M=M, th=theta)
e5 = Equation(rw_geom, rw_geom=rw, rm=rm, l=l,theta=theta)
e6 = Equation(L_yforce, L_yforce=lift, M=M, m=m,lam=roll)
e7 = Equation(tension_centrip_top, tension_centrip_top=t, theta=theta,M=M,omega=omega, rm=rm)
e8 = Equation(tension_centrip_bot, tension_centrip_bot=t, theta=theta, M=M, m=m, L=lift, lam=roll, omega=omega, rw=rw)

equations = [e1,e2,e3,e4,e5,e6,e7,e8]


# Specify 6: the two from the other solver, plus M, m, Sref, CL!

# guessDict = {'theta':10,'roll':10,'Vm':3, 'rm':3, 'rw':20,'Vw':20, 'l':100,'omega':1,'lift':12000,'t':12000,'rw':35,'CL':1.1}
guessDict = {theta:10,roll:10,Vm:3, rm:3, rw:20,Vw:20, l:100,omega:1,lift:12000,t:12000,rw:35,CL:1.1}


def printSoln(x):
    maxlen = max(len(k) for k in x.keys())
    fmt = "{:"+str(maxlen)+"s} = {:.2f}"
    print("--- x Solution ---")

    for k,v in x.items():
        globals()[k].print(v)
        # print(fmt.format(k,v))
    print("------------------")

sys1 = System(equations, [m, M, Sref, CL, rw, l])
sys1.setConstVals(dict(m=200, M=1000, Sref=40,CL=1.1,rw=35,l=200))
# lbDict = {"roll":-90,"Vw":0,"rw":0,"l":0,"CL":0,"Sref":0,"m":0,"M":0} #,"theta":-90}
# ubDict = {"roll":90} #,"theta":90}

lbDict = {"Vw":0,"rw":0,"l":0,"CL":0,"Sref":0,"m":0,"M":0} #,"theta":-90}
ubDict = {} #,"theta":90}

# x = sys1.solve(minBoundDict=lbDict, maxBoundDict=ubDict)
# printSoln(x)
# sys.exit()

x = solve(equations,guessDict, dict(m=200, M=1000, Sref=40,CL=1.1,rw=35,l=200))
printSoln(x)

x = solve(equations,guessDict, dict(m=200, M=1000, Sref=40,Vw=21.6,rw=35,l=200))
printSoln(x)

# This combination is NOT solvable, and not obviously not so!
# We get a conflict for lift via:
#   equation 1 -- Lift(Vw,CL,Sref)
#   equation 6 -- Lift(M,m,lambda)
#
#   By specifying both lambda and CL/Vw, we have a conflict!
try:
    x = solve(equations,guessDict,dict(m=200, M=1000, Sref=40,Vw=21.6,roll=10,CL=1.1))
    printSoln(x)
except Exception as e:
    print("CORRECT FAILURE:\n ",e)

lineLengths = [20*i for i in range(2,10)]
rws = [10*j for j in range(1,8)]

designMass = 375
wingLoading = 5 # (cruise wingloading -- is higher when banking)

# x = solve(equations,guessDict, dict(m=200, M=1000, Sref=40,CL=1.1,rw=35,l=200))

nullVal = 100
sys1 = System(equations, [m, M, Sref, CL, rw, l])

# sys1.analyze()
sys.exit()

guessDict = {theta:10,roll:20,Vm:5, rm:5, Vw:25, omega:1,lift:designMass*10,t:designMass*10}

outputs = np.zeros((len(lineLengths),len(rws)))
for i,this_L in enumerate(lineLengths):
    for j,this_rW in enumerate(rws):

        sys1.setConstVals(dict(m=125, M=250, Sref=designMass/wingLoading,CL=0.9,rw=this_rW,l=this_L))
        # sys1.setConstVals(dict(m=200, M=1000, Sref=40,CL=1.1,rw=this_rW,l=this_L))
        try:
            x = sys1.solve(guessDict,minBoundDict=lbDict, maxBoundDict=ubDict)
        except Exception as e:
            print(this_L, this_rW,"failed: "+str(e))
            output = nullVal
        else:
            guessDict = x # Next guess is previous solution
            # output = x["rm"]
            output = x["theta"]
            # printSoln(x)
            print(this_L, this_rW)
            theta.print(output)

        outputs[i,j] = output



label = "Mass radius"

outputs = outputs.transpose() # Get line lengths on x-axis
extent = (min(lineLengths), max(lineLengths), min(rws), max(rws))

# from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import


cm = mpl.cm.get_cmap('gray')

invertColors = True
from lib_Plot import reverse_colourmap
if invertColors:
    cm = reverse_colourmap(cm)

def getContours(data, Ncontours, maxLevelVal=1000):
    minout = np.min(data)
    maxout = np.max(data)
    print(minout, maxout)
    # Create auto contour levels
    if minout > 0:
        # minout = 0
        start_levels = [0]
    else:
        start_levels = []
    dx = (maxout*0.99-minout)/(Ncontours-1)
    assert dx > 0, dx
    levels = start_levels + [minout+dx*_ for _ in range(Ncontours)]
    levels = [_ for _ in levels if _ < maxLevelVal]
    levels.sort()
    return levels

levels = getContours(outputs, Ncontours=20)
fig, ax = plt.subplots()
plt.title(label)
plt.xlabel('Length (m)')
plt.ylabel('Wing radius (m)')

im = ax.imshow(outputs, interpolation='bicubic', origin='lower',
            cmap=cm, extent=extent, aspect='auto')

CS = ax.contour(outputs, levels, origin='lower', linewidths=1, extent=extent)

plt.clabel(CS, levels[::2],  # label every second level
       inline=1, fmt='%1.1f', fontsize=12)
plt.show()

