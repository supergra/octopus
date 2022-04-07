#!/usr/bin/env python

import os
from math import radians, degrees, pi
import sys

from octoBootstrap import setupSolver
from cleanJson import ReadJson

performanceParams = ReadJson("perfModel.json")
aircraftParams = ReadJson("aircraft.json")

inputs = ReadJson("inputs.json")
specifiedVarNames = list(inputs.keys())

system, _fixedVals, guessDict, lbDict, ubDict = setupSolver(aircraftParams, performanceParams, specifiedVarNames)

# Add in the variables we are specifying
fixedVals = _fixedVals.copy()
for k,v in inputs.items():
    if k in ('gamma','bank','gammaX','propOmega'): # DUMMY DUMMY HACK
        v = radians(v)
    fixedVals[k] = v

# modify the guess dict, lb dict, ub dict, if desired
guessDict["L"] = inputs["mass"]*9.81

system.setConstVals(fixedVals)

def printSoln(x):
    print("--- Specified ---")
    for k,v in fixedVals.items():
        print(k,"=",v)
    system.printSoln(x)

x = system.solve(guessDict=guessDict,minBoundDict=lbDict,maxBoundDict=ubDict)

printSoln(x)

Vlo = 10
Vhi = 30
# Vlo = x["Vm"] - 1
# Vhi = x["VM"] + 1
N = 15

Vs = []
solns = []

for i in range(N+1):
    V = Vlo + (Vhi-Vlo)*float(i)/N
    fixedVals["V"] = V


    system.setConstVals(fixedVals)
    try:
        x = system.solve(guessDict=guessDict,minBoundDict=lbDict,maxBoundDict=ubDict)
    except Exception as e:
        print("Failed: ",e)
    else:
        Vs.append(V)
        solns.append(x)
        print("V: {:.1f} m/s --> Thrust {:.1f}N".format(V,x["thrust"]))
        # system.printSoln(x)


import matplotlib.pyplot as plt
from lib_Plot import NewSetOfPlots, AddLegend
axes = NewSetOfPlots(n=4)
hdots = [x["hdot"] for x in solns]
drags = [x["drag"] for x in solns]
thrusts = [x["thrust"] for x in solns]
power_reqs = [x["drag"]*V for (x,V) in zip(solns,Vs)]
power_app = [x["Power"] for x in solns]
power_avail = [x["thrust"]*V for x,V in zip(solns,Vs)]


axes[-1].set_xlabel("Airspeed (m/s)")

axes[0].set_ylabel("Thrust/Drag (N)")
axes[0].plot(Vs, thrusts, label='Thrust available')
axes[0].plot(Vs, drags, label='Drag')
axes[0].plot(Vs, [x["torque"] for x in solns], label='Shaft torque applied')

axes[1].set_ylabel("Power (W)")
axes[1].plot(Vs, power_avail, label='Power available')
axes[1].plot(Vs, power_reqs, label='Power required')
axes[1].plot(Vs, power_app, label='Shaft Power applied')
axes[1].plot(Vs, [x["propOmega"]*60/(2*pi) for x in solns], label='RPM')

axes[2].set_ylabel("Climb rate (m/s)/ angle (deg)")
axes[2].plot(Vs, [x["hdot"] for x in solns], label='hdot')
axes[2].plot(Vs, [degrees(x["gamma"]) for x in solns],label='gamma')

axes[3].set_ylabel("Efficiency")
axes[3].plot(Vs, [x["J"] for x in solns], label='Advance ratio')
axes[3].plot(Vs, [x["propEff"] for x in solns], label='prop eff')


for ax in axes:
    AddLegend(ax)

plt.show()


