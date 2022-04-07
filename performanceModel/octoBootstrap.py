from __future__ import print_function, division
import sys
assert sys.version_info >= (3,0)
from math import radians, degrees, sin, cos, atan, sqrt, pi, asin

import os

from pathlib import Path

octopus_dir = os.path.abspath(os.path.join(Path(__file__).parent.resolve(),".."))
# print(octopus_dir)
sys.path.append(octopus_dir)
from octopus import Variable, Equation, System, getAssignedName
import numpy as np

# Constants
g = 9.81
RHO0 = 1.225

def checkValidity(aircraftParams, performanceParams):
    Sref = aircraftParams["Sref"]
    d = aircraftParams["propDiam"]

    AR = aircraftParams["AR"]
    P0 = aircraftParams["P0"]
    N0 = aircraftParams["N0"]
    M0 = P0/(2*pi*N0/60)
    e = performanceParams["e"]
    CD0 = performanceParams["CD0"]
    m = performanceParams["m"]
    b = performanceParams["b"]

    arg = 2*pi*M0*m/(mass*g*d) - 2*sqrt((Sref*CD0-2*d*d*b)/(Sref*pi*e*AR))

    if not -1 <= arg <= 1:
        raise Exception("Invalid inputs! "+str(arg))




# Define functions in the natural way, with a single output
def LiftDef(rho, V, CL, Sref):
    return 0.5*rho*V*V*CL*Sref

def Lbank(bankAngleRad, mass):
    assert -pi/2 < bankAngleRad < pi/2, bankAngleRad
    return g*mass/cos(bankAngleRad)

def CD_def(CD0, e, CL, AR):
    assert e > 0, e
    assert AR > 0, AR
    return CD0 + CL*CL/(pi*e*AR)

def DragDef(rho, V, CD, Sref):
    return 0.5*rho*V*V*CD*Sref

def hdotEq(V, gamma):
    assert -pi/2 < gamma < pi/2, gamma
    return V*sin(gamma)

def bs_PHI(rho, Cdropoff, hack):
    # For C=0.12, this is only good up to about 15000 ft, rho = ~ 0.15
    assert rho > 0
    # if hack:
    #     # return 0.955 # must be lower when density is higher.
    #     return 0.968 # 0.969+ fails

    if Cdropoff < 0:
        return 1.0

    sigma = rho/RHO0
    # This can go negative at high altitude....
    # if sigma <= Cdropoff:
    #     return 0
    return (sigma-Cdropoff)/(1-Cdropoff)

def bs_E(rho, Cdropoff, M0, m, propDiam):
    # eq. 7.53
    phi = bs_PHI(rho, Cdropoff, hack=False)
    return phi*2*pi*M0*m/propDiam

def bs_F(rho, propDiam, b):
    # eq. 7.54
    assert b < 0, b
    assert propDiam > 0
    assert rho > 0
    return rho*propDiam**2*b

def bs_G(rho, Sref, CD0):
    # eq. 7.55
    return 0.5*rho*Sref*CD0

def bs_H(mass, rho, Sref, e, AR):
    # eq. 7.56
    return 2*(mass*g)**2/(rho*Sref*pi*e*AR)

# other params are combos of E-H:
#  K = F-G  (< 0)
#  Q = E/K  (< 0)
#  R = H/K  (< 0)
#  U = H/G  (> 0)



def bs_all(rho, Cdropoff, M0, m, propDiam, b, Sref, CD0, mass, e, AR):
    # Compute all the bootstrap helper parameters at once.
    E = bs_E(rho, Cdropoff, M0, m, propDiam)
    F = bs_F(rho, propDiam, b)
    G = bs_G(rho, Sref, CD0)
    H = bs_H(mass, rho, Sref, e, AR)
    K = F-G
    Q = E/K
    R = H/K
    U = H/G
    assert F < 0
    assert K < 0
    assert Q < 0, (Q, E, K, F, G)
    assert R < 0
    assert U > 0
    # Check arguments:
    if Q*Q/4+R < 0 or Q*Q/36-R/3 < 0:
        print("Invalid Q,R for sqrt args:")
        print("E=",E)
        print("F=",F)
        print("G=",G)
        print("H=",H)
        raise Exception("asdf")

    return E,F,G,H,K,Q,R,U

def bs_no_EQ(rho, propDiam, b, Sref, CD0, mass, e, AR):
    # Compute all the bootstrap helper parameters at once.
    F = bs_F(rho, propDiam, b)
    G = bs_G(rho, Sref, CD0)
    H = bs_H(mass, rho, Sref, e, AR)
    K = F-G
    R = H/K
    U = H/G
    assert F < 0
    assert K < 0
    assert R < 0
    assert U > 0
    return F,G,H,K,R,U

def bs_VM(rho, Cdropoff, M0, m, propDiam, b, Sref, CD0, mass, e, AR):
    # eq. 7.19
    E,F,G,H,K,Q,R,U = bs_all(rho, Cdropoff, M0, m, propDiam, b, Sref, CD0, mass, e, AR)
    # print(Q,R)

    arg = Q*Q/4+R
    if arg < 0:
        raise Exception("Invalid Q,R: "+str(Q)+" "+str(R))

    arg = -Q/2 + sqrt(arg)
    if arg < 0:
        raise Exception("Invalid Q,R: "+str(Q)+" "+str(R))
    return sqrt(arg)

def bs_Vm(rho, Cdropoff, M0, m, propDiam, b, Sref, CD0, mass, e, AR):
    # eq. 7.21 [different sign on sqrt from 7.19]
    E,F,G,H,K,Q,R,U = bs_all(rho, Cdropoff, M0, m, propDiam, b, Sref, CD0, mass, e, AR)
    return sqrt(-Q/2-sqrt(Q*Q/4+R))

def bs_Vy(rho, Cdropoff, M0, m, propDiam, b, Sref, CD0, mass, e, AR):
    # eq. 7.21 [different sign on sqrt from 7.19]
    E,F,G,H,K,Q,R,U = bs_all(rho, Cdropoff, M0, m, propDiam, b, Sref, CD0, mass, e, AR)
    return sqrt(-Q/6+sqrt(Q*Q/36-R/3))

def bs_Vx(rho, propDiam, b, Sref, CD0, mass, e, AR):
    # Lowry eq. 7.27
    F,G,H,K,R,U = bs_no_EQ(rho, propDiam, b, Sref, CD0, mass, e, AR)
    return (-R)**0.25


def bs_Vbg(rho, propDiam, b, Sref, CD0, mass, e, AR):
    # Lowry eq. 7.27
    F,G,H,K,R,U = bs_no_EQ(rho, propDiam, b, Sref, CD0, mass, e, AR)
    return (U)**0.25

def bs_Vmd(rho, propDiam, b, Sref, CD0, mass, e, AR):
    # Lowry eq. 7.27
    F,G,H,K,R,U = bs_no_EQ(rho, propDiam, b, Sref, CD0, mass, e, AR)
    return (U/3)**0.25

def bs_gammaX(rho, Cdropoff, M0, m, propDiam, b, Sref, CD0, mass, e, AR):
    # eq. 7.44
    # Note: This is the maximum climb angle for full throttle.
    # For a heavily oversized motor, it could be 90 degrees, and
    # even exceed it (i.e. be accelerating vertically).
    E,F,G,H,K,Q,R,U = bs_all(rho, Cdropoff, M0, m, propDiam, b, Sref, CD0, mass, e, AR)

    if (-K*H) < 0:
        raise Exception("Invalid K and H: "+str(K)+" "+str(H))

    arg = (E-2*sqrt(-K*H))/(mass*g)
    if abs(arg) > 1:
        # print(E,2*sqrt(-K*H))
        if arg > 1: # Super oversized motor, can climb vertically
            return pi/2
        # else:
        #     return -pi/2
        raise Exception("Argument to arcsin outside -1,1: "+str(arg)+" K={}, H={}".format(K,H))
    return asin(arg)

def bs_hdot_md(rho, propDiam, b, Sref, CD0, mass, e, AR):
    # eq. 7.46, using Vmd
    F,G,H,K,R,U = bs_no_EQ(rho, propDiam, b, Sref, CD0, mass, e, AR)

    Vmd = bs_Vmd(rho, propDiam, b, Sref, CD0, mass, e, AR)

    return (-G*Vmd**3-H/Vmd)/(mass*g)



def reqTorque(V,gamma,mass,rho,propDiam,m,b,Sref,CD0,e,AR,bank):
    # Typically, think of this as (V, bank) -> torque required
    assert -pi/2 < gamma < pi/2, gamma
    assert V > 0, V
    F = bs_F(rho,propDiam, b)
    G = bs_G(rho, Sref, CD0)
    H = bs_H(mass,rho,Sref,e,AR)
    K = F-G
    R = H/K
    # print(F,G,H,K,R)
    assert F < 0 and K < 0 and R < 0, (F,K,R) # (because b < 0)
    # Vx = bs_Vx_R(R)
    # print("VX: ",Vx)
    sigma = rho/RHO0

    # Now use Lowry eq. 10.6 to get required Torque
    # (note: Kb is K derived at full throttle)
    # Treq = (mass*g*sin(gamma) + sigma*K*(R-V**4)/V**2)*propDiam/(2*pi*m)

    # Now using eq. 10.35 to account for increased induced drag due to
    # increased lift required in bank.
    Treq = (mass*g*sin(gamma) - K*V**2 + H/(V**2 * (cos(bank))**2))*propDiam/(2*pi*m)

    return Treq

def CpJsqTorque(torque, rho, propDiam, V):
    ''' returns CP/J^2 (power coeff over advance ratio squared)'''
    return 2*pi*torque/(rho*propDiam**3*V**2)

def advanceRatio(V, omega, propDiam):
    # omega in rad/s
    return V/(propDiam*omega/(2*pi))

def powerApplied(torque, omega):
    return torque*omega

def torqueFactor2Torque(torqueFactor, M0, rho, Cdropoff):
    ''' Returns torque at this altitude and torque factor '''
    phi = bs_PHI(rho, Cdropoff, hack=True)

    torque = torqueFactor*phi*M0
    # print(phi,M0,torqueFactor,torque)
    return torque

def CP_def(power, rho, propOmega, propDiam):
    # Definition of power coefficient
    n = propOmega/(2*pi)
    assert n > 0, n
    return power/(rho*n**3*propDiam**5)

def CT_def(thrust, rho, propOmega, propDiam):
    # Definition of thrust coefficient
    n = propOmega/(2*pi)
    assert n > 0, n
    return thrust/(rho*n**2*propDiam**4)

def propEfficiency(J, CT, CP):
    if CP == 0:
        return 0
    # assert CP != 0, CP
    return J*CT/CP

def CT_linearPolar(CP, J, m, b):
    # Lowry eq. 6.53
    return m*CP + b*J**2


def setupSolver(aircraftParams, performanceParams, specifiedVarNames):

    # try:
    #     checkValidity(aircraftParams, performanceParams)
    # except:
    #     remove_gammaX = True
    # else:
    #     remove_gammaX = False

    isElectric = (aircraftParams["Cdropoff"] < 0)

    throttle2torquefactor_array = performanceParams["throttle2torquefactor"]
    CPJsq2J_array = performanceParams["CPJsq2J"]


    # define functions that require closure (e.g. interpolation-style ones)
    # and can't be handled as residual-style equations
    # [Not sure how to describe this.]

    _cps = tuple([_[0] for _ in CPJsq2J_array])
    _js = tuple([_[1] for _ in CPJsq2J_array])
    def CpJsq2J(CpJsq):
        # Lowry p. 326
        if CpJsq < _cps[0]:
            raise Exception("CpJsq too low "+str(CpJsq))
        if CpJsq > _cps[-1]:
            raise Exception("CpJsq too high "+str(CpJsq))
        J = float(np.interp(CpJsq,_cps,_js))
        return J

    _ts = tuple([_[0] for _ in throttle2torquefactor_array])
    _tfs = tuple([_[1] for _ in throttle2torquefactor_array])
    def throttle2TorqueFactor(throttle):
        assert 0 <= throttle <= 1, throttle
        torqueFactor = float(np.interp(throttle,_ts,_tfs))
        return torqueFactor

    if not isElectric:
        power2bsfc_array = performanceParams["power2bsfc"]
        _ps = tuple([_[0] for _ in power2bsfc_array])
        _bsfcs = tuple([_[1] for _ in power2bsfc_array])
        def power2Bsfc(power):
            assert 0 <= power, power
            bsfc = float(np.interp(power,_ps,_bsfcs))
            return bsfc


    # Wrappers to register and track all POSSIBLE variables and equations
    # Which ones are used will depend on what we're trying to solve
    variablesAvailable = []
    equationsAvailable = []
    def addvar(**kwargs):
        name = getAssignedName()
        v = Variable(name=name, **kwargs)
        assert v.name == name
        variablesAvailable.append(v)
        # print("Registered variable",v.name)
        return v

    def addeqn(*args, **kwargs):
        e = Equation(*args, **kwargs)
        equationsAvailable.append(e)
        # print("Registered",e.name)
        return e

    # Make variables on one-line each, as follows.
    # Internal names automatically become the variable names given here
    CL = addvar()

    if isElectric:
        rho = addvar(unitlabel="kg/m^3",minVal=0,maxVal=1.5)
    else:
        # Have to limit gas engines to where the dropoff factor is valid
        rho = addvar(unitlabel="kg/m^3",minVal=0.15,maxVal=1.5)

    Sref = addvar(unitlabel="m^2",minVal=0)
    Lift = addvar(unitlabel="N")
    drag = addvar(unitlabel="N",minVal=0)
    V = addvar(unitlabel="m/s",minVal=0)
    mass = addvar(unitlabel="kg",minVal=0)
    torque = addvar(unitlabel="N-m",minVal=0)
    bank = addvar(unitlabel="rad",cyclicBounds=(-pi,pi),minVal=-pi/2,maxVal=pi/2)

    CD = addvar(minVal=0)
    CD0 = addvar(minVal=0,maxVal=1.0)
    e = addvar(minVal=0, maxVal=2.0)
    AR = addvar(minVal=0, maxVal=50)

    m = addvar(minVal=0, maxVal=1.5)
    b = addvar(maxVal=0,minVal=-0.5)

    propDiam = addvar(unitlabel="m",minVal=0)
    gamma = addvar(unitlabel="rad",cyclicBounds=(-pi,pi),minVal=-pi/2,maxVal=pi/2)

    hdot = addvar(unitlabel="m/s")

    CpJsq = addvar(minVal=0)
    J = addvar() # advance ratio
    propOmega = addvar(unitlabel="rad/s",minVal=0)

    Power = addvar(unitlabel="W",minVal=0)

    throttle = addvar(minVal=0,maxVal=1)

    TorqueFactor = addvar(minVal=0,maxVal=1)
    M0 = addvar(unitlabel="N-m",minVal=0,maxVal=1e5)

    if isElectric:
        Cdropoff = addvar(minVal=-1,maxVal=-1) # Force -1 to signal electric
    else:
        Cdropoff = addvar(minVal=0,maxVal=0.3) # usually 0.12
        bsfc = addvar(unitlabel="g/(kW-hr)",minVal=100,maxVal=600)

    # V-speeds and best rates/angles

    Vm = addvar(unitlabel="m/s",minVal=0)
    Vx = addvar(unitlabel="m/s",minVal=0)
    Vbg = addvar(unitlabel="m/s",minVal=0)
    Vmd = addvar(unitlabel="m/s",minVal=0)
    Vy = addvar(unitlabel="m/s",minVal=0)
    VM = addvar(unitlabel="m/s",minVal=0)


    gammaX = addvar(unitlabel="rad",minVal=-90,maxVal=90)
    hdot_md = addvar(unitlabel="m/s", maxVal=0)

    thrust = addvar(unitlabel="N")

    CP = addvar()
    CT = addvar()
    propEff = addvar(maxVal=1) # Can go negative, and must let it!

    # Create equations, which are functions with a map between
    # variable names and the function argument names.
    addeqn(LiftDef, CL=CL, rho=rho, Sref=Sref,V=V,LiftDef=Lift)
    addeqn(reqTorque,V=V,gamma=gamma,mass=mass,rho=rho,propDiam=propDiam,m=m,b=b,Sref=Sref,CD0=CD0,e=e,AR=AR,bank=bank,reqTorque=torque)

    addeqn(CpJsqTorque, torque=torque, rho=rho, propDiam=propDiam, V=V,CpJsqTorque=CpJsq)

    addeqn(advanceRatio, V=V, omega=propOmega, propDiam=propDiam,advanceRatio=J)
    addeqn(powerApplied, torque=torque, omega=propOmega, powerApplied=Power)

    # Data interpolators:
    addeqn(CpJsq2J, CpJsq=CpJsq, CpJsq2J=J)
    addeqn(throttle2TorqueFactor, throttle=throttle, throttle2TorqueFactor=TorqueFactor)

    addeqn(torqueFactor2Torque, torqueFactor=TorqueFactor, M0=M0, rho=rho, Cdropoff=Cdropoff, torqueFactor2Torque=torque)

    addeqn(Lbank, bankAngleRad=bank, mass=mass, Lbank=Lift)

    addeqn(bs_VM, rho=rho, Cdropoff=Cdropoff, M0=M0, m=m, propDiam=propDiam, b=b, Sref=Sref, CD0=CD0, mass=mass, e=e, AR=AR, bs_VM=VM)
    addeqn(bs_Vm, rho=rho, Cdropoff=Cdropoff, M0=M0, m=m, propDiam=propDiam, b=b, Sref=Sref, CD0=CD0, mass=mass, e=e, AR=AR, bs_Vm=Vm)
    addeqn(bs_Vy, rho=rho, Cdropoff=Cdropoff, M0=M0, m=m, propDiam=propDiam, b=b, Sref=Sref, CD0=CD0, mass=mass, e=e, AR=AR, bs_Vy=Vy)
    addeqn(bs_Vx, rho=rho, propDiam=propDiam, b=b, Sref=Sref, CD0=CD0, mass=mass, e=e, AR=AR, bs_Vx=Vx)
    addeqn(bs_Vbg, rho=rho, propDiam=propDiam, b=b, Sref=Sref, CD0=CD0, mass=mass, e=e, AR=AR, bs_Vbg=Vbg)
    addeqn(bs_Vmd, rho=rho, propDiam=propDiam, b=b, Sref=Sref, CD0=CD0, mass=mass, e=e, AR=AR, bs_Vmd=Vmd)

    addeqn(bs_gammaX, rho=rho, Cdropoff=Cdropoff, M0=M0, m=m, propDiam=propDiam, b=b, Sref=Sref, CD0=CD0, mass=mass, e=e, AR=AR, bs_gammaX=gammaX)
    addeqn(bs_hdot_md, rho=rho, propDiam=propDiam, b=b, Sref=Sref, CD0=CD0, mass=mass, e=e, AR=AR, bs_hdot_md=hdot_md)

    addeqn(hdotEq, V=V, gamma=gamma, hdotEq=hdot)


    addeqn(CP_def, power=Power, rho=rho, propOmega=propOmega, propDiam=propDiam, CP_def=CP)
    addeqn(CT_def, thrust=thrust, rho=rho, propOmega=propOmega, propDiam=propDiam, CT_def=CT)
    addeqn(propEfficiency, J=J, CT=CT, CP=CP, propEfficiency=propEff)
    addeqn(CT_linearPolar, CP=CP, J=J, m=m, b=b, CT_linearPolar=CT)
    addeqn(CD_def, CD0=CD0, e=e, CL=CL, AR=AR, CD_def=CD)
    addeqn(DragDef, rho=rho,V=V,CD=CD,Sref=Sref,DragDef=drag)

    if not isElectric:
        addeqn(power2Bsfc, power=Power, power2Bsfc=bsfc)


    # -- DONE REGISTERING ALL VARIABLES AND EQUATIONS! --

    varName2Var = {v.name: v for v in variablesAvailable}


    # Trying to auto-detect which variables are required to get
    # a set of outputs... WIP
    # outputsRequested = ["J","propEff"]
    # variableNamesInvolved = specifiedVarNames + outputsRequested

    # for param in list(aircraftParams.keys()) + list(performanceParams.keys()):
    #     for v in variablesAvailable:
    #         if param == v.name:
    #             variableNamesInvolved.append(param)

    # For now just use everything
    variableNamesInvolved = list(varName2Var.keys())
    # if remove_gammaX:
    #     variableNamesInvolved.remove("gammaX")
    # print(variableNamesInvolved)

    varNames = list(varName2Var.keys())

    # Determine which equations are actually necessary to solve for all
    # the variables.
    equations = []
    for eq in equationsAvailable:

        if all(_ in variableNamesInvolved for _ in eq.allVarNames):
            # print("Using",eq.name)
            equations.append(eq)
        # else:
            # print("Don't need",eq.name)
            # print(tuple(n for n in eq.allVarNames if n not in variableNamesInvolved))
            # print(set(eq.allVarNames).remove(set(variableNamesInvolved)))


    fixedVals = {
        Sref:aircraftParams["Sref"],
        AR:aircraftParams["AR"],
        Cdropoff:aircraftParams["Cdropoff"] if aircraftParams["Cdropoff"] >= 0 else -1.0,
        propDiam:aircraftParams["propDiam"],
        CD0:performanceParams["CD0"],
        e:performanceParams["e"],
        m:performanceParams["m"],
        b:performanceParams["b"],
        M0:aircraftParams["P0"]/(2*pi*aircraftParams["N0"]/60)
    }



    fixedVars = list(fixedVals.keys())

    specifiedVars = fixedVars
    variablesInvolved = [varName2Var[name] for name in variableNamesInvolved]
    for var in variablesInvolved:
        if var.name in specifiedVarNames:
            specifiedVars.append(var)

    system = System(equations, specifiedVars)

    def printSoln(x):
        print("--- Solution ---")
        for k,v in x.items():
            for variable in variablesInvolved:
                if variable.name == k:
                    variable.print(v)
                    break
            else:
                print("Couldn't find",k)
        print("------------------")

    # Hack-town: register solution printer...
    system.printSoln = printSoln

    # Convert to name:val format
    fixedVals = {var.name:val for var,val in fixedVals.items()}


    defaultGuessDict = {

        rho: 1.225,
        V:   20,
        mass:    15,
        gamma:   radians(0),
        bank:    radians(0),

        # possibly variables:
        J:    0.5,
        propOmega:   200,
        Power:   5000,
        throttle:    0.5,

        # Not intended to be used as inputs:
        Lift:   150,
        CpJsq:   0.1,
        TorqueFactor:    1.0,
        torque:  5,

        CD:  0.05,
        CP:  0.1,
        CT:  0.1,
        CL:  0.5,

        Vx:15,
        Vy:20,
        VM:35,
        Vm:10,
        Vbg:25,
        Vmd:15,
        hdot:4,
        propEff:0.5,
        thrust:100,

        hdot_md: -2,
        gammaX: radians(20)
    }
    if not isElectric:
        defaultGuessDict[bsfc] = 300


    newGuessDict = {}
    for var,guess in defaultGuessDict.items():
        if type(var) is Variable:
            newGuessDict[var.name] = guess
        else:
            newGuessDict[var] = guess
    defaultGuessDict = newGuessDict



    for var in system.unknownVarsOrder:
        assert type(var) is str, type(var)
        if var not in defaultGuessDict:
            print("Default guess for ",var)
            defaultGuessDict[var] = 0.1

    defaultLbDict = {v.name:v.minVal for v in variablesInvolved}
    defaultUbDict = {v.name:v.maxVal for v in variablesInvolved}

    # print(defaultUbDict)
    return system, fixedVals, defaultGuessDict, defaultLbDict, defaultUbDict




