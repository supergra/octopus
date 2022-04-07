import random
import inspect # getargspec for converting function signatures
import string # For random name maker
import traceback # For magic variable maker

import numpy as np
import scipy.optimize

import math

# NOLIMIT = 1e10
NOLIMIT = float('inf')


def randomname(length):
   letters = string.ascii_lowercase
   return ''.join(random.choice(letters) for i in range(length))


def packDict(keyOrder, dictToPack, defaultval):
    ''' Order values from dictionary into a list '''
    if dictToPack is None:
        dictToPack = {}

    out = []
    for v in keyOrder:
        if v in dictToPack:
            out.append(dictToPack[v])
        else:
            # print("Using default: {}={}".format(v,defaultval))
            out.append(defaultval)
    return out

def packDictStrict(keyOrder, dictToPack):
    return [dictToPack[v] for v in keyOrder]

def unpackToDict(varOrder, packedList):
    ''' Pack an ordered list into a dictionary '''
    return {var:val for var,val in zip(varOrder, packedList)}

def getAssignedName(depth=2):
    # Should be called by an __init__ block of a class.
    # Returns the local name assigned to the instance in the code calling
    # the instantiator!
    stack = traceback.extract_stack()
    filename, lineno, function_name, code = stack[-depth-1]
    name = code.split('=')[0].strip()
    return name

class Variable:
    # This is mostly a way to help map variable names to function
    # argument names.

    def __init__(self, name="", unitlabel="",minVal=-NOLIMIT, maxVal=NOLIMIT, cyclicBounds=None):

        if not name:
            try:
                name = getAssignedName()
            except:
                name = randomname(length=3)
                print("Failed to make good name -- giving random name: "+name)

        # print("Made variable",name)

        self.name = name
        self.cyclicBounds = cyclicBounds
        self.minVal = minVal
        self.maxVal = maxVal
        self.unitlabel = unitlabel

        if self.cyclicBounds:
            self.minVal = max(self.minVal, self.cyclicBounds[0])
            self.maxVal = min(self.maxVal, self.cyclicBounds[1])


    def print(self, value):
        fmt = "{:10s} : {:8.4g} {}"
        if self.cyclicBounds:

            if not self.cyclicBounds[0] <= value <= self.cyclicBounds[1]:
                width = self.cyclicBounds[1] - self.cyclicBounds[0]
                realvalue = ((value-self.cyclicBounds[0]) % width)+self.cyclicBounds[0]
                print(fmt+" [{}]".format(self.name,realvalue,self.unitlabel,value))
                return

        print(fmt.format(self.name,value,self.unitlabel))



class Equation:
    # An equation is a function with a mapping of variables to parameters.

    # Example: take a function defined as:
    #   def L(rho, V, CL, Sref):
    #       return 0.5*rho*V*V*CL*Sref

    # we could use this to compute lift relations at stall.
    # (1) Create 5 variables (could be constants):
    #   rho1000 = Variable()
    #   Vstall = Variable()
    #   CLmax = Variable()
    #   Sref = Variable()
    #   Lmax = Variable() # Note we need to create the output as a var too!

    # (2) Map these to the lift variables when instantiating Equation():
    #   e = Equation(L, rho=rho1000, V=Vstall, CL=CLmax, Sref=Sref, L=Lmax)



    # Some instance variables:
    #   allParamNames: All kwargs in func + output param name
    #   allVarNames:   Variable names that have been mapped to the params
    #   resFunc
    #   param2var, var2param [maps variables to/from equation parameters]

    def __init__(self, func, funcIsRes=False, **kwargs):
        '''
        funcIsRes: function is already in residual form [returns 0 when satisified]
        '''
        try:
            func.__name__
        except:
            raise Exception("func passed in is not a function")

        argSpec = inspect.getargspec(func)
        inputParamNames = argSpec.args
        self.name = "eqn_"+func.__name__
        self.allParamNames = inputParamNames[:] # copy!

        if not funcIsRes:
            outputParamName = func.__name__
            # self.allParamNames += [outputParamName]
            self.allParamNames.append(outputParamName)


        if set(kwargs.keys()) != set(self.allParamNames):
            raise Exception("kwargs = "+str(kwargs.keys())+ " but should be "+str(self.allParamNames))

        if any(type(_) is not Variable for _ in kwargs.values()):
            raise Exception("Input types must all be var")

        # Store the map of variables to function inputs
        self.var2param = {}
        self.param2var = {}
        self.allVarNames = []

        minValidValues = []
        maxValidValues = []
        for k,v in kwargs.items():
            self.param2var[k] = v.name
            self.var2param[v.name] = k
            self.allVarNames.append(v.name)

            minValidValues.append(-10 if v.minVal == -NOLIMIT else v.minVal)
            maxValidValues.append( 10 if v.maxVal == NOLIMIT else v.maxVal)

        checkValidFunction(func)

        if funcIsRes:
            # Function is already in residual form
            self.resFunc = func
        else:
            # Create a new function in residual form
            self.resFunc = makeResidualFunc(func)

            # Test that everything works properly

            # Make randomly generated inputs and try calling function.
            # If variable validity bounds are set, respect those to avoid
            # predictable function failures
            testKwargs = {}
            for k in inputParamNames:
                pname = self.param2var[k]
                i = self.allVarNames.index(pname)
                lo = minValidValues[i]
                hi = maxValidValues[i]

                testKwargs[k] = lo + random.random()*(hi-lo)

            try:
                out = func(**testKwargs)
            except:
                print("Warning: Residual function failed with random inputs"+str(testKwargs))
                return
                # raise Exception("Residual function failed with random inputs"+str(testKwargs))

            testKwargs[outputParamName] = out
            try:
                res = self.resFunc(**testKwargs)
            except:
                raise Exception("Residual function failed with random inputs: res = ",res)
            else:
                if res != 0:
                    raise Exception("Weird failure. Function failed to be converted to residual form")
                # else:
                #     print("Residual is happy with random inputs!")
                #



    def makeCompressedResFunc(self, globalSetVarNameList):
        # Tell the equation which variables will be set, and get
        # back a function with the correctly reduced signature

        # globalSetVarValueDict contains ALL set values, even those
        # that don't apply to this equation. This function processes
        # those that it understands.

        # Closure using var-to-param mapping
        # Returns a new residual function with constant values baked in

        # (1) First, extract only the variables we recognize
        setVarNames = []
        setParamNames = []
        for name in globalSetVarNameList:
            if name in self.allVarNames:
                setVarNames.append(name)
                setParamNames.append(self.var2param[name])

        self.setParamValues = {} # Will be populated later

        # (2) Determine what parameters will have to be passed in every call,
        #     thus determining signature of new compressed residual func.
        unsetParams = list(set(self.allParamNames) - set(setParamNames))
        unsetVarNames = [self.param2var[p] for p in unsetParams]

        self.setVarNames = setVarNames

        firstCall = True
        def compressedResFunc(**inputs):
            # Inputs are the unset variable values as kwargs
            # References: unsetVarNames to error check,
            #             setParamValues dict to close on specified values.
            # Invokes the original residual function,
            # which invokes the original function


            # Error-checking only on first time through
            nonlocal firstCall # Python3-only syntax!
            if firstCall:
                if len(self.setParamValues) != len(self.setVarNames):
                    raise Exception(self.name+" resFunc called before constant values were populated")

                if set(inputs.keys()) != set(unsetVarNames):
                    raise Exception(self.name+" resFunc called incorrectly: got "+str(sorted(inputs.keys()))+" expected "+str(sorted(unsetVarNames)))
                firstCall = False

            fullResFuncDict = self.setParamValues.copy()
            for var,val in inputs.items():
                # print(var,":",val)
                p = self.var2param[var]
                fullResFuncDict[p] = val

            return self.resFunc(**fullResFuncDict)

        compressedResFunc.__name__ = self.resFunc.__name__+"_compressed"

        return compressedResFunc, unsetVarNames


    def getSolverReadyFunc(self, setVarNames):

        # (1) Reduce from the full residual function to one with the
        #     specified variables "baked in".
        #     varList is the list of variables that must be set.
        #     [The values themselves are set later. This just generates
        #     a function with the right signature, and with access to
        #     a dictionary that will eventually hold the constant values]
        compressedResFunc, varList = self.makeCompressedResFunc(setVarNames)

        # (2) Convert from kwarg-only function to arg-only function,
        #     so that scipy solvers work. This then requires input to
        #     be passed in in a consistent order, so we return that too,
        #     so that the outer solver can reorder things properly.
        solverReadyFunc, funcVarOrder = mkKwarglessResFunc(compressedResFunc, varList)

        return solverReadyFunc, funcVarOrder

    def setSpecifiedValues(self, setVarValuesDict):
        # Call right before solving to set the values of the variables
        # that are being specified.

        # Convert {var:value} into {param:value}
        self.setParamValues = {}
        for varName, value in setVarValuesDict.items():
            # Get the equation parameter name that this variable has been
            # mapped to
            eqnParamName = self.var2param[varName]
            self.setParamValues[eqnParamName] = value


    # Four levels of indirection from original function
    #   L = f(...)
    #   resL = f(** L,...) --> makes into residual form
    #   compressedRes = f(**unknowns) --> remaps *AND* simplifies passing 1 unknown
    #   scipyRes = f(*unknowns) [in order] --> allows scipy to call

def mkKwarglessResFunc(kwargFunc, expectedArguments):
    # Fourth level of indirection from original function
    #
    order = expectedArguments[:] # Could reorder if we wanted, but no reason
    def kwargless(inputs):
        kwargs = {o:i for o,i in zip(order, inputs)}
        return kwargFunc(**kwargs)
    kwargless.__name__ = kwargFunc.__name__ + "_kwargless"

    return kwargless, order

class System:
    # A system of equations, with some variables set and others to be
    # solved for. This class manages the "closure" of the equations, given
    # the prescribed variables, and creates a more compact system of
    # residual functions only in terms of the unknown variables.
    #
    # Meanwhile, it does basic checks for solvability of the system:
    #   [e.g. are there the right number of variables and equations, and is
    #    the dependency matrix invertible (full rank).]
    #
    # Finally, it invokes a solver -- various ones from scipy can be used.


    def __init__(self, equations, constVariables):
        '''
        constVariables can be either the Variables themselves or
        their names
        '''
        if not all(isinstance(eq, Equation) for eq in equations):
            raise Exception("Pass in equations")

        if all(isinstance(v, Variable) for v in constVariables):
            specifiedVarNames = [v.name for v in constVariables]
        elif all(type(v) == str for v in constVariables):
            specifiedVarNames = constVariables[:]
        else:
            raise Exception("Must pass in Variables or variable string names. Got: "+str(constVariables))

        # Given the specified values (constVals), create new solver-ready
        # functions with those values baked in.
        # The functions have to be called with parameters in the correct
        # order, so the equation maker returns the order in which the functions
        # will expect the parameters to be passed in.
        resFuncs = []
        funcVarOrders = [] # order functions expect vars to be passed in
        for eq in equations:
            resFunc, var_order = eq.getSolverReadyFunc(specifiedVarNames)
            resFuncs.append(resFunc)
            funcVarOrders.append(var_order)


        # Determine the total set of variables to be solved for.
        allVars = set()
        for theseVars in funcVarOrders:
            allVars.update(theseVars)

        nVar = len(allVars)
        nEq = len(equations)

        # Check that all specified values were recognized by at least one eqn.
        # and determine the total number of "participating" variables,
        # including those set by the user.
        nParams = nVar # initialize, increment later

        constValReferenced = {var:False for var in specifiedVarNames}
        for eq in equations:
            processedconsts = eq.setVarNames
            for var in processedconsts:
                constValReferenced[var] = True

        for constVar, refd in constValReferenced.items():
            if refd:
                nParams += 1
            else:
                print("Warning: Ignoring set variable '{}', which is not in any equation".format(constVar))

        # And an arbitrary order for the specified variables
        setVarOrder = sorted(specifiedVarNames)

        self.nParams = nParams
        self.equations = equations
        self.setVarOrder = setVarOrder

        # Set an arbitrary order of variables for the solver.


        self.resFuncs = resFuncs
        self.funcVarOrders = funcVarOrders

        unknownVarsOrder = sorted(list(allVars))
        self.analyze(unknownVarsOrder)

        self.setUnknownVarOrder(unknownVarsOrder)

        self.minBoundDict = {}
        self.maxBoundDict = {}

        self.verbose = False

        # Disable solver until the constant values have been set!
        if len(setVarOrder) > 0:
            self.constValuesSet = False


    def setUnknownVarOrder(self, unknownVarsOrder):
        self.unknownVarsOrder = unknownVarsOrder

        def combinedResFunc(x):
            '''
            The residual function for the multi-equation problem,
            returning one residual per equation
            '''
            valDict = unpackToDict(unknownVarsOrder, x)
            # valDict = {var:val for var,val in zip(unknownVarsOrder,x)}

            all_residuals = []
            for spResFunc, varOrder in zip(self.resFuncs, self.funcVarOrders):
                try:
                    xEq = packDictStrict(varOrder, valDict)
                    # xEq = [valDict[v] for v in varOrder]
                except KeyError:
                    print("Var in x:",x)
                    print("doesn't belong in valDict:", valDict.keys())
                    raise
                res = spResFunc(xEq)
                all_residuals.append(res)

            return all_residuals

        # Save some stuff
        self.resFunc = combinedResFunc


    def getAdjacencyMatrix(self, unknownVarsOrder):
        nParams = self.nParams
        equations = self.equations
        setVarOrder = self.setVarOrder
        nEq = len(equations)
        nVar = len(unknownVarsOrder)


        M = np.zeros((nParams,nParams)) # (eq, variable)
        for i,eq in enumerate(equations):
            eqVars = list(eq.var2param.keys()) # variables in this equation
            for j,v in enumerate(unknownVarsOrder): # unknown?
                if v in eqVars:
                    M[i,j] = 1
            for j,v in enumerate(setVarOrder): # or fixed?
                if v in eqVars:
                    M[i,j+nVar] = 1
        for i,varName in enumerate(setVarOrder):
            j = nVar + i
            M[i+nEq,j] = 1
        return M

    def getAdjacencyMatrixSmall(self,unknownVarsOrder):

        # N = np.zeros((nVar,nVar)) # (eq, variable)
        # for i,eq in enumerate(equations):
        #     eqVars = list(eq.var2param.keys())
        #     for j,v in enumerate(unknownVarsOrder):
        #         if v in eqVars:
        #             N[i,j] = 1

        # if not np.array_equal(M[:nVar,:nVar], N):
        #     raise Exception("Huh?")
        nVar = len(unknownVarsOrder)
        M = self.getAdjacencyMatrix(unknownVarsOrder)
        return M[:nVar,:nVar]

    def analyze(self, unknownVarsOrder):

        nParams = self.nParams
        equations = self.equations
        setVarOrder = self.setVarOrder

        nVar = len(unknownVarsOrder)

        nEq = len(equations)

        print("System has {} equations involving {} variables".format(nEq,nVar))
        print("Equation names: ")
        for i,eq in enumerate(equations,start=1):
            print("{:2d}) {}".format(i,eq.name.replace("eqn_","")))
        print("Variable names: ")
        for i,v in enumerate(unknownVarsOrder,start=1):
           print("{:2d}) {}".format(i,v.lstrip('v')))

        # Check for obvious unsolvability
        if nVar > nEq:
            raise Exception("Underspecified: Must set {} more variable(s) (or provide {} more equation(s))".format(nVar-nEq,nVar-nEq))
        elif nVar < nEq:
            raise Exception("Overspecified: Must set {} fewer variable(s) (or remove {} equation(s))".format(nEq-nVar,nEq-nVar))

        M = self.getAdjacencyMatrix(unknownVarsOrder)

        N = self.getAdjacencyMatrixSmall(unknownVarsOrder)

        # Try replacing 1's with random numbers, because sometimes a large set
        # of valid equations gives a non-invertible matrix when it's filled
        # with ones. Not sure why yet.
        Mfloat = np.zeros((nParams,nParams)) # (eq, variable)
        for i in range(nParams):
            for j in range(nParams):
                if M[i,j] != 0:
                    Mfloat[i,j] = random.random()

        Nfloat = np.zeros((nVar,nVar)) # (eq, variable)
        for i in range(nVar):
            for j in range(nVar):
                if N[i,j] != 0:
                    Nfloat[i,j] = random.random()

        def printN(N):
            print("--> ",unknownVarsOrder)
            print(N)

        def printM(M):
            print("Top-> ",unknownVarsOrder, setVarOrder)
            print(M[:nEq,:])
            print("Side: equations")

        # print("DET: Mfloat ",np.linalg.det(Mfloat))
        # print("DET: Nfloat ",np.linalg.det(Nfloat))

        valid = False
        def passedDeterminant(matrix):
            return abs(np.linalg.det(matrix)) > 1e-15

        mpassed = passedDeterminant(Mfloat)
        npassed = passedDeterminant(Nfloat)
        if (mpassed ^ npassed):
            print("Conflicting N and M results:")
            print("DET: Mfloat ",np.linalg.det(Mfloat))
            print("DET: Nfloat ",np.linalg.det(Nfloat))
            valid = True
        elif not mpassed:
            valid = False
        else:
            valid = True


        printN(N)
        if not valid:
            # printM(M)
            printN(N)

            diagnosed = False
            # Try to diagnose:
            for i in range(nEq):
                for j in range(i+1,nEq):
                    if np.array_equal(N[i], N[j]):
                        print("Eqn {} is the same as eqn {}".format(i+1,j+1))
                        diagnosed = True
            if not diagnosed:
                r =  np.linalg.matrix_rank(Nfloat)
                print("RANK Nfloat = {}/{}".format(r,len(Nfloat)))

            raise Exception("This is an invalid set of variables to specify. Equations are not solvable")


    def setConstVals(self, constVals):

        if all(isinstance(v, Variable) for v in constVals.keys()):
            constVals = {k.name:v for k,v in constVals.items()}
        elif all(type(v) is str for v in constVals.keys()):
            pass
        else:
            raise Exception("Must pass in Variables or variable string names. Got: "+str(constVals))

        # Finally, register the specific constant values
        # Later this will be moved farther out when I decompose this
        # function into many calls with different inputs!
        for eq in self.equations:
            # And tell this equation what the constants are
            theseConstValues = {}
            for varName in eq.setVarNames:
                if varName not in constVals:
                    raise Exception("Equation '{}' expected constant variable '{}' to be set. Got variables: {}".format(eq.name,varName,list(constVals.keys())))
                theseConstValues[varName] = constVals[varName]

            eq.setSpecifiedValues(theseConstValues)

        self.constValuesSet = True

    def randomGuess(self, initGuess, lbs, ubs):
        guess = []
        for xInit, lb, ub in zip(initGuess,lbs,ubs):

            if ub < 1e8 and lb > -1e8:
                pass
            elif lb > -1e8:
                ub = xInit + (xInit-lb)*2
            elif ub < 1e8:
                lb = xInit - (ub-xInit)*2
            else:
                ub = xInit*2
                lb = xInit*0.5
                if ub < lb:
                    lb, ub = ub, lb
            guess.append(random.random()*(ub-lb)+lb)


        return guess

    def solve(self, guessDict={}, minBoundDict={}, maxBoundDict={}, verbose=False):
        if not self.constValuesSet:
            assert len(self.setVarOrder) > 0
            raise Exception("Must set the constant values: "+str(self.setVarOrder))

        self.verbose = verbose

        # If guessDict was passed in as {Variable: guess}, which is more
        # convenient, convert it to {string: guess}
        newGuessDict = {}
        for k,v in guessDict.items():
            if isinstance(k, Variable):
                newGuessDict[k.name] = v
            elif type(k) is str:
                newGuessDict[k] = v
            else:
                raise Exception("Invalid guessDict format: type="+type(k))
        guessDict = newGuessDict

        # Defaults of 1 seem safer sometimes, but it's really sketchy
        guess = packDict(self.unknownVarsOrder, guessDict, defaultval=1)
        if len(guess) != len(self.unknownVarsOrder):
            raise Exception("Bad guess formatting!")

        lb = packDict(self.unknownVarsOrder, minBoundDict, defaultval=float("-inf"))
        ub = packDict(self.unknownVarsOrder, maxBoundDict, defaultval=float("inf"))

        bounds = (lb, ub)

        def solver(x0, bnds=None):
            ''' pick your poison '''
            # Only least_squares supports bounds naturally, but you have
            # to be careful that it actually solved the problem!
            xSoln, failed, msg = solve_scipy_lstsq(self.resFunc, x0, bounds=bnds, maxIter=2000)
            # xSoln, failed, msg = solve_scipy_root(self.resFunc, x0)
            # xSoln, failed, msg = solve_scipy_fsolve(self.resFunc, x0)
            # xSoln, failed, msg = solve_broyden(self.resFunc, x0)

            return xSoln, failed, msg

        # Optional, retry with randomize guesses
        maxTries = 1
        x0 = guess
        for i in range(maxTries):
            xSoln, failed, msg = solver(x0, bnds=bounds)
            if not failed:
                if i > 0:
                    print("Succeeded in {} tries".format(i+1))
                break
            if i < maxTries - 1:
                x0 = self.randomGuess(guess, lb, ub)
        else:
            raise Exception("Failed {} times: msg={}".format(maxTries,msg))


        # Remap to a dict for easy parsing
        xSolnDict = unpackToDict(self.unknownVarsOrder, xSoln)

        return xSolnDict


def solve(equations, guessDict={}, constVals={}):

    sys = System(equations, list(constVals.keys()))
    sys.setConstVals(constVals)
    x = sys.solve(guessDict)
    return x


def solve_scipy_lstsq(resFunc, initGuess, bounds=None, maxIter=1000, quiet=False):
    # This is the gold standard root solver in scipy.
    # It allows bounds, unlike other solvers, and has excellent
    # automatic rescaling of the system to improve convergence.

    if bounds:
        if len(bounds) != 2:
            raise Exception("Bad bounds spec")
        N = len(initGuess)
        if len(bounds[0]) != N:
            raise Exception("Bad LB bounds spec")
        if len(bounds[1]) != N:
            raise Exception("Bad UB bounds spec")
        # Autofix guess, if bounds violated
        for i in range(N):
            if not (bounds[0][i] <= initGuess[i] <= bounds[1][i]):
                if bounds[0][i] == -NOLIMIT:
                    initGuess[i] = bounds[1][i] - 1
                elif bounds[1][i] == NOLIMIT:
                    initGuess[i] = bounds[0][i] + 1
                else:
                    initGuess[i] = 0.5*(bounds[0][i]+bounds[1][i])
                print("Warning: initial guess for var index {} outside bounds, resetting to {}".format(i,initGuess[i]))

    vLevel = 1 # 0/1/2

    # 'jac' does automatic variable scaling. Seems to really help
    # convergence for even these moderately badly scaled aero problems
    # optRes = scipy.optimize.least_squares(resFunc, initGuess, bounds=bounds,  max_nfev=maxIter,verbose=vLevel,ftol=1e-8,xtol=1e-8,x_scale='jac')
    optRes = scipy.optimize.least_squares(resFunc, initGuess, bounds=bounds,  max_nfev=maxIter,verbose=vLevel,ftol=1e-10,gtol=1e-10,xtol=1e-10,x_scale='jac')
    message = str(optRes.status)+" "+str(optRes.message)
    failed = False
    if not optRes.success:
        failed = True

    COST_TOL_WARN = 1e-7
    COST_TOL_ERR = 1e-4
    if optRes.status > 0 and optRes.cost > COST_TOL_ERR:
        failed = True
        message += " Converged, but no solution found. Cost = {}".format(optRes.cost)
    elif optRes.status > 0 and optRes.cost > COST_TOL_WARN:
        print(" Converged, but to poor solution. Cost = {}".format(optRes.cost))
    elif optRes.status == 0: # maxfev reached
        failed = True

    # print(optRes.status,"COST:",optRes.cost)
    # if (optRes.cost > 1e-6):
    #     raise Exception("Failed to converge")
    return optRes.x, failed, message

    # fun, x0, jac='2-point', bounds=(-inf, inf), method='trf', ftol=1e-08, xtol=1e-08, gtol=1e-08, x_scale=1.0, loss='linear', f_scale=1.0, diff_step=None, tr_solver=None, tr_options={}, jac_sparsity=None, max_nfev=None, verbose=0, args=(), kwargs={})

def solve_scipy_fsolve(resFunc, initGuess, maxIter=1000, quiet=False):

    xopt, infodict, ierr, msg = scipy.optimize.fsolve(resFunc, initGuess, full_output=True, maxfev=maxIter)
    # print(infodict["fvec"])
    print("RES",math.sqrt(sum(_*_ for _ in infodict["fvec"])))

    failed = False
    message = str(ierr)+" "+msg
    if ierr != 1:
        failed = True

    return xopt, failed, message

def solve_scipy_root(resFunc, initGuess, maxIter=1000, quiet=False):

    options = {'maxfev':maxIter}
    # scipy.optimize.show_options('root', method='hybr',disp=True)
    # import sys
    # sys.exit()
    optRes = scipy.optimize.root(resFunc, initGuess, options=options)

    # (fun, x0, args=(), method='hybr', jac=None, tol=None, callback=None, options=None)
    failed = False
    message = str(optRes.status)+" "+str(optRes.message)
    if not optRes.success:
        failed = True

    if not failed:
        res = math.sqrt(sum(_*_ for _ in infodict["fvec"]))
        if res > 1e-6:
            message += " converged, but not to solution: res = "+str(res)
            failed = True

    # if optRes.status > 0 and optRes.cost > 1e-6:
    #     raise Exception("Converged, but no solution found: "+optRes.message)
    # elif optRes.status == 0:
    #     raise Exception("hit max nfev "+optRes.message)

    # print(optRes.status,"COST:",optRes.cost)

    return optRes.x, failed, message



def solve_broyden(resFunc, initGuess, maxIter=1000, quiet=False):
    optRes = scipy.optimize.broyden2(resFunc, initGuess, verbose=True, maxiter=maxIter)

    message = str(optRes.status)+" "+str(optRes.message)
    failed = False
    if not optRes.success:
        failed = True

    if optRes.status > 0 and optRes.cost > 1e-6:
        failed = True
        message += " Converged, but no solution found."
    elif optRes.status == 0: # maxfev reached
        failed = True

    # print(optRes.status,"COST:",optRes.cost)
    # if (optRes.cost > 1e-6):
    #     raise Exception("Failed to converge")
    return optRes.x, failed, message

    # return xopt


def checkValidFunction(func):
    '''
    Currently the input function must have a fixed number of arguments,
    and all must be required. (No *args, **kwargs, or default values)
    '''
    argSpec = inspect.getargspec(func)

    if argSpec.varargs or argSpec.keywords:
        raise Exception("Can't make a residual function from a function accepting variable argument lists")
    if argSpec.defaults:
        # This seems to work, but it's worrisome, so I will raise an
        # exception until I'm sure...
        raise Exception("Can't make a residual function from a function with argument defaults yet")


def makeResidualFunc(func):
    '''
    Convert a function z = f(x,y) into res = z - f(x,y)

    Currently the input function must have a fixed number of arguments,
    and all must be required. (No *args, **kwargs, or default values)
    '''
    argSpec = inspect.getargspec(func)

    inputParamNames = argSpec.args
    N = len(inputParamNames)
    outputParamName = func.__name__
    resFuncName = "res_"+func.__name__

    errorCheck = True
    verbose = False
    # Notes:
    #   1) Want to forbid args, and only allow kwargs
    #      Could remove *args parameter to enforce this, but for now,
    #      I want to explicitly catch and warn about this.
    #   2) Can eventually disable error checking at every call for speed
    #      -- Or could disable after the first valid call.
    def resFunc(*args, **kwargs):
        if errorCheck:
            if args:
                raise Exception("Can't call residual function "+resFuncName+" without kwargs")
            if outputParamName not in kwargs:
                raise Exception("Residual function "+resFuncName+" missing parameter "+outputParamName)
            if len(kwargs) != N + 1:
                raise Exception("Got wrong number of ")

        outputValue = kwargs[outputParamName]
        del kwargs[outputParamName]

        # Invoke the original function with the original argument signature
        # and subtract from the predicted output to get the residual
        res = outputValue - func(**kwargs)

        if verbose:
            print("{:.2e}".format(res),resFunc.__name__,kwargs,outputValue)
        return res

    resFunc.__name__ = resFuncName

    return resFunc

