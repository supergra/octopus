#!/usr/bin/env python

# Stub for using a combination of data sources to build the
# most compact possible performance model, via the bootstrap method.
# Ultimately, we can combine:
#   1) Flight test data
#   2) Analysis, CFD, etc.
#   3) Hypothetical data

import os
from cleanJson import WriteCleanJson

model = {}

model['CD0'] = 0.04122
model['e'] = 0.848
model['m'] = 0.6831
model['b'] = -0.04104
model['throttle2torquefactor'] = [[0,0],[0.5,0.5],[1,1]] # nondim:nondim
model['CPJsq2J'] = [[0.14,0.77],[0.65,0.32]] # nondim: nondim
model['power2bsfc'] = [[0,300],[100,300]] # kW: g/hr/kW

WriteCleanJson('perfModel.json',model)
