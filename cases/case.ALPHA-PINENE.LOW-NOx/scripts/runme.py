'''==========================================================================
                THIS IS THE INTERFACE FOR THE PYTHON WRAPPER
=========================================================================='''

import sys
import numpy as np
import scipy.optimize as sio

from functions.f_comp_model import comp_model
from functions.f_fit_model  import fit_model
from functions.f_feed4ward  import feed4ward

import variables as ns

# SOM PARAMETERS:
# =======================================================================================
xx = np.array([1.588,0.112,0.431,0.002,0.552,0.015])

# FIT CONFIGURATIONS:
# =======================================================================================
# CONTROL SWITCH
# 0 FOR FORWARD
# 1 FOR FITTING:
switch = 0

# UPPER AND LOWER BOUNDS:
bounds = ((1.0,0.0,1e-5,1e-5,1e-5,1e-5),(3.0,10.0,1e4,1e4,1e4,1e4))

# RELATIVE DIFFERENTIAL STEP FOR FITTING:
ds = ([0.001,0.001,0.001,0.001,0.001,0.001])

# RUN THE MODEL (FIT/FORWARD/GLOBAL SEARCH):
# =======================================================================================
if switch == 1:
    # COMPILE MODEL:
    comp_model('SOM-TOMAS')
    
    # REFRESH RECORD FILE: 
    f1 = open('outputs/records.dat','w').close()
    
    # INITIAL GUESS:
    x0 = np.array([1.788,0.112,431.,2.,552.,15.])
    
    # RUN FITTING:
    results = lsq(fit_model,x0,diff_step=ds,bounds=bounds,loss='linear',\
                  ftol=1e-6,xtol=1e-6,verbose=0,max_nfev=10)
    
    # PRINT FITS:
    fits = results.x
    print('xx = np.array([%.5f,%.5f,%.5f,%.5f,%.5f,%.5f])'%tuple(list(fits)))
    
    # SAVE FITS:
    np.save('outputs/out.FITS.npy',results.x)

elif switch == 0:
    # COMPILE MODEL:
    comp_model('SOM-TOMAS')
    
    # VISUALIZE:
    feed4ward(ns,FREEPARAMS=xx); print('FINISHED')
