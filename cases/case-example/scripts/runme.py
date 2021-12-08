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

# FIT CONFIGURATIONS:
# =======================================================================================
# CONTROL SWITCH
# 0 FOR FORWARD
# 1 FOR FITTING:
switch = int(sys.argv[1])

# UPPER AND LOWER BOUNDS:
bounds = ((-7.,12.))

# RELATIVE DIFFERENTIAL STEP FOR FITTING:
ds = ([0.001])

# RUN THE MODEL (FIT/FORWARD/GLOBAL SEARCH):
# =======================================================================================
if switch == 1:
    # CLEAR RECORD FILE:
    f1 = open('outputs/records.dat','w').close()
    
    # COMPILE MODEL:
    comp_model(sys.argv[2])
    
    # INITIAL GUESS:
    x0 = np.array([7.0])

    # RUN FITTING:
    results = sio.least_squares(fit_model,x0,diff_step=ds,bounds=bounds,loss='linear',\
                                ftol=1e-8,xtol=1e-8,verbose=0,max_nfev=10,kwargs={'ns':ns})
    
    # PRINT FITS:
    fits = results.x
    print('xx = %.2e'%tuple(list(fits))); np.save('outputs/fits.npy',fits)
    
elif switch == 0:
    # COMPILE MODEL:
    comp_model(sys.argv[2])
    
    # VISUALIZE:
    feed4ward(ns); print('FINISHED')
