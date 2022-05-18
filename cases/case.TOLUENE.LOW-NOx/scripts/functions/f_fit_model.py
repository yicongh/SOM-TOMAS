'''===========================================================
This is the objective function for fitting the SOM-TOMAS model
==========================================================='''

import os
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt

#from data.counter import ctr

from functions.f_get_root  import get_root
from functions.f_run_model import run_model

import data.mdata as mdata

# ROOT DIRECTORY:
root = get_root()

# ITERATION COUNTER:
ctr = 0

def fit_model(xx,ns=None):
    
    # UNPACK THE PARAMETERS:
    # ========================================================    
    UU = ns.UU
    SS = ns.SS
    NN = 10.**(xx)

    ns.NUC = np.array([NN/SS/np.sqrt(2.*np.pi)*np.exp(-0.5*((i - UU)/SS)**2.) for i in np.arange(0.,ns.t_end,ns.mdt)])

    # RECORD:
    # ========================================================
    global ctr
    
    f1 = open('outputs/records.dat','a')
    f1.write('%i... \n'%(ctr))
    f1.write('NN = %.2e \n'%NN)
    f1.close()

    ctr += 1
    
    # RUN THE MODEL:
    # ========================================================
    os.system('rm ../outputs/* 2>/dev/null')
    
    outputs = run_model(ns)
    
    df_NDIST = outputs['df_NDIST']
    
    tt = df_NDIST.values[:,0]*3600.
    yy = df_NDIST.values[:,1:]

    # RESIDENCE TIME DIST.:
    df_RT = pd.read_csv('%s/data/db.OFR_RESIDENCE_TIME/d0.RT_DIST.csv'%root,header=None,sep=',')

    ttt = df_RT[0].values
    yyy = df_RT[1].values

    yyy = np.interp(tt,ttt,yyy,left=yyy[0],right=yyy[-1])
    yyy = yyy/np.sum(yyy)
    
    # AVERAGED NUMBER:
    index = np.where(ns.bDIAM*1e9 >= 10.)[0][0]
    
    NUM = np.sum(np.sum(yy[:,index:]*ns.bWLOG[index:],axis=1)*yyy)
    
    # CALCULATE OBJECTIVE FUNCTION:
    # ========================================================
    LL = np.array(['A','B','C','D','E','F'])

    index = np.where(LL == ns.runname.split('_')[1])[0][0]
    y_obs = mdata.df_mNUM['NUM'].values[index]

    diff = abs((NUM - y_obs)/y_obs)*100.
    diff = diff[~np.isnan(diff)]
    
    return diff
