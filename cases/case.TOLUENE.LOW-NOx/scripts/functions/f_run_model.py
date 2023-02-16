'''==================================================================================
                        THIS FUNCTION RUNS THE SOM-TOMAS MODEL
=================================================================================='''

import os
import time   as tm
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt

from functions.f_write_input import write_input

def run_model(ns):
    
    # WRITE THE INPUT FILE:
    # ===============================================================================
    write_input(ns)
    
    # RUN THE MODEL:
    # ===============================================================================
    print('MODEL IS RUNNING...')

    # && REMOVE PREVIOUS OUTPUTS:
    os.system('rm ../outputs/* 2>/dev/null')

    # && DEFAULT INPUT AND OUTPUT PATHS:
    PATHi = './inputs/'
    PATHo = './outputs/'

    line = PATHi + '\n' + PATHo

    # >> RUN MODEL:
    os.system('cd ../ ; echo -e "%s" | ./box.exe'%line)
    
    # READ OUTPUT FILES:
    # ===============================================================================
    # DATAFRAME FOR SOA:
    df_SOA = pd.read_csv('../outputs/out.SOA.%s'%ns.runname,delim_whitespace=True)
    
    # DATAFRAME FOR NUMBER CONC.:
    df_NCONC = pd.read_csv('../outputs/out.NCONC.%s'%ns.runname,delim_whitespace=True)

    # DATAFRAME FOR NUMBER DIST.:
    df_NDIST = pd.read_csv('../outputs/out.NDIST.%s'%ns.runname,delim_whitespace=True)
    
    df_GCMASS = pd.read_csv('../outputs/out.GCMASS.%s'%ns.runname,delim_whitespace=True)
    
    # OUTPUTS:
    # ===============================================================================
    OUT = {}
     
    OUT['df_SOA']   = df_SOA
    OUT['df_NCONC'] = df_NCONC
    OUT['df_NDIST'] = df_NDIST
    
    return OUT
    
