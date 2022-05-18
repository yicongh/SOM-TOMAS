'''========================================================================================
                      THIS IS THE LIST OF INPUT VARIABLES FOR SOM-TOMAS
========================================================================================'''

import sys,json
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt

import data.mdata as mdata

from functions.f_get_root import get_root

# VARIABLES - GENERAL:
#==========================================================================================
# RUN NAME:
runname = 'example'

# MODEL TIME STEP [s] AND
# OUTPUT FREQUENCY [# of steps]:  
mdt  = 10.
freq = 60

# TOTAL SIMULATION TIME [s]:
t_end = 3600.*10.

# PRESSURE [Pa],
# TEMPERATURE [K],
# RELATIVE HUMIDITY AND
# REACTOR VOLUME [cm3]:
Pres    = 101325.0
Temp    = 298.0
RH      = 0.2
BoxVol  = 24e6

# ACCOMMODATION COEFFICIENT,
# BULK DIFFUSIVITY [m2 s-1] AND
# OH UPTAKE COEFFICIENT:,
Alpha = 1.0        
Db    = 1e-10
yOH   = 0.0

# SWITCHES FOR COAGULATION,
# VAPOR WALL LOSS,
# PARTICLE WALL LOSS,
# HETEROGENEOUS CHEMISTRY,
# ENDOGENOUS DB,
# FIRST GENERATION ONLY, AND
# ABSORBING WALLS [0 or 1]:
mCOAG  = 1
mVWL   = 1
mPWL   = 1
mV2PWL = 1
mHET   = 0
mOLIG  = 0
mENDB  = 0
mGEN1  = 0
mWSORB = 1

# VAPOR WALL LOSS RATE [s-1]:
kvap_on = 4e-4

# OLIGOMER FORMATION [cm3 s-1] AND
# DISSOCIATION [s-1] RATES:
k_f = 1e-24*0.
k_r = 30.0/3600.*0.

# TIME PROFILE - OXIDANT CONCENTRATION:
#==========================================================================================
# EXPONENTIAL DECAY - CONC. [cm-3]** AND 
# DECAY RATE [s-1]**:
A1 = 3.30e6
B1 = 0.0
A2 = 0.0
B2 = 0.0

# PROFILE:
OH = np.array([A1*np.exp(-B1*i) + A2*np.exp(-B2*i) for i in np.arange(0.,t_end,mdt)])

# NUCLEATION RATES:
#==========================================================================================
# MODE:
mNUCL = 1

# FOR NUCL. BASED ON 
# ELVOC [MODE 0]:
NUCL_A = 3.0
NUCL_B = 1.5
NUCL_C = 100.

# FOR NUCL. BASED ON 
# PRESCRIBED PROFILE [MODE 1]:
UU = 1.0
SS = 1.0
NN = 0.0

NUC = np.array([NN/SS/np.sqrt(2.*np.pi)*np.exp(-0.5*((i - UU)/SS)**2.) for i in np.arange(0.,t_end,mdt)])

# SIZE PROFILE - BIN DIAMETERS AND WIDTH:
#==========================================================================================
# NUMBER OF BINS**:
nBINS = 36

# OVERALL LOWER** AND UPPER** MASS BOUNDS [kg]:
bLOWER = 2.5023e-23
bUPPER = 7.4142e-15

# BIN BOUNDS:
bBOUND = np.logspace(np.log10(bLOWER),np.log10(bUPPER),num=(nBINS+1))

# INITIAL SEED DENSITY [kg m-3]**:
RhoSeed = 1770.

# BIN DIAMETERS [m]:
bDIAM = np.array([np.sqrt(bBOUND[i]*bBOUND[i+1]) for i in range(nBINS)])
bDIAM = np.array([(i/RhoSeed*6./np.pi)**(1./3.)  for i in bDIAM])

# LOG BIN WIDTH:
bWLOG = np.log10(bBOUND[1:]/bBOUND[:-1])/3.

# SIZE PROFILE - INITIAL SIZE DISTRIBUTION:
#==========================================================================================
# MEDIAN DIAMETER [um]**,
# SPREAD** AND
# NUMBER CONCENTRATION [# cm-3]**:
N1   = 35000.
CMD1 = 0.093
Sig1 = 1.7

N2   = 0.
CMD2 = 0.09
Sig2 = 1.7

# ABSORBING SEED [0 or 1]**:
mABSORB = 0

# BI-MODAL DISTRIBUTIONS:
MODE1 = 2.303*N1/np.sqrt(2.*np.pi)/np.log(Sig1)* \
        np.exp(-np.log(bDIAM*1e6/CMD1)**2./2./np.log(Sig1)**2.)

MODE2 = 2.303*N2/np.sqrt(2.*np.pi)/np.log(Sig2)* \
        np.exp(-np.log(bDIAM*1e6/CMD2)**2./2./np.log(Sig2)**2.)

# INITIAL DIST. [# box-1]:
bDIST0 = (MODE1 + MODE2)*bWLOG*BoxVol

# SIZE PROFILE - PARTICLE WALL-LOSS RATES:
#==========================================================================================
df = mdata.df_PWL
ss = df['size'].values
bb = df['beta'].values

bPWL = np.interp(bDIAM,ss,bb,left=bb[0],right=bb[-1])

# PRECURSORS - READ INPUT TABLES:
#==========================================================================================
df_prec = pd.read_excel('data/precursors.xlsx',sheet_name='precs').sort_values(by=['Code'])
df_grid = pd.read_excel('data/precursors.xlsx',sheet_name='grids').sort_values(by=['Grid Code'])

# MODIFY THE DATAFRAME:
df_prec.at[0,'kOH']        = 5.63e-12
df_prec.at[0,'Emiss']      = 0.053
df_prec.at[0,'MW']         = 92.14
df_prec.at[0,'Carbon Num'] = 7
df_prec.at[0,'f_HOMS']     = 0.001

# PRECURSORS - MAKE INPUT ARRAYS:
#==========================================================================================
# NAMES FOR SOM GRIDS:
grid_names = df_grid['Grid Code'] + 'SOMG'
grid_names = grid_names.str.cat(sep=' ')

# DLVP FOR SOM GRIDS:
grid_dlvps = df_grid['dlvp'].astype(str).str.cat(sep=' ')

# NUMBER OF SOM GRIDS:
ngrids = len(grid_names.split())

# NAMES FOR SOM PRECURSORS:
prec_names = df_prec['Code'].str.cat(sep=' ')

# INITIAL CONC. OF PRECURSORS (ppm):
prec_emiss = df_prec['Emiss'].astype(str).str.cat(sep=' ')

# NUMBER OF PRECURSORS:
nprecs = len(prec_names.split())

# NUMBER OF PRECURSORS IN EACH GRID:
void,uprecs = np.unique(df_prec['Grid Code'].values,return_counts=True)
    
# APPEND TO PRECURSORS:
df_prec = pd.merge(df_prec,df_grid,on='Grid Code',how='outer')

# GET SOM PARAMETERS:
params = np.empty(shape=(0))

df_aggr = df_prec[['Grid Code','mfrag','p1','p2','p3','p4']]
df_aggr = df_aggr.drop_duplicates(subset='Grid Code',keep='first')

for i,grid_name in enumerate(df_aggr['Grid Code']):
    
    line   = df_grid.loc[df_grid['Grid Code']==grid_name][['mfrag','p1','p2','p3','p4']].values[0]
    params = np.append(params,line,axis=0)

# SOM PRODUCTS:
#==========================================================================================
df_aggr  = df_prec.groupby('Grid Code').aggregate({'Carbon Num':'max'})
cmaxgrid = df_aggr['Carbon Num'] + 1
omaxgrid = np.ones(len(cmaxgrid))*7
iorggrid = [i*7 - 9 for i in cmaxgrid]

# NUMBER OF ORGANIC SPECIES,
# NUMBER OF DIAGNOSTIC SPECIES AND
# NUMBER OF ALL SPECIES:
iORG  = sum(iorggrid)*2 + 1
iDIAG = 2
iCOMP = iORG + iDIAG + 1
