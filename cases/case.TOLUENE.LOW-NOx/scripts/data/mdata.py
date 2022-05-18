'''=====================================================================
                  THIS FILE READS THE MEASURED DATA
====================================================================='''

import os
import numpy  as np
import pandas as pd

from functions.f_get_root import get_root

# DATABASE DIRECTORY:
# ======================================================================
# ROOT DIRECTORY:
root = get_root()

# DATABASE:
path = '%s/data/db.CALTECH/'%root

# READ SOA:
# ======================================================================
filename = 'SOA.TOLUENE.LOW-NOx.xlsx'

# READ DATA:
df = pd.read_excel('%s/%s'%(path,filename))

x_obs = df['time'].values
y_obs = df['SOA'].values

index = x_obs <= 10.

x_obs = x_obs[index]
y_obs = y_obs[index]

# READ PWL RATES:
# ======================================================================
filename = 'PWL_RATES.xlsx'

# READ DATA:
df = pd.read_excel('%s/%s'%(path,filename))

ss = df['size'].values*1e-9
bb = df['beta'].values/60.

df_PWL = pd.DataFrame({'size':ss,'beta':bb})

# READ INITIAL SIZE DISTRIBUTION:
# ======================================================================
filename = 'PSD0.xlsx'

# READ DATA:
df = pd.read_excel('data/%s'%filename)

ss = df['Diameter'].values*1e-9
dd = df['dNdLogDp'].values

df_PSD0 = pd.DataFrame({'size':ss,'dist':dd})

# INTERFACE:
# ======================================================================
class namespace:
    def __init__(self,**kwargs):
        self.__dict__.update(kwargs)

mdata = namespace(t_soa=x_obs,y_soa=y_obs)
