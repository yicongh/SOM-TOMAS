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
path = '%s/data/data.NG.CALTECH/'%root

# READ SOA:
# ======================================================================
filename = 'SOA.NG.CALTECH.TOLUENE.xlsx'

# READ DATA:
df = pd.read_excel('%s/%s'%(path,filename))

x_obs = df['time'].values
y_obs = df['SOA'].values

yy1 = df['SOA_w0'].values
yy2 = df['SOA_w1'].values

index = x_obs <= 10.

x_obs = x_obs[index]
y_obs = y_obs[index]

yy1 = yy1[index]
yy2 = yy2[index]

# READ PWL RATES:
# ======================================================================
filename = 'd.PWL.50nm.csv'

# READ DATA:
df = pd.read_csv('data/%s'%filename,sep=',')

ss = df['size'].values*1e-9
bb = df['beta'].values

df_PWL = pd.DataFrame({'size':ss,'beta':bb})

# READ INITIAL SIZE DISTRIBUTION:
# ======================================================================
filename = 'PSD0.xlsx'

# READ DATA:
df = pd.read_excel('data/%s'%filename)

ss = df['Diameter'].values*1e-9
dd = df['dNdLogDp'].values

df_PSD0 = pd.DataFrame({'size':ss,'dist':dd})

# READ SOA:
# ======================================================================
path = '%s/data/db.LAMBE/'%root; filename = 'DATA.NUM_MEASURED.xlsx'

# READ DATA:
df_mNUM = pd.read_excel('%s/%s'%(path,filename))

# INTERFACE:
# ======================================================================
class namespace:
    def __init__(self,**kwargs):
        self.__dict__.update(kwargs)

mdata = namespace(t_soa=x_obs,y_soa=y_obs,yy1=yy1,yy2=yy2)
