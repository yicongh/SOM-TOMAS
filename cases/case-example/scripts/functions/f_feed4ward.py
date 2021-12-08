'''===============================================================================
                      THIS FUNCTION RUNS THE MODEL FORWARD
================================================================================'''

import os,sys
import time   as tm
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt

import matplotlib #; matplotlib.use('Agg')

from functions.f_run_model import run_model
from functions.f_save_data import save_data

from data.mdata import mdata

def feed4ward(ns):
    
    # RUN MODEL AND EXTRACT RESULTS:
    # =============================================================================
    OUT = run_model(ns)
    
    df_SOA   = OUT['df_SOA']
    df_NCONC = OUT['df_NCONC']
    df_NDIST = OUT['df_NDIST']

    t_obs = mdata.t_soa
    y_obs = mdata.y_soa
    
    # SAVE DATA:
    # =============================================================================
    os.system('cd outputs/ ; rm *')

    df1 = df_SOA
    df2 = df_NDIST
    df3 = pd.DataFrame({'size':ns.bDIAM})
    df5 = df_NCONC
    
    tt = df_NDIST.values[:,0]
    yy = df_NDIST.values[:,1:]

    index = np.where(ns.bDIAM*1e9 >= 10.)[0][0]
    
    NUM = np.sum(yy[:,index:]*ns.bWLOG[index:],axis=1)
    
    df5 = pd.DataFrame({'time':tt,'NUM_SUSP':NUM})

    with pd.ExcelWriter('outputs/MODEL_OUTS.xlsx') as w:
        df1.to_excel(w,sheet_name='df_SOA')
        df2.to_excel(w,sheet_name='df_NDIST')
        df3.to_excel(w,sheet_name='df_SIZE')
        df5.to_excel(w,sheet_name='df_NCONC'); w.close()
    
    # MAKE PLOTS:
    # =============================================================================
    plt.figure()
    ax = plt.gca()
        
    tt = df_SOA['time'].values*3600.
    y1 = df_SOA['SOA']

    ax.plot(tt,y1,color='red',label='Susp.')
    
    ax.set_xlim((0.,None))
    ax.set_ylim((0.,None))
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('SOA ($\mu$g m$^{-3}$)')
    ax.legend(loc=2)
    plt.savefig('outputs/soa.png',bbox_inches='tight')

    # MAKE PLOTS:
    # =============================================================================
    plt.figure()
    ax = plt.gca()
        
    tt = df_NDIST.values[:,0]
    yy = df_NDIST.values[:,1:]
    
    NUM = np.sum(yy[:,10:]*ns.bWLOG[10:],axis=1)
    
    ax.plot(tt,NUM,color='r')
    
    print('NUM = %.0f'%NUM[-1])
    
    ax.set_xlim((0.,None))
    ax.set_ylim((0.,None))
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Number Conc. (cm$^{-3}$)')
    ax.legend(loc=2)
    #plt.savefig('outputs/soa.png',bbox_inches='tight')

    # MAKE PLOTS:
    # =============================================================================
    plt.figure()
    ax = plt.gca()
        
    tt = df_SOA['time']
    y1 = df_SOA['O2C']

    ax.plot(tt,y1,color='red',label='Susp.')
    
    ax.set_xlim((0.,None))
    ax.set_ylim((0.,None))
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('O:C')
    ax.legend(loc=2)
    plt.savefig('outputs/o2c.png',bbox_inches='tight')

    # MAKE PLOTS:
    # =============================================================================
    plt.figure()
    ax = plt.gca()
        
    tt = df_NDIST.values[:,0]
    y1 = df_NDIST.values[10,1:]
    y2 = df_NDIST.values[40,1:]
    y3 = df_NDIST.values[70,1:]
    y4 = df_NDIST.values[99,1:]

    #ax.plot(ns.bDIAM*1e9,y1/max(y1),marker='o')
    #ax.plot(ns.bDIAM*1e9,y2/max(y2),marker='o')
    #ax.plot(ns.bDIAM*1e9,y3/max(y3),marker='o')
    #ax.plot(ns.bDIAM*1e9,y4/max(y4),marker='o')

    #ax.plot(ns.bDIAM*1e9,y1,marker='o')
    #ax.plot(ns.bDIAM*1e9,y2,marker='o')
    #ax.plot(ns.bDIAM*1e9,y3,marker='o')
    ax.plot(ns.bDIAM*1e9,y4,marker='o',label='Modeled')

    #df = pd.read_excel('../../../../../data/DATA.LAMBE_RAW/DATA.DIST_A.xlsx')
    
    #ss = df['size'].values
    #dd = df['dist'].values
    
    #ax.plot(ss,dd/max(dd),color='grey',marker='o',linestyle='None',zorder=0)
    #ax.plot(ss,dd,color='grey',marker='o',linestyle='None',label='Measured',zorder=0)
    
    ax.set_xscale('log')
    ax.set_ylim((0.,None))
    
    ax.ticklabel_format(axis='y',style='sci',scilimits=(0.,0.),useMathText=True)

    ax.set_xlabel('Size (nm)')
    ax.set_ylabel('dN/dLogD$_{p}$ (cm$^{-3}$)')

    ax.legend(loc=2)
    
    plt.savefig('outputs/fig.SIZE_DIST.png',bbox_inches='tight')
    plt.show()
