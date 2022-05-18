'''============================================================
     THIS FUNCTION WRITES THE INPUT FILE FOR THE BOX MODEL
============================================================'''

import os

def write_input(ns):

    # WRITE THE INPUT VARIABLES:
    # =========================================================
    f1 = open('../inputs/in.vars','w')
    
    # RUN NAME:
    f1.write('%s\n'%ns.runname)

    # MODEL TIMESTEP [s] AND
    # OUTPUT FREQUENCY:
    f1.write('%f\n'%ns.mdt)
    f1.write('%i\n'%ns.freq)
    
    # TOTAL SIMULATION TIME [s]:
    f1.write('%.3f\n'%ns.t_end)

    # SWITCHES [0 or 1]:
    f1.write('%i\n'%ns.mCOAG)
    f1.write('%i\n'%ns.mVWL)
    f1.write('%i\n'%ns.mPWL)
    f1.write('%i\n'%ns.mHET)
    f1.write('%i\n'%ns.mENDB)
    f1.write('%i\n'%ns.mGEN1)
    f1.write('%i\n'%ns.mV2PWL)
    f1.write('%i\n'%ns.mABSORB)
    f1.write('%i\n'%ns.mNUCL)
    f1.write('%i\n'%ns.mWSORB)
    f1.write('%i\n'%ns.mOLIG)    
    
    # REACTOR CONDITIONS - PRESSURE [Pa],
    # TEMPERATURE [K],
    # RELATIVE HUMIDITY AND
    # REACTOR VOLUME [cm3]:
    f1.write('%.5f\n'%ns.Pres)
    f1.write('%.5f\n'%ns.Temp)
    f1.write('%.5f\n'%ns.RH)
    f1.write('%.5f\n'%ns.BoxVol)
    
    # VAPOR WALL LOSS RATE [s-1]:
    f1.write('%.5e\n'%ns.kvap_on)
    
    # AEROSOL PROPERTIES - ACCOMMODATION COEFFICIENT,
    # BULK DIFFUSIVITY [m2 s-1]:
    f1.write('%.5f\n'%ns.Alpha)
    f1.write('%.5e\n'%ns.Db)
    
    # OLIGOMER FORMATION [cm3 s-1] AND
    # DISSOCIATION [s-1] RATES:
    f1.write('%.3e\n'%ns.k_f)
    f1.write('%.3e\n'%ns.k_r)

    # NUCLEATION RATE PARAMETERIZATIONS:
    f1.write('%.3e\n'%ns.NUCL_A)
    f1.write('%.3e\n'%ns.NUCL_B)
    f1.write('%.3e\n'%ns.NUCL_C)    
    
    # SOM GRID NAMES AND DLVPS:
    f1.write('%s\n'%ns.grid_names)
    f1.write('%s\n'%ns.grid_dlvps)
    
    # PRECURSORS AND EMISSIONS [ppm]:
    f1.write('%s\n'%ns.prec_names)
    f1.write('%s\n'%ns.prec_emiss)
    
    # SOM PARAMETERS:
    line = ['%.4e' for i in range(len(ns.params))]
    line = ' '.join(line) + '\n'
    f1.write(line%tuple(ns.params))
    
    # HOMS MOLAR YIELD:
    f_homs = ns.df_prec['f_HOMS'].astype(str).str.cat(sep=' ')
    f1.write('%s\n'%f_homs)

    # PRECURSOR kOH:
    koh = ns.df_prec['kOH'].astype(str).str.cat(sep=' ')
    f1.write('%s\n'%koh)
    
    # PRECURSOR CARBON NUMBER:
    cno = ns.df_prec['Carbon Num'].astype(str).str.cat(sep=' ')
    f1.write('%s\n'%cno)

    # PRECURSOR OXYGEN NUMBER:
    ono = ns.df_prec['Oxygen Num'].astype(str).str.cat(sep=' ')
    f1.write('%s\n'%ono)
    
    # BIN MASS BOUNDS:
    f1.write('%.3e\n'%ns.bLOWER)
    f1.write('%.3e\n'%ns.bUPPER)

    # OH UPTAKE COEFFICIENT:
    f1.write('%.3e\n'%ns.yOH)
    
    f1.close()
    
    # WRITE THE OH CONCENTRATION PROFILE:
    # =========================================================
    with open('../inputs/in.OH_PROFILE','w') as ff:
        for i in ns.OH:
            ff.write('%.4e\n'%i)
     
    # WRITE THE INITIAL SIZE DISTRIBUTIONS:
    # =========================================================
    with open('../inputs/in.PSD0','w') as ff:
        for i in ns.bDIST0:
            ff.write('%.3f\n'%float(i))

    # WRITE PWL RATES:
    # =========================================================
    with open('../inputs/in.PWL','w') as ff:
        for i in ns.bPWL:
            ff.write('%.4e\n'%float(i))
    
    # WRITE THE NUCLEATION RATES:
    # =========================================================
    with open('../inputs/in.JNUC','w') as ff:
        for i in ns.NUC:
            ff.write('%.4e\n'%float(i))
