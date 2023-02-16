'''============================================================================
               THIS FUNCTION WRITES THE MODULES FOR THE MODEL
============================================================================'''

import numpy as np

def write_modules(path,ns):
    
    # WRITE THE HEADER FOR TOMAS ARRAYS:
    # =========================================================================
    with open('%s/source.fort/raw_files/raw.mod.tomas.f90'%path,'r') as f1:
        lines = np.array(f1.read().splitlines(),dtype='object')
    
    # WRITE NUMBER OF BINS,
    # COMPONENTS, DIAGNOSTIC SPECIES AND
    # ORGANIC SPECIES:
    index = np.where(lines == '!FLAG1')[0][0]
    
    lines[index+1] = \
    '      PARAMETER(iBINS=%i, iCOMP=%i, iDIAG=%i, iORG=%i)'%(ns.nBINS,ns.iCOMP,ns.iDIAG,ns.iORG)
    
    # WRITE TO MODEL AND
    # DIAGNOSTIC FILES:
    f1 = open('%s/source.fort/modules/mod.tomas.f90'%path,'w')
    f2 = open('../outputs/diag.mod.tomas.f90','w')
    
    for line in lines:
        f1.write(line + '\n')
        f2.write(line + '\n')
        
    f1.close()
    f2.close()
    
    # WRITE THE HEADER FOR SOM ARRAYS:
    # =========================================================================
    with open('%s/source.fort/raw_files/raw.mod.som.f90'%path,'r') as f1:
        lines = np.array(f1.read().splitlines(),dtype='object')

    # NUMBER OF GRIDS AND PRECURSORS:
    index = np.where(lines=='!FLAG1')[0][0]
    lines[index+1] = \
    '      INTEGER,PARAMETER :: nSOMCLASS = %i'%ns.ngrids
    lines[index+2] = \
    '      INTEGER,PARAMETER :: nSOMPRECS = %i'%ns.nprecs
    
    # NUMBER OF PRECURSORS IN EACH GRID:
    index = np.where(lines=='!FLAG2')[0][0]
    lines[index+1] = \
    '      INTEGER,PARAMETER :: uSOMPRECS(nSOMCLASS) = (/%s/)'%(','.join(['%i']*ns.ngrids))%tuple(ns.uprecs)
    lines[index+2] = \
    '      INTEGER,PARAMETER :: uSOMCELLS(nSOMCLASS) = (/%s/)'%(','.join(['%i']*ns.ngrids))%tuple(ns.iorggrid)
    
    # FIRST INDICES:
    PREC_1ST_ALL = np.zeros(ns.ngrids)
    GRID_1ST_ALL = np.zeros(ns.ngrids)
    GRID_LST_ALL = np.zeros(ns.ngrids)
    
    PREC_1ST_ALL[0] = 0
    GRID_1ST_ALL[0] = 0 + ns.nprecs

    for i in range(1,ns.ngrids):
        PREC_1ST_ALL[i] = PREC_1ST_ALL[i-1] + ns.uprecs[i-1]
        GRID_1ST_ALL[i] = GRID_1ST_ALL[i-1] + ns.iorggrid[i-1]
    
    for i in range(ns.ngrids):
        GRID_LST_ALL[i] = GRID_1ST_ALL[i] + ns.iorggrid[i] - 1
    
    index = np.where(lines=='!FLAG3')[0][0]
    lines[index+1] = \
    '      INTEGER,PARAMETER :: xPREC_1ST(nSOMCLASS) = (/%s/)'%(','.join(['%i']*ns.ngrids))%tuple(PREC_1ST_ALL+1)
    lines[index+2] = \
    '      INTEGER,PARAMETER :: xGRID_1ST(nSOMCLASS) = (/%s/)'%(','.join(['%i']*ns.ngrids))%tuple(GRID_1ST_ALL+1)
    lines[index+3] = \
    '      INTEGER,PARAMETER :: xGRID_LST(nSOMCLASS) = (/%s/)'%(','.join(['%i']*ns.ngrids))%tuple(GRID_LST_ALL+1)
    
    # SIZE OF THE SOM ARRAY:
    index = np.where(lines=='!FLAG4')[0][0]
    lines[index+1] = \
    '      INTEGER,PARAMETER :: nSOM = %i'%(ns.nprecs+sum(ns.iorggrid))
    
    # NUMBER OF SOM PARAMETERS:
    index = np.where(lines=='!FLAG5')[0][0]
    lines[index+1] = \
    '      INTEGER,PARAMETER :: nSOMPARAM = %i'%(ns.ngrids*5)
    
    # GRID MAX CARBON AND OXYGEN NUMBERS:
    index = np.where(lines=='!FLAG6')[0][0]
    lines[index+1] = \
    '      INTEGER,PARAMETER :: GRID_CMAX(nSOMCLASS) = (/%s/)'%(','.join(['%i']*ns.ngrids))%tuple(ns.cmaxgrid)
    lines[index+2] = \
    '      INTEGER,PARAMETER :: GRID_OMAX(nSOMCLASS) = (/%s/)'%(','.join(['%i']*ns.ngrids))%tuple(ns.omaxgrid)
    
    # WRITE TO FILE:
    # =========================================================================
    # WRITE TO MODEL AND
    # DIAGNOSTIC FILES:
    f1 = open('%s/source.fort/modules/mod.som.f90'%path,'w')
    f2 = open('../outputs/diag.mod.som.f90','w')
    
    for line in lines:
        f1.write(line + '\n')
        f2.write(line + '\n')
        
    f1.close()
    f2.close()
