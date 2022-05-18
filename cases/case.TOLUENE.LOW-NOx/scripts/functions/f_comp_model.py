''' ========================================================================
       THIS FUNCTION COMPILES THE SOM-TOMAS MODEL (FORTRAN AND C++)
========================================================================='''

import os
import numpy  as np
import pandas as pd
import time   as tm

import variables as ns

from functions.f_get_root      import get_root
from functions.f_write_modules import write_modules

def comp_model(model):
    
    # MAKE A COPY OF THE SOURCE FOLDER:
    # ======================================================================
    # GET ROOT DIRECTORY FOR EXPERIMENT:
    root = get_root()
    
    # GET COMPILE PATH:
    path = '%s/models/%s.comp.%s'%(root,model,ns.runname)
    
    # MAKE COPY:
    os.system("cd %s/models/ ; cp -a %s %s.comp.%s"%(root,model,model,ns.runname))
    
    # WRITE THE HEADERFILES:
    # ======================================================================
    write_modules(path,ns)
        
    # COMPILE THE SOM-TOMAS MODEL:
    # ======================================================================
    os.system('cd %s ; make clean >/dev/null; make'%path)
    print('SOM-TOMAS Compiled in %s'%path)
    
    os.system('cp %s/box.exe ../box.exe'%(path))
    os.system('rm -rf %s'%path)
    print('SOM-TOMAS Copied to Current Case.')
