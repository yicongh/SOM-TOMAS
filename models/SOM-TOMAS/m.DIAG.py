
import os,sys

import numpy  as np
import pandas as pd

# INPUT AND OUTPUT FOLDERS:
PATHi = './diag.inputs/'
PATHo = './diag.outputs/'

# CLEAN OUTPUTS:
os.system('cd ./%s/ ; rm *'%PATHo)

# RUN MODEL:
line = \
PATHi + '\n' + PATHo

os.system('echo -e "%s" | ./box.exe'%line)
