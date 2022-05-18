'''===========================================================
THIS FUNCTION SAVES THE RESULTS FROM SIMULATIONS
==========================================================='''

import pandas as pd

def save_data(outputs,mdata):

    # EXTRACT SIMULATION RESULTS:
    # ========================================================
    time  = outputs['time']
    SOA   = outputs['SOA']
    O2C   = outputs['O2C']
    Yield = outputs['Yield']
    size  = outputs['bin_size']
    
    # NUMBER DISTRIBUTION AND
    # NUMBER DENSITY DISTRIBUTION:
    df_aenum      = outputs['df_aenum']
    df_aenum_norm = outputs['df_aenum_norm']

    # SAVE DATA:
    # ========================================================
    df1 = pd.DataFrame({'time':time,'SOA':SOA,'O2C':O2C,'Yield':Yield})
    df2 = pd.DataFrame({'time':mdata.t_soa,'SOA':mdata.y_soa})
    #df3 = pd.DataFrame({'size':size})
    #df4 = pd.DataFrame({'size':mdata.s_dist})
    #df5 = df_aenum_norm
    #df6 = mdata.y_dist
    #df7 = pd.DataFrame({'time':mdata.t_o2c,'O2C':mdata.y_o2c})
    
    with pd.ExcelWriter('outputs/model_out.xlsx') as w:
        df1.to_excel(w,sheet_name='SOA_mod')
        df2.to_excel(w,sheet_name='SOA_mea')
        #df3.to_excel(w,sheet_name='size_mod')
        #df4.to_excel(w,sheet_name='size_mea')
        #df5.to_excel(w,sheet_name='dist_mod')
        #df6.to_excel(w,sheet_name='dist_mea')
        #df7.to_excel(w,sheet_name='O2C_mea')
        w.close()
