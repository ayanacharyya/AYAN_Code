#-----utility routine to extract OIII2320 and OIII1666 fluxes from all spectra; for David Nicholl's paper----------
#-----by Ayan, Jan 2017------------
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
import sys
sys.path.append('../')
import ayan.mage as m
#----------------------------------------------------------
path = '/Users/acharyya/Dropbox/MagE_atlas/Contrib/EWs/'
fn_arr = ['allspec_fitted_emission_linelist.txt', 'all_byneb_stack_fitted_emission_linelist.txt', 'allesi_fitted_emission_linelist.txt']
line_arr = ['[OIII]2320', 'OIII]1666']
fout = 'o3_linelist.csv'
f_SNR_thresh = 1.
f_signi_thresh = 3.
EW_signi_thresh = 3.
#----------------------------------------------------------
df_final = pd.DataFrame()
for line in line_arr:
    df = pd.DataFrame()
    for fn in fn_arr:
        line_table = pd.read_table(path+fn, delim_whitespace=True, comment="#")
        o3_table = line_table[line_table.line_lab == line]
        #o3_table['f_SNR'] = np.abs(o3_table['f_line'])/o3_table['f_line_u']
        #o3_table = o3_table[o3_table['f_SNR'] > f_SNR_thresh] #eliminating those with flux uncertainty greater than flux
        #o3_table = o3_table[o3_table['f_signi'] > f_signi_thresh] #eliminating those below detection limit i.e. f_signi_thresh
        o3_table = o3_table[o3_table['EW_signi'] > EW_signi_thresh] #eliminating those below detection limit i.e. EW_signi_thresh
        #o3_table['f_redcor'],o3_table['f_redcor_u'] = m.extinct(o3_table.rest_wave, o3_table.f_line, o3_table.f_line_u) #calling de-redenning routine
        #o3_table['f_Suplim_redcor'],dummy = m.extinct(o3_table.rest_wave, o3_table.f_Suplim, o3_table.f_line_u)
        o3_short_table = o3_table[['label','line_lab','rest_wave','f_line','f_line_u','f_Suplim','f_signi','EWr_fit', 'EWr_fit_u','EW_signi']]
        df = df.append(o3_short_table, ignore_index=True)
    df.f_signi = df.f_signi.astype(np.float64)
    df.sort_values('f_signi',ascending=False, inplace=True) #sorting
    df_final = df_final.append(df, ignore_index=True)
df_final.to_csv(path+fout, sep='\t',mode ='a', index=None)
print df_final
print '\n'

