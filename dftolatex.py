##---to import line_table as pandas dataframe, apply detection criteria, and pump out latex tables----##
##----by Ayan-------##
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
'''
input_path = '/Users/acharyya/Work/astro/mageproject/AYAN_Code/' #'/Users/acharyya/Dropbox/MagE_atlas/Contrib/EWs/' #where is the dataframe resides
infile = 'fitted_line_list.txt'
output_path = '/Users/acharyya/Documents/writings/papers/abundance_pap/'
outfile = 'stacked_lineEW'
galaxy = 'magestackbyneb'
'''
input_path = '/Users/acharyya/Documents/writings/papers/abundance_pap/' #'/Users/acharyya/Dropbox/MagE_atlas/Contrib/EWs/' #where is the dataframe resides
infile = 'lineflux.txt'
output_path = '/Users/acharyya/Documents/writings/papers/abundance_pap/'
outfile = 'lineflux_detected'
galaxy = 'rcs0327-E'
line_table = pd.read_table(input_path+infile, delim_whitespace=True, comment="#") #input dataframe file
line_table = line_table[line_table['label'].eq(galaxy)]
thresh = 3

tab=line_table[['line_lab','rest_wave','EWr_fit','EWr_fit_u','EW_signi','EWr_Suplim','f_redcor','f_line_u','f_Suplim_redcor']]  
tab.EW_signi=tab.EW_signi.astype(np.float64) 
tab.EWr_Suplim=tab.EWr_Suplim.astype(np.str) 
tab.f_Suplim_redcor=tab.f_Suplim_redcor.map('{:,.3e}'.format).astype(np.str) 
tab.EWr_fit=tab.EWr_fit.astype(np.str) 
tab.f_redcor=tab.f_redcor.map('{:,.3e}'.format).astype(np.str)
tab.EWr_fit = np.where(tab.EW_signi<thresh,'< '+tab.EWr_Suplim,tab.EWr_fit) 
tab.f_redcor = np.where(tab.EW_signi<thresh,'< '+tab.f_Suplim_redcor,tab.f_redcor) 
tab.EWr_fit_u=tab.EWr_fit_u.astype(np.str) 
tab.f_line_u=tab.f_line_u.astype(np.str) 
tab.EWr_fit_u[tab.EW_signi<thresh] = None 
tab.f_line_u[tab.EW_signi<thresh] = None 
#tab.EWr_Suplim[tab.EW_signi<thresh] = None
tab.EW_signi[tab.EW_signi<thresh] = None 
tab.f_line_u[tab.f_line_u<thresh] = None 
tab=tab[['line_lab','rest_wave','EWr_fit','EWr_fit_u','EW_signi','f_redcor','f_line_u']]  
tab['f_line_u']=tab['f_line_u'].astype(np.float64).map('{:,.3e}'.format)
tab=tab.replace('nan','-') 
tab['EWr_fit_u'].fillna(value='-', inplace = True)
tab['f_line_u'].fillna(value='-', inplace = True)
print tab #
tab.to_latex(output_path+outfile+'.tex', index=None)
tab.to_csv(output_path+outfile+'.txt', sep='\t',mode ='w', index=None)