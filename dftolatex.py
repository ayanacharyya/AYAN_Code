##---to import line_table as pandas dataframe, apply detection criteria, and pump out latex tables----##
##----by Ayan-------##
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
import argparse as ap
import os
HOME = os.getenv('HOME')+'/'
#-----------------------------------------------
parser = ap.ArgumentParser(description="Mage spectra fitting tool")
parser.add_argument("--infile")
parser.add_argument("--outfile")
parser.add_argument("--shortlabel")
parser.add_argument("--thresh")
args, leftovers = parser.parse_known_args()
if args.infile is not None:
    infile = args.infile
else:
    infile = HOME+'Documents/writings/papers/abundance_pap/lineflux_restUV.txt' #ffull-path-name of input file
if args.outfile is not None:
    outfile = args.outfile
else:
    outfile = infile[:-4]+'_detected'
if args.shortlabel is not None:
    shortlabel = args.shortlabel
else:
    shortlabel = 'rcs0327-E'
if args.thresh is not None:
    thresh = float(args.thresh)
else:
    thresh = 3.
#-----------------------------------------------
line_table = pd.read_table(infile, delim_whitespace=True, comment="#") #input dataframe file
line_table = line_table[line_table['label'].eq(shortlabel)]

tab=line_table[['line_lab','rest_wave','EWr_fit','EWr_fit_u','EW_signi','EWr_Suplim','f_redcor','f_line_u','f_Suplim_redcor']]  
tab.EW_signi=tab.EW_signi.astype(np.float64) 
tab.EWr_Suplim=tab.EWr_Suplim.astype(np.str) 
tab.f_Suplim_redcor=tab.f_Suplim_redcor.map('{:,.3e}'.format).astype(np.str) 
tab.EWr_fit=tab.EWr_fit.astype(np.str) 
tab.f_redcor=tab.f_redcor.map('{:,.3e}'.format).astype(np.str)
tab.EWr_fit = np.where(tab.EW_signi<thresh,'<'+tab.EWr_Suplim,tab.EWr_fit) 
tab.f_redcor = np.where(tab.EW_signi<thresh,'<'+tab.f_Suplim_redcor,tab.f_redcor) 
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
tab.to_latex(outfile+'.tex', index=None)
tab.to_csv(outfile+'.txt', sep='\t',mode ='w', index=None)
print 'Written files '+ outfile+'.tex'
print 'and '+ outfile+'.txt'
print 'Finished!'