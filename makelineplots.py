'''
Makes different line ratio plots, using information from the dataframe resulting from line fitting using EW_fitter.py
For now its not automated, hence you have to manually comment out/modify this each time to get a different line ratio
plot. Still a the crude version.
Started July 2016, Ayan Acharyya
'''
import sys
sys.path.append('../')
import jrr
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import argparse as ap

parser = ap.ArgumentParser(description="Mage spectra fitting tool")
parser.add_argument('--n')
args, leftovers = parser.parse_known_args()
if args.n is not None:
    n = args.n
else:
    n=50
#---------to see if detected-----------
def detect(table, line):
    if table[table['line_lab'].eq(line)].EW_signi.values[0] > ew_thresh and (table[table['line_lab'].eq(line)].f_line.values[0]/table[table['line_lab'].eq(line)].f_line_u.values[0] > fSNR_thresh):
        return True
    else:
        return False
#--------------------------------------
output_path = '/Users/acharyya/Desktop/mage_plot/' #where is the output to be put
input_path = '/Users/acharyya/Dropbox/MagE_atlas/Contrib/EWs/' #where is the dataframe resides
table1 = pd.read_table(input_path+'allspec_fitted_emission_linelist.txt', delim_whitespace=True, comment="#") #input dataframe file
table2 = pd.read_table(input_path+'all_byneb_stack_fitted_emission_linelist.txt', delim_whitespace=True, comment="#") #input dataframe file
table3 = pd.read_table(input_path+'allesi_fitted_emission_linelist.txt', delim_whitespace=True, comment="#") #input dataframe file
frames = [table1, table2, table3]
fulltable = pd.concat(frames)
lines = pd.read_table('labframe.shortlinelist_emission', delim_whitespace=True, comment="#") #input list of lines fitted
excludelabel = ['magestack_byneb_highZ', 'magestack_byneb_lowZ', 'magestack_byneb_midage8to16Myr', \
'magestack_byneb_oldgt16Myr', 'magestack_byneb_younglt8Myr']#'S0004-0103']#'rcs0327-E']#'S0957+0509']#,'S1050+0017',]# 'S1429+1202'] #which spectra to exclude
fulltable = fulltable[~fulltable['label'].isin(excludelabel)] #
labels = np.unique(fulltable['label'])
colors = 'rcmybgrmcybgrmcybg' #color schemes
#------------------------------------------------
#lines_num = ['OIII]1660', 'OIII1666']
lines_num = ['OIII2320', '[OIII]2331']
#lines_den = ['OII2470mid']

#lines_num = ['CIII977']
#lines_num = ['CIII]1906', 'CIII]1908']#, 'CIII2323', 'CII2325b', 'CII2325c', 'CII2325d', 'CII2328']
lines_den = ['CII2323', 'CII2325b', 'CII2325c', 'CII2325d', 'CII2328']
#lines_den = ['CII1335b', 'CII1335c']

#lines_num = ['SiIII1882', 'SiIII1892']
#lines_den = ['SiII2335a', 'SiII2335b']
#lines_den = ['SiII1533']

#lines_num = ['NIV]1486']
#lines_den = ['NII1084', 'NII1085']
#lines_den = ['NII]2140']

#lines_den = ['HeII1640']
#------------------------------------------------
#for line in lines.LineID.values:
fig, ax1 = plt.subplots(1,1)
fig.subplots_adjust(hspace=0.7, top=0.95, bottom=0.1, left=0.1, right=0.8)
h = np.zeros(len(lines)) #array to store line indices, to plot on the x-axis
g = np.zeros(len(labels)) #array to store galaxy indices, to plot on the x-axis
ew_thresh = 4.0
fSNR_thresh = 1.0
quantity = 'f_line'#'EWr_fit' #
quantity_u = quantity + '_u'
lim = 'f_Suplim' #EWr_Suplim' #
#------------------------------------------------

for ii in range(0,len(labels)):
    color = colors[ii%len(colors)]
    if 'stack' in labels[ii]: color ='k' #black color for stack
    ew = np.arange(len(lines)+2)*np.nan
    ewu = np.arange(len(lines)+2)*np.nan
    #filtering criteria as follows
    table = fulltable[(~np.isnan(fulltable['f_line'])) & (fulltable['label'].eq(labels[ii]))]
    
    
    a, avar, b, bvar, c, cvar, d, dvar, z, zu = 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
    for line in np.array(lines_num):
        try:
            if detect(table,line) :
                a += table[table['line_lab'].eq(line)][quantity].values[0]
                avar += table[table['line_lab'].eq(line)][quantity_u].values[0]**2
            if lim in locals(): c += table[table['line_lab'].eq(line)][lim].values[0]
        except:
            #print 'passing' #
            pass
    
    for line in lines_den:
        try:
            if detect(table,line) :
                b += table[table['line_lab'].eq(line)][quantity].values[0]
                bvar += table[table['line_lab'].eq(line)][quantity_u].values[0]**2
            if lim in locals(): d += table[table['line_lab'].eq(line)][lim].values[0]
        except:
            pass
    '''
    #----------------------for A/B vs C/D plots------------------
    for line in lines_num2:
        try:
            c += table[table['line_lab'].eq(line)][quantity].values[0]
            cvar += table[table['line_lab'].eq(line)][quantity_u].values[0]**2
            #print 'added c' #
        except:
            pass

    for line in lines_den2:
        try:
            d += table[table['line_lab'].eq(line)][quantity].values[0]
            dvar += table[table['line_lab'].eq(line)][quantity_u].values[0]**2
            #print 'added d' #
        except:
            pass
    #--------------------------------------------------------------
    '''
    print ii, labels[ii], a, b, c, d #
    '''
    #-----------------For A/B vs object index plot---------------------------------------------
    if a > 0 and b > 0:
        err = jrr.util.sigma_adivb(a, np.sqrt(avar), b, np.sqrt(bvar))
        #err2 = jrr.util.sigma_adivb(c, np.sqrt(cvar), d, np.sqrt(dvar))
        #pl=ax1.errorbar(np.log10(np.divide(c,d)), np.log10(np.divide(a,b)), fmt='o', lw=0.5, xerr=np.log10(err2), \
        #yerr=np.log10(err))
        pl=ax1.errorbar(ii+1, a/b, fmt='o', lw=0.5, yerr=err, color=color)
    
    elif a > 0 and d > 0:
        err = jrr.util.sigma_adivb(a, np.sqrt(avar), d, np.sqrt(0.))
        pl=ax1.errorbar(ii+1, a/d, capsize=5, elinewidth=1, yerr=err, lolims=True, color=color)
    elif b > 0 and c > 0:
        err = jrr.util.sigma_adivb(c, np.sqrt(0.), b, np.sqrt(bvar))
        pl=ax1.errorbar(ii+1, c/b, capsize=5, elinewidth=1, yerr=err, uplims=True, color=color)
    elif c > 0 and d > 0:
        err = jrr.util.sigma_adivb(c, np.sqrt(0.), d, np.sqrt(0.))
        pl=ax1.errorbar(ii+1, c/d, fmt='*', lw=2, yerr=np.log10(err), color=color)
    #--------------------------------------------------------------
    '''
    '''
    #----------------------for A/B vs C/D plots------------------
        err2 = jrr.util.sigma_adivb(c, np.sqrt(cvar), d, np.sqrt(dvar))
        pl=ax1.errorbar(np.log10(np.divide(c,d)), np.log10(np.divide(a,b)), fmt='o', lw=0.5, xerr=np.log10(err2), \
        yerr=np.log10(err))
    #--------------------------------------------------------------
    '''
    #-----------------For A vs B plot---------------------------------------------
    if np.abs(a) > 0 and np.abs(b) > 0:
        pl=ax1.errorbar(b, a, fmt='o', lw=0.5, yerr=avar, xerr=bvar, color=color)
        height = ax1.get_ylim()[1] - ax1.get_ylim()[0]
        ax1.annotate(labels[ii], xy=(b,a), xytext=(b,a-0.05*height),color=color,fontsize=10)
    elif np.abs(a) > 0 and np.abs(d) > 0:
        pl=ax1.errorbar(d, a, capsize=5, elinewidth=1, yerr=avar, lolims=True, color=color)
    elif np.abs(b) > 0 and np.abs(c) > 0:
        pl=ax1.errorbar(b, c, capsize=5, elinewidth=1, xerr=bvar, uplims=True, color=color)
    elif np.abs(c) > 0 and np.abs(d) > 0:
        pl=ax1.errorbar(d, c, fmt='*', lw=2, color=color)
    #--------------------------------------------------------------
    '''
    #------------------For A (line) vs object index plot--------------------------------------------
    if a < 0:
        #pl=ax1.errorbar(z/len(lines_num), -1*a, fmt='o', lw=0.5, xerr=zu/len(lines_num), yerr=np.sqrt(avar))
        pl=ax1.errorbar(z, -1*a, fmt='o', lw=0.5, xerr=zu, yerr=np.sqrt(avar))
    
    for jj, line in enumerate(lines.LineID):
        if table['line_lab'].isin([line]).any():
            #h[jj] += 1
            ew[jj+1] = table[table['line_lab'].isin([line])].EWr_fit.values[0]
            ewu[jj+1] = table[table['line_lab'].isin([line])].EWr_fit_u.values[0]

    try:
        #plt.bar(range(len(lines)), h, lw=0, align = 'center', color=color])
        #pl=plt.errorbar(range(len(ew)), ew, fmt='o', lw=0.5, yerr=ewu)
        pl=plt.errorbar(range(len(table)), table.EWr_fit.values, fmt='o', lw=0.5, yerr=table.EWr_fit_u.values)
        plt.xticks(range(len(lines)+2),np.concatenate(([' '],lines.LineID.values,[' '])), rotation = 90, fontsize='small')
        print pl[0].get_color(), table.label.values[0], len(table)    
    except:
        pass
    #--------------------------------------------------------------
    '''

    t = 'Comparing_'+quantity
    ylab = '('+'+'.join(lines_num)+')'#'('+'+'.join(lines_num)+')_by_('+'+'.join(lines_den)+')'
    xlab = '('+'+'.join(lines_den)+')'#'obj'
    
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    plt.title(t)
    if 'f' in quantity:
        plt.yscale('log')
        plt.xscale('log')
    #plt.pause(0.1)

'''
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(lines.restwave)
ax2.set_xticklabels(lines.LineID, rotation = 45, ha='left', fontsize='small')
'''
'''
#----------For anything vs object index plot----------------------------------------------------
ax1.set_xticks(np.arange(len(labels))+1)
ax1.set_xticklabels(labels, rotation = 45, ha='right', fontsize='small')
plt.xlim(0,len(labels)+1)
#--------------------------------------------------------------
'''
out = output_path+t+'_'+ylab+'_vs_'+xlab+'.png'
fig.savefig(out)
print 'Saved as', out
plt.show(block=False)