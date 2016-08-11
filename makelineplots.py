'''
Makes different line ratio plots, using information from the dataframe resulting from line fitting using EW_fitter.py
For now its not automated, hence you have to manually comment out/modify this each time to get a different line ratio
plot. Still a the crude version.
Started July 2016, Ayan Acharyya
Usage options:
--n N ; N = number to be appended at the end of the output .png file to avoid overwriting of files (yes thats a stupid
        way of doing it)
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

plt.close('all')
fig = plt.figure()
output_path = '/Users/acharyya/Desktop/mage_plot/' #where is the output to be put
fulltable = pd.read_table('fitted_emission_list_allspec.txt', delim_whitespace=True, comment="#") #input dataframe file
lines = pd.read_table('labframe.shortlinelist_emission', delim_whitespace=True, comment="#") #input list of lines fitted
excludelabel = ['S0957+0509']#,'S1050+0017',]# 'S1429+1202'] #which spectra to exclude
fulltable = fulltable[~fulltable['label'].isin(excludelabel)] #
labels = np.unique(fulltable['label'])
colors = 'rcmykbgrmcykbgrmcykbg' #color schemes
#------------------------------------------------
#lines_num = ['OIII]1660', 'OIII1666']
#lines_num = ['OIII2320', '[OIII]2331']
lines_den = ['OII2470mid']

#lines_num = ['CIII977']
lines_num2 = ['CIII]1906', 'CIII]1908']#, 'CIII2323', 'CII2325b', 'CII2325c', 'CII2325d', 'CII2328']
#lines_den = ['CIII2323', 'CII2325b', 'CII2325c', 'CII2325d', 'CII2328']
lines_den2 = ['CII1335b', 'CII1335c']

#lines_num2 = ['SiIII1882', 'SiIII1892']
#lines_den2 = ['SiII2335a', 'SiII2335b']
#lines_den2 = ['SiII1533']

#lines_num = ['NIV]1486']
#lines_num = ['NII1084', 'NII1085']
lines_num = ['NII]2140']

#lines_den = ['HeII1640']
#------------------------------------------------
ax1 = fig.add_subplot(111)
h = np.zeros(len(lines))
for ii in range(0,len(labels)):
    ew = np.arange(len(lines)+2)*np.nan
    ewu = np.arange(len(lines)+2)*np.nan
    #filtering criteria as follows
    table = fulltable[(~np.isnan(fulltable['f_line'])) & (fulltable['label'].eq(labels[ii])) & \
    (~fulltable['type'].eq('ISM')) & \
    (fulltable['EW_significance'] > 4.)]# & \
    #(fulltable['MAD_significance'] > 2.) & \
    #(fulltable['f_line']/fulltable['f_line_u'] > 2.)]
    table = table[table['EWr_fit_u']/np.abs(table['EWr_fit']) < 3.]
    
    a, avar, b, bvar, c, cvar, d, dvar = 0.,0.,0.,0.,0.,0.,0.,0.
    #print 'for', labels[ii] #
    for line in lines_num:
        try:
            a += table[table['line_lab'].eq(line)].f_line.values[0]
            avar += table[table['line_lab'].eq(line)].f_line_u.values[0]**2
            #print 'added a' #
        except:
            pass
    for line in lines_den:
        try:
            b += table[table['line_lab'].eq(line)].f_line.values[0]
            bvar += table[table['line_lab'].eq(line)].f_line_u.values[0]**2
            #print 'added b' #
        except:
            pass
    
    for line in lines_num2:
        try:
            c += table[table['line_lab'].eq(line)].f_line.values[0]
            cvar += table[table['line_lab'].eq(line)].f_line_u.values[0]**2
            #print 'added c' #
        except:
            pass
    
    for line in lines_den2:
        try:
            d += table[table['line_lab'].eq(line)].f_line.values[0]
            dvar += table[table['line_lab'].eq(line)].f_line_u.values[0]**2
            #print 'added d' #
        except:
            pass
    
    print ii, labels[ii], a, b, c, d #
    if a > 0 and c > 0:
        err = jrr.util.sigma_adivb(a, np.sqrt(avar), b, np.sqrt(bvar))
        err2 = jrr.util.sigma_adivb(c, np.sqrt(cvar), d, np.sqrt(dvar))
        pl=ax1.errorbar(np.log10(np.divide(c,d)), np.log10(np.divide(a,b)), fmt='o', lw=0.5, xerr=np.log10(err2), \
        yerr=np.log10(err))
    
    '''
    for jj, line in enumerate(lines.LineID):
        if table['line_lab'].isin([line]).any():
            #h[jj] += 1
            ew[jj+1] = table[table['line_lab'].isin([line])].EWr_fit.values[0]
            ewu[jj+1] = table[table['line_lab'].isin([line])].EWr_fit_u.values[0]
    
    try:
        #plt.bar(range(len(lines)), h, lw=0, align = 'center', color=colors[ii])
        #pl=plt.errorbar(range(len(ew)), ew, fmt='o', lw=0.5, yerr=ewu)
        pl=plt.errorbar(range(len(table)), table.EWr_fit.values, fmt='o', lw=0.5, yerr=table.EWr_fit_u.values)
        plt.xticks(range(len(lines)+2),np.concatenate(([' '],lines.LineID.values,[' '])), rotation = 90, fontsize='small')
        print pl[0].get_color(), table.label.values[0], len(table)    
    except:
        pass
    '''
    
    plt.xlabel(lines_num2+['/']+lines_den2)
    #plt.xlabel(lines_den)
    plt.ylabel(lines_num+['/']+lines_den)
    plt.title('Flux ratios')
    #plt.title('EW ratios')
    '''
    #plt.xlabel('Rest wavelength (A)')
    plt.ylabel('Measured EW (A)')
    #plt.ylabel('Number of galaxies line is detected in')
    '''
    plt.draw()
    #plt.pause(1)

'''
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(lines.restwave)
ax2.set_xticklabels(lines.LineID, rotation = 45, ha='left', fontsize='small')
'''
fig.savefig(output_path+'plot'+str(n)+'.png')
plt.show(block=False)