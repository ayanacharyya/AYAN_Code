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
#---------to see if detected-----------
def detect(table, line):
    if table[table['line_lab'].eq(line)].EW_signi.values[0] > ewthresh and (table[table['line_lab'].eq(line)].f_line.values[0]/table[table['line_lab'].eq(line)].f_line_u.values[0] > 1.):
        return True
    else:
        return False
#--------------------------------------
output_path = '/Users/acharyya/Desktop/mage_plot/' #where is the output to be put
input_path = '/Users/acharyya/Dropbox/MagE_atlas/Contrib/EWs/' #where is the dataframe resides
fulltable = pd.read_table(input_path+'allspec_fitted_emission_linelist.txt', delim_whitespace=True, comment="#") #input dataframe file
lines = pd.read_table('labframe.shortlinelist_emission', delim_whitespace=True, comment="#") #input list of lines fitted
excludelabel = []#'S0004-0103']#'rcs0327-E']#'S0957+0509']#,'S1050+0017',]# 'S1429+1202'] #which spectra to exclude
fulltable = fulltable[~fulltable['label'].isin(excludelabel)] #
labels = np.unique(fulltable['label'])
colors = 'rcmykbgrmcykbgrmcykbg' #color schemes
#------------------------------------------------
#lines_num = ['OIII]1660', 'OIII1666']
lines_num = ['OIII2320', '[OIII]2331']
lines_den = ['OII2470mid']

#lines_num = ['CIII977']
#lines_num = ['CIII]1906', 'CIII]1908']#, 'CIII2323', 'CII2325b', 'CII2325c', 'CII2325d', 'CII2328']
#lines_den = ['CIII2323', 'CII2325b', 'CII2325c', 'CII2325d', 'CII2328']
#lines_den2 = ['CII1335b', 'CII1335c']

#lines_num = ['SiIII1882', 'SiIII1892']
#lines_den = ['SiII2335a', 'SiII2335b']
#lines_den2 = ['SiII1533']

#lines_num = ['NIV]1486']
#lines_num = ['NII1084', 'NII1085']
#lines_num = ['NII]2140']

#lines_den = ['HeII1640']
#------------------------------------------------
#for line in lines.LineID.values:
fig, ax1 = plt.subplots(1,1)
h = np.zeros(len(lines)) #array to store line indices, to plot on the x-axis
g = np.zeros(len(labels)) #array to store galaxy indices, to plot on the x-axis
ewthresh = 4.0

for ii in range(0,len(labels)):
    ew = np.arange(len(lines)+2)*np.nan
    ewu = np.arange(len(lines)+2)*np.nan
    #filtering criteria as follows
    table = fulltable[(~np.isnan(fulltable['f_line'])) & (fulltable['label'].eq(labels[ii]))]# & \
    #(fulltable['EW_signi'] > 4.) & \
    #(fulltable['f_line']/fulltable['f_line_u'] > 1.)]
    
    
    a, avar, b, bvar, c, cvar, d, dvar, z, zu = 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
    #print 'for', labels[ii] #
    for line in np.array(lines_num):
        try:
            if detect(table,line) :
                a += table[table['line_lab'].eq(line)].f_line.values[0]
                avar += table[table['line_lab'].eq(line)].f_line_u.values[0]**2
            #a += table[table['line_lab'].eq(line)].EWr_fit.values[0]
            #avar += table[table['line_lab'].eq(line)].EWr_fit_u.values[0]**2
            #z += table[table['line_lab'].eq(line)].zz.values[0]
            #zu += table[table['line_lab'].eq(line)].zz_u.values[0]
            #print 'added a' #
            c += table[table['line_lab'].eq(line)].f_Suplim.values[0]
        except:
            #print 'passing' #
            pass
    
    for line in lines_den:
        try:
            if detect(table,line) :
                b += table[table['line_lab'].eq(line)].f_line.values[0]
                bvar += table[table['line_lab'].eq(line)].f_line_u.values[0]**2
            #print 'added b' #
            d += table[table['line_lab'].eq(line)].f_Suplim.values[0]
        except:
            pass
    '''
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
    '''
    print ii, labels[ii], a, b, c, d #
    if a > 0 and b > 0:
        err = jrr.util.sigma_adivb(a, np.sqrt(avar), b, np.sqrt(bvar))
        #err2 = jrr.util.sigma_adivb(c, np.sqrt(cvar), d, np.sqrt(dvar))
        #pl=ax1.errorbar(np.log10(np.divide(c,d)), np.log10(np.divide(a,b)), fmt='o', lw=0.5, xerr=np.log10(err2), \
        #yerr=np.log10(err))
        pl=ax1.errorbar(ii+1, a/b, fmt='o', lw=0.5, yerr=err)
    
    elif a > 0 and d > 0:
        err = jrr.util.sigma_adivb(a, np.sqrt(avar), d, np.sqrt(0.))
        pl=ax1.errorbar(ii+1, a/d, capsize=5, elinewidth=1, yerr=err, lolims=True)
    elif b > 0 and c > 0:
        err = jrr.util.sigma_adivb(c, np.sqrt(0.), b, np.sqrt(bvar))
        pl=ax1.errorbar(ii+1, c/b, capsize=5, elinewidth=1, yerr=err, uplims=True)
    elif c > 0 and d > 0:
        err = jrr.util.sigma_adivb(c, np.sqrt(0.), d, np.sqrt(0.))
        pl=ax1.errorbar(ii+1, c/d, fmt='*', lw=2, yerr=np.log10(err))
    '''
    if a < 0:
        #pl=ax1.errorbar(z/len(lines_num), -1*a, fmt='o', lw=0.5, xerr=zu/len(lines_num), yerr=np.sqrt(avar))
        pl=ax1.errorbar(z, -1*a, fmt='o', lw=0.5, xerr=zu, yerr=np.sqrt(avar))
    
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

    #plt.xlabel('zz')
    #plt.xlabel(lines_num2+['/']+lines_den2)
    #plt.xlabel(lines_den)
    #plt.ylabel(lines_num)
    plt.ylabel('('+' + '.join(lines_num)+')/('+' + '.join(lines_den)+')')
    plt.title('Flux ratios')
    #plt.title('EW ratios')
    '''
    #plt.xlabel('Rest wavelength (A)')
    plt.ylabel('Measured EW (A)')
    #plt.ylabel('Number of galaxies line is detected in')
    '''
    plt.yscale('log')
    #plt.draw()
    #plt.pause(1)

'''
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(lines.restwave)
ax2.set_xticklabels(lines.LineID, rotation = 45, ha='left', fontsize='small')
'''
fig.savefig(output_path+'plot'+str(n)+'.png')
ax1.set_xticks(np.arange(len(labels))+1)
ax1.set_xticklabels(labels, rotation = 45, ha='right', fontsize='small')
plt.xlim(0,len(labels)+1)
plt.show(block=False)