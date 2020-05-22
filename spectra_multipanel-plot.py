'''
Based on JRR_Code/mage_multipanel-redo_as_function.py, but modified to not-mage-specific i.e. to work for any generic spectra file with any given linelist
by Ayan, May 2020
'''

import os
HOME = os.getenv('HOME')
import sys
sys.path.append(HOME + '/Work/astro/ayan_codes/mageproject/')
from ayan import mage as amage
from jrr import mage as jmage
from jrr import plot as jplot
import matplotlib.pyplot as plt
mage_mode = "reduction"

# --------------to mo=ake new linelist---------------
def create_linelist(linelist_template, zz_dic, zz_err_dic):
    color_dict = {'EMISSION':'red', 'FINESTR':'purple', 'PHOTOSPHERE':'blue', 'ISM':'green', 'WIND':'cyan'}
    LL = amage.getlist(linelist_template, zz_dic, zz_err_dic)
    LL.rename(columns={'wave':'obswav', 'label':'lab1'}, inplace=True)
    LL['color'] = LL['type'].map(lambda x: color_dict[x])
    return LL

# -----------------to grab the linelist from wherever available-----------
def getlinelist(thisgal, linelist_path, my_linelist=None, backup_linelist=None):
    linelist = jmage.get_linelist_name(thisgal['filename'], linelist_path)
    if os.path.exists(linelist):
        print 'Found and using corresponding linelist', linelist
        LL, z_systemic = jmage.get_linelist(linelist)
    elif os.path.exists(my_linelist):
        print 'Using provided linelist', my_linelist
        LL, z_systemic = jmage.get_linelist(my_linelist)
    else:
        print 'Creating new linelist from', backup_linelist
        zz_dic = {'EMISSION': thisgal['z_neb'], 'FINESTR': thisgal['z_neb'],
                  'PHOTOSPHERE': thisgal['z_stars'] if thisgal['fl_st'] == 0 else thisgal['z_neb'],
                  'ISM': thisgal['z_ISM'], 'WIND': thisgal['z_ISM']}
        zz_err_dic = {'EMISSION': thisgal['sig_neb'] if thisgal['fl_neb'] == 0 else thisgal['sig_ISM'],
                      'FINESTR': thisgal['sig_neb'] if thisgal['fl_neb'] == 0 else thisgal['sig_ISM'],
                      'PHOTOSPHERE': thisgal['sig_st'] if thisgal['fl_st'] == 0 else thisgal['sig_neb'],
                      'ISM': thisgal['sig_ISM'], 'WIND': thisgal['sig_ISM']}
        LL, z_systemic = create_linelist(backup_linelist, zz_dic, zz_err_dic), thisgal['z_syst']
    return LL, z_systemic

# --------modify the following lines as required-----------
my_spectra_list = HOME + '/Desktop/mage_plot/Spectra/esi-spectra-filenames-redshifts.txt'
my_linelist = HOME + '/Dropbox/MagE_atlas/Linelists/rcs0327-knotE-lores-comb1'
output_path = HOME + '/Desktop/mage_plot/Spectra/'
my_labels = ['s1723_center_a_esi'] # include as many galaxies as you want from my_spectra_list file
ylim = (0, 10)

# ------------main code: do not modify-----------------------------------------------
linelist_path = HOME + '/Dropbox/MagE_atlas/Linelists/'
backup_linelist = HOME + '/Dropbox/MagE_atlas/Contrib/EWs/linelists/labframe.shortlinelist'
plt.close('all')

specs = jmage.wrap_getlist(mage_mode, which_list="labels", optional_file=my_spectra_list, labels=my_labels, MWdr=False)
for index in range(len(specs)):
    thisgal = specs.iloc[index]
    sp = amage.open_esi_spectrum(thisgal['filename'], getclean=False)
    LL, z_systemic = getlinelist(thisgal, linelist_path, my_linelist=my_linelist, backup_linelist=backup_linelist)
    amage.fit_autocont(sp, z_systemic, linelist_path, thisgal['filename'], boxcar=1001)
    the_dfs = [sp]
    the_zzs = [z_systemic]
    the_pdf = output_path + 'multipanel_' + thisgal['short_label'] + '.pdf'
    zzstr = ','.join(['{0:0.8f}'.format(i) for i in the_zzs])
    jplot.echelle_spectrum(the_dfs, the_zzs, LL, outfile=the_pdf, plot_cont=True, norm_by_cont=True, apply_bad=False, title=(thisgal['short_label']+" "+zzstr), topfid=(1.1,1), ylim=ylim)#, waverange=(1000,3000))
    plt.show(block=False)
