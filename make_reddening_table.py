# python routine to compute reddening using several Balmer lines for S1723 grism spectra
# just to find out which extinction measurements are reliable and the mean E(B-V) that should be used for S1723
# by Ayan, Sep 2017

import numpy as np
from matplotlib import pyplot as plt
from uncertainties import unumpy as unp
import pandas as pd
import os

HOME = os.getenv('HOME')
import sys

sys.path.append('../')
import ayan.mage as m
import ayan.splot_util as u
import calc_reddening as cr
from uncertainties import ufloat
import subprocess

# ------------function to get weighted average in pandas groups----------
def wavg(group, avg_name, uncert_name):
    """ http://stackoverflow.com/questions/10951341/pandas-dataframe-aggregate-function-using-multiple-columns
    In rare instance, we may not have weights, so just return the mean. Customize this if your business case
    should return otherwise.
    """
    d = group[avg_name]
    u = group[uncert_name]
    try:
        return ufloat((d / u**2).sum() / (1./u**2).sum(), np.sqrt(1./(1./u**2).sum()))
    except ZeroDivisionError:
        return ufloat(d.mean(), 0.)
# ----------------------------------------------------------------------

input_path = HOME + '/Dropbox/Grism_S1723/WFC3_fit_1Dspec/1Dsum/'
fout = HOME + '/Dropbox/Grism_S1723/Latex/Tabs/s1723_red_table.tex'
filenames = ['bothroll', 'roll308', 'roll139']
edf = pd.DataFrame()

labels = np.array(['Halpha', 'Hbeta', 'Hgamma', 'Hdelta', 'Hepsilon', 'Hzeta', 'Heta'])
wave = np.array([6564.61, 4862.683, 4341.684, 4102.892, 3971.195, 3890.158, 3836.479])  # rest-frame vaccuum wavelengths of Balemr lines
k = m.getfullkappa(wave)
ratio_labels = [r'H$\alpha$/H$\beta$', r'H$\alpha$/H$\beta$', r'H$\beta$/H$\gamma$', r'H$\beta$/H$\delta$']  # , r'H$\beta$/H$\epsilon$', r'H$\beta$/H$\zeta$', r'H$\zeta$/H$\epsilon$']
spec_labels = ['G141', 'G141/G102', 'G102', 'G102']  # , 'G102', 'G102/ESI', 'ESI']
delta_k = np.array([k[1] - k[0], k[1] - k[0], -k[1] + k[2], -k[1] + k[3]])  # , -k[1]+k[4], -k[1]+k[5], -k[5]+k[6]])
ratio_theoretical = np.array([2.8785 / 1, 2.8785 / 1, 1 / 0.46756, 1 / .25825])  # , 1/.15856, 1/0.10472, 0.10472/.072916]) #from /Users/acharyya/Desktop/Lisa_UV_diag/P_spherical/spherical/sp_P60_a05modelfiles/Q800/spec0003.csv
# MAPPINGS V model for lpok=6, logq=8, Z=8.23
g102_hlist = ['Hbeta', 'Hgamma', 'Hdelta', 'Hepsilon']
g141_hlist = ['Halpha', 'Hbeta']

plt.close('all')
fig, ax = plt.subplots()
fig.subplots_adjust(hspace=0.1, wspace=0.03, top=0.97, bottom=0.2, left=0.1, right=0.98)
col_ar = ['r', 'g', 'b', 'orange']

for (index,item) in enumerate(filenames): # for each roll
    wfc_g102_file = input_path + 'sgas1723_1Dsum_'+item+'_G102_wcontMWdr_meth2.fitdf'
    wfc_g102_file = cr.deredden(wfc_g102_file, 0.3, 0.02, niter=int(1e5), constant=1e-17, dered=True, change_ID=True, change_errors=True,
                        SNR_thresh=0.014, readcsv=True)  # to deredden WFC3 spectra and re-write the file

    wfc_g141_file = input_path + 'sgas1723_1Dsum_'+item+'_G141_wcontMWdr_meth2.fitdf'
    wfc_g141_file = cr.deredden(wfc_g141_file, 0.3, 0.02, niter=int(1e5), constant=1e-17, dered=True, change_ID=True, change_errors=True,
                        SNR_thresh=0.014, readcsv=True)  # to deredden WFC3 spectra and re-write the file

    g102 = pd.read_table(wfc_g102_file, delim_whitespace=True, comment='#')
    g141 = pd.read_table(wfc_g141_file, delim_whitespace=True, comment='#')
    H_g102 = unp.uarray(cr.getf(g102, g102_hlist), cr.gete(g102, g102_hlist))
    H_g141 = unp.uarray(cr.getf(g141, g141_hlist), cr.gete(g141, g141_hlist))

    ratio_obs = [H_g141[0] / H_g141[1], H_g141[0] / H_g102[0], H_g102[0] / H_g102[1],
                 H_g102[0] / H_g102[2]]  # , H_g102[0]/H_g102[3], H_g102[0]/(H_esi[0]*esi_to_mmt_scale), H_esi[0]/H_esi[1]]
    Eb_v = (2.5 / delta_k) * unp.log10(ratio_obs / ratio_theoretical)
    print '\nMean E(B-V) =', np.mean(Eb_v)
    edf = pd.concat([edf, pd.DataFrame(np.transpose(np.vstack([ratio_labels, spec_labels, ['%.3F' % x.n for x in Eb_v], ['%.2F' % x.s for x in Eb_v], ['%.2F' % x.n for x in ratio_obs], ['%.2F' % x.s for x in ratio_obs]])),
        columns=['Line ratio', 'Grism(s)', 'E(B-V)', r'E(B-V)$_u$', 'Observed ratio', 'Observed ratio uncert'])])
    # ------to plot the E(B-V) values for visual comparison-----
    for ii in range(len(Eb_v)):
        ax.scatter(index + ii*0.1, Eb_v[ii].n, s=40, color=col_ar[ii], lw=0, label=ratio_labels[ii]+' '+spec_labels[ii] if not index else None)
        ax.errorbar(index + ii*0.1, Eb_v[ii].n, yerr=Eb_v[ii].s, color=col_ar[ii])

# --------to plot weighted means of each extinction---------------
for column in edf.columns[2:]: edf[column] = edf[column].astype(np.float64)
wavg_arr = edf.groupby(['Line ratio', 'Grism(s)'], sort=False).apply(wavg, 'E(B-V)', 'E(B-V)$_u$').values
for (ii,mean) in enumerate(wavg_arr):
    plt.axhline(mean.n, color=col_ar[ii], ls='dashed', label='Weighted means' if not ii else None)
print wavg_arr

fs = 15 # plot label fontsize
plt.legend(loc='lower left', fontsize=fs)
ax.set_xticks(np.arange(len(filenames)))
ax.set_xticklabels(filenames, rotation=50, fontsize=fs, ha='right', va='top')
plt.ylim(-0.7,0.8)
plt.ylabel('E(B-V)', fontsize=fs)
plt.xlabel('Rolls', fontsize=fs)
plt.show(block=False)
fig.savefig(os.path.splitext(fout)[0]+'.eps')

edf.to_latex(fout, index=None, escape=False)
for i in range(len(filenames)): u.insert_line_in_file('\hline\n\multicolumn{4}{c}{'+filenames[i]+'} \\\ \n\hline \n', 4 + i*(3+len(ratio_labels)), fout)
u.insert_line_in_file(r'%Generated by source code '+__file__+'\n', 0, fout)
with open(fout, 'a') as f:
    f.write(r'\caption{E(B-V) measurements using different Balmer line pairs (column "Line ratio"). The spectrographs corresponding to \
the line pairs are listed in "Grism(s)". The reddening values and corresponding uncertainties are quoted in the "E(B-V)"\
and "E(B-V)$_u$" columns respectively. Intrinsic ratios used to calculate the reddening are from an isobaric, MAPPINGS-V photoionisation \
model assuming spherical geometry and $\log{(\mathrm P/k)}=6$, $\log{(\mathrm q)}=8$ and $\log{(\mathrm O/H)} + 12=8.23$}')

print '\nList of E(B-V) values\n', edf
print 'Written file', fout
