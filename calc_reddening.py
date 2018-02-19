# python routine to compute E(B-V) using several Balmer lines for S1723
# by Ayan, Sep 2017

import numpy as np
from uncertainties import unumpy as unp
import pandas as pd
import os

HOME = os.getenv('HOME')
import sys

sys.path.append('../')
import ayan.mage as m
import ayan.splot_util as u
import subprocess
import jrr as jrr
import re


# ------------------------------------------------------
def getf(df, lines, fluxcol='integrated_flux', labelcol='ID'):  # function to extract flux values from dataframe
    return [df[df[labelcol] == line][fluxcol].values[0] if len(df[df[labelcol] == line][fluxcol].values) == 1 else
            df[df[labelcol] == line][fluxcol].values for line in lines]


def gete(df, lines, errorcol='uncert_flux', labelcol='ID'):  # function to extract flux uncertainties from dataframe
    return [df[df[labelcol] == line][errorcol].values[0] if len(df[df[labelcol] == line][errorcol].values) == 1 else
            df[df[labelcol] == line][errorcol].values for line in lines]


# ------------------------------------------------------
def getheader(filename, comment='#'):
    header = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(comment):
                header.append(line[1:])
            else:
                break
    return ''.join(header)


# ------------------------------------------------------
def scale_df(df, scale_factor, columns_to_scale=['integrated_flux', 'uncert_flux', 'dered_flux', 'uncert_dered_flux']):
    for column in columns_to_scale:
        if column in df:
            df[column] *= scale_factor
        else:
            continue
    return df


# ------------------------------------------------------
def format_df(df, columns_to_scale=['integrated_flux', 'uncert_flux', 'dered_flux', 'uncert_dered_flux']):
    for column in columns_to_scale:
        if column in df:
            df[column] = df[column].map('{:.3f}'.format)
        else:
            continue
    return df


# ------------------------------------------------------
def deredden(filename, E, E_u, niter=int(1e5), constant=1e-17, dered=True, change_ID=False, change_errors=False,
             SNR_thresh=0.):
    header = getheader(filename)
    df = pd.read_table(filename, delim_whitespace=True, comment='#')

    if change_errors:
        df['SNR'] = df['uncert_flux'] / df['integrated_flux']
        concerned = (df['uncert_flux'] > 0.) & (df['SNR'] < SNR_thresh)
        print 'Before changing uncertainties for\n', df[concerned]
        df.ix[concerned, 'uncert_flux'] = df.ix[
                                              concerned, 'integrated_flux'] * SNR_thresh  # modifying uncertainties by hand
        df['SNR'] = df['uncert_flux'] / df['integrated_flux']
        print 'After changing uncertainties for\n', df[concerned]
        df.drop('SNR', axis=1, inplace=True)
        filename = filename[:-4] + '_umod.txt'
        header += '\
    "uncert_flux" has been modified as follows:\n\
    If any line had a "integrated_flux"/"uncert_flux" < ' + str(
            SNR_thresh) + ', "uncert_flux" was replaced such that this ratio becomes equal to ' + str(SNR_thresh) + '\n\
'
    if change_ID:
        df['ID'] = ['OII3727+9', 'NeIII3869', 'Hzeta+HeI', 'NeIII3968', 'Hepsilon', 'Hdelta', 'Hgamma', 'OIII4363',
                    'Hbeta', 'Hbeta', 'OIII4959', 'OIII5007', 'HeI5875', 'SIII6312', 'Halpha', 'SII6717', 'ArIII7136',
                    'NII6549', 'NII6584']
        filename = filename[:-4] + '_lmod.txt'
        header += '\
    "ID" has been modified as follows:\n\
    Line labels have been changed to match a generalised format.\n\
'
    if dered:
        df.replace('3727+3729', '3728.5',
                   inplace=True)  # replacing doublet rest wavelength to its mean value for dereddening
        df.replace('6718+6732', '6725.',
                   inplace=True)  # replacing doublet rest wavelength to its mean value for dereddening
        df.replace('6302,6313', '6307.5',
                   inplace=True)  # replacing doublet rest wavelength to its mean value for dereddening

        df['rest_wave'] = df['rest_wave'].astype(np.float64)
        dered_flux, uncert_dered_flux = m.extinct(df.rest_wave, df.integrated_flux * constant,
                                                  df.uncert_flux * constant, E, E_u, doMC=False, size=niter)
        df.insert(4, 'dered_flux', dered_flux)
        df.insert(5, 'uncert_dered_flux', uncert_dered_flux)
        df['dered_flux'] /= constant
        df['uncert_dered_flux'] /= constant

        filename = filename[:-4] + '_dered.txt'
        header += '\
    dered_flux:    dereddened flux with E(B-V) = ' + str(E) + ' +/- ' + str(E_u) + ' (units of ' + str(constant) + ' erg/s/cm^2)\n\
    uncert_dered_flux:    error in above qty. (units of ' + str(constant) + ' erg/s/cm^2)\n\
'
    df['integrated_flux'] = df['integrated_flux'].astype(np.float64).map('{:,.2f}'.format)
    df['uncert_flux'] = df['uncert_flux'].astype(np.float64).map('{:,.2f}'.format)
    if 'dered_flux' in df:
        df['dered_flux'] = df['dered_flux'].astype(np.float64).map('{:,.2f}'.format)
        df['uncert_dered_flux'] = df['uncert_dered_flux'].astype(np.float64).map('{:,.2f}'.format)

    for i in range(len(df)):
        if df.loc[i]['uncert_flux'] == '%.2F' % (-99.):
            df.ix[i, 'integrated_flux'] = '<' + str(
                df.loc[i]['integrated_flux'])  # to add '<' to fluxes that are upper limits
            df.ix[i, 'uncert_flux'] = '-'
            if 'dered_flux' in df:
                df.ix[i, 'dered_flux'] = '<' + str(
                    df.loc[i]['dered_flux'])  # to add '<' to fluxes that are upper limits
                df.ix[i, 'uncert_dered_flux'] = '-'

    np.savetxt(filename, [], header=header, comments='#')
    df.to_csv(filename, index=None, mode='a', sep='\t')
    print 'Written modified file as', filename
    return filename


# ------------------------------------------------------
def combine_tables(mmt, esi, wfc, esi_to_mmt_scale, esi_to_wfc_scale, master_filename_txt, master_filename_tex):
    print '\nCombining 3 files into one master file...'
    # esi = scale_df(esi, 1./esi_to_wfc_scale)
    # mmt = scale_df(mmt, esi_to_mmt_scale/esi_to_wfc_scale)
    master = pd.concat([mmt, esi, wfc], ignore_index=True).sort_values('uncert_flux').drop_duplicates(subset='ID',
                                                                                                      keep='first').sort_values(
        'rest_wave').reset_index(drop=True)
    specdict = {'MMT': 0, 'ESI': 1, 'WFC3,~G102': 2, 'WFC3,~G141': 3,
                'GEMINI/GNIRS': 4}  # what order of spectrograph we want the final table in
    master['speclabel'] = master['spectrograph'].map(specdict)
    master = master.sort_values('rest_wave').sort_values('speclabel')
    master.drop('speclabel', axis=1, inplace=True)
    master = m.get_flux_from_atomic(master, labelcol='ID', fluxcol='integrated_flux', fluxucol='uncert_flux', \
                                    dered_fluxcol='dered_flux', dered_fluxucol='uncert_dered_flux',
                                    notescol='Notes')  # modify some undetected flux values by tying them to relevant atomic ratios

    # -------for writing txt file-----------
    master_txt = master.copy(True)  # copy df to make changes specifically for txt file (for running PIZI on it)
    master_txt = master_txt.sort_values('rest_wave')
    header = getheader(esi_file)
    np.savetxt(master_filename_txt, [], header=header, comments='#')
    master_txt.to_csv(master_filename_txt, index=None, sep='\t', mode='a')
    print '\nWritten master file', master_filename_txt

    # -------for modifying txt file for PIZI-----------
    master_txt.replace(['OII2470', 'SII6717', 'SIII6312'], ['OII2470a;OII2470b', 'SII6717;SII6731', 'OI6300;SIII6312'],
                       regex=True, inplace=True)
    no_PIZI_lines = ['Fe', 'Mg', 'Hzeta+HeI',
                     'OII3727+9']  # lines certainly not to be included for PIZI; OII doublet from WFC isn't icncluded bcz we have them resolved & measured from ESI
    for l in no_PIZI_lines: master_txt['ID'] = master_txt['ID'].str.replace(re.sub('\+', '\+', l), '#' + l)
    np.savetxt(master_filename_txt[:-4] + '_forPIZI.txt', [], header=header, comments='#')
    master_txt.to_csv(master_filename_txt[:-4] + '_forPIZI.txt', index=None, sep='\t', mode='a')
    print '\nWritten master file', master_filename_txt[:-4] + '_forPIZI.txt'
    print 'PIZI\n', master_txt  #

    # -------for writing tex file-----------
    master.replace(['Halpha', 'Hbeta', 'Hgamma', 'Hdelta', 'Hepsilon', 'Hzeta', 'Heta'], \
                   [r'H$\alpha$', r'H$\beta$', r'H$\gamma$', r'H$\delta$', r'H$\epsilon$', r'H$\zeta$', r'H$\eta$'],
                   inplace=True)  # replacing names of Balmer lines
    try:
        master = format_df(master)
    except:
        pass

    master.replace(float('%.4F' % 3728.5), '3727+3729',
                   inplace=True)  # replacing back doublet rest wavelength after dereddening
    master.replace(float('%.4F' % 6725.), '6718+6732',
                   inplace=True)  # replacing back doublet rest wavelength after dereddening
    master.replace(float('%.4F' % 6307.5), '6302,6313',
                   inplace=True)  # replacing back doublet rest wavelength after dereddening

    master.rename(columns={'ID': 'Line ID', 'rest_wave': '$\lambda_{\mathrm{rest}}$ (\AA)', 'integrated_flux': 'flux', \
                           'uncert_flux': '$\delta$ flux', 'dered_flux': 'dereddened flux',
                           'uncert_dered_flux': '$\delta$ dereddened flux'}, inplace=True)
    # master.columns = master.columns.str.replace('_','\_')
    master.to_latex(master_filename_tex, index=None, escape=False)
    with open(master_filename_tex, 'a') as f:
        f.write(r'\caption{Measured flux values ($flux_{\mathrm{integrated}}$) and corresponding uncertainties ($\Delta flux_{\mathrm{integrated}}$) for nebular emission lines \
in the re-scaled S1723 spectra. The "spectrograph" column denotes which instrument the concerned line was detected with. \
We quote the flux values after accounting for the intrinsic reddening in $flux_{\mathrm{dereddened}}$ and $\Delta flux_{\mathrm{dereddened}}$ respectively. \
All the flux and corresponding uncertainties are in 10$^{-17}$ ergs/s/$\mathrm{cm^2}$ units.}')
    u.insert_line_in_file('\scriptsize\n', 0, master_filename_tex)
    print 'Written master file', master_filename_tex, '\n'

    return master_txt


# ------------------------------------------------------
labels = np.array(['Halpha', 'Hbeta', 'Hgamma', 'Hdelta', 'Hepsilon', 'Hzeta', 'Heta'])
wave = np.array([6564.61, 4862.683, 4341.684, 4102.892, 3971.195, 3890.158,
                 3836.479])  # rest-frame vaccuum wavelengths of Balemr lines
k = m.getfullkappa(wave)
ratio_labels = [r'H$\alpha$/H$\beta$', r'H$\alpha$/H$\beta$', r'H$\beta$/H$\gamma$',
                r'H$\beta$/H$\delta$']  # , r'H$\beta$/H$\epsilon$', r'H$\beta$/H$\zeta$', r'H$\zeta$/H$\epsilon$']
spec_labels = ['G141', 'G141/G102', 'G102', 'G102']  # , 'G102', 'G102/ESI', 'ESI']
delta_k = np.array([k[1] - k[0], k[1] - k[0], -k[1] + k[2], -k[1] + k[3]])  # , -k[1]+k[4], -k[1]+k[5], -k[5]+k[6]])
ratio_theoretical = np.array([2.8785 / 1, 2.8785 / 1, 1 / 0.46756,
                              1 / .25825])  # , 1/.15856, 1/0.10472, 0.10472/.072916]) #from /Users/acharyya/Desktop/Lisa_UV_diag/P_spherical/spherical/sp_P60_a05modelfiles/Q800/spec0003.csv
# MAPPINGS V model for lpok=6, logq=8, Z=8.23
# ------------------------------------------------------
mmt_file = HOME + '/Dropbox/mage_atlas/Contrib/EWs/emission/s1723_MMT_emission_measured_lines.out.txt'
# subprocess.call(['python ../../abundance_pap/dftolatex.py --short s1723_MMT_wcont_new-format --infile '+mmt_file+' --EW_thresh 3 --SNR_thresh 0 --const 1e-17 --nopizi --notex'], shell=True)
mmt_file = mmt_file[:-4] + '_detected.txt'

esi_file = HOME + '/Dropbox/mage_atlas/Contrib/EWs/emission/s1723_ESI_emission_measured_lines.out.txt'
# subprocess.call(['python ../../abundance_pap/dftolatex.py --short s1723_ESI_wcont_new-format --infile '+esi_file+' --EW_thresh 3 --SNR_thresh 0 --const 1e-17 --nopizi --notex'], shell=True)
esi_file = esi_file[:-4] + '_detected.txt'

wfc_file = HOME + '/Dropbox/grism_s1723/WFC3_fit_1Dspec/WFC3_measuredlines.txt'
wfc_file = deredden(wfc_file, 0.3, 0.02, niter=int(1e5), constant=1e-17, dered=True, change_ID=True, change_errors=True,
                    SNR_thresh=0.014)  # to deredden WFC3 spectra and re-write the file

mmt = pd.read_table(mmt_file, delim_whitespace=True, comment='#')
esi = pd.read_table(esi_file, delim_whitespace=True, comment='#')
wfc = pd.read_table(wfc_file, delim_whitespace=True, comment='#')

g102 = wfc[wfc['spectrograph'].str.contains('G102')]
g141 = wfc[wfc['spectrograph'].str.contains('G141')]
# ------------------------------------------------------
esi_hlist = ['Hzeta', 'Heta']
g102_hlist = ['Hbeta', 'Hgamma', 'Hdelta', 'Hepsilon']
g141_hlist = ['Halpha', 'Hbeta']
c3list = ['CIII1906', 'CIII1908']
esi_o2list = ['OII3727', 'OII3729']
g102_o2list = ['OII3727+9']  # ['[O~II]']
# ---------------scaling---------------------------------------
C3_mmt = unp.uarray(getf(mmt, c3list), gete(mmt, c3list))
C3_esi = unp.uarray(getf(esi, c3list), gete(esi, c3list))
O2_esi = unp.uarray(getf(esi, esi_o2list), gete(esi, esi_o2list))
O2_g102 = unp.uarray(getf(g102, g102_o2list), gete(g102, g102_o2list))

mmt_CIII, esi_CIII = np.sum(C3_mmt), np.sum(C3_esi)
wfc_OII, esi_OII = np.sum(O2_g102), np.sum(O2_esi)
esi_to_mmt_scale = esi_CIII / mmt_CIII
esi_to_wfc_scale = esi_OII / wfc_OII

print '\nCIII1907-9 ESI/MMT ratio =', esi_to_mmt_scale
print 'OII3727-9 ESI/WFC ratio =', esi_to_wfc_scale
# ---------------reddening---------------------------------------
H_esi = unp.uarray(getf(esi, esi_hlist), gete(esi, esi_hlist))
H_g102 = unp.uarray(getf(g102, g102_hlist), gete(g102, g102_hlist))
H_g141 = unp.uarray(getf(g141, g141_hlist), gete(g141, g141_hlist))

ratio_obs = [H_g141[0] / H_g141[1], H_g141[0] / H_g102[0], H_g102[0] / H_g102[1],
             H_g102[0] / H_g102[2]]  # , H_g102[0]/H_g102[3], H_g102[0]/(H_esi[0]*esi_to_mmt_scale), H_esi[0]/H_esi[1]]
Eb_v = (2.5 / delta_k) * unp.log10(ratio_obs / ratio_theoretical)
print '\nMean E(B-V) =', np.mean(Eb_v)
edf = pd.DataFrame(
    np.transpose(np.vstack([ratio_labels, spec_labels, ['%.2F' % x.n for x in Eb_v], ['%.2F' % x.s for x in Eb_v]])),
    columns=['Lines', 'Spectrograph', 'E(B-V)', 'E(B-V)_u'])
fout = HOME + '/Dropbox/Grism_S1723/Latex/Tabs/s1723_red_table.tex'
edf.to_latex(fout, index=None, escape=False)
with open(fout, 'a') as f:
    f.write(r'\caption{E(B-V) measurements using different Balmer line pairs (column "Lines"). The spectrographs corresponding to \
the line pairs are listed in "Spectrograph". The reddening values and corresponding uncertainties are quoted in the "E(B-V)"\
and "E(B-V)\_u" columns respectively.}')

print '\nList of E(B-V) values\n', edf
print 'Written file', fout
# ------------------------------------------------------
master = combine_tables(mmt, esi, wfc, esi_to_mmt_scale.n, esi_to_wfc_scale.n,
                        HOME + '/Dropbox/Mage_atlas/Contrib/EWs/emission/s1723_measured_emission.all.txt',
                        HOME + '/Dropbox/Grism_S1723/Latex/Tabs/s1723_measured_emission.all.tex')  # to combine all 3 files into one master file, by scaling fluxes at OII3727-9 and CIII1907-9
