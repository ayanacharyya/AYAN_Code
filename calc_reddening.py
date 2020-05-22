# python routine to compute E(B-V) using several Balmer lines for S1723
# and then scale the dereddened fluxes across different spectrographs to create a master linelist
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
import dftolatex as d2l


# ------------------------------------------------------
def getf(df, lines, fluxcol='integrated_flux', labelcol='ID'):  # function to extract flux values from dataframe
    return [df[df[labelcol] == line][fluxcol].values[0] for line in lines]


def gete(df, lines, errorcol='uncert_flux', labelcol='ID'):  # function to extract flux uncertainties from dataframe
    return [df[df[labelcol] == line][errorcol].values[0] for line in lines]


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
            for i in range(len(df)):
                item = str(df.iloc[i][column])
                if item[0].isdigit(): item = '%.2F'%(float(item)*scale_factor)
                elif len(item)>1: item = item[0] + '%.2F'%(float(item[1:])*scale_factor)
                df = df.set_value(i, column, item)
        else:
            continue
    return df


# ------------------------------------------------------
def format_df(df, columns_to_scale=['integrated_flux', 'uncert_flux', 'dered_flux', 'uncert_dered_flux']):
    to_replace = [('<', r'$<$'), ('>', r'$>$'), ('nan', '..'), ('Ly-alpha', r'Ly$\\alpha$')]  # [(replace_this, with_this),..]
    to_replace = np.array(to_replace)
    df2 = df.replace(to_replace[:, 0].tolist(), to_replace[:, 1].tolist(), regex=True)
    for column in columns_to_scale:
        if column in df2:
            try: df2[column] = df2[column].map('{:.3f}'.format)
            except: pass
        else:
            continue
    return df2


# ------------------------------------------------------
def deredden(filename, E, E_u, niter=int(1e5), constant=1e-17, dered=True, change_ID=False, change_errors=False,
             SNR_thresh=0., readcsv=False):
    header = getheader(filename)
    if readcsv:
        df = pd.read_csv(filename, comment='#')
        df = df.rename(columns={'linename':'ID', 'restwave':'rest_wave', 'flux':'integrated_flux', 'flux_u':'uncert_flux'})
        df['rest_wave'] = df['rest_wave'].astype(np.float64)
        if 'spectrograph' not in df:
            if 'G102' in filename: df['spectrograph'] = 'G102'
            elif 'G141' in filename: df['spectrograph'] = 'G141'
        if 'Notes' not in df: df['Notes'] = '--'
    else:
        df = pd.read_table(filename, delim_whitespace=True, comment='#')

    if change_errors:
        df['SNR'] = df['uncert_flux'] / df['integrated_flux']
        concerned = (df['uncert_flux'] > 0.) & (df['SNR'] < SNR_thresh)
        if len(df[concerned]) > 0:
            print 'Before changing uncertainties for\n', df[concerned]
            df.ix[concerned, 'uncert_flux'] = df.ix[concerned, 'integrated_flux'] * SNR_thresh  # modifying uncertainties by hand
            df['SNR'] = df['uncert_flux'] / df['integrated_flux']
            print 'After changing uncertainties for\n', df[concerned]
        df.drop('SNR', axis=1, inplace=True)
        filename = os.path.splitext(filename)[0] + '_umod.txt'
        header += '\
    "uncert_flux" has been modified as follows:\n\
    If any line had a "integrated_flux"/"uncert_flux" < ' + str(
            SNR_thresh) + ', "uncert_flux" was replaced such that this ratio becomes equal to ' + str(SNR_thresh) + '\n\
'
    if change_ID:
        label_dict = {3727.092:'OII3727', 3729.875:'OII3729', 3869.860:'NeIII3869', 3890.151:'Hzeta+HeI', 3968.593:'NeIII3968', \
                      3971.195:'Hepsilon', 4025.000:'HeI4025;HeII', 4102.892:'Hdelta', 4341.684:'Hgamma', 4364.436:'OIII4363', \
                      4472.700:'HeI4471', 4687.020:'HeII4685', 4741.449:'ArIV4741', 4862.683:'Hbeta', 4960.295:'OIII4959', 5008.240:'OIII5007', \
                      5877.590:'HeI5875', 6310.000:'OI6300;SIII6312', 6549.850:'NII6549', 6564.610:'Halpha', 6585.280:'NII6584', \
                      6725.000:'SII6717;SII6731'}

        df['ID'] = df['rest_wave'].apply(lambda x: label_dict[float('%.3F'%x)])
        filename = os.path.splitext(filename)[0] + '_lmod.txt'
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

        filename = os.path.splitext(filename)[0] + '_dered.txt'
        header += '\
    dered_flux:    dereddened flux with E(B-V) = ' + str(E) + ' +/- ' + str(E_u) + ' (units of ' + str(constant) + ' erg/s/cm^2)\n\
    uncert_dered_flux:    error in above qty. (units of ' + str(constant) + ' erg/s/cm^2)\n\
'
    df['integrated_flux'] = df['integrated_flux'].astype(np.float64).map('{:,.2f}'.format)
    df['uncert_flux'] = df['uncert_flux'].astype(np.float64).map('{:,.2f}'.format)
    if 'dered_flux' in df:
        df['dered_flux'] = df['dered_flux'].astype(np.float64).map('{:,.2f}'.format)
        df['uncert_dered_flux'] = df['uncert_dered_flux'].astype(np.float64).map('{:,.2f}'.format)

    df = df[['ID', 'rest_wave', 'integrated_flux', 'uncert_flux', 'dered_flux', 'uncert_dered_flux', 'spectrograph', 'Notes']]

    for i in range(len(df)):
        if df.loc[i]['uncert_flux'] == '%.2F' % (-99.):
            df.ix[i, 'integrated_flux'] = '<' + str(
                df.loc[i]['integrated_flux'])  # to add '<' to fluxes that hiscolare upper limits
            df.ix[i, 'uncert_flux'] = '-'
            if 'dered_flux' in df:
                df.ix[i, 'dered_flux'] = '<' + str(
                    df.loc[i]['dered_flux'])  # to add '<' to fluxes that are upper limits
                df.ix[i, 'uncert_dered_flux'] = '-'

    cols_to_add = ['EWr_fit', 'EWr_fit_u', 'EW_signi'] # these columns to be added as dummy, if they don't already exist
    for thiscol in cols_to_add:
        if thiscol not in df.columns: df[thiscol] = '-'

    np.savetxt(filename, [], header=header, comments='#')
    df.to_csv(filename, index=None, mode='a', sep='\t')
    print 'Written modified file as', filename
    return filename


# ------------------------------------------------------
def combine_tables(mmt, esi, wfc, esi_to_mmt_scale, esi_to_wfc_scale, master_filename_txt, master_filename_tex):
    print '\nCombining 3 files into one master file...'
    esi = scale_df(esi, 1./esi_to_wfc_scale)
    mmt = scale_df(mmt, esi_to_mmt_scale/esi_to_wfc_scale)
    master = pd.concat([mmt, esi, wfc], ignore_index=True).sort_values('uncert_flux').sort_values('rest_wave').reset_index(drop=True)
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
    master_txt.replace(['OII2470'], ['OII2470a;OII2470b'],
                       regex=True, inplace=True)
    no_PIZI_lines = ['Fe', 'Mg', 'Hzeta+HeI',
                     'OII3727+9']  # lines certainly not to be included for PIZI; OII doublet from WFC isn't icncluded bcz we have them resolved & measured from ESI
    for l in no_PIZI_lines: master_txt['ID'] = master_txt['ID'].str.replace(re.sub('\+', '\+', l), '#' + l)
    pizi_filename = os.path.splitext(master_filename_txt.replace('txt_tables' ,'PIZI_tables'))[0] + '_forPIZI.txt'
    np.savetxt(pizi_filename, [], header=header, comments='#')
    master_txt.to_csv(pizi_filename, index=None, sep='\t', mode='a')
    print '\nWritten master file', pizi_filename

    # -------for writing tex file-----------
    master.replace(['Halpha', 'Hbeta', 'Hgamma', 'Hdelta', 'Hepsilon', 'Hzeta', 'Heta', 'Hzeta+HeI'], \
                   [r'H$\alpha$', r'H$\beta$', r'H$\gamma$', r'H$\delta$', r'H$\epsilon$', r'H$\zeta$', r'H$\eta$', r'H$\zeta$;HeI'],
                   inplace=True)  # replacing names of Balmer lines
    master = format_df(master)
    master.drop(['Notes'], axis=1, inplace=True)
    column_order = ['ID', 'rest_wave', 'spectrograph', 'EWr_fit', 'EWr_fit_u', 'EW_signi', 'integrated_flux', 'uncert_flux', 'dered_flux', 'uncert_dered_flux']
    master = master[column_order]

    master.replace(float('%.4F' % 3728.5), '3727+3729',
                   inplace=True)  # replacing back doublet rest wavelength after dereddening
    master.replace(float('%.4F' % 6725.), '6718+6732',
                   inplace=True)  # replacing back doublet rest wavelength after dereddening
    master.replace(float('%.4F' % 6307.5), '6302,6313',
                   inplace=True)  # replacing back doublet rest wavelength after dereddening

    master.rename(columns={'ID':'Line ID', 'rest_wave':r'$\lambda_{\mathrm{rest}}$', 'EWr_fit':r'W$_{\mathrm{r,fit}}$', 'EWr_fit_u':r'$\Delta$ W$_{\mathrm{r,fit}}$', \
     'EW_signi':r'W$_{\mathrm{r,signi}}$', 'integrated_flux':r'flux', 'uncert_flux':r'$\Delta$ flux', \
     'dered_flux':r'flux$_{\mathrm{deredenned}}$', 'uncert_dered_flux':r'$\Delta$ flux$_{\mathrm{deredenned}}$'}, inplace=True)

    # master.columns = master.columns.str.replace('_','\_')
    master['Line ID'] = master['Line ID'].apply(lambda x: d2l.format_forbidden(x))  # to format line labels with space and brackets

    master.to_latex(master_filename_tex, index=None, escape=False)
    with open(master_filename_tex, 'a') as f:
        f.write(r'\caption{Measured flux values and corresponding uncertainties ($\Delta$ flux) for nebular emission lines \
in the re-scaled S1723 spectra. Each emission line is labeled with "Line ID" and its rest-frame wavelength is indicated as $\lambda_{\mathrm{rest}}$. \
The "spectrograph" column denotes which instrument the concerned line was detected with. \
W$_{\mathrm{r,fit}}$ denotes the fitted, rest-frame equivalent width. The corresponding uncertainty and signifcance (see Section XYZ) are quoted as $\Delta$ W$_{\mathrm{r,fit}}$ and W$_{\mathrm{r,signi}}$ respectively.\
We also quote the flux values after accounting for the intrinsic reddening and the corresponsing uncertainty as $flux_{\mathrm{dereddened}}$ and $\Delta flux_{\mathrm{dereddened}}$ respectively. \
All the flux and corresponding uncertainties are in 10$^{-17}$ ergs/s/$\mathrm{cm^2}$ units.}')
    u.insert_line_in_file('\scriptsize\n', 0, master_filename_tex)
    
    print 'Written master file', master_filename_tex, '\n'
    print 'tex file\n', master  #

    return master_txt


# ------------------------------------------------------
if __name__ == '__main__':
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
    txt_table_path = '/Dropbox/Mage_atlas/Contrib/EWs/emission/txt_tables/'

    mmt_file = HOME + txt_table_path + 's1723_MMT_emission_measured.txt'
    subprocess.call(['python dftolatex.py --short s1723_MMT_wcont_new-format --infile '+mmt_file+' --outpath '+'/'.join(os.path.split(mmt_file)[0].split('/')[:-1])+'/'+' --EW_thresh 3 --SNR_thresh 0 --const 1e-17 --nopizi --notex'], shell=True)
    mmt_file = os.path.splitext(mmt_file)[0]+'_detected.txt'

    esi_file = HOME + txt_table_path + 's1723_ESI_emission_measured.txt'
    subprocess.call(['python dftolatex.py --short s1723_ESI_wcont_new-format --infile '+esi_file+' --outpath '+'/'.join(os.path.split(esi_file)[0].split('/')[:-1])+'/'+' --EW_thresh 3 --SNR_thresh 0 --const 1e-17 --nopizi --notex'], shell=True)
    esi_file = os.path.splitext(esi_file)[0]+'_detected.txt'

    wfc_g102_file = HOME + '/Dropbox/Grism_S1723/WFC3_fit_1Dspec/1Dsum/sgas1723_1Dsum_bothroll_G102_wcontMWdr_meth2.fitdf'
    wfc_g102_file = deredden(wfc_g102_file, 0.028, 0.04, niter=int(1e5), constant=1e-17, dered=True, change_ID=True, change_errors=True,
                        SNR_thresh=0.014, readcsv=True)  # to deredden WFC3 spectra and re-write the file
    # E and E_u are value Ha/Hb from G102/G141 new reduction of bothrolls, computed by AYAN_Code/make_reddening_table.py

    wfc_g141_file = HOME + '/Dropbox/Grism_S1723/WFC3_fit_1Dspec/1Dsum/sgas1723_1Dsum_bothroll_G141_wcontMWdr_meth2.fitdf'
    wfc_g141_file = deredden(wfc_g141_file, 0.028, 0.04, niter=int(1e5), constant=1e-17, dered=True, change_ID=True, change_errors=True,
                        SNR_thresh=0.014, readcsv=True)  # to deredden WFC3 spectra and re-write the file

    mmt = pd.read_table(mmt_file, delim_whitespace=True, comment='#')
    esi = pd.read_table(esi_file, delim_whitespace=True, comment='#')
    g102 = pd.read_table(wfc_g102_file, delim_whitespace=True, comment='#')
    g141 = pd.read_table(wfc_g141_file, delim_whitespace=True, comment='#')
    wfc = pd.concat([g102,g141])
    # ------------------------------------------------------
    esi_hlist = ['Heta']
    c3list = ['CIII1906', 'CIII1908']
    esi_o2list = ['OII3727', 'OII3729']
    g102_o2list = ['OII3727', 'OII3729']  # ['[O~II]']
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
    # ------------------------------------------------------
    master = combine_tables(mmt, esi, wfc, esi_to_mmt_scale.n, esi_to_wfc_scale.n,
                            HOME + txt_table_path + 's1723_measured_emission_all.txt',
                            HOME + '/Dropbox/Grism_S1723/Latex/Tabs/s1723_measured_emission_all.tex')  # to combine all 3 files into one master file, by scaling fluxes at OII3727-9 and CIII1907-9
