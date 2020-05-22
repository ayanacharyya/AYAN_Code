##---to import line_table as pandas dataframe, apply detection criteria, and pump out latex tables----##
##----by Ayan-------##
import numpy as np
import pandas as pd

pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
import argparse as ap
import os
HOME = os.getenv('HOME') + '/'
import sys

sys.path.append(HOME + 'Work/astro/ayan_codes/mageproject/')
import ayan.splot_util as u
import ayan.mage as m
import re
import subprocess

# -----------to check if line is detected given our etection criteria------------------------------------
def detect(table, EW_thresh=0, SNR_thresh=0, EW_signi='EW_signi', flux='f_line', flux_u='f_line_u'):
    return ((np.abs(table[EW_signi]) >= EW_thresh) & (table[flux].values / table[flux_u].values >= SNR_thresh))


# ------------to check if line is (semi/)forbidden and format Id accordingly-----------------------------------
def format_forbidden(label1):
    forbidden = ['CIII1906', 'SiII1808', 'SiII1816', 'SiIII1882', 'OIII2320', 'OII2470', 'OII3727', 'OII3729', 'OII3727,9', 'NeIII3869',
                 'NeIII3968', \
                 'OIII4362', 'OIII4959', 'OIII5007', 'ArIV4741', 'NII6584', 'SII6717', 'SII6730', 'ArIII7136']
    semiforbidden = ['CIII1908', 'NII1430', 'NII1431', 'NIV1486', 'OIII1660', 'OIII1666', 'NIII1750', 'SiIII1892',
                     'NII2140', \
                     'CII2323', 'CII2325c', 'CII2325d', 'CII2328', 'SiII2335a', 'SiII2335b']
    new_label = []
    for label in label1.split(';'):
        if label in forbidden:
            brackets = ('{[}', ']')  # curly braces otherwise latex table does not compile
        elif label in semiforbidden:
            brackets = ('', ']')
        else:
            brackets = ('', '')
        label = brackets[0] + label
        label = re.sub(r"(\d+)", brackets[1] + r"\1", label, count=1)
        label = re.sub(r"(I+)", r" \1", label, count=1)
        label = re.sub(r"(\d+)", r" \1", label, count=1)  # to format line labels with space
        new_label.append(label)

    return '; '.join(new_label)


# -----to format latex column headers------------------------------------------
def column_header_format(header, string_replace_dict, substring_replace_dict, sep='_'):
    n = len(sep)
    index1 = header.find(sep)
    if index1 == -1:
        string = header
        substring, subsubstring = '', ''
    else:
        string = header[:index1]
        index2 = header[index1 + n:].find(sep)
        if index2 == -1:
            substring = header[index1 + n:]
            subsubstring = ''
        else:
            substring = header[index1 + n:index1 + n + index2]
            subsubstring = header[index1 + index2 + n + n:]
    isdelta = r'$\delta$ ' if subsubstring == 'u' else r''
    string_replace = string_replace_dict[string] if string in string_replace_dict else string
    substring_replace = substring_replace_dict[substring] if substring in substring_replace_dict else substring

    if string_replace == 'rest' and substring_replace == '$\lambda$': string_replace, substring_replace = '$\lambda$', 'rest'

    header_string = isdelta + string_replace
    header_substring = ' ID' if substring_replace == 'ID' else '$_{\mathrm{' + substring_replace + '}}$'

    new_header = header_string + header_substring
    new_header = new_header.replace('$$', '')
    return new_header


# -----------------------------------------------
if __name__ == '__main__':
    parser = ap.ArgumentParser(description="Tool to transform pandas dataframe to tex and ASCII table")
    parser.add_argument('--onlyem', dest='onlyem', action='store_true')
    parser.set_defaults(onlyem=False)
    parser.add_argument('--onlyinterv', dest='onlyinterv', action='store_true')
    parser.set_defaults(onlyinterv=False)
    parser.add_argument('--nopizi', dest='nopizi', action='store_true')
    parser.set_defaults(nopizi=False)
    parser.add_argument('--notex', dest='notex', action='store_true')
    parser.set_defaults(notex=False)

    parser.add_argument("--infile")
    parser.add_argument("--outpath")
    parser.add_argument("--outfile")
    parser.add_argument("--shortlabel")
    parser.add_argument("--EW_thresh")
    parser.add_argument("--SNR_thresh")
    parser.add_argument("--const")
    args, leftovers = parser.parse_known_args()
    if args.infile is not None:
        infile = args.infile
    else:
        infile = HOME + 'Dropbox/papers/abundance_pap/AA_working/lineflux_restUV.txt'  # full-path-name of input file
    inpath = os.path.split(infile)[0] + '/'
    if args.shortlabel is None: args.shortlabel = 'rcs0327-E'

    if args.outpath is not None:
        outpath = args.outpath
    else:
        outpath = inpath  # full-path of output directory
    subprocess.call(['mkdir -p '+outpath+'PIZI_tables/'], shell=True)
    subprocess.call(['mkdir -p '+outpath+'txt_tables/'], shell=True)
    subprocess.call(['mkdir -p '+outpath+'tex_tables/'], shell=True)

    if args.outfile is not None:
        outfile = args.outfile
    elif 'allspec' in infile:
        outfile = os.path.splitext(os.path.basename(infile))[0].replace('allspec', args.shortlabel) + '_detected'
    else:
        outfile = os.path.splitext(os.path.basename(infile))[0] + '_detected'

    if args.EW_thresh is not None:
        EW_thresh = float(args.EW_thresh)
    else:
        EW_thresh = 3.
    if args.SNR_thresh is not None:
        SNR_thresh = float(args.SNR_thresh)
    else:
        SNR_thresh = 1.
    if args.const is not None:
        const = float(args.const)
    else:
        const = 1e-17
    # -----------------------------------------------
    line_table = pd.read_table(infile, delim_whitespace=True, comment="#")  # input dataframe file
    print 'Deb136:', args.shortlabel, infile #
    if args.shortlabel is not None: line_table = line_table[line_table['label'].eq(args.shortlabel)]
    if args.onlyem: line_table = line_table[line_table.EWr_fit <= 0.].reset_index(drop=True)
    if args.onlyinterv:
        line_table = line_table[line_table.type == 'INTERVE'].reset_index(drop=True)
        quantities_to_extract = ['line_lab', 'rest_wave', 'EWr_fit', 'EWr_fit_u', 'zz', 'zz_u']
        outfile += '_interv'
        tab = line_table[quantities_to_extract]

    else:
        quantities_to_extract = ['line_lab', 'rest_wave', 'EWr_fit', 'EWr_fit_u', 'EW_signi', 'EWr_Suplim', 'f_line',
                                 'f_line_u', 'f_Suplim']
        if 'f_redcor' in line_table:
            quantities_to_extract += ['f_redcor', 'f_redcor_u', 'f_Suplim_redcor']

        tab = line_table[quantities_to_extract]
        tab.loc[:, 'f_line'] /= const
        tab.loc[:, 'f_line_u'] /= const
        tab.loc[:, 'f_Suplim'] /= const
        # print tab#
        if 'f_redcor' in tab:
            tab.loc[:, 'f_redcor'] /= const
            tab.loc[:, 'f_redcor_u'] /= const
            tab.loc[:, 'f_Suplim_redcor'] /= const

        bad = ~detect(tab, EW_thresh=EW_thresh, SNR_thresh=SNR_thresh)

        tab.EW_signi = tab.EW_signi.astype(np.float64)
        tab.EWr_Suplim = tab.EWr_Suplim.astype(np.str)
        if 'f_redcor' in tab: tab.loc[:, 'f_Suplim_redcor'] = tab.loc[:, 'f_Suplim_redcor'].map('{:,.2f}'.format).astype(
            np.str)
        tab.loc[:, 'f_Suplim'] = tab.loc[:, 'f_Suplim'].map('{:,.2f}'.format).astype(np.str)
        tab.EWr_fit = tab.EWr_fit.astype(np.str)
        if 'f_redcor' in tab: tab.f_redcor = tab.f_redcor.map('{:,.2f}'.format).astype(np.str)
        tab.loc[:, 'f_line'] = tab.loc[:, 'f_line'].map('{:,.2f}'.format).astype(np.str)
        tab.loc[:, 'EWr_fit'] = np.where(bad, '>' + tab.EWr_Suplim, tab.EWr_fit) # sign is other way round because EW is -ve
        if 'f_redcor' in tab: tab.loc[:, 'f_redcor'] = np.where(bad, '<' + tab.f_Suplim_redcor, tab.f_redcor)
        tab.loc[:, 'f_line'] = np.where(bad, '<' + tab.f_Suplim, tab.f_line)
        tab.EWr_fit_u = tab.EWr_fit_u.astype(np.str)
        tab.f_line_u = tab.f_line_u.astype(np.str)
        tab.loc[bad, 'EWr_fit_u'] = None
        tab.loc[bad, 'f_line_u'] = None
        if 'f_redcor' in tab: tab.loc[bad, 'f_redcor_u'] = None
        tab.loc[bad, 'EW_signi'] = None

        quantities_to_show = ['line_lab', 'rest_wave', 'EWr_fit', 'EWr_fit_u', 'EW_signi', 'f_line', 'f_line_u']
        if 'f_redcor' in tab: quantities_to_show += ['f_redcor', 'f_redcor_u']
        tab = tab[quantities_to_show]

        tab['f_line_u'] = tab['f_line_u'].astype(np.float64).map('{:,.2f}'.format)
        if 'f_redcor' in tab: tab['f_redcor_u'] = tab['f_redcor_u'].astype(np.float64).map('{:,.2f}'.format)
        tab['EW_signi'] = tab['EW_signi'].astype(np.float64).map('{:,.2f}'.format)
        tab['EWr_fit_u'] = tab['EWr_fit_u'].astype(np.float64).map('{:,.2f}'.format)

    tab.line_lab = tab.line_lab.str.replace('_', '').str.replace('*', '')
    if not args.onlyinterv: tab = m.get_flux_from_atomic(tab, labelcol='line_lab', fluxcol='f_line', fluxucol='f_line_u', \
                                                         dered_fluxcol='f_redcor', dered_fluxucol='f_redcor_u',
                                                         bad_value='nan')  # modify some undetected flux values by tying them to relevant atomic ratios
    # ------modifying column names for S1723----------
    if 's1723' in args.shortlabel:
        # tab = tab.drop(['EWr_fit', 'EWr_fit_u', 'EW_signi'], axis=1) # commented out by AA on 4th Mar '19 because JRR needs EW values in final files
        tab.rename(columns={'line_lab': 'ID', 'f_line': 'integrated_flux', 'f_line_u': 'uncert_flux'}, inplace=True)
        if 'f_redcor' in tab: tab.rename(columns={'f_redcor': 'dered_flux', 'f_redcor_u': 'uncert_dered_flux'},
                                         inplace=True)
        if 'ESI' in args.shortlabel:
            spectroscope = 'ESI'
        elif 'MMT' in args.shortlabel:
            spectroscope = 'MMT'
        tab['spectrograph'] = spectroscope
        tab['Notes'] = '--'

        header = 'This file contains the measurements of lines in the MagE sample. Generated by dftolatex.py.\n\
        From file ' + infile + '\n\
        Columns are:\n\
         line_lab:    label of the line the code was asked to fit\n\
         rest_wave:   rest frame vacuum wavelength of the line (A)\n\
         integrated_flux:      flux i.e. area under Gaussian fit (units of ' + str(const) + ' erg/s/cm^2)\n\
         uncert_flux:    error in above qty. (units of ' + str(const) + ' erg/s/cm^2)\n\
        '
        if 'dered_flux' in tab:
            header += '\
            dered_flux:    dereddened flux (units of ' + str(const) + ' erg/s/cm^2)\n\
            uncert_dered_flux:    error in above qty. (units of ' + str(const) + ' erg/s/cm^2)\n\
            '
    if 'stack' in args.shortlabel:
        tab = tab.drop(['f_line', 'f_line_u'], axis=1)
        if 'f_redcor' in tab: tab.drop(['f_redcor', 'f_redcor_u'], axis=1)
        tab.rename(columns={'line_lab': 'ID'}, inplace=True)

        header = 'This file contains the measurements of lines in the MagE stacked spectra. Generated by dftolatex.py.\n\
        From file ' + infile + '\n\
        Columns are:\n\
         line_lab:    label of the line the code was asked to fit\n\
         rest_wave:   rest frame vacuum wavelength of the line (A)\n\
         EWr_fit:      Equivalent width of Gaussian fit (units of A)\n\
         EWr_fit_u:    error in above qty. (units of A)\n\
         EW_signi:      significance of the detection\n\
        '
    else:
        header = '#Fluxes are in %.2E flambda units\n' % (const)
    # --replacing nans and line summed IDs in txt file--
    tab_txt = tab.replace(['nan'], ['-'], regex=True)

    # --adding header to txt file--
    fullfilename = outpath + 'txt_tables/' + outfile + '.txt'
    np.savetxt(fullfilename, [], header=header, comments='#')
    tab_txt.to_csv(fullfilename, sep='\t', mode='a', index=None)
    print 'Written files ' + fullfilename

    if not args.nopizi and not args.onlyinterv:
        # --writing to separate txt file for use by PIZI--
        fullfilename = outpath + 'PIZI_tables/' + outfile + '_forPIZI.txt'
        tab_txt = tab_txt.replace(['OII2470'], ['OII2470a;OII2470b'], regex=True)
        no_PIZI_lines = ['Fe', 'Mg', 'Blend']  # lines certainly not to be included for PIZI
        labelcol = 'ID' if 's1723' in args.shortlabel else 'line_lab'
        for l in no_PIZI_lines: tab_txt[labelcol] = tab_txt[labelcol].str.replace(l, '#' + l)
        with open(fullfilename, 'w') as file:
            file.write(header)
        tab_txt.to_csv(fullfilename, sep='\t', mode='a', index=None)

        # --adding few lines by hand to txt file for use by PIZI--
        if args.shortlabel == 'rcs0327-E' and not args.onlyinterv:
            with open(fullfilename, 'a') as f:
                f.write('\
    #Added OII lines to use for PIZI input for UV+O2 case\n\
    #OII3727	3727.092		999     999     509.228	17.400	1281.769	28.999\n\
    #OII3729	3729.875		999     999     443.935	23.182	1512.625	34.773\n\
    ')
        print 'Written ' + fullfilename

    if not args.notex:
        # --replacing stuff in tex file--
        to_replace = [('<', r'$<$'), ('nan', '..'), ('Ly-alpha', r'Ly$\\alpha$')]  # [(replace_this, with_this),..]
        to_replace = np.array(to_replace)
        tab_tex = tab.replace(to_replace[:, 0].tolist(), to_replace[:, 1].tolist(), regex=True)
        if 'f_redcor_u' in tab_tex: tab_tex.drop('f_redcor_u', axis=1,
                                                 inplace=True)  # to exclude column of uncertainty in dereddened flux

        # --formatting column names in tex file--
        string_replace_dict = {'line': 'Line', 'f': 'flux'}
        substring_replace_dict = {'lab': 'ID', 'wave': '$\lambda$', 'redcor': 'dereddened', 'line': 'obs'}
        if 'stack' in args.shortlabel:
            tab_tex.rename(columns={'EWr_fit': r'W$_{r}$', 'EWr_fit_u': r'$\delta$ W$_{r}$', 'EW_signi': r'significance'},
                           inplace=True)
        else:
            for i in range(len(tab_tex.columns)):
                tab_tex.columns.values[i] = column_header_format(tab_tex.columns.values[i], string_replace_dict,
                                                                 substring_replace_dict)
            tab_tex.rename(columns={r'EWr$_{\mathrm{fit}}$': r'W$_{\mathrm{r,fit}}$',
                                    r'$\delta$ EWr$_{\mathrm{fit}}$': r'$\delta$ W$_{\mathrm{r,fit}}$',
                                    r'EW$_{\mathrm{signi}}$': r'W$_{\mathrm{r,signi}}$'}, inplace=True)

        tab_tex['Line ID'] = tab_tex['Line ID'].apply(lambda x: format_forbidden(x))  # to format line labels with space and brackets
        # --writing the tex file--
        fullfilename = outpath + 'tex_tables/' + outfile + '.tex'
        tab_tex.to_latex(fullfilename, index=None, escape=False)

        # --assigning captions--
        if args.shortlabel == 'rcs0327-E':
            if args.onlyinterv:
                caption = 'Intervening absorption lines in \knotE Magellan/MagE spectrum. W$_{\mathrm{r,fit}}$ denotes the rest-frame \
        equivalent width measured, in \AA. $z$ an $\Delta z$ are the redshift and corresponding uncertainty respectively, as \
        measured from our line fitting code.'
                subcaption = r'& (\AA) & (\AA) & (\AA) & & \\'
            elif args.onlyem:
                caption = 'MagE/Magellan line flux measurements. flux$_{\mathrm{obs}}$ and $\delta$ flux$_{\mathrm{obs}}$ denote \
    the observed flux and uncertainty respectively. flux$_{\mathrm{dereddened}}$ is the dereddened flux using \
    E(B-V) = 0.4 $\pm$ 0.07. W$_\mathrm{r,fit}$ and $\delta$ W$_{\mathrm{r,fit}}$ denote the rest-frame \
    equivalent width measured and the corresponding uncertainty in \AA\, respectively. For cases of \
    non-detection (i.e. < ' + str(int(EW_thresh)) + ' $\sigma$ detection), the 3 $\sigma$ upper limit on equivalent widths and fluxes are quoted. \
    Uncertainty estimates for these entries are not quoted because they do not provide any meaningful information.'
                subcaption = '& (\AA) & (\AA) & (\AA) & (\AA) & (10$^{' + str(int(np.log10(const))) + '}$ ergs/s/cm$^2$) & \
                (10$^{' + str(int(np.log10(const))) + '}$ ergs/s/cm$^2$) & (10$^{' + str(int(np.log10(const))) + r'}$ ergs/s/cm$^2$) \\'
        else:
            caption = ''
            subcaption = ''
        # --adding caption to beginning of table--
        u.insert_line_in_file(r'\caption{' + caption + '}\n', 0, fullfilename)  # use -1 instead of 0 to append caption to end of table
        u.insert_line_in_file(subcaption + '\n', 4, fullfilename)  # use -1 instead of 0 to append caption to end of table
        print 'Written ' + fullfilename

    if not args.notex:
        print tab_tex
    else:
        print tab_txt

    print 'Finished!'
