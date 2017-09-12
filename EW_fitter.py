''' Equivalent width and flux fitter.  Given a spectrum, a continuum, a linelist, and a redshift,
fit a bunch of emission or absorption lines.
Started july 2016, Jane Rigby and Ayan acharyya

README
Usage: python EW_fitter.py --<options>
<options> are:
--fcen FCEN ; FCEN = 1 if initial guess for gaussian centers are to be fixed, default 0(False)
--fcon FCON ; FCON = 1 if initial guess for continuum is to be fixed, default 1(True)
--dx DX ; DX = how many angstroms of wavelength should be plotted in one panel in the final plot, default = 300A
--resoln RESOLUTION ; RESOLUTION = instrument spectral resolution R, needs to be specified only for 
                    non-esi, non-mage and non-stack spectra (i.e. any spectrum file which does not have 'esi' or
                    'stack' in its name and also doesn't belong to the MagE sample), default R = 4000
--dresoln RESOLUTION_U ; RESOLUTION_U = uncertainty in instrument spectral resolution R, needs to be specified only for 
                    non-esi, non-mage and non-stack spectra (i.e. any spectrum file which does not have 'esi' or
                    'stack' in its name and also doesn't belong to the MagE sample), default RESOLUTION_U = 40
--only ONLY ; ONLY = 1 if final plot should have only those panels(patches of wavelength axis) where a line was fitted,
                else 0 to display all frames, default = 1
--vmax VMAX ; VMAX = in km/s, to set the maximum FWHM that can be fit as a line, default = 300km/s
--frame FRAME ; FRAME = index of the panel if you want to see specific panels, the panel indices can be found in top 
                left of each panel in the final plot, default plots all frames
--nbin NBIN ; NBIN = # of wavelength points to be binned together to calculate median and MAD binned fluxes and errors; NOT USED anymore
--ndlambda NDLAMBDA ; NDLAMBDA = # of spectral resolution elements to the left & right of the group of lines being fitted, for which the chunk of spectra will be fit
--lines LINES ; LINES = 'emission' OR 'photospheric' depending on which lines you want to be fit, linelists can be 
                found in files labframe.shortlinelist_emission and labframe.shortlinelist_photospheric, respectively.
                default = emission
--fout FOUT ; FOUT = filename (only file name, excluding path, path can be set by --path option (see below)) you want the output \
                ASCII file to have, you may or may not include .txt, the code automatically adds the extension if required, \
                if FOUT is explicitly specified, output files are saved in PATH+FOUT
                else, by default saved in current working directory as fitted_line_list.txt
--spec_list_file FILENAME ; FILENAME = filename of the spectra list file you want to be used, instead of the usual spectra_redshift.txt
--see LABEL ; LABEL = line label you want to check (string) if you want to see only the frame containing a relevant line \
                    and not the whole spectrum, default = None (shows all lines)
--keepprev ; boolean option, if present then doesn't kill the previous matplotlib plots
--silent ; boolean option, if present then doesn't print whole bunch of info on console while running
--mymask ; boolean option, if present then MODIFIES the masking around skylines (as opposed to jrr.mage.flag_skylines), 
            using the function flag_skylines in ayan.mage
--check ; boolean option, if present then prints at the end the root mean square deviation of 1 sigma error between EWs
        computed via fitting and via summing
--allspec ; boolean option, if present then run the code over all the individual spectra files
--stackbyneb ; boolean option, if present then run the code over all the byneb stacked spectra files, automatically
            knows only to fit the emission line list and suitable fout (output file name)
--stackbystars ; boolean option, if present then run the code over all the bystars stacked spectra files, automatically
            knows only to fit the photospheric line list and suitable fout (output file name)
--savepdf ; boolean option, if present then save the plots as pdfs
--noplot ; boolean option, if present then does not create any plot
--hide ; boolean option, if present then does not show plots at the end
--showbin ; boolean option, if present then show the binning (on either side of line/s used to calculate median and
            MAD errors) on the resultant plot
--fullmad ; boolean option, if present then calculate the MAD at every point on spectrum and add a column to dataframe
--showerr ; boolean option, if present then plots the Schneider precription EWs with errors
--extract LABEL; LABEL = line label you want to extract (string), to plot for papers
--spec_list_file FILENAME ; FILENAME = name of a file analogous to spectra-filenames-redshifts file, if you don't want to use \
                the latter. Used only for spectrum files which have 'new-format' in their names, i.e. that are not the usual mage/esi spectra \
                for usual mage spectra the code knows to use spectra-filenames-redshifts in Dropbox and \
                for esi spectra (reduced by Ayan) the code knows to use the special version of spectra-filenames-redshifts file in /Users/acharyya/Desktop/mage/.
--plotfnu ; boolean option, if present then plots theflux axis in fnu units instead of flambda
--nrow NROW ; NROW = maximum number of rows to be sub-plotted in the the resulting plots/PDFs, default=4
--ncol NCOL ; NCOL = maximum number of columns to be sub-plotted in the the resulting plots/PDFs, default=1
--usefnucont COLFNUCONT ; COLFNUCONT = name of the column carrying continuum (in fnu units) to be used for line fitting, instead of default jrr.mage.auto_fit_cont
--useflamcont COLFLAMCONT ; COLFLAMCONT = name of the column carrying continuum (in flam units) to be used for line fitting, instead of default jrr.mage.auto_fit_cont
--path PATH ; PATH = full path of directory where resulting dataframe and PDFs would be saved, default is ~/Dropbox/MagE_atlas/Contrib/EWs/.
--nodered ; boolean option, if present then does not perform dereddening corrections
--nofit ; boolean option, if present then plots only the spectra, does not fit the lines
--linelistpath LINELISTPATH ; LINELISTPATH=path to labframe.shortlinelist_emission or labframe.shortlinelist_photospheric
--linelistfile LINELISTFILE; LINELISFILE=only filename of the line-list-to-be-fit file, if the name is not of the format labframe.shortlinelist_
--debug ; boolean option, if present prints out debugging statements
--fitinterven ; boolean option, if present tries to fit intervenning absorption lines from labframe.shortlinelist.interven
--fixgroupz ; boolean option, if present tries to fit all lines in a group of neighboring lines with same redshift
'''
import sys
sys.path.append('../')
import ayan.mage as m
import jrr
import ayan.splot_util as s
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
from  astropy.io import ascii
from matplotlib import pyplot as plt
mage_mode = "released"
import argparse as ap
from matplotlib.backends.backend_pdf import PdfPages
import os
HOME = os.getenv('HOME')+'/'
import subprocess
#-----------Main function starts------------------
parser = ap.ArgumentParser(description="Mage spectra fitting tool")
parser.add_argument("--path")
parser.add_argument("--shortlabel")
parser.add_argument("--fcen")
parser.add_argument("--fcon")
parser.add_argument("--dx")
parser.add_argument("--resoln")
parser.add_argument("--dresoln")
parser.add_argument("--only")
parser.add_argument("--vmax")
parser.add_argument("--frame")
parser.add_argument("--nbin")
parser.add_argument("--ndlambda")
parser.add_argument("--lines")
parser.add_argument("--fout")
parser.add_argument("--see")
parser.add_argument("--extract")
parser.add_argument("--spec_list_file")
parser.add_argument("--nrow")
parser.add_argument("--ncol")
parser.add_argument("--usefnucont")
parser.add_argument("--useflamcont")
parser.add_argument("--linelistpath")
parser.add_argument("--linelistfile")
parser.add_argument('--keepprev', dest='keepprev', action='store_true')
parser.set_defaults(keepprev=False)
parser.add_argument('--silent', dest='silent', action='store_true')
parser.set_defaults(silent=False)
parser.add_argument('--mymask', dest='mymask', action='store_true')
parser.set_defaults(mymask=False)
parser.add_argument('--check', dest='check', action='store_true')
parser.set_defaults(check=False)
parser.add_argument('--allspec', dest='allspec', action='store_true')
parser.set_defaults(allspec=False)
parser.add_argument('--allesi', dest='allesi', action='store_true')
parser.set_defaults(allesi=False)
parser.add_argument('--stackbyneb', dest='stackbyneb', action='store_true')
parser.set_defaults(stackbyneb=False)
parser.add_argument('--stackbystars', dest='stackbystars', action='store_true')
parser.set_defaults(stackbystars=False)
parser.add_argument('--savepdf', dest='savepdf', action='store_true')
parser.set_defaults(savepdf=False)
parser.add_argument('--saveeps', dest='saveeps', action='store_true')
parser.set_defaults(saveeps=False)
parser.add_argument('--hide', dest='hide', action='store_true')
parser.set_defaults(hide=False)
parser.add_argument('--noplot', dest='noplot', action='store_true')
parser.set_defaults(noplot=False)
parser.add_argument('--showbin', dest='showbin', action='store_true')
parser.set_defaults(showbin=False)
parser.add_argument('--fullmad', dest='fullmad', action='store_true')
parser.set_defaults(fullmad=False)
parser.add_argument('--showerr', dest='showerr', action='store_true')
parser.set_defaults(showerr=False)
parser.add_argument('--plotfnu', dest='plotfnu', action='store_true')
parser.set_defaults(plotfnu=False)
parser.add_argument('--makelatex', dest='makelatex', action='store_true')
parser.set_defaults(makelatex=False)
parser.add_argument('--nodered', dest='nodered', action='store_true')
parser.set_defaults(nodered=False)
parser.add_argument('--nofit', dest='nofit', action='store_true')
parser.set_defaults(nofit=False)
parser.add_argument('--debug', dest='debug', action='store_true')
parser.set_defaults(debug=False)
parser.add_argument('--fitinterven', dest='fitinterven', action='store_true')
parser.set_defaults(fitinterven=False)
parser.add_argument('--fixgroupz', dest='fixgroupz', action='store_true')
parser.set_defaults(fixgroupz=False)
args, leftovers = parser.parse_known_args()
if args.path is not None:
    path = args.path
else:
    path = HOME+'Dropbox/MagE_atlas/Contrib/EWs/emission/' #directory to store the resulting output files
if args.dx is not None:
    dx = float(args.dx)
else:
    dx = 300.
if args.nrow is not None:
    nrow = int(args.nrow)
else:
    nrow = 4
if args.ncol is not None:
    ncol = int(args.ncol)
else:
    ncol = 1
if args.resoln is not None:
    resoln = float(args.resoln)
else:
    resoln = 4000.
if args.dresoln is not None:
    dresoln = float(args.dresoln)
else:
    dresoln = 40.
if args.only is not None:
   display_only_success = int(args.only)
else:
    display_only_success = 1
if args.frame is not None:
    frame = args.frame
else:
    frame = None
if args.usefnucont is not None:
    colfnucont = args.usefnucont
else:
    colfnucont = ''
if args.useflamcont is not None:
    colflamcont = args.useflamcont
else:
    colflamcont = ''
if args.linelistpath is not None:
    linelistpath = args.linelistpath
else:
    linelistpath = ''
if args.shortlabel is not None:
    labels = [item for item in args.shortlabel.split(',')]
else:
    labels = ['rcs0327-E']
if args.stackbyneb:
    labels = [
    'magestack_byneb_standard',\
    ]
    listname = 'emission'
    fout = path + 'all_byneb_stack_fitted_emission_linelist.txt'
if args.stackbystars:
    labels = [
    'magestack_bystars_standard',\
    ]
    listname = 'photospheric'
    fout = path + 'all_bystars_stack_fitted_photospheric_linelist.txt'
if args.lines is not None:
    listname = str(args.lines)
elif 'listname' not in locals():
    listname = 'emission'
if args.allesi:
    labels=[    
    's1723_center_a_esi',\
    's1723_center_b_esi',\
    's1723_arc_a_esi',\
    's1723_arc_b_esi',\
    's1723_counter_a_esi',\
    's1723_counter_b_esi',\
    's1723_side_a_esi',\
    's1723_side_b_esi',\
    's2340_a1_a_esi',\
    's2340_a2_a_esi',\
    's2340_a3_a_esi',\
    's2340_a4_a_esi',\
    's2340_a1_b_esi',\
    's2340_a2_b_esi',\
    's2340_a3_b_esi',\
    's2340_a4_b_esi',\
    'J1050_2_b_esi',\
    'J1050_3_b_esi',\
    'J1050_arc_a_esi',\
    'J1458_arc_a_esi',\
    'J1458_arc_b_esi',\
    ]
    fout = path + 'allesi_fitted_'+listname+'_linelist.txt'
if args.allspec:
    labels = [
    #'rcs0327-B',\
    'rcs0327-E',\
    #'rcs0327-Ehires',\
    #'rcs0327-Elores',\
    'rcs0327-G',\
    'rcs0327-U',\
    #'rcs0327-BDEFim1',\
    #'rcs0327-counterarc',\
    'S0004-0103',\
    #'S0004-0103alongslit',\
    #'S0004-0103otherPA',\
    'S0033+0242',\
    'S0108+0624',\
    'S0900+2234',\
    'S0957+0509',\
    'S1050+0017',\
    'Horseshoe',\
    'S1226+2152',\
    #'S1226+2152hires',\
    #'S1226+2152lores',\
    'S1429+1202',\
    'S1458-0023',\
    'S1527+0652',\
    #'S1527+0652-fnt',\
    'S2111-0114',\
    'Cosmic~Eye',\
    'S2243-0935',\
    ]
    fout = path + 'allspec_fitted_'+listname+'_linelist.txt'
if args.fout is not None:
    fout = path + str(args.fout)
    if not fout[-4:] == '.txt': fout += '.txt'
elif 'fout' not in locals():
    fout = 'fitted_line_list.txt'
if not args.keepprev:
    plt.close('all')
if args.linelistfile is not None:
    linelistfile = args.linelistfile
    listname = ''
else:
    linelistfile = 'labframe.shortlinelist_'+listname
#-------------------------------------------------------------------------
specs = jrr.mage.wrap_getlist(mage_mode, which_list="wcont", drop_s2243=True, optional_file=False, labels=(), zchoice='neb', MWdr=True) #loading default spectra_redshift.txt file
loaded_esi, loaded_stack, loaded_other = False, False, False
for thislabel in labels:
    if 'esi' in thislabel and not loaded_esi:
        optional_file = HOME+'Desktop/mage_plot/Spectra/esi-spectra-filenames-redshifts.txt'
        loaded_esi = True
        specs = specs.append(jrr.mage.wrap_getlist(mage_mode, which_list="wcont", drop_s2243=True, optional_file=optional_file, labels=(), zchoice='neb', MWdr=True)) #appending spectra_redshift.txt file for esi if there is any esi spectra to be fit
    if 'new-format' in thislabel and not loaded_other:
        optional_file = args.spec_list_file
        loaded_other = True
        specs = specs.append(jrr.mage.wrap_getlist(mage_mode, which_list="labels", drop_s2243=True, optional_file=optional_file, labels=[thislabel], zchoice='neb', MWdr=False)) #appending spectra_redshift.txt file for other spectra if there is any
    elif 'stack' in thislabel and not loaded_stack:
        optional_file = HOME+'Dropbox/mage_atlas/Spectra/stacked-spectra-filenames-redshifts.txt'
        loaded_stack = True
        specs = specs.append(jrr.mage.wrap_getlist(mage_mode, which_list="wcont", drop_s2243=True, optional_file=optional_file, labels=(), zchoice='neb', MWdr=True)) #appending spectra_redshift.txt file for stack if there is any stacked spectra to be fit

specs = specs[specs['short_label'].isin(labels)] #curtailing specs to only those spectra that are to be fit
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
line_table = pd.DataFrame(columns=['label', 'line_lab', 'obs_wav', 'rest_wave', 'type','EWr_fit','EWr_fit_u', 'EWr_sum', \
'EWr_sum_u', 'f_line','f_line_u', \
'EWr_Suplim', 'EW_signi', 'f_Suplim', 'f_signi', 'fit_cont','fit_f','fit_cen', 'fit_cen_u', \
'fit_sig','zz','zz_u'])

#m.makelist(line_path+'stacked.linelist') #required if you need to make a new labframe.shortlinelist file
#m.make_interven_list('/Users/acharyya/Dropbox/mage_atlas/Linelists/MINE/interven.lst', zz_interv=0.983) #required if you need to make a new labframe.shortlinelist_interven file
#sys.exit() #

for ii in range(0, len(specs)) :                  
    #try:
    shortlabel = specs['short_label'][ii]
    print 'Spectrum', (ii+1), 'of', len(specs),':', shortlabel #Debugging
    filename  = specs['filename'][ii]
    zz_sys = specs['z_syst'][ii] # from the new z_syst column in spectra_filename file
    zz_dic = {'EMISSION':specs['z_neb'][ii], 'FINESTR':specs['z_neb'][ii], 'PHOTOSPHERE': specs['z_stars'][ii] if specs['fl_st'][ii]==0 else specs['z_neb'][ii], 'ISM':specs['z_ISM'][ii], 'WIND':specs['z_ISM'][ii]}
    zz_err_dic = {'EMISSION':specs['sig_neb'][ii] if specs['fl_neb'][ii]==0 else specs['sig_ISM'][ii], 'FINESTR':specs['sig_neb'][ii] if specs['fl_neb'][ii]==0 else specs['sig_ISM'][ii], 'PHOTOSPHERE': specs['sig_st'][ii] if specs['fl_st'][ii]==0 else specs['sig_neb'][ii], 'ISM':specs['sig_ISM'][ii], 'WIND':specs['sig_ISM'][ii]}    
    #-----------reading spec in different formats----------------------------
    if 'esi' in shortlabel:
        specdir = HOME+'Documents/esi_2016b/2016aug27_2x1/IRAF/reduced/'
        sp_orig = m.open_esi_spectrum(specdir+filename, getclean=True)
        resoln = 4000.   # ESI spectral resoln for 1" slit
        dresoln = 40.    #         
    elif 'new-format' in shortlabel:
        specdir = os.path.expanduser(specs['origdir'][ii])
        sp_orig = m.open_esi_spectrum(specdir+filename, getclean=True) # open_esi_spectrum() can open any other spectra as well, if the other spectra has been converted to desired format
    elif 'stack'in shortlabel:
        (sp_orig, LL_dummy) = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=filename) # alt_infile= Put the filename of the stacked spectrum file here
        resoln = 3e5/200.   # vel resol of 200km/s from file /Users/acharyya/Dropbox/mage_atlas/Contrib/S99/stack-A-sb99-fit.txt
        dresoln = 40.       # 
    else:
        (sp_orig, resoln, dresoln, LL_dummy, zz_sys)  = jrr.mage.wrap_open_spectrum(shortlabel, mage_mode, addS99=False, zchoice='neb', MWdr=True)
    #-----------fitting continuum unless asked to use existing continuum column in spectrum dataframe----------------------------
    if args.usefnucont is not None:
        if colfnucont in sp_orig:
            sp_orig.drop('fnu_autocont',1,inplace=True) #removing any pre-existing fnu_autocont column, as now we're going to create a new fnu_autocont column
            sp_orig.rename(columns={colfnucont:'fnu_autocont'},inplace=True)
            sp_orig['flam_autocont'] = jrr.spec.fnu2flam(sp_orig.wave, sp_orig.fnu_autocont)
            print 'Using continuum values from', colfnucont, 'and NOT automatically fitting continuum.'
        else:
            print 'Column', colfnucont, 'does not exist in', filename, '. Exiting..'
            sys.exit()
    elif args.useflamcont is not None:
        if colflamcont in sp_orig:
            if 'flam_autocont' in sp_orig: sp_orig.drop('flam_autocont',1,inplace=True) #removing any pre-existing flam_autocont column, as now we're going to create a new flam_autocont column
            sp_orig.rename(columns={colflamcont:'flam_autocont'},inplace=True)
            sp_orig.fnu_autocont = jrr.spec.flam2fnu(sp_orig.wave, sp_orig.flam_autocont)   
            print 'Using continuum values from', colflamcont, 'and NOT automatically fitting continuum.'
        else:
            print 'Column', colflamcont, 'does not exist in', filename, '. Exiting..'
            sys.exit()
    else:
        m.fit_autocont(sp_orig, zz_sys, line_path,filename)
        if 'flam_autocont' not in sp_orig: sp_orig['flam_autocont'] = jrr.spec.fnu2flam(sp_orig['wave'], sp_orig.fnu_autocont) #added this because now I use the new jrr.spec.fit_autocont() which, unlike jrr.mage.fit_autocont(), does not create flam_autocont and rest_***_cont columns
        (dummy, rest_fnu_autocont, dummy) =  jrr.spec.convert2restframe(sp_orig['wave'], sp_orig.fnu_autocont, sp_orig.fnu_autocont,  zz_sys, 'fnu')
        sp_orig['rest_fnu_autocont'] = pd.Series(rest_fnu_autocont)
        (dummy, rest_flam, dummy)  = jrr.spec.convert2restframe(sp_orig['wave'], sp_orig.flam_autocont,  sp_orig.flam_u,  zz_sys, 'flam')    
        sp_orig['rest_flam_autocont'] = pd.Series(rest_flam)

    #-------masking sky lines-----------------
    if 'stack' not in shortlabel:
        if args.mymask:
            m.flag_skylines(sp_orig) #modified masking for skylines, as compared to jrr.mage.flag_skylines
        elif 'esi' not in shortlabel:
            sp_orig = sp_orig[~sp_orig['badmask']].copy(deep=True)
    #-----calculating MAD error over entire spectrum--------
    if args.fullmad:
        m.calc_mad(sp_orig, resoln, 5)
        continue
    #------calculating the EW limits at every point following Schneider et al. 1993---------
    m.calc_schneider_EW(sp_orig, resoln, plotit=args.showerr)
    print 'Loading line list to be fitted from '+linelistpath+linelistfile
    line_full = m.getlist(linelistpath+linelistfile, zz_dic, zz_err_dic) #reading the linelist to be used for fittting
    if os.path.exists(linelistpath+'labframe.shortlinelist_interven') and args.fitinterven:
        line_interven = m.get_interven_list(linelistpath+'labframe.shortlinelist_interven', zz_err = 0.0004) #reading the intervenning linelist
        print 'Including intervening linelist from', linelistpath+'labframe.shortlinelist_interven'
        line_full = pd.concat([line_full, line_interven], ignore_index=True) #appending the intervenning linelist to emission linelist
    line_full.sort_values('wave', inplace=True)
    #------------Preparing to plot----------------------------------------

    if args.extract is not None:
        lines_to_extract = [item for item in args.extract.split(',')]
        xmid = [line_full[line_full.label == item].wave.values[0] for item in lines_to_extract]
        xstart = xmid[0] -dx/2
        xlast = xmid[-1] + dx/2
    else:
        xstart = max(np.min(line_full.wave) - 50.,np.min(sp_orig.wave))
        xlast = min(np.max(line_full.wave) + 50.,np.max(sp_orig.wave))
    if args.extract is not None:
        n_arr = np.arange(len(lines_to_extract))
    elif frame is None:
        n_arr = np.arange(int(np.ceil((xlast-xstart)/dx))).tolist()
    else:
        n_arr = [int(ar)-1 for ar in frame.split(',')] #Use this to display selected frame/s
    if args.extract: name = path + shortlabel+'-individual-lines_fit-'+args.extract
    else: name = path + shortlabel+'-'+listname+'_fit'
    if args.savepdf:
        pdf = PdfPages(name+'.pdf')
    #---------pre check in which frames lines are available if display_only_success = 1---------
    #---------------------------just a visualisation thing-----------------------------------
    if display_only_success and args.extract is None:
        for jj in range(len(n_arr)):
            xmin = xstart + n_arr[jj]*dx
            xmax = min(xmin + dx, xlast)
            sp = sp_orig[sp_orig['wave'].between(xmin,xmax)]
            try:
                line = line_full[line_full['wave'].between(xmin*(1.+5./resoln), xmax*(1.-5./resoln))]
            except IndexError:
                continue
            if not len(line) > 0 or not line['wave'].between(np.min(sp.wave),np.max(sp.wave)).all():
                n_arr[jj] = np.ma.masked
            if args.see is not None and not any(args.see in x for x in line.label.values):
                n_arr[jj] = np.ma.masked
        n_arr = np.ma.compressed(n_arr)
        if len(n_arr) < 1:
            print 'None of the requested frames have any line in them. Try with a different frame number.'
            continue
    #------------------------------------------------------------
    n_subplot = nrow*ncol

    n_subarr = np.split(n_arr,np.arange(n_subplot,n_subplot*len(n_arr)/n_subplot+1,n_subplot)) #np.array_split(n_arr, int(np.ceil(len(n_arr)/float(nrow*ncol))))
    if len(n_subarr[-1]) == 0: n_subarr = n_subarr[:-1] #trimming last sub array if empty
    for ss in range(len(n_subarr)):
        n_arr = n_subarr[ss]
        n = len(n_arr)
        nrow_actual = int(np.ceil(n/float(ncol)))
        ncol_actual = min(n,ncol)
        if not args.noplot:
            if args.extract:
                fig = plt.figure(figsize=(8+1*ncol_actual,5+1.5*nrow_actual))
            elif args.saveeps:
                fig = plt.figure(figsize=(18+10/(n+1),(12 if n > 2 else n*6)))
            else:
                fig = plt.figure(figsize=(18+10/(n+1),(12 if n > 2 else n*3)))
            #fig = plt.figure(figsize=(14+8/(n+1),(9 if n > 2 else n*3)))
            if not args.see and not args.extract: plt.title(shortlabel + "  z=" + str(zz_sys)+'.\n Vertical lines legend: Blue=initial guess of center,'+\
            ' Red=fitted center, Black=no detection(upper limit)', y=1.02)
        for fc, jj in enumerate(n_arr):
            if args.extract is not None:
                xmid = line_full[line_full.label == lines_to_extract[fc]].wave.values[0]
                xmin = xmid - dx/2
                xmax = xmid + dx/2
            else:
                xmin = xstart + jj*dx
                xmax = min(xmin + dx, xlast)
            if not args.noplot:
                if args.extract is None: ax1 = fig.add_subplot(n,1,fc+1)
                else: ax1 = fig.add_subplot(nrow_actual,ncol_actual,fc+1)
            sp = sp_orig[sp_orig['wave'].between(xmin,xmax)]
            if not args.plotfnu:
                ymin = np.min(sp.flam_u)*0.98 #setting ylimits for plotting, to little lower than minimum value of the error
                ymax = np.max(sp.flam)*1.5 #little higher than maximum flux value
            else:
                ymin = np.min(sp.fnu_u)*0.98 #setting ylimits for plotting, to little lower than minimum value of the error
                ymax = np.max(sp.fnu)*1.5 #little higher than maximum flux value
            try:
                line = line_full[line_full['wave'].between(xmin*(1.+5./resoln), xmax*(1.-5./resoln))]
            except IndexError:
                continue
            #------------Plot the results------------
            if not args.noplot:
                try:
                    max_xticks = 3
                    tick_size = 10
                    if not args.plotfnu:
                        spec_color = 'k'
                        plt.step(sp.wave, sp.flam, color=spec_color)
                        plt.step(sp.wave, sp.flam_u, color='gray')
                        plt.plot(sp.wave, sp.flam_autocont, color='y')
                        if ('stack' not in shortlabel and not args.saveeps) and ('esi' not in shortlabel and 'new-format' not in shortlabel):
                            if 'flam_cont' in sp: plt.step(sp.wave, sp.flam_cont, color='b')
                            plt.ylim(0, 1.2E-17)
                        elif 'stack' in shortlabel:
                            plt.ylim(ymin,ymax)
                        else:
                            try:
                                plt.ylim(ymin,ymax)
                            except:
                                pass
                        if args.extract: 
                            plt.ylim(0,ymax)
                    else:
                        spec_color = 'k'
                        plt.step(sp.wave, sp.fnu, color=spec_color)
                        plt.step(sp.wave, sp.fnu_u, color='gray')
                        plt.plot(sp.wave, sp.fnu_autocont, color='y')
                        if ('stack' not in shortlabel and not args.saveeps) and ('esi' not in shortlabel and 'new-format' not in shortlabel):
                            if 'fnu_cont' in sp: plt.step(sp.wave, sp.fnu_cont, color='b')
                            plt.ylim(0, 1.2E-28)
                        else:
                            try:
                                plt.ylim(min(0,ymin),min(3,ymax))
                            except:
                                pass
                        if args.extract:
                            plt.ylim(0,ymax)
                    plt.xlim(xmin, xmax)
                    if args.extract:
                        ax1.set_xticks(np.round(np.arange(xmin+dx/(max_xticks+1),xmax,dx/(max_xticks+1))))
                        ax1.tick_params(axis='x', labelsize=tick_size)
                except:
                    print 'failed at', shortlabel
                    break
                if not args.extract and not args.saveeps: plt.text(xmin+dx*0.005, ax1.get_ylim()[1]*0.9, 'Frame '+str(int(jj)+1))
            if not args.fullmad and not args.nofit:
                m.fit_some_EWs(line, sp, resoln, shortlabel, line_table, dresoln, sp_orig, args=args) #calling line fitter
            if not args.noplot:
                ax2 = ax1.twiny()
                ax2.set_xlim(ax1.get_xlim())
                if args.extract: 
                    ax2.set_xticks(ax1.get_xticks()[1:])
                    ax2.tick_params(axis='x', labelsize=tick_size)
                elif args.saveeps:
                    ax2.set_xticks(ax2.get_xticks()[2:])
                ax2.set_xticklabels(np.round(ax2.get_xticks()/(1.+zz_sys)))
        if not args.noplot:
            if args.see:
                fig.subplots_adjust(hspace=0.7, top=0.8, bottom=0.15, left=0.05, right=0.98)
            elif args.extract:
                fig.subplots_adjust(hspace=0.4, top=0.90, bottom=0.10, left=0.10, right=0.95)
            else:
                fig.subplots_adjust(hspace=0.7, top=0.94, bottom=0.05, left=0.06, right=0.95)
            AA = '$\mathrm{\AA}$'
            fig.text(0.5, 0.02, r'Observed Wavelength ('+AA+')', ha='center')
            fig.text(0.5, 0.96, r'Rest-frame Wavelength ('+AA+')', ha='center')
            if not args.plotfnu: fig.text(0.02, 0.5, r'$f_\lambda$ (ergs/s/cm$^2$/'+AA+')', va='center', rotation='vertical')
            else: fig.text(0.02, 0.5, r'$f_\nu$ (ergs/s/cm$^2$/Hz)', va='center', rotation='vertical')
        if args.savepdf: pdf.savefig(fig)
        if args.saveeps: fig.savefig(name+'.eps')
        if not args.hide: 
            print 'Debug524: Attempting to plot..' #
            plt.show(block=False)
            print 'Debug526: Done plotting.' #
    if args.savepdf:
        pdf.close()
    '''
    except Exception, e:
        print 'Could not successfully complete', shortlabel, 'due to:'
        print e, 'in line', sys.exc_info()[-1].tb_lineno
        continue
    '''
#------------changing data types------------------------------
line_table.obs_wav = line_table.obs_wav.astype(np.float64)
line_table.rest_wave = line_table.rest_wave.astype(np.float64)
line_table.EWr_fit = line_table.EWr_fit.astype(np.float64)
line_table.EWr_fit_u = line_table.EWr_fit_u.astype(np.float64)
line_table.EWr_sum = line_table.EWr_sum.astype(np.float64)
line_table.EWr_sum_u = line_table.EWr_sum_u.astype(np.float64)
line_table.f_line = line_table.f_line.astype(np.float64)
line_table.f_line_u = line_table.f_line_u.astype(np.float64)
line_table.EWr_Suplim = line_table.EWr_Suplim.astype(np.float64)
line_table.f_Suplim = line_table.f_Suplim.astype(np.float64)
line_table.fit_cont = line_table.fit_cont.astype(np.float64)
line_table.fit_f = line_table.fit_f.astype(np.float64)
line_table.fit_cen = line_table.fit_cen.astype(np.float64)
line_table.fit_cen_u = line_table.fit_cen_u.astype(np.float64)
line_table.fit_sig = line_table.fit_sig.astype(np.float64)
line_table.zz = line_table.zz.astype(np.float64)
line_table.zz_u = line_table.zz_u.astype(np.float64)
line_table.EW_signi=line_table.EW_signi.astype(np.float64)
line_table.f_signi=line_table.f_signi.astype(np.float64)
if any(lab in shortlabel for lab in ['rcs0327-E','s1723']) and not args.nodered and args.see is None and args.extract is None:
    niter = int(1e5) #no. of Monter Carlo realisations to be done while dereddening
    print 'Extinction available for '+shortlabel+'. Performing redenning correction...'
    if 'rcs0327-E' in shortlabel:
        E, E_u = 0.4, 0.07 #extinction from Whitaker et al. 2014
    elif 's1723' in shortlabel:
        E, E_u = 0.3, 0.02 #mean value from those computed by AYAN_Code/calc_reddening.py  
    line_table['f_redcor'],line_table['f_redcor_u']=m.extinct(line_table.rest_wave, line_table.f_line, line_table.f_line_u, E, E_u, doMC=False, size=niter)
    line_table['f_Suplim_redcor'],dummy=m.extinct(line_table.rest_wave, line_table.f_Suplim, np.zeros(len(line_table)), E, E_u, doMC=False, size=niter)
    line_table['f_redcor']=line_table['f_redcor'].map('{:.3e}'.format)
    line_table['f_redcor_u']=line_table['f_redcor_u'].map('{:.3e}'.format)
    line_table['f_Suplim_redcor']=line_table['f_Suplim_redcor'].map('{:.3e}'.format)
#----------------Saving dataframe to ASCII file--------------------------------------------
head = 'This file contains the measurements of lines in the MagE sample. Generated by EW_fitter.py.\n\
Columns are:\n\
 label:       shortlabel of the galaxy/knot\n\
 line_lab:    label of the line the code was asked to fit\n\
 obs_wav:     observed frame vacuum wavelength of the line (A)\n\
 rest_wave:   rest frame vacuum wavelength of the line (A)\n\
 type:        is it emission or ism etc.\n\
 EWr_fit:     eqv width as calculated from the Gaussian fit to the line (A)\n\
 EWr_fit_u:  error in above qty. (A)\n\
 EWr_sum:     eqv width as calculated by summing the flux (A)\n\
 EWr_sum_u:   error in above qty. (A)\n\
 f_line:      flux i.e. area under Gaussian fit (erg/s/cm^2)\n\
 f_line_u:    error in above qty. (erg/s/cm^2)\n\
 EWr_Suplim:  3sigma upper limit for unresolved OR detection criteria for resolved EWs, as determined using \
              Schneider et al. 1993 prescription\n\
 EW_signi:    = 3*EWr_fit/EWr_Suplim. Probably use this for SIGNIFICANCE\n\
 f_Suplim:    3sigma upper limit for unresolved OR detection criteria for resolved FLUXES following above prescription\n\
 f_signi:     = 3*f_line/f_Suplim\n\
 fit_cont:    continuum from the Gaussian fit (observed frame, continuum normalised fit, hence dimensionless)\n\
 fit_f:       amplitude (i.e. height of Gaussian above the continuum) from the Gaussian fit (observed frame, continuum normalised fit, hence dimensionless)\n\
 fit_cen:     center from the Gaussian fit (observed frame, continuum normalised fit, units=Angstrom)\n\
 fit_cen_u:   error in above qty. (units of Angstrom)\n\
 fit_sig:     linewidth of the Gaussian fit (observed frame, continuum normalised fit, units=Angstrom)\n\
 zz:          Corrected redshift of this line, from the fitted center\n\
 zz_u:        error in above qty.\n\
 NaN means the code was unable to fit the line.\n\
'
np.savetxt(fout, [], header=head, comments='#')
line_table.to_csv(fout, sep='\t',mode ='a', index=None)
print 'Full table saved to', fout
#----------Displaying and saving part of dataframe if asked to---------------------------
line_table['f_SNR']=np.abs(line_table['f_line'])/line_table['f_line_u']
if shortlabel == 'rcs0327-E' and not args.nodered:
    short_table = line_table[['line_lab','rest_wave','EWr_fit','EWr_fit_u','EWr_Suplim','EW_signi','f_line','f_line_u','f_SNR','f_Suplim','f_signi','f_redcor','f_redcor_u','f_Suplim_redcor']]
elif listname == 'trial': #For correcting zz_sys: to check redshifts using only strong emission lines
    short_table = line_table[['label', 'line_lab','rest_wave','EWr_fit','EWr_fit_u','EW_signi','zz','zz_u']]
    short_table.EW_signi = short_table.EW_signi.astype(np.float64)
    short_table = short_table[short_table.EW_signi > 3.] #taking out undetected lines
    fout = HOME+'Dropbox/MagE_atlas/Contrib/EWs/zz_list_few_galx.txt'
    head = 'This file contains the measurements of lines in the MagE sample. Generated by EW_fitter.py.\n\
    Its for a few galaxies, using only bright emission lines, for the purpose of investigating the nebular redshifts.\n\
    Columns are:\n\
     label:       shortlabel of the galaxy/knot\n\
     line_lab:    label of the line the code was asked to fit\n\
     rest_wave:   rest frame vacuum wavelength of the line (A)\n\
     EWr_fit:     eqv width as calculated from the Gaussian fit to the line (A)\n\
     EWr_fit_u:  error in above qty. (A)\n\
     EW_signi:    significance of detection\n\
     zz:          Corrected redshift of this line, from the fitted center\n\
     zz_u:        error in above qty.\n\
     Only lines that had a EW_signi>3 (i.e. were detected) are included in this table and the analysis that follows.\n\
    '
    np.savetxt(fout, [], header=head, comments='#')
    short_table.to_csv(fout, sep='\t',mode ='a', index=None)
    with open(fout,'a') as myfile:
        myfile.write('label\t\tzz_neb\t\tmean zz\t\tmedian zz\t\tmedian zz_u\n')
        for ii in range(0, len(specs)) :                  
            shortlabel = specs['short_label'][ii]
            temp_table = short_table[short_table['label'].eq(shortlabel)]
            myfile.write(shortlabel+'\t'+ str(specs['z_neb'][ii])+'\t'+str(np.mean(temp_table['zz']))+'\t'+str(np.median(temp_table['zz']))+'\t'+str(np.median(temp_table['zz_u']))+'\n')
    print 'Created file', fout
elif 's1723' in shortlabel and args.see is None and args.extract is None:
    EW_signi_thresh, constant = 1., 1e-17
    quantities_to_extract = ['line_lab','rest_wave','f_line','f_line_u','EW_signi', 'f_SNR']
    if 'f_redcor' in line_table : quantities_to_extract += ['f_redcor', 'f_redcor_u']
    short_table = line_table[quantities_to_extract]
    short_table = short_table[(short_table['EW_signi'] >= EW_signi_thresh) & (short_table['f_SNR'] > 0.)]
    short_table.drop('EW_signi', axis=1, inplace=True)
    short_table.drop('f_SNR', axis=1, inplace=True)
    short_table.rename(columns={'line_lab':'ID'}, inplace=True)
    short_table.rename(columns={'f_line':'integrated_flux'}, inplace=True)
    short_table.rename(columns={'f_line_u':'uncert_flux'}, inplace=True)
    if 'f_redcor' in short_table:
        short_table.rename(columns={'f_redcor':'dered_flux'}, inplace=True)
        short_table.rename(columns={'f_redcor_u':'uncert_dered_flux'}, inplace=True)
    
    
    short_table['integrated_flux'] /= constant
    short_table['uncert_flux'] /= constant
    if 'dered_flux' in short_table:
        short_table['dered_flux'] = short_table['dered_flux'].astype(np.float64)/constant
        short_table['uncert_dered_flux'] = short_table['uncert_dered_flux'].astype(np.float64)/constant
    
    if 'ESI' in shortlabel: spectroscope = 'ESI'
    elif 'MMT' in shortlabel: spectroscope = 'MMT'
    short_table['spectrograph'] = spectroscope 
    short_table['Notes'] = '--'

    head = 'This file contains the measurements of lines in the MagE sample. Generated by EW_fitter.py.\n\
    Columns are:\n\
     line_lab:    label of the line the code was asked to fit\n\
     rest_wave:   rest frame vacuum wavelength of the line (A)\n\
     integrated_flux:      flux i.e. area under Gaussian fit (units of '+str(constant)+' erg/s/cm^2)\n\
     uncert_flux:    error in above qty. (units of '+str(constant)+' erg/s/cm^2)\n\
    '
    if 'dered_flux' in short_table:
        head += '\
        dered_flux:    dereddened flux with E(B-V) = '+str(E)+' +/- '+str(E_u)+' (units of '+str(constant)+' erg/s/cm^2)\n\
        uncert_dered_flux:    error in above qty. (units of '+str(constant)+' erg/s/cm^2)\n\
        '

    if 'MMT' in shortlabel: suffix = '.mmt' 
    elif 'ESI' in shortlabel: suffix = '.esi'
    fout = path+'s1723_measured_emission' + suffix
    np.savetxt(fout+'.txt', [], header=head, comments='#')
    short_table.to_csv(fout+'.txt', sep='\t',mode ='a', index=None)
    short_table.to_latex(fout+'.tex', index=None)
    print 'Short table saved to '+fout+'.txt and .tex'
    print 'Detections below EW_signi = '+str(EW_signi_thresh)+' have not been included.'
else:
    short_table = line_table[['line_lab','rest_wave','EWr_fit','EWr_fit_u','EWr_Suplim','EW_signi','f_line','f_line_u','f_SNR','f_Suplim','f_signi']]

print 'Some columns of the table are:'
print short_table
#-------------For correcting zz_ism------------------
if listname is 'ism':
    print 'label', 'median zz_ism', 'median zz_ism_u'
    for ii in range(0, len(specs)) :                  
        shortlabel = specs['short_label'][ii]
        short_table = line_table[line_table['label'].eq(shortlabel)]
        print shortlabel, np.median(short_table['zz']), np.median(short_table['zz_u'])
#----------------Sanity check: comparing 2 differently computed EWs------------------
if args.check:
    err_sum, es, n = 0., 0., 0
    for p in range(0,len(line_table)):
        EWr_fit = float(line_table.iloc[p].EWr_fit)
        EWr_sum = float(line_table.iloc[p].EWr_sum)
        EWr_sum_u = float(line_table.iloc[p].EWr_sum_u)
        #print line_table.iloc[p].line_lab, EWr_fit, EWr_sum, EWr_sum_u #Debugging
        if np.abs(EWr_sum) > EWr_sum_u and EWr_fit > 0.:
            err_sum += EWr_fit - EWr_sum
            es += (EWr_fit - EWr_sum)**2.
            n += 1
    if n > 0:
        print err_sum/n, np.sqrt(es/n)
    else:
        print 'No lines detected.'
#------------------------------------------End of main function------------------------------------------------
if args.makelatex:
    print 'List of spec', labels
    for lab in labels:
        print 'Converting df to tex for spec', lab
        subprocess.call(['python '+HOME+'Documents/writings/papers/abundance_pap/dftolatex.py --infile '+fout+\
    ' --outfile '+HOME+'Documents/writings/papers/magesample/fluxes/'+lab+'_fitted_detected --shortlabel '+lab],shell=True)