''' Equivalent width and flux fitter.  Given a spectrum, a continuum, a linelist, and a redshift,
fit a bunch of emission or absorption lines.
Started july 2016, Jane Rigby and Ayan Acharyya

README
Usage: python EW_fitter.py --<options>
<options> are:
--fcen FCEN ; FCEN = 1 if initial guess for gaussian centers are to be fixed, default 0(False)
--fcon FCON ; FCON = 1 if initial guess for continuum is to be fixed, default 1(True)
--dx DX ; DX = how many angstroms of wavelength should be plotted in one panel in the final plot, default = 300A
--only ONLY ; ONLY = 1 if final plot should have only those panels(patches of wavelength axis) where a line was fitted,
                else 0 to display all frames, default = 1
--vmax VMAX ; VMAX = in km/s, to set the maximum FWHM that can be fit as a line, default = 300km/s
--frame FRAME ; FRAME = index of the panel if you want to see specific panels, the panel indices can be found in top 
                left of each panel in the final plot, default plots all frames
--nbin NBIN ; NBIN = # of wavelength points to be binned together to calculate median and MAD binned fluxes and errors
--lines LINES ; LINES = 'emission' OR 'photospheric' depending on which lines you want to be fit, linelists can be 
                found in files labframe.shortlinelist_emission and labframe.shortlinelist_photospheric, respectively.
                default = emission
--fout FOUT ; FOUT = filename you want the output ASCII file to have, default = fitted_line_list.txt
--see LABEL ; LABEL = line label you want to check (string) if you want to see only the frame containing a relevant line \
                    and not the whole spectrum, default = None (shows all lines)
--keepprev ; boolean option, if present then doesn't kill the previous matplotlib plots
--silent ; boolean option, if present then doesn't print whole bunch of info on console while running
--mymask ; boolean option, if present then MODIFIES the masking around skylines (as opposed to jrr.mage.flag_skylines), 
            using the function flag_skylines in ayan.mage
--check ; boolean option, if present then prints at the end the root mean square deviation of 1 sigma error between EWs
        computed via fitting and via summing
--allspec ; boolean option, if present then run the code over all the individual spectra files
--allbyneb ; boolean option, if present then run the code over all the byneb stacked spectra files, automatically
            knows only to fit the emission line list and suitable fout (output file name)
--allbystars ; boolean option, if present then run the code over all the bystars stacked spectra files, automatically
            knows only to fit the photospheric line list and suitable fout (output file name)
--savepdf ; boolean option, if present then save the plots as pdfs
--noplot ; boolean option, if present then does not create any plot
--hide ; boolean option, if present then does not show plots at the end
--showbin ; boolean option, if present then show the binning (on either side of line/s used to calculate median and
            MAD errors) on the resultant plot
--fullmad ; boolean option, if present then calculate the MAD at every point on spectrum and add a column to dataframe
--showerr ; boolean option, if present then plots the Schneider precription EWs with errors
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

#-----------Main function starts------------------
path = '/Users/acharyya/Dropbox/MagE_atlas/Contrib/EWs/' #directory to store the resulting output files
parser = ap.ArgumentParser(description="Mage spectra fitting tool")
parser.add_argument("--shortlabel")
parser.add_argument("--fcen")
parser.add_argument("--fcon")
parser.add_argument("--dx")
parser.add_argument("--only")
parser.add_argument("--vmax")
parser.add_argument("--frame")
parser.add_argument("--nbin")
parser.add_argument("--lines")
parser.add_argument("--fout")
parser.add_argument("--see")
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
parser.add_argument('--allbyneb', dest='allbyneb', action='store_true')
parser.set_defaults(allbyneb=False)
parser.add_argument('--allbystars', dest='allbystars', action='store_true')
parser.set_defaults(allbystars=False)
parser.add_argument('--savepdf', dest='savepdf', action='store_true')
parser.set_defaults(savepdf=False)
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
args, leftovers = parser.parse_known_args()
if args.dx is not None:
    dx = float(args.dx)
else:
    dx = 310.
if args.only is not None:
   display_only_success = args.only
else:
    display_only_success = 1
if args.frame is not None:
    frame = args.frame
else:
    frame = None
if args.shortlabel is not None:
    labels = [item for item in args.shortlabel.split(',')]
else:
    labels = ['rcs0327-E']
if args.allbyneb:
    labels = [
    'magestack_byneb_standard',\
    'magestack_byneb_highZ',\
    'magestack_byneb_lowZ',\
    'magestack_byneb_midage8to16Myr',\
    'magestack_byneb_oldgt16Myr',\
    'magestack_byneb_younglt8Myr',\
    ]
    listname = 'emission'
    fout = path + 'all_byneb_stack_fitted_emission_linelist.txt'
if args.allbystars:
    labels = [
    'magestack_bystars_standard',\
    'magestack_bystars_highZ',\
    'magestack_bystars_lowZ',\
    'magestack_bystars_midage8to16Myr',\
    'magestack_bystars_oldgt16Myr',\
    'magestack_bystars_younglt8Myr',\
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
    fout = str(args.fout)+'.txt'
elif 'fout' not in locals():
    fout = 'fitted_line_list.txt'
if not args.keepprev:
    plt.close('all')

#-------------------------------------------------------------------------
(specs) = jrr.mage.getlist_labels(mage_mode, labels, optional_file='/Users/acharyya/Desktop/mage_plot/Spectra/spectra-filenames-redshifts.txt')
(spec_path, line_path) = jrr.mage.getpath(mage_mode)
line_table = pd.DataFrame(columns=['label', 'line_lab', 'obs_wav', 'rest_wave', 'type','EWr_fit','EWr_fit_u', 'EWr_sum', \
'EWr_sum_u', 'f_line','f_line_u', \
'EWr_Suplim', 'f_Suplim', 'fit_cont','fit_f','fit_cen', 'fit_cen_u', \
'fit_sig','zz','zz_u'])

for ii in range(0, len(specs)) :                  
    try:
        shortlabel     = specs['short_label'][ii]
        print 'Spectrum', (ii+1), 'of', len(specs),':', shortlabel #Debugging
        filename  = specs['filename'][ii]
        zz_sys = specs['z_syst'][ii] # from the new z_syst column in spectra_filename file
        zz_dic = {'EMISSION':specs['z_neb'][ii], 'FINESTR':specs['z_neb'][ii], 'PHOTOSPHERE': specs['z_stars'][ii] if specs['fl_st'][ii]==0 else specs['z_neb'][ii], 'ISM':specs['z_ISM'][ii], 'WIND':specs['z_ISM'][ii]}
        zz_err_dic = {'EMISSION':specs['sig_neb'][ii] if specs['fl_neb'][ii]==0 else specs['sig_ISM'][ii], 'FINESTR':specs['sig_neb'][ii] if specs['fl_neb'][ii]==0 else specs['sig_ISM'][ii], 'PHOTOSPHERE': specs['sig_st'][ii] if specs['fl_st'][ii]==0 else specs['sig_neb'][ii], 'ISM':specs['sig_ISM'][ii], 'WIND':specs['sig_ISM'][ii]}    
        #-----------reading spec in different formats----------------------------
        if 'esi' in shortlabel:
            specdir = '/Users/acharyya/Documents/esi_2016b/2016aug27_2x1/IRAF/reduced/'
            sp_orig = m.open_esi_spectrum(specdir+filename, getclean=True) # alt_infile= Put the filename of the stacked spectrum file here
            resoln = 4000.   # ESI spectral resoln for 1" slit
            dresoln = 40.       #         
        elif 'stack'in shortlabel:
            sp_orig = jrr.mage.open_stacked_spectrum(mage_mode, alt_infile=filename) # alt_infile= Put the filename of the stacked spectrum file here
            resoln = 3e5/200.   # vel resol of 200km/s from file /Users/acharyya/Dropbox/mage_atlas/Contrib/S99/stack-A-sb99-fit.txt
            dresoln = 40.       # 
        else:
            (sp_orig, resoln, dresoln)  = jrr.mage.open_spectrum(filename, zz_sys, mage_mode)
        #-----------fitting continuum----------------------------
        m.fit_autocont(sp_orig, zz_sys, line_path,filename)
        #-------masking sky lines-----------------
        if 'stack' not in shortlabel:
            if args.mymask:
                m.flag_skylines(sp_orig) #modified masking for skylines, as compared to jrr.mage.flag_skylines
            if 'esi' not in shortlabel:
                sp_orig = sp_orig[~sp_orig['badmask']].copy(deep=True)
        #-----calculating MAD error over entire spectrum--------
        if args.fullmad:
            m.calc_mad(sp_orig, resoln, 5)
            continue
        #------calculating the EW limits at every point following Schneider et al. 1993---------
        m.calc_schneider_EW(sp_orig, resoln, plotit=args.showerr)
        #m.makelist(line_path+'stacked.linelist') #required if you need to make a new labframe.shortlinelist file
        line_full = m.getlist('labframe.shortlinelist_'+listname, zz_dic, zz_err_dic)
        #------------Preparing to plot----------------------------------------
        xstart = max(np.min(line_full.wave) - 50.,np.min(sp_orig.wave))
        xlast = min(np.max(line_full.wave) + 50.,np.max(sp_orig.wave))
        if frame is None:
            n_arr = np.arange(int(np.ceil((xlast-xstart)/dx))).tolist()
        else:
            n_arr = [int(ar)-1 for ar in frame.split(',')] #Use this to display selected frame/s
        name = path + listname+'/'+shortlabel+'-'+listname+'_fit'
        if args.savepdf:
            pdf = PdfPages(name+'.pdf')
        #---------pre check in which frames lines are available if display_only_success = 1---------
        #---------------------------just a visualisation thing-----------------------------------
        if display_only_success:
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
        nrow = 4 #frames per page
        n_subarr = np.array_split(n_arr, int(np.ceil(float(len(n_arr)/nrow))))
        for ss in range(len(n_subarr)):
            n_arr = n_subarr[ss]
            n = len(n_arr)
            if not args.noplot:
                fig = plt.figure(figsize=(18+10/(n+1),(12 if n > 2 else n*3)))
                #fig = plt.figure(figsize=(14+8/(n+1),(9 if n > 2 else n*3)))
                plt.title(shortlabel + "  z=" + str(zz_sys)+'.\n Vertical lines legend: Blue=initial guess of center,'+\
                ' Red=fitted center, Green=no detection(< 3sigma), Black=unable to fit gaussian', y=1.02)
            for fc, jj in enumerate(n_arr):
                xmin = xstart + jj*dx
                xmax = min(xmin + dx, xlast)
                if not args.noplot:
                    ax1 = fig.add_subplot(n,1,fc+1)
                sp = sp_orig[sp_orig['wave'].between(xmin,xmax)]
                ymin = min(0,np.min(sp.fnu_u)*0.98) #setting ylimits for plotting, to little lower than minimum value of the error
                ymax = min(3,np.max(sp.fnu)*1.01) #little higher than maximum flux value
                try:
                    line = line_full[line_full['wave'].between(xmin*(1.+5./resoln), xmax*(1.-5./resoln))]
                except IndexError:
                    continue
                #------------Plot the results------------
                if not args.noplot:
                    try:
                        plt.step(sp.wave, sp.fnu, color='b')
                        plt.step(sp.wave, sp.fnu_u, color='gray')
                        plt.plot(sp.wave, sp.fnu_autocont, color='k')
                        if 'stack' not in shortlabel and 'esi' not in shortlabel:
                            plt.step(sp.wave, sp.fnu_cont, color='y')
                            plt.ylim(0, 1.2E-28)
                        else:
                            try:
                                plt.ylim(ymin,ymax)
                            except:
                                pass
                        plt.xlim(xmin, xmax)
                    except:
                        print 'failed at', shortlabel
                        break
                    plt.text(xmin+dx*0.05, ax1.get_ylim()[1]*0.8, 'Frame '+str(int(jj)+1))
                if not args.fullmad:
                    m.fit_some_EWs(line, sp, resoln, shortlabel, line_table, dresoln, sp_orig, args=args) #calling line fitter
                if not args.noplot:
                    ax2 = ax1.twiny()
                    ax2.set_xlim(ax1.get_xlim())       
                    ax2.set_xticklabels(np.round(np.divide(ax1.get_xticks(),(1.+zz_sys)),decimals=0))
                    labels2 = [item.get_text() for item in ax2.get_xticklabels()]
                    ax2.set_xticks(np.concatenate([ax2.get_xticks(), line.wave]))
                    ax2.set_xlim(ax1.get_xlim())       
                    ax2.set_xticklabels(np.concatenate([labels2,np.array(line.label.values).tolist()]), rotation = 45, ha='left', fontsize='small')
                    fig.subplots_adjust(hspace=0.7, top=0.94, bottom=0.05, left=0.05, right=0.95)
            if not args.hide:
                plt.show(block=False)
            if args.savepdf:
                pdf.savefig(fig)
        if args.savepdf:
            pdf.close()
        #fig.savefig(name+'.png')
    except:
        print 'Could not successfully complete', shortlabel
        continue
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
line_table['EW_signi']=3.*line_table['EWr_fit']/line_table['EWr_Suplim']
line_table['EW_signi']=line_table['EW_signi'].map('{:,.3f}'.format)
line_table['f_signi']=3.*line_table['f_line']/line_table['f_Suplim']
line_table['f_signi']=line_table['f_signi'].map('{:,.3f}'.format)
#----------------Saving dataframe to ASCII file--------------------------------------------
head = 'This file contains the measurements of lines in the MagE sample. Generated by EW_fitter.py.\n\
Columns are:\n\
 label:       shortlabel of the galaxy/knot\n\
 line_lab:    label of the line the code was asked to fit\n\
 obs_wav:     observed frame wavelength of the line (A)\n\
 rest_wave:   rest frame wavelength of the line (A)\n\
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
#----------Displaying part of dataframe if asked to---------------------------
line_table['f_SNR']=np.abs(line_table['f_line'])/line_table['f_line_u']
short_table = line_table[['line_lab','EWr_fit','EWr_fit_u','EWr_Suplim','EW_signi','f_line','f_line_u','f_SNR','f_Suplim','f_signi']]
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
