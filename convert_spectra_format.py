'''-----python routine to convert format of a given spectrum file to one that can be read by EW_fitter.py----
--------also creates a other-spectra-filenames-redshifts.txt file (analogous to the usual spectra-filenames-redshifts.txt)----
--------which is then fed into EW_fitter.py along with the spectrum file in the new-format-------------
--------by Ayan, March 2017----------'''
import numpy as np
import pandas as pd
import sys
import argparse as ap
import os
HOME = os.getenv('HOME')+'/'
#---------------------------------------------------------------- 
def writespec_txt(spec, fout, z='', z_u='', filename=''):
    head = '#1D spectrum for object file '+filename+'\n\
    #made using python convert_spectra_format.py by Ayan, Mar 2017\n\
    redshift z = '+str(z)+'\n\
    redshift uncertainty z_u = '+str(z_u)+'\n'
    np.savetxt(fout, [], header=head, comments='#')
    spec.to_csv(fout, sep='\t',mode ='a', index=None)
    print 'Written dataframe to file ', fout
#---------------------------------------------------------------- 
def writespeclist_txt(fout, specdir, shortlabel, z, z_u):
    if not os.path.exists(fout):
        head = '#Spectra: listed in the same format as individual ones for ease of being read by the code jrr.mage.getlist_labels()\n\
#Columns:\n\
#1)    spectrum filename  (path is /Users/acharyya/Dropbox/MagE_atlas/Contrib/ESI_Spectra/Temp/  )\n\
#2)   short object name\n\
#3-5)   Redshift of stars,  uncertainty, flag: same as nebular redshift\n\
#6-8) Redshift of NEBULAR GAS,  uncertainty, flag\n\
#9-11)Redshift of ISM, uncertainty, flag.  May be centroid since broad: same as nebular redshift\n\
#12  Notes\n\
origdir\t\tfilename\t\tshort_label\t\tz_stars\t\tsig_st\t\tfl_st\t\tz_neb\t\tsig_neb\t\tfl_neb\t\tz_ISM\t\tsig_ISM\t\tfl_ISM\t\tNOTES\n'
        np.savetxt(fout, [], header=head, comments='')
    with open(fout, 'a') as f:
        np.savetxt(f, np.column_stack([specdir, shortlabel+'.txt', shortlabel, z, z_u, 0., z, z_u, 0., z, z_u, 0., 'converted-format']),"%s")
    print 'Written speclist to file ', fout
#---------------------------------------------------------------- 
#---------------------------------------------------------------- 
parser = ap.ArgumentParser(description="Spectra format conversion tool")
parser.add_argument("--inpath")
parser.add_argument("--infile")
parser.add_argument("--wavecol")
parser.add_argument("--flamcol")
parser.add_argument("--flamucol")
parser.add_argument("--flamconst")
parser.add_argument("--flamcontcol")
parser.add_argument("--z")
parser.add_argument("--zu")
args, leftovers = parser.parse_known_args()

if args.inpath is not None:
    inpath = args.inpath
else:
    inpath = HOME+'Dropbox/MagE_atlas/Contrib/Temp/'
if args.infile is not None:
    infile = args.infile
else:
    infile = ''
    print 'Must provide an input filename as --infile <FILENAME-WITHOUT-PATH>. Exiting.'
    sys.exit()
if args.wavecol is not None:
    wavecol = args.wavecol
else:
    wavecol = 'wave'
if args.flamcol is not None:
    flamcol = args.flamcol
else:
    flamcol = 'flam'
if args.flamucol is not None:
    flamucol = args.flamucol
else:
    flamucol = 'flam_u'
if args.flamcontcol is not None:
    flamcontcol = args.flamcontcol
else:
    flamcontcol = 'flam_cont'
if args.flamconst is not None:
    flamconst = np.float(args.flamconst)
else:
    flamconst = 1.0
if args.z is not None:
    z = np.float(args.z)
else:
    z = 0.0
if args.zu is not None:
    zu = np.float(args.zu)
else:
    zu = 0.0

sp =  pd.read_table(inpath+infile, delim_whitespace=True, comment="#", header=0)
sp.rename(columns= {wavecol  : 'obswave'}, inplace=True)
sp.rename(columns= {flamcol  : 'flam'}, inplace=True)
sp.rename(columns= {flamucol  : 'flam_u'}, inplace=True)

sp['flam'] = sp['flam']*flamconst
sp['flam_u'] = sp['flam_u']*flamconst
sp['restwave'] = sp['obswave']/(1+z)
sp['badmask'] = False

if flamcontcol in sp:
    sp.rename(columns= {flamcontcol  : 'flam_cont'}, inplace=True)
    sp['flam_cont'] = sp['flam_cont']*flamconst

outfile = infile[:-4] + '_new-format.txt'
writespec_txt(sp, inpath+outfile, z=z, z_u=zu, filename=inpath+infile)
writespeclist_txt(inpath+'other-spectra-filenames-redshifts.txt', inpath, outfile[:-4], z, zu)
print 'All Done!'
