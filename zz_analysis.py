import numpy as np
import sys

sys.path.append('../')
import ayan.mage as m
import jrr

mage_mode = "released"
import pandas as pd
from matplotlib import pyplot as plt
import argparse as ap

parser = ap.ArgumentParser(description="Mage z analysis tool")
parser.add_argument('--allspec', dest='allspec', action='store_true')
parser.set_defaults(allspec=False)
parser.add_argument('--allesi', dest='allesi', action='store_true')
parser.set_defaults(allesi=False)
parser.add_argument('--stack', dest='stack', action='store_true')
parser.set_defaults(stack=False)
parser.add_argument('--hide', dest='hide', action='store_true')
parser.set_defaults(hide=False)
parser.add_argument('--keep', dest='keep', action='store_true')
parser.set_defaults(keep=False)
parser.add_argument('--silent', dest='silent', action='store_true')
parser.set_defaults(silent=False)
parser.add_argument('--noweight', dest='noweight', action='store_true')
parser.set_defaults(noweight=False)
parser.add_argument("--path")
parser.add_argument("--fname")
parser.add_argument("--lines")
parser.add_argument("--exclude")
parser.add_argument("--shortlabels")
parser.add_argument("--ew_thresh")
parser.add_argument("--nbins")
parser.add_argument("--weightby")
parser.set_defaults(showerr=False)
args, leftovers = parser.parse_known_args()
if args.lines is not None:
    lines = args.lines
else:
    lines = 'emission'
if args.fname is not None:
    fname = args.fname
else:
    fname = ''
if args.path is not None:
    path = args.path
else:
    path = '/Users/acharyya/Dropbox/MagE_atlas/Contrib/EWs/'  # directory to store the resulting output files
if args.stack:
    fname = 'all_byneb_stack_fitted_' + lines + '_linelist.txt'
if args.allesi:
    fname = 'allesi_fitted_' + lines + '_linelist.txt'
if args.allspec:
    fname = 'allspec_fitted_' + lines + '_linelist.txt'
if args.shortlabels is not None:
    shortlabels = [item for item in args.shortlabels.split(',')]
else:
    shortlabels = None
if args.exclude is not None:
    excludelabel = [item for item in args.exclude.split(',')]
else:
    excludelabel = []
if args.ew_thresh is not None:
    ew_thresh = float(args.ew_thresh)
else:
    ew_thresh = 3.
if args.nbins is not None:
    nbins = int(args.nbins)
else:
    nbins = 10
if args.weightby is not None:
    weightby = args.weightby
else:
    weightby = 'EWr_fit'

c = 3e5  # km/s
if not args.keep: plt.close('all')
fig = plt.figure(figsize=(10, 6))
ftemp = path + fname
line_table = pd.read_table(ftemp, delim_whitespace=True, comment="#")
line_table = line_table[~line_table['label'].isin(excludelabel)]  # take out 'bad' galaxies
if args.shortlabels is not None: line_table = line_table[line_table['label'].isin(shortlabels)]  # select few galaxies
line_table = line_table[line_table.EW_signi > ew_thresh]  # taking out undetected lines
labels = pd.unique(line_table['label'])
(specs) = jrr.mage.getlist_labels(mage_mode, labels,
                                  optional_file='/Users/acharyya/Desktop/mage_plot/Spectra/spectra-filenames-redshifts.txt')
for ii in range(0, len(labels)):
    short_table = line_table[line_table['label'].eq(labels[ii])]
    vel_offset = (short_table['zz'] - specs['z_neb'][ii]) * c
    if not args.noweight: vel_offset = np.multiply(vel_offset, short_table[weightby]) / np.sum(
        short_table[weightby])  # weighting by EWr_fit
    y, bin_edges = np.histogram(vel_offset, bins=nbins)
    y = np.pad(y, pad_width=1, mode='constant', constant_values=0)  # just to make histogram look nicer
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    bin_centers = np.pad(bin_centers, pad_width=1, mode='constant', constant_values=(
    bin_centers[0] - np.diff(bin_edges)[0],
    bin_centers[-1] + np.diff(bin_edges)[0]))  # just to amke histogram look nicer
    plt.plot(bin_centers, y, label=labels[ii])
    if not args.silent: print labels[ii], np.min(vel_offset), np.max(vel_offset), np.median(vel_offset)

plt.legend(bbox_to_anchor=(0.35, 0.9), bbox_transform=plt.gcf().transFigure)
plt.xlim(-100, 100)
# plt.ylim(0,5)
plt.ylabel('# of lines')
plt.xlabel('Velocity offset from zz_neb (km/s)')
t = 'Velocity offset of ' + lines + ' lines from fitted z and z_neb for ' + str(fname)
plt.title(t)
if not args.hide:
    plt.show(block=False)
else:
    plt.close()
