import pandas as pd
from matplotlib import pyplot as plt

gal = 'rcs0327-knotE'
fig = plt.figure()
# ftemp = '/Users/acharyya/Documents/writings/papers/abundance_pap/lineflux_restUV.txt'
ftemp = '/Users/acharyya/Work/astro/mageproject/AYAN_Code/fitted_line_list.txt'
line_table = pd.read_table(ftemp, delim_whitespace=True, comment="#")

plt.errorbar(line_table.rest_wave, line_table.f_line, yerr=line_table.f_line_u, xerr=0, c='b', fmt='o', lw=0.5,
             label='obs')
plt.errorbar(line_table.rest_wave, line_table.f_redcor, yerr=line_table.f_redcor_u, xerr=0, c='g', fmt='o', lw=0.5,
             label='redcor')
plt.plot(line_table.rest_wave, line_table.f_Suplim, c='b', linestyle='--', label='Suplim')
plt.plot(line_table.rest_wave, line_table.f_Suplim_redcor, c='g', linestyle='--', label='Suplim_redcor')
'''
plt.plot(line_table.rest_wave,line_table.f_line/line_table.f_line_u,c='b',linestyle='--',label='pre-redcor-SNR')
plt.plot(line_table.rest_wave,line_table.f_redcor/line_table.f_redcor_u,c='g',linestyle='--',label='post-redcor-SNR')
'''
plt.legend()
plt.ylim(-1e-16, 1e-16)
plt.ylabel('flux')
plt.xlabel('rest wave')
t = 'Redenning_correction_effects_for_' + gal
plt.title(t)
fig.savefig(t + '.png')
plt.show(block=False)
