"""
Plot the Lomb-Scargle peak period of your janky light curves vs KC19 reported
age.
"""

###########
# imports #
###########
import os, socket, requests
from glob import glob
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from numpy import array as nparr

from scipy.stats import iqr

from stringcheese import pipeline_utils as pu
from plot_iqr_vs_age import get_data

host = socket.gethostname()
if 'brik' in host:
    basedir = '/home/luke/local/stringcheese/'
else:
    raise NotImplementedError

###############
# main driver #
###############

def plot_lspeakperiod_vs_age(addnoise=0, showpctile=0, onlynamedgroups=0,
                             onlybluestars=0):

    if showpctile:
        assert addnoise
    if onlynamedgroups:
        assert addnoise and showpctile
    if onlybluestars:
        assert addnoise and showpctile

    ykeys = ['ls_period']

    for ykey in ykeys:

        plt.close('all')
        f,ax = plt.subplots(figsize=(4,3))

        df, sdf = get_data()

        if onlynamedgroups:
            df = df[~pd.isnull(df['name'])]
            sdf = sdf[~pd.isnull(sdf['name'])]

        if onlybluestars:
            mag_cutoff = 1.8
            sel = (df['phot_bp_mean_mag'] - df['phot_rp_mean_mag']) < mag_cutoff
            df = df[sel]
            sel = (sdf['phot_bp_mean_mag'] - sdf['phot_rp_mean_mag']) < mag_cutoff
            sdf = sdf[sel]

        xval = df['age']
        yval = df[ykey]

        if addnoise:
            np.random.seed(42)
            xval += np.random.uniform(-0.1, 0.1, size=len(df))

        cutoff_days = 10.5
        ax.scatter(
            xval[yval<cutoff_days], yval[yval<cutoff_days],
            color='black', edgecolors='none',
            alpha=1, linewidths=0.4, zorder=2, s=3, marker='o',
        )
        ax.scatter(
            xval[yval>cutoff_days], yval[yval>cutoff_days],
            color='gray', edgecolors='none',
            alpha=1, linewidths=0.4, zorder=2, s=3, marker='o',
        )

        if showpctile:

            bin_diff = 0.25
            bin_min, bin_max = 7.0, 9.5
            bin_left = np.arange(bin_min, bin_max-bin_diff, bin_diff)
            bin_right = bin_left+bin_diff

            pctiles = {50:[], 25:[], 75:[], 'agebin_middle':[]}

            for bl, br in zip(bin_left, bin_right):

                sel_y = yval[(xval>=bl) & (xval<br) & (yval<cutoff_days)]

                for pct in [25,50,75]:
                    pctiles[pct].append(np.nanpercentile(sel_y, pct))

                pctiles['agebin_middle'].append(np.mean([bl, br]))

            for ix, pct in enumerate([25,50,75]):
                ax.plot(pctiles['agebin_middle'], pctiles[pct],
                        label='{}%'.format(pct), color='C{}'.format(ix),
                        marker='s', markersize=4,
                        zorder=3)

        ax.set_xlim((7,9.8))
        ax.set_ylim((0.1, 17.5))

        #ax.set_yscale('log')

        ax.set_xlabel('log10(age [yr])')
        ax.set_ylabel('LS peak period [day]')

        if showpctile:
            ax.legend(loc='lower right', fontsize='x-small')

        titlestr = (
            'T<14, plx>5mas, {} LCs ({} allsky)'.
            format(len(df), len(sdf))
        )
        if onlynamedgroups:
            titlestr = (
                'T<14, plx>5mas, ONLYNAMED {} LCs ({} allsky)'.
                format(len(df), len(sdf))
            )
        if onlybluestars:
            titlestr = (
                'T<14, plx>5mas, Bp-Rp<{} {} LCs ({} allsky)'.
                format(mag_cutoff, len(df), len(sdf))
            )
        ax.set_title(titlestr, fontsize='x-small')

        ax.get_yaxis().set_tick_params(which='both', direction='in',
                                       labelsize='small', top=True, right=True)
        ax.get_xaxis().set_tick_params(which='both', direction='in',
                                       labelsize='small', top=True, right=True)

        outpath = (
            '../results/'
            'lspeakperiod_vs_age.png'
        )
        if addnoise:
            outpath = (
                '../results/'
                'lspeakperiod_vs_age_addnoise.png'
            )
        if showpctile:
            outpath = (
                '../results/'
                'lspeakperiod_vs_age_addnoise_showpctile.png'
            )
        if onlynamedgroups:
            outpath = (
                '../results/'
                'lspeakperiod_vs_age_addnoise_showpctile_onlynamedgroups.png'
            )
        if onlybluestars:
            outpath = (
                '../results/'
                'lspeakperiod_vs_age_addnoise_showpctile_onlybluestars.png'
            )

        f.savefig(outpath, dpi=300, bbox_inches='tight')
        print('made {}'.format(outpath))



if __name__=="__main__":

    plot_lspeakperiod_vs_age(addnoise=1, showpctile=1, onlybluestars=1)

    plot_lspeakperiod_vs_age(addnoise=1, showpctile=1)
    plot_lspeakperiod_vs_age(addnoise=0, showpctile=0)
    plot_lspeakperiod_vs_age(addnoise=1, showpctile=0)
    plot_lspeakperiod_vs_age(addnoise=1, showpctile=1)
    plot_lspeakperiod_vs_age(addnoise=1, showpctile=1, onlynamedgroups=1)
