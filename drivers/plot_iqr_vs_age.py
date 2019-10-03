"""
Plot the IQR of your janky light curves vs KC19 reported age.
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

host = socket.gethostname()
if 'brik' in host:
    basedir = '/home/luke/local/stringcheese/'
else:
    basedir = np.nan

###############
# main driver #
###############

def get_data():

    source_df = pd.read_csv('../data/kounkel_table1_sourceinfo.csv')

    sel = (
        (source_df['Tmag_pred'] < 14)
        &
        (source_df['parallax'] > 5)
    )

    sdf = source_df[sel]

    outpath = '../data/Tmag_lt_14_plx_gt_5_concatenated_lc_statistics.csv'

    if os.path.exists(outpath):
        return pd.read_csv(outpath), sdf

    df2 = pd.read_csv('../data/string_table2.csv')
    sdf2 = df2[df2['parallax'] > 5]

    starinfos, varinfos = [], []
    for ix, r in sdf.iterrows():

        source_id = np.int64(r['source_id'])
        ra, dec = float(r['ra_x']), float(r['dec_x'])
        name = str(r['name'])
        group_id = str(r['group_id'])

        workingdir = os.path.join(basedir,
                                  'fits_pkls_results_pngs',
                                  'group{}_name{}'.format(group_id, name))
        workingdir = os.path.join(workingdir, str(source_id))

        resultfile = os.path.join(workingdir, 'GLS_rotation_period.results')

        if os.path.exists(resultfile):
            d = pu.load_status(resultfile)

            starinfos.append(dict(d['starinfo']))
            varinfos.append(dict(d['variability_info']))

        else:
            continue

    stardf = pd.DataFrame(starinfos)
    vardf = pd.DataFrame(varinfos)

    mdf = stardf.merge(
        vardf,
        on=list(np.intersect1d(stardf.columns, vardf.columns)),
        how="outer"
    )

    assert len(stardf) == len(mdf)

    mdf.to_csv(outpath, index=False)

    return mdf, sdf


def plot_iqr_vs_age(addnoise=0, showpctile=0, onlynamedgroups=0):

    if showpctile:
        assert addnoise
    if onlynamedgroups:
        assert addnoise and showpctile

    ykeys = ['rel_flux_iqr', 'rel_flux_15_to_85', 'rel_flux_5_to_95',
             'rel_flux_mad']

    for ykey in ykeys:

        plt.close('all')
        f,ax = plt.subplots(figsize=(6,3))

        df, sdf = get_data()

        if onlynamedgroups:
            df = df[~pd.isnull(df['name'])]
            sdf = sdf[~pd.isnull(sdf['name'])]

        xval = df['age']
        yval = df[ykey]

        if addnoise:
            np.random.seed(42)
            xval += np.random.uniform(-0.1, 0.1, size=len(df))

        ax.scatter(
            xval, yval,
            color='black', edgecolors='k',
            alpha=1, linewidths=0.4, zorder=2, s=3, marker='o',
        )

        if showpctile:

            bin_diff = 0.20
            bin_min, bin_max = 7.0, 9.6
            bin_left = np.arange(bin_min, bin_max-bin_diff, bin_diff)
            bin_right = bin_left+bin_diff

            pctiles = {50:[], 25:[], 75:[], 'agebin_middle':[]}

            for bl, br in zip(bin_left, bin_right):

                sel_y = yval[(xval>=bl) & (xval<br)]

                for pct in [25,50,75]:
                    pctiles[pct].append(np.nanpercentile(sel_y, pct))

                pctiles['agebin_middle'].append(np.mean([bl, br]))

            for ix, pct in enumerate([25,50,75]):
                ax.plot(pctiles['agebin_middle'], pctiles[pct],
                        label='{}%'.format(pct), color='C{}'.format(ix),
                        marker='s', markersize=4,
                        zorder=3)

        ax.set_xlim((7,9.5))
        ax.set_ylim((1e-4,2))

        ax.set_yscale('log')

        ax.set_xlabel('log10(age [yr])')
        ax.set_ylabel('{}'.format(ykey))

        if showpctile:
            ax.legend(loc='upper right')

        titlestr = (
            'T<14, plx>5mas, {} LCs ({} allsky)'.
            format(len(df), len(sdf))
        )
        if onlynamedgroups:
            titlestr = (
                'T<14, plx>5mas, ONLYNAMED {} LCs ({} allsky)'.
                format(len(df), len(sdf))
            )
        ax.set_title(titlestr, fontsize='x-small')

        ax.get_yaxis().set_tick_params(which='both', direction='in',
                                       labelsize='small', top=True, right=True)
        ax.get_xaxis().set_tick_params(which='both', direction='in',
                                       labelsize='small', top=True, right=True)

        outpath = (
            '../results/'
            'iqr_vs_age_yval_{}.png'.
            format(ykey)
        )
        if addnoise:
            outpath = (
                '../results/'
                'iqr_vs_age_yval_{}_addnoise.png'.
                format(ykey)
            )
        if showpctile:
            outpath = (
                '../results/'
                'iqr_vs_age_yval_{}_addnoise_showpctile.png'.
                format(ykey)
            )
        if onlynamedgroups:
            outpath = (
                '../results/'
                'iqr_vs_age_yval_{}_addnoise_showpctile_onlynamedgroups.png'.
                format(ykey)
            )

        f.savefig(outpath, dpi=300, bbox_inches='tight')
        print('made {}'.format(outpath))



if __name__=="__main__":

    plot_iqr_vs_age(addnoise=0, showpctile=0)
    plot_iqr_vs_age(addnoise=1, showpctile=0)
    plot_iqr_vs_age(addnoise=1, showpctile=1)
    plot_iqr_vs_age(addnoise=1, showpctile=1, onlynamedgroups=1)
