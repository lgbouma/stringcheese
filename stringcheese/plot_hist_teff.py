import os
from glob import glob
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.io import ascii as ap_ascii
from numpy import array as nparr

from stringcheese import pipeline_utils as pu

def get_my_data(groupid=113, groupname='nan', classifxndate=20190907,
                is_field_star_comparison=False):
    #
    # for anything flagged manually as good (in other words, the rotation
    # period found just from the LS peak was OK), get the rotation period and
    # teff from the .results file.
    #

    if is_field_star_comparison:
        fs_str = 'field_star_comparison_'
    else:
        fs_str = ''

    classifixndir = (
        '../results/manual_classification/'
        '{}_{}group{}_name{}_classification/'.
        format(classifxndate, fs_str, groupid, groupname)
    )

    all_paths = glob(os.path.join(classifixndir,'*.png'))
    n_paths = len(all_paths)

    # get Teffs from all paths (irrespective of Prot measurement)
    if n_paths==0:
        raise AssertionError('expected some good sourceids')

    gd_sourceids = []
    for p in all_paths:
        if 'good' in p:
            gd_sourceids.append(np.int64(
                os.path.basename(p).split("_")[-1].replace('[good].png','')
            ))
        else:
            gd_sourceids.append(np.int64(
                os.path.basename(p).split("_")[-1].replace('.png','')
            ))

    # now get the Teff results
    datadir = (
        '../results/pkls_statuses_pages/{}group{}_name{}'.
        format(fs_str, groupid, groupname)
    )

    teffs = []
    for sourceid in gd_sourceids:

        status_file = os.path.join(datadir, str(sourceid),
                                   'GLS_rotation_period.results')

        if not os.path.exists(status_file):
            raise AssertionError('expected {} to exist'.format(status_file))

        d = pu.load_status(status_file)

        teffs.append(d['lomb-scargle']['teff'])

    df = pd.DataFrame(
        {'teff': nparr(teffs).astype(float)}
    )

    return df, n_paths


def plot_hist_teff(classifxndate=20190907, groupid=113, groupname='nan',
                   is_field_star_comparison=False):

    df, n_paths = get_my_data(
        groupid=groupid,
        groupname=groupname,
        classifxndate=classifxndate,
        is_field_star_comparison=is_field_star_comparison
    )


    ##########################################

    plt.close('all')
    f,ax = plt.subplots(figsize=(4,3))

    bins = np.arange(2500,
                     7500+500,
                     500)

    if is_field_star_comparison:
        label = 'Group {} field neighbors'.format(groupid)
    else:
        label = 'Group {}'.format(groupid)

    ax.hist(df['teff'], bins=bins, cumulative=False, color='black', fill=False,
            linewidth=0.5, label=label)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize('small')
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize('small')
    ax.set_xlabel('Effective temperature [K]')
    ax.set_ylabel('Number per bin')
    #ax.set_yscale('log')

    ax.legend(loc='best', fontsize='x-small')

    ax.set_xlim((7000,3000))
    ax.set_ylim((0,55))

    ax.get_yaxis().set_tick_params(which='both', direction='in',
                                   labelsize='small', top=True, right=True)
    ax.get_xaxis().set_tick_params(which='both', direction='in',
                                   labelsize='small', top=True, right=True)

    if is_field_star_comparison:
        fs_str = 'field_star_comparison_'
    else:
        fs_str = ''

    outpath = (
        '../results/hist_teff_{}group{}_name{}.png'.
        format(fs_str, groupid, groupname)
    )
    f.savefig(outpath, dpi=300, bbox_inches='tight')
    print('made {}'.format(outpath))


if __name__ == "__main__":

    plot_hist_teff(groupid=113, groupname='nan', classifxndate=20190910)
    plot_hist_teff(groupid=113, groupname='nan', classifxndate=20190910,
                   is_field_star_comparison=True)
