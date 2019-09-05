import os
from glob import glob
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.io import ascii as ap_ascii
from numpy import array as nparr

from stringcheese import pipeline_utils as pu

def get_reference_data():
    #
    # Table 4 of Douglas, Curtis et al (2019) Praesepe rotation periods, teffs.
    #
    # Quoting Curtis+2019: Prot for 743 members were amassed from the literature
    # and measured from K2 Campaign 5 light curves by Douglas et al. (2017).
    # Douglas et al. (2019) crossmatched this list with DR2 and filtered out stars
    # that failed membership, multiplicity, and data quality criteria, leaving us
    # with 359 single star members.
    #
    # And following Douglas+2019 table 4 caption for the flags...
    #
    praesepe_tab = ap_ascii.read("../data/apjab2468t4_mrt.txt")

    sel = (
        (praesepe_tab['SFlag'] == 'YYYYY')
        |
        (praesepe_tab['SFlag'] == 'YYY-Y')
    )

    praesepe_tab = praesepe_tab[sel]

    assert len(praesepe_tab) == 359

    praesepe_df = praesepe_tab.to_pandas()

    #
    # Figure 4 of Curtis+2019. Has a "gold sample" of Pleiades members that
    # involved some crossmatching, and removal of binaries. I used WebPlotDigitizer
    # to measure the rotation periods from that figure (rather than reproduce the
    # actual procedure Jason discusses in the text.)
    #
    pleiades_df = pd.read_csv('../data/pleaides_prot_vs_teff.csv')

    return praesepe_df, pleiades_df


def get_my_data(groupid=113, groupname='nan'):
    #
    # for anything flagged manually as good (in other words, the rotation
    # period found just from the LS peak was OK), get the rotation period and
    # teff from the .results file.
    #

    # NOTE: might need to fix this date-based naming convention...
    classifixndir = (
        '../results/manual_classification/'
        '20190905_group{}_name{}_classification/'.format(groupid, groupname)
    )

    gd_paths = glob(os.path.join(classifixndir,'*good*.png'))

    gd_sourceids = [
        np.int64(os.path.basename(p).split("_")[-1].replace('[good].png',''))
        for p in gd_paths
    ]

    if len(gd_sourceids)==0:
        raise AssertionError('expected some good sourceids')

    # now get the LS results
    datadir = (
        '../results/pkls_statuses_pages/group{}_name{}'.
        format(groupid, groupname)
    )

    prots, teffs = [], []
    for sourceid in gd_sourceids:

        status_file = os.path.join(datadir, str(sourceid),
                                   'GLS_rotation_period.results')

        if not os.path.exists(status_file):
            raise AssertionError('expected {} to exist'.format(status_file))

        d = pu.load_status(status_file)

        teffs.append(d['lomb-scargle']['teff'])
        prots.append(d['lomb-scargle']['ls_period'])

    df = pd.DataFrame(
        {'teff': teffs, 'prot': prots}
    )

    return df


def plot_prot_vs_teff(groupid=113, groupname='nan'):

    praesepe_df, pleiades_df = get_reference_data()
    group113_df = get_my_data(groupid=groupid, groupname=groupname)

    plt.close('all')
    f,ax = plt.subplots(figsize=(4,3))

    ax.scatter(
        nparr(praesepe_df['Teff']), nparr(praesepe_df['Prot']),
        color='C0', edgecolors='k',
        alpha=1, linewidths=0.4, zorder=2, s=9, marker='o',
        label='Praesepe 670 Myr'
    )
    ax.scatter(
        nparr(pleiades_df['teff']), nparr(pleiades_df['prot']),
        color='C1', edgecolors='k',
        alpha=1, linewidths=0.4, zorder=1, s=9, marker='o',
        label='Pleiades 120 Myr'
    )
    ax.scatter(
        nparr(group113_df['teff']).astype(float), nparr(group113_df['prot']).astype(float),
        color='C2', edgecolors='k',
        alpha=1, linewidths=0.4, zorder=3, s=9, marker='o',
        label='Group {}'.format(groupid)
    )
    ax.legend(loc='best', fontsize='x-small')

    ax.set_xlim((7000,3000))
    ax.set_ylim((0,16.2))

    ax.set_xlabel('Effective temperature [K]')
    ax.set_ylabel('Rotation period [days]')

    ax.get_yaxis().set_tick_params(which='both', direction='in',
                                   labelsize='small', top=True, right=True)
    ax.get_xaxis().set_tick_params(which='both', direction='in',
                                   labelsize='small', top=True, right=True)

    outpath = (
        '../results/prot_vs_teff_group{}_name{}.png'.
        format(groupid, groupname)
    )
    f.savefig(outpath, dpi=300, bbox_inches='tight')
    print('made {}'.format(outpath))



if __name__ == "__main__":
    plot_prot_vs_teff()
