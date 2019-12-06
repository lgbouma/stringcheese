import os
from glob import glob
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.io import ascii as ap_ascii
from numpy import array as nparr
from astrobase.services.gaia import objectid_search
from mpl_toolkits.axes_grid1 import make_axes_locatable
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

    gd_paths = glob(os.path.join(classifixndir,'*good*.png'))

    gd_sourceids = [
        np.int64(os.path.basename(p).split("_")[-1].replace('[good].png',''))
        for p in gd_paths
    ]

    if len(gd_sourceids)==0:
        raise AssertionError('expected some good sourceids')

    # now get the LS results
    datadir = (
        '../results/pkls_statuses_pages/{}group{}_name{}'.
        format(fs_str, groupid, groupname)
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
        {'teff': teffs, 'prot': prots, 'source_id':gd_sourceids}
    )

    return df, n_paths


def plot_prot_vs_teff_singlegroup(classifxndate=20190907, groupid=113,
                                  groupname='nan',
                                  is_field_star_comparison=False,
                                  remove_outliers=False):

    praesepe_df, pleiades_df = get_reference_data()
    group_df, n_paths = get_my_data(
        groupid=groupid,
        groupname=groupname,
        classifxndate=classifxndate,
        is_field_star_comparison=is_field_star_comparison
    )
    kc19_df = pd.read_csv('../data/string_table2.csv')

    if remove_outliers:
        # remove outliers manually selected from glue (RVs or HR diagram
        # offset)
        _hr = pd.read_csv(
            '../data/kc19_group{}_table1_hr_diagram_weirdos.csv'.
            format(groupid)
        )
        _rv = pd.read_csv(
            '../data/kc19_group{}_table1_rv_weirdos.csv'.
            format(groupid)
        )
        outlier_df = pd.concat((_hr, _rv))

        common = group_df.merge(outlier_df, on='source_id', how='inner')

        print('before pruning RV and HR diagram outliers, had {} Prots'.
              format(len(group_df)))

        group_df = group_df[~group_df.source_id.isin(common.source_id)]

        print('after pruning RV and HR diagram outliers, had {} Prots'.
              format(len(group_df)))

    row = kc19_df[kc19_df['group_id'] == groupid]
    age = 10**(float(row['age']))
    age_gyr = age/(1e9)
    age_myr = age_gyr*1e3

    ##########################################

    plt.close('all')
    f,ax = plt.subplots(figsize=(4,3))

    ax.scatter(
        nparr(praesepe_df['Teff']), nparr(praesepe_df['Prot']),
        color='gray', edgecolors='k',
        alpha=1, linewidths=0.4, zorder=2, s=6, marker='s',
        label='Praesepe 670 Myr'
    )
    ax.scatter(
        nparr(pleiades_df['teff']), nparr(pleiades_df['prot']),
        color='whitesmoke', edgecolors='gray',
        alpha=1, linewidths=0.4, zorder=1, s=6, marker='X',
        label='Pleiades 120 Myr'
    )
    if is_field_star_comparison:
        label = 'Group {} field neighbors'.format(groupid)
    else:
        label = 'Group {}'.format(groupid)
    ax.scatter(
        nparr(group_df['teff']).astype(float), nparr(group_df['prot']).astype(float),
        color='darkorange', edgecolors='k',
        alpha=1, linewidths=0.4, zorder=3, s=9, marker='o',
        label=label
    )
    ax.legend(loc='best', fontsize='x-small')

    ax.set_xlim((7000,3000))
    ax.set_ylim((0,16.2))

    ax.set_xlabel('Effective temperature [K]')
    ax.set_ylabel('Rotation period [days]')

    titlestr = (
        'Name: {}. KC19 isochrone age: {:d} Myr.\n{}/{} ({:d}%) with Prot.'.
        format(groupname, int(age_myr), len(group_df), n_paths,
               int(100*len(group_df)/n_paths))
    )
    ax.set_title(titlestr, fontsize='x-small')

    ax.get_yaxis().set_tick_params(which='both', direction='in',
                                   labelsize='small', top=True, right=True)
    ax.get_xaxis().set_tick_params(which='both', direction='in',
                                   labelsize='small', top=True, right=True)

    if is_field_star_comparison:
        fs_str = '_field_star_comparison'
    else:
        fs_str = ''

    outpath = (
        '../results/prot_vs_teff/prot_vs_teff_group{}_name{}{}.png'.
        format(groupid, groupname, fs_str)
    )
    if remove_outliers:
        outpath = (
            '../results/prot_vs_teff/prot_vs_teff_{}group{}_name{}_outliers_removed.png'.
            format(fs_str, groupid, groupname)
        )
    f.savefig(outpath, dpi=300, bbox_inches='tight')
    print('made {}'.format(outpath))

    if groupid==113 and not is_field_star_comparison:
        sel = (
            (group_df['teff'].astype(float)>4000) &
            (group_df['prot'].astype(float)>8)
        )
        print('group 113 prot vs teff outliers are...')
        print(group_df[sel].source_id)


def plot_prot_vs_teff_allgroups(classifxndate=20190907, groupids=None,
                                groupnames=None, colorisBpmRp=0, xisBpmRp=None,
                                xisTeff=None):

    assert isinstance(groupids, list)
    assert isinstance(groupnames, list)
    if xisBpmRp:
        assert not xisTeff
    if xisTeff:
        assert not xisBpmRp

    praesepe_df, pleiades_df = get_reference_data()

    outpath = '../results/prot_vs_teff/prot_vs_teff_allgroups.csv'

    if not os.path.exists(outpath):
        group_df_list = []
        n_paths = []
        for groupid, groupname in zip(groupids, groupnames):

            i_group_df, i_n_paths = get_my_data(
                groupid=groupid,
                groupname=groupname,
                classifxndate=classifxndate,
                is_field_star_comparison=False
            )
            i_group_df['groupid'] = groupid
            i_group_df['groupname'] = groupname

            group_df_list.append(i_group_df)
            n_paths.append(i_n_paths)

        group_df = pd.concat(group_df_list)

        kc19_df = pd.read_csv('../data/string_table2.csv')

        group_df = group_df.merge(kc19_df[['group_id','age','parallax']],
                                  left_on='groupid', right_on='group_id',
                                  how='left')
        group_df = group_df.drop('group_id', axis=1)

        gaia_info_list = []
        for sid in nparr(group_df['source_id']):
            res = objectid_search(sid)
            gaia_info_list.append(pd.read_csv(res['result']))
        gaia_info_df = pd.concat(gaia_info_list)

        group_df = group_df.merge(gaia_info_df, on='source_id')
        group_df.to_csv(outpath, index=False)
        print('made {}'.format(outpath))

    else:
        group_df = pd.read_csv(outpath)

    ##########################################

    plt.close('all')
    f,ax = plt.subplots(figsize=(4,3))

    if xisTeff:
        ax.scatter(
            nparr(praesepe_df['Teff']), nparr(praesepe_df['Prot']),
            color='gray', edgecolors='k',
            alpha=1, linewidths=0.4, zorder=2, s=6, marker='s',
            label='Praesepe 670 Myr'
        )
        ax.scatter(
            nparr(pleiades_df['teff']), nparr(pleiades_df['prot']),
            color='whitesmoke', edgecolors='gray',
            alpha=1, linewidths=0.4, zorder=1, s=6, marker='X',
            label='Pleiades 120 Myr'
        )

    cval = nparr(group_df['age'])
    if xisTeff:
        xval = nparr(group_df['teff']).astype(float)
    if xisBpmRp:
        xval = nparr(group_df['phot_bp_mean_mag'] - group_df['phot_rp_mean_mag'])
    yval = nparr(group_df['prot']).astype(float)

    cm = ax.scatter(
        xval, yval,
        c=cval, edgecolors='k',
        cmap='plasma',
        alpha=1, linewidths=0.4, zorder=3, s=9, marker='o',
    )

    divider0 = make_axes_locatable(ax)
    cax0 = divider0.append_axes('right', size='5%', pad=0.05)
    cbar = f.colorbar(cm, ax=ax, cax=cax0)
    cbar.set_label('age')

    ax.legend(loc='best', fontsize='x-small')

    ax.set_ylim((0,16.2))

    if xisBpmRp:
        ax.set_xlabel('Bp-Rp')
    if xisTeff:
        ax.set_xlabel('Teff [K]')
        ax.set_xlim((8000,3000))
    ax.set_ylabel('Rotation period [days]')

    ax.get_yaxis().set_tick_params(which='both', direction='in',
                                   labelsize='small', top=True, right=True)
    ax.get_xaxis().set_tick_params(which='both', direction='in',
                                   labelsize='small', top=True, right=True)

    if xisBpmRp:
        outpath = '../results/prot_vs_teff/prot_vs_teff_allgroups_xisBpmRp.png'
    if xisTeff:
        outpath = '../results/prot_vs_teff/prot_vs_teff_allgroups_xisTeff.png'
    f.savefig(outpath, dpi=300, bbox_inches='tight')
    print('made {}'.format(outpath))



if __name__ == "__main__":

    group113_special = 0
    individual_groups = 0
    all_groups = 1

    if all_groups:

        allgroupids = [113, 208, 424,
                       676, 508, 509,
                       676, 786, 1089]
        groupnames = ['nan', 'Columba','nan',
                      'nan', 'nan', 'nan',
                      'nan', 'nan', 'nan']

        plot_prot_vs_teff_allgroups(groupids=allgroupids,
                                    groupnames=groupnames, xisBpmRp=1)

        plot_prot_vs_teff_allgroups(groupids=allgroupids,
                                    groupnames=groupnames, xisTeff=1)


    #
    # rotation period plots for each group with rotation period measurements
    #
    if individual_groups:
        newgroupids = [676, 424, 113, 1089, 508, 509, 786]

        for newgroupid in newgroupids:
            plot_prot_vs_teff_singlegroup(groupid=newgroupid, groupname='nan')

        plot_prot_vs_teff_singlegroup(groupid=208, groupname='Columba')


    #
    # group 113 specific processing
    #
    if group113_special:

        plot_prot_vs_teff_singlegroup(groupid=113, groupname='nan',
                                      classifxndate=20190910)

        plot_prot_vs_teff_singlegroup(groupid=113, groupname='nan',
                                      classifxndate=20190910,
                                      remove_outliers=True)

        plot_prot_vs_teff_singlegroup(groupid=113, groupname='nan',
                                      classifxndate=20190910,
                                      is_field_star_comparison=True)
