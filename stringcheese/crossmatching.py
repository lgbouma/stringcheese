import os
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.io import fits

from numpy import array as nparr
from collections import Counter

from astropy.table import Table
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import units as u, constants as const
from astropy.time import Time

from stringcheese.plotutils import savefig, format_ax

def match_kc19_to_cs16(max_sep=3*u.arcsec):
    """
    Match Kounkel & Covey 2019 to Cotten & Song 2016, which claims to be the
    most comprehensive list of IR excess stars in existence.
    """

    csvpath = '../data/kc19_to_cs16_within_{}_arcsec.csv'.format(max_sep.value)

    if not os.path.exists(csvpath):

        #
        # Get Cotten & Song 2016 IR excess sources.
        #
        raise NotImplementedError
        t = Vizier.query
        # J/ApJS/225/15/table5
        #TODO

        #FIXME #TODO
        #FIXME #TODO
        #FIXME #TODO
        #FIXME #TODO


def plot_kc19_to_2rxs_age_vs_parallax(max_sep=12*u.arcsec, overwrite=0):
    #
    # Compare KC2019 isochrone ages for things with matches to things without.
    #

    outpath = (
        '../results/crossmatching/kc19_to_2rxs_age_vs_parallax_sep{}.png'.
        format(max_sep.value)
    )

    if not os.path.exists(outpath) or overwrite:

        kc19_match_df, df = match_kc19_to_2rxs(max_sep=max_sep)

        plt.close('all')
        f, ax = plt.subplots(figsize=(4,3))

        label = 'KC2019: {} sources'.format(len(df))
        ax.scatter(df.age, 1/(1e-3*df.parallax), label=label, c='gray',
                   marker='o', zorder=2, rasterized=True, s=1, linewidths=0)

        label = (
            'KC2019 & ROSAT match < {}": {} sources'.
            format(int(max_sep.value), len(kc19_match_df))
        )

        ax.scatter(kc19_match_df.age, 1/(1e-3*kc19_match_df.parallax), label=label,
                   c='k', marker='*', zorder=3, rasterized=True, s=5, linewidths=0)

        ax.legend(loc='best', fontsize='xx-small')

        format_ax(ax)

        ax.set_xlabel('KC2019 log10(isochrone age)')
        ax.set_ylabel('dist [pc] = 1/(parallax [arcsec])')
        ax.set_yscale('log')

        f.tight_layout(pad=0.2)
        savefig(f, outpath)

    else:
        print('found {}'.format(outpath))





def match_kc19_to_2rxs(max_sep=12*u.arcsec):
    """
    2RXS has RA/dec of ROSAT sources in J2000.

    Cut on the ROSAT detection likelihood to be >=9. Following the
    Boller et al (2016) abstract, this is "conservative".
    Cut on the "Sflag' screening flag (0=good, not equal to 0 = screening flag
    set).

    https://ui.adsabs.harvard.edu/abs/2016A%26A...588A.103B/abstract
    """

    csvpath = '../data/kc19_to_2rxs_within_{}_arcsec.csv'.format(max_sep.value)

    if not os.path.exists(csvpath):

        #
        # Get 2RXS sources. From Boller+2016, they were taken with the PSPC between
        # June 1990 and Aug 1991.
        #
        with fits.open("../data/Boller_2016_2RXS_ROSAT_sources.fits") as hdul:
            t = Table(hdul[1].data)

        sel = (t['ExiML'] >= 9) & (t['Sflag'] == 0)

        t = t[sel]

        rosat_coord = SkyCoord(nparr(t['_RAJ2000']), nparr(t['_DEJ2000']),
                               frame='icrs', unit=(u.deg, u.deg),
                               equinox=Time(2000.0, format='jyear'),
                               obstime=Time(1991.0, format='jyear')
                               )

        #
        # Get Kounkel & Covey 2019 sources with gaia info crossmatched. Set an
        # annoyingly high 1" absolute tolerance on the ra/dec match between what I
        # got from gaiadr2.source table, and what was listed in KC19 table 1.  Use
        # the gaiadr2.source ra/dec (in this case, "_x", not "_y"). rtol is set as
        # what it needed to be to pass. Set the correct equinox and estimate of
        # observing time by the Gaia satellite, and include the proper motions.
        #
        df = pd.read_csv('../data/kounkel_table1_sourceinfo.csv')

        for k_x, k_y in zip(['ra_x', 'dec_x'],['ra_y', 'dec_y']):
            np.testing.assert_allclose(nparr(df[k_x]), nparr(df[k_y]),
                                       rtol=3e-2, atol=1/3600)

        kc19_coord = SkyCoord(nparr(df.ra_x), nparr(df.dec_x), frame='icrs',
                              unit=(u.deg, u.deg),
                              pm_ra_cosdec=nparr(df.pmra)*u.mas/u.yr,
                              pm_dec=nparr(df.pmdec)*u.mas/u.yr,
                              equinox=Time(2015.5, format='jyear'),
                              obstime=Time(2015.5, format='jyear')
                             )

        kc19_coord_rosat_epoch = kc19_coord.apply_space_motion(
            new_obstime=Time(1991.0, format='jyear')
        )

        #
        # For each ROSAT coordinate, get the nearest KC2019 match. Ensure they are
        # each in the same equinox system (J2000.0).  Examine the separation
        # distribution of the match.
        #

        idx, d2d, _ = rosat_coord.match_to_catalog_sky(
            kc19_coord_rosat_epoch.transform_to(rosat_coord)
        )

        sel = (d2d < 1000*u.arcsec) # ~1/3rd of a degree

        bins = np.logspace(-2, 3, 11)

        plt.close('all')
        f, ax = plt.subplots(figsize=(4,3))
        ax.hist(d2d[sel].to(u.arcsec).value, bins=bins, cumulative=False,
                color='black', fill=False, linewidth=0.5)

        format_ax(ax)

        ax.set_xlabel('2RXS to nearest KC19 source [arcsec]')
        ax.set_ylabel('Number per bin')
        ax.set_xscale('log')
        ax.set_yscale('log')

        f.tight_layout(pad=0.2)
        outpath = '../results/crossmatching/kc19_to_2rxs_separation.png'
        savefig(f, outpath)

        #
        # The above shows a hint of a bump at like ~10 arcsec. This is roughly the
        # level of positional uncertainty expected from ROSAT (Ayres+2004,
        # https://ui.adsabs.harvard.edu/abs/2004ApJ...608..957A/abstract).  6" was
        # the design specification accuracy. Seems like 12" (2x the expected
        # positional uncertainty) is a decent place to set the cut. It will give
        # ~1.5k matches.
        #

        sep_constraint = (d2d < max_sep)

        rosat_matches = t[sep_constraint]
        kc19_match_df = df.iloc[idx[sep_constraint]]

        kc19_match_df['has_2RXS_match'] = 1

        kc19_match_df['2RXS_match_dist_arcsec'] = (
            d2d[sep_constraint].to(u.arcsec).value
        )

        kc19_match_df['2RXS_name'] = rosat_matches['_2RXS']

        kc19_match_df['2RXS_ExiML'] = rosat_matches['ExiML']

        kc19_match_df.to_csv(csvpath, index=False)
        print('made {}'.format(csvpath))

    else:
        kc19_match_df = pd.read_csv(csvpath)
        df = pd.read_csv('../data/kounkel_table1_sourceinfo.csv')

    return kc19_match_df, df


def plot_kc19_to_2rxs_age_vs_parallax(max_sep=12*u.arcsec, overwrite=0):
    #
    # Compare KC2019 isochrone ages for things with matches to things without.
    #

    outpath = (
        '../results/crossmatching/kc19_to_2rxs_age_vs_parallax_sep{}.png'.
        format(max_sep.value)
    )

    if not os.path.exists(outpath) or overwrite:

        kc19_match_df, df = match_kc19_to_2rxs(max_sep=max_sep)

        plt.close('all')
        f, ax = plt.subplots(figsize=(4,3))

        label = 'KC2019: {} sources'.format(len(df))
        ax.scatter(df.age, 1/(1e-3*df.parallax), label=label, c='gray',
                   marker='o', zorder=2, rasterized=True, s=1, linewidths=0)

        label = (
            'KC2019 & ROSAT match < {}": {} sources'.
            format(int(max_sep.value), len(kc19_match_df))
        )

        ax.scatter(kc19_match_df.age, 1/(1e-3*kc19_match_df.parallax), label=label,
                   c='k', marker='*', zorder=3, rasterized=True, s=5, linewidths=0)

        ax.legend(loc='best', fontsize='xx-small')

        format_ax(ax)

        ax.set_xlabel('KC2019 log10(isochrone age)')
        ax.set_ylabel('dist [pc] = 1/(parallax [arcsec])')
        ax.set_yscale('log')

        f.tight_layout(pad=0.2)
        savefig(f, outpath)

    else:
        print('found {}'.format(outpath))


def get_kc19_to_2rxs_match_frac(max_sep=12*u.arcsec):

    kc19_match_df, df = match_kc19_to_2rxs(max_sep=max_sep)

    rxs_group_id_cnt = Counter(nparr(kc19_match_df.group_id))
    kc19_group_id_cnt = Counter(nparr(df.group_id))

    count_df = pd.DataFrame({
        'group_id': nparr(rxs_group_id_cnt.most_common())[:,0],
        'N_RXS_matches': nparr(rxs_group_id_cnt.most_common())[:,1]
    })

    N_KC19 = []
    for ix, r in count_df.iterrows():
        N_KC19.append(
            len(df[df.group_id == r.group_id])
        )

    count_df['N_KC19_members'] = N_KC19
    count_df['frac_of_KC19_with_RXS_matches'] = (
        count_df.N_RXS_matches / count_df.N_KC19_members
    )

    t2 = Table.read('../data/ajab339at2_mrt.txt', format='ascii.cds')
    t2_df = t2.to_pandas()
    sel_cols = ['Theia', 'OName', 'logage', 'Avmag', 'string', 'GLON', 'GLAT',
               'Par']
    t2_df = t2_df[sel_cols]
    t2_df['d_pc'] = 1/(t2_df.Par*1e-3)

    out_df = count_df.merge(t2_df, left_on='group_id', right_on='Theia',
                            how='left')

    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)
    print(42*'=')
    print('MAXIMUM SEPARATION: {}'.format(max_sep))
    print('SORTED BY NUMBER OF ROSAT SOURCES')
    print(out_df.head(n=40))

    out_df = out_df.sort_values(by='frac_of_KC19_with_RXS_matches',
                                ascending=False)

    outpath = (
        '../data/kc19_to_2rxs_match_frac_within_{}_arcsec.csv'.
        format(max_sep.value)
    )
    out_df.to_csv(outpath, index=False)
    print('SORTED BY FRACTION OF KC2019 WITH ROSAT MATCHES')
    print(out_df.head(n=40))
    print('saved as {}'.format(outpath))
