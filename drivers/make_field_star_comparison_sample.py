"""
For a given Kounkel+Covey 2019 groupid, construct a sample of stars in a
similar direct, at similar parallaxes from Gaia. Ensure none are in the KC19
list. This is your reference star sample. Measure the rotation periods of these
stars.

usage:
    python -u make_field_star_comparison_sample.py &> logs/comparison.log &
"""

###########
# imports #
###########
import os, socket
from glob import glob
import numpy as np, pandas as pd
from numpy import array as nparr

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import LombScargle
from astropy.io.votable import from_table, writeto, parse

from astroquery.gaia import Gaia
from astroquery.mast import Catalogs

from stringcheese.wrangling import get_fficutout as gfc
from stringcheese.wrangling import get_lc_given_fficutout as glgf
from stringcheese import verification_page as vp
from stringcheese import pipeline_utils as pu

host = socket.gethostname()
if 'phtess2' in host:
    basedir = '/home/lbouma/stringcheese'
    homedir = '/home/lbouma/'
elif 'brik' in host:
    basedir = '/home/luke/local/stringcheese'
    homedir = '/home/luke/'

###############
# main driver #
###############

def main(kc19_groupid=113, Tmag_cutoff=14, clean_gaia_cache=False):

    #
    # get info needed to query gaia for comparison stars
    #
    source_df = pd.read_csv('../data/kounkel_table1_sourceinfo.csv')
    sdf = source_df[(source_df['Tmag_pred'] < Tmag_cutoff) &
                    (source_df['group_id']==kc19_groupid)]
    n_sel_sources_in_group = len(sdf)

    df2 = pd.read_csv('../data/string_table2.csv')

    gdf = df2[df2['group_id']==kc19_groupid]

    group_coord = SkyCoord(float(gdf['l'])*u.deg,
                           float(gdf['b'])*u.deg, frame='galactic')
    ra = group_coord.icrs.ra
    dec = group_coord.icrs.dec
    plx_mas = float(gdf['parallax'])

    #
    # define relevant directories / paths
    #
    gaiadir = os.path.join(basedir, 'gaia_queries')
    if not os.path.exists(gaiadir):
        os.mkdir(gaiadir)

    outfile = os.path.join(
        gaiadir,
        'group{}_comparison_sample.xml.gz'.format(kc19_groupid)
    )

    #
    # run the gaia query. require the same cuts imposed by Kounkel & Covey 2019
    # on stellar quality. also require close on-sky (within 5 degrees of KC19
    # group position), and close in parallax space (within +/-20% of KC19
    # parallax).
    #
    if clean_gaia_cache and os.path.exists(outfile):
        os.remove(outfile)

    if not os.path.exists(outfile):

        Gaia.login(credentials_file=os.path.join(homedir, '.gaia_credentials'))

        jobstr = (
        '''
        SELECT TOP 2000 *
        FROM gaiadr2.gaia_source
        WHERE 1=CONTAINS(
          POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {ra:.8f}, {dec:.8f}, {sep_deg:.1f}))
        AND parallax < {plx_upper:.2f} AND parallax > {plx_lower:.2f}
        AND parallax > 1
        AND parallax_error < 0.1
        AND 1.0857/phot_g_mean_flux_over_error < 0.03
        AND astrometric_sigma5d_max < 0.3
        AND visibility_periods_used > 8
        AND (
                (astrometric_excess_noise < 1)
                OR
                (astrometric_excess_noise > 1 AND astrometric_excess_noise_sig < 2)
        )
        '''
        )

        query = jobstr.format(sep_deg=5.0, ra=ra.value, dec=dec.value,
                              plx_upper=1.2*plx_mas, plx_lower=0.8*plx_mas)


        if not os.path.exists(outfile):
            print(42*'-')
            print('launching\n{}'.format(query))
            print(42*'-')
            j = Gaia.launch_job(query=query,
                                verbose=True,
                                dump_to_file=True, output_file=outfile)

        Gaia.logout()

    vot = parse(outfile)
    tab = vot.get_first_table().to_table()
    field_df = tab.to_pandas()

    #
    # require the same Tmag cutoff for the nbhd stars. then randomly sample the
    # collection of stars.
    #

    Tmag_pred = (
        field_df['phot_g_mean_mag']
            - 0.00522555 * (field_df['phot_bp_mean_mag'] - field_df['phot_rp_mean_mag'])**3
            + 0.0891337 * (field_df['phot_bp_mean_mag'] - field_df['phot_rp_mean_mag'])**2
            - 0.633923 * (field_df['phot_bp_mean_mag'] - field_df['phot_rp_mean_mag'])
            + 0.0324473
    )

    field_df['Tmag_pred'] = Tmag_pred

    sfield_df = field_df[field_df['Tmag_pred'] < Tmag_cutoff]

    n_field = len(sfield_df)

    if n_sel_sources_in_group < n_field:
        errmsg = (
            'ngroup: {}. nfield: {}. plz tune gaia query to get more stars'.
            format(n_sel_sources_in_group, n_field)
        )
        raise AssertionError(errmsg)

    srfield_df = sfield_df.sample(n=n_sel_sources_in_group)

    #
    # now given the gaia ids, get the rotation periods
    #
    for ix, r in srfield_df.iterrows():

        source_id = np.int64(r['source_id'])
        ra, dec = float(r['ra']), float(r['dec'])

        #FIXME FIXME: 
        #FIXME FIXME: all the paths below need to be redefined
        #FIXME FIXME: 
        name = str(r['name'])
        group_id = str(r['group_id'])

        ##########################################
        # NOTE: change often

        #if int(group_id) not in subset4:
        #    continue
        #if int(group_id) not in close_middle_aged:
        #    continue

        #if int(group_id) != 113:
        #    continue
        #if source_id != 5220404075366707584:
        #    continue
        #if int(group_id) not in np.array(sdf2_str['group_id']).astype(int):
        #    # require that we only look at things Kounkel labelled as strings
        #    continue
        #if name != 'AB_Dor':
        #    continue
        #if source_id != 5579169050153502976:
        #    continue
        ##########################################

        c_obj = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')

        #
        # require that we are on-silicon. for year 1, this roughly means --
        # were are in southern ecliptic hemisphere
        #
        if c_obj.barycentrictrueecliptic.lat > 0*u.deg:
            print('group{}, {}: found in northern hemisphere. skip!'.
                  format(group_id, name))
            continue

        workingdir = os.path.join(basedir,
                                  'group{}_name{}'.format(group_id, name))
        if not os.path.exists(workingdir):
            os.mkdir(workingdir)
        workingdir = os.path.join(workingdir, str(source_id))
        if not os.path.exists(workingdir):
            os.mkdir(workingdir)

        outvppath = os.path.join(
            workingdir, 'verification_page_{}.png'.format(source_id)
        )
        if os.path.exists(outvppath):
            print('found {}, continue'.format(outvppath))
            continue

        #
        # if you already downloaded ffi cutouts for this object, dont get any
        # more. otherwise, get them
        #
        cutouts = glob(os.path.join(workingdir,'*.fits'))
        if len(cutouts)>=1:
            print('found {} cutouts in {}, skip'.format(
                len(cutouts), workingdir)
            )
        else:
            gfc.get_fficutout(c_obj, cutoutdir=workingdir)

        #
        # given the FFI cutouts, make simple light curves.
        #
        cutouts = glob(os.path.join(workingdir,'*.fits'))
        if len(cutouts)>=1:
            d = glgf.get_lc_given_fficutout(workingdir, cutouts, c_obj,
                                            return_pkl=False)
        else:
            d = np.nan
            print('WRN! did not find fficutout for {}'.format(workingdir))

        if not isinstance(d, dict):
            print('WRN! got bad light curve for {}. skipping.'.
                  format(workingdir))
            continue

        outpath = os.path.join(workingdir, 'GLS_rotation_period.results')

        #
        # do Lomb scargle w/ uniformly weighted points.
        #
        ls = LombScargle(d['time'], d['rel_flux'])
        period_min = 0.1
        period_max = np.min([
            0.9*(np.max(d['time']) - np.min(d['time'])), 16
        ])
        freq, power = ls.autopower(minimum_frequency=1/period_max,
                                   maximum_frequency=1/period_min)
        try:
            _ = power.max()
        except ValueError:
            print('WRN! got bad Lomb-Scargle for {}. skipping.'.
                  format(workingdir))
            continue

        ls_fap = ls.false_alarm_probability(power.max(), method='baluev')
        ls_period = 1/freq[np.argmax(power)]

        d['ls_fap'] = ls_fap
        d['ls_period'] = ls_period

        #
        # try to get TIC Teff. search TIC within 5 arcseconds, then take the
        # Gaia-ID match.  (removing sources with no gaia ID, which do exist in
        # TICv8.
        #
        radius = 5.0*u.arcsecond

        stars = Catalogs.query_region(
            "{} {}".format(float(c_obj.ra.value), float(c_obj.dec.value)),
            catalog="TIC",
            radius=radius
        )

        nbhr_source_ids = np.array(stars['GAIA'])

        stars = stars[ nbhr_source_ids != '' ]
        nbhr_source_ids = nbhr_source_ids[ nbhr_source_ids != '' ]

        sel = nbhr_source_ids.astype(int) == source_id

        if len(sel[sel])==1:
            star = stars[sel]
        else:
            raise NotImplementedError('did not get any TIC match. why?')

        teff = float(star['Teff'])
        if not isinstance(teff,float) and np.isfinite(teff):
            raise NotImplementedError('got nan TIC teff. what do?')


        #
        # make "check plot" analog for visual inspection
        #
        outd = {
            'ls_fap': d['ls_fap'],
            'ls_period': d['ls_period'],
            'source_id': source_id,
            'ra': ra,
            'dec': dec,
            'name': name,
            'group_id': group_id,
            'teff':teff
        }
        pu.save_status(outpath, 'lomb-scargle', outd)

        vp.generate_verification_page(d, ls, freq, power, cutouts, c_obj,
                                      outvppath, outd)

        #TODO: maybe only make plot if significant FAP?
        #TODO: and only if teff is within some range...


        #FIXME: need to verify your photometry implementation got the 0.5 pixel
        #convention right (it probably did not!)


if __name__ == "__main__":
    main()
