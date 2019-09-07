"""
Measure rotation periods for bright, young stars in strings. In other words:
T < 13 7.5 < log(age) < 8.5   (like 30-300 Myr. no younger, b/c reddening.)

usage:
    python -u generate_reports.py &> logs/make_ABDor_pages.log &
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

from astroquery.gaia import Gaia
from astroquery.mast import Catalogs

from stringcheese.wrangling import get_fficutout as gfc
from stringcheese.wrangling import get_lc_given_fficutout as glgf
from stringcheese import verification_page as vp
from stringcheese import pipeline_utils as pu

host = socket.gethostname()
if 'phtess2' in host:
    basedir = '/home/lbouma/stringcheese_cutouts/'
elif 'brik' in host:
    basedir = '/home/luke/stringcheese_cutouts/'

###############
# main driver #
###############

def main():
    source_df = pd.read_csv('../data/kounkel_table1_sourceinfo.csv')

    sel = (
        (source_df['Tmag_pred'] < 13)
        &
        (source_df['age'] < 9.0)
        &
        (source_df['age'] > 7.5)
    )

    sdf = source_df[sel]
    print(
        'after making cuts on T<13, log(age) from 7.5-8.5, got {} stars, {} groups'.
        format(len(sdf), len(np.unique(sdf['group_id'])))
    )

    df2 = pd.read_csv('../data/string_table2.csv')
    sdf2 = df2[(df2['age']>7.5) & (df2['age']<8.5)]
    sdf2_str = sdf2[sdf2['string']=='y']

    close_middle_aged = [
        1005, 208, 506, 424, 676, 507, 594, 209, 425, 595, 677, 905, 45, 7, 63
    ]

    # now given the gaia ids, get the rotation periods
    for ix, r in sdf.iterrows():

        source_id = np.int64(r['source_id'])
        ra, dec = float(r['ra_x']), float(r['dec_x'])
        name = str(r['name'])
        group_id = str(r['group_id'])

        ##########################################
        # NOTE: change often
        # require that we ony look at close, middle-aged objects as flagged
        # from glue visualizations

        # if int(group_id) not in close_middle_aged:
        #     continue

        if int(group_id) != 113:
            continue
        #if source_id != 5220404075366707584:
        #    continue
        ##########################################

        ## require that we only look at things Kounkel labelled as strings
        #if int(group_id) not in np.array(sdf2_str['group_id']).astype(int):
        #    continue
        #if name != 'AB_Dor':
        #    continue
        # if source_id != 5579169050153502976:
        #     continue

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
