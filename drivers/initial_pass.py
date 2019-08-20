"""
measure rotation periods for bright, young candidates.
(in strings, and not in strings).

as an initial pass, by this I mean:
    T < 13
    7.5 < log(age) < 8.5   (like 30-300 Myr. no younger, b/c reddening.)
"""

import os
from glob import glob
import numpy as np, pandas as pd
from numpy import array as nparr

from astropy.coordinates import SkyCoord
from astropy import units as u

from stringcheese.wrangling import get_fficutout as gfc

def main():
    source_df = pd.read_csv('../data/kounkel_table1_sourceinfo.csv')

    sel = (
        (source_df['Tmag_pred'] < 13)
        &
        (source_df['age'] < 8.5)
        &
        (source_df['age'] > 7.5)
    )

    sdf = source_df[sel]
    print(
        'after making cuts on T<13, log(age) from 7.5-8.5, got {} stars, {} groups'.
        format(len(sdf), len(np.unique(sdf['group_id'])))
    )

    # now given the gaia ids, get the rotation periods
    for ix, r in sdf.iterrows():

        source_id = np.int64(r['source_id'])
        ra, dec = float(r['ra_x']), float(r['dec_x'])
        name = str(r['name'])
        group_id = str(r['group_id'])

        c_obj = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')

        basedir = '/home/luke/stringcheese_cutouts/'
        workingdir = os.path.join(basedir,
                                  'group{}_name{}'.format(group_id, name))
        if not os.path.exists(workingdir):
            os.mkdir(workingdir)
        workingdir = os.path.join(workingdir, str(source_id))
        if not os.path.exists(workingdir):
            os.mkdir(workingdir)

        # if you already downloaded ffi cutouts for this object, dont get any
        # more. otherwise, do.
        if len(glob(os.path.join(workingdir,'*.fits')))>=1:
            pass
        else:
            gfc.get_fficutout(c_obj, cutoutdir=workingdir)

        # get_lc_given_fficutout(source_id)

        # measure_rotation_period(source_id)

if __name__ == "__main__":
    main()
