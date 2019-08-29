from astropy.io import fits
import astrobase.imageutils as iu
from astropy import wcs
import numpy as np

import pickle, os

from photutils import aperture_photometry, CircularAperture, CircularAnnulus

def get_lc_given_fficutout(workingdir, cutouts, c_obj, return_pkl=False):
    """
    Do simple aperture photometry on FFI cutouts. Uses world's simplest
    background subtraction -- the cutout median. Imposes an aperture radius of
    3 pixels. Invents the error bars as 1/sqrt(n_counts).

    Saves the lightcurve and related data to a pickle file, in workingdir. If
    the pickle is found to already exist, and return_pkl is True, it is loaded
    and returned

    args:
        workingdir (str): directory to which the light curve is written

        cutouts (list of length at least 1): paths to fficutouts. assumed to be
        from different sectors.

        c_obj (astropy.skycoord): coordinate of object, used to project wcs to
        image.
    """

    outpath = os.path.join(workingdir, 'multisector_lc.pkl')
    if os.path.exists(outpath) and not return_pkl:
        print('WRN! found {}, returning without load'.format(outpath))
        return
    elif os.path.exists(outpath) and return_pkl:
        print('found {}, returning with load'.format(outpath))
        with open(outpath, 'rb') as f:
            out_dict = pickle.load(f)
        return out_dict

    if len(cutouts) == 0:
        raise AssertionError('something wrong in tesscut! FIXME')

    # img_flux and img_flux_err are image cubes of (time x spatial x spatial).
    # make lists of them for each sector.
    img_fluxs = [iu.get_data_keyword(f, 'FLUX') for f in cutouts]

    # for the background, just take the median of the image. very simple.
    bkgd_fluxs = [np.array([np.median(img_flux[ix, :, :]) for ix in
                            range(len(img_flux))]) for img_flux in img_fluxs]

    img_flux_errs = [iu.get_data_keyword(f, 'FLUX_ERR') for f in cutouts]

    times = [iu.get_data_keyword(f, 'TIME')+2457000 for f in cutouts]
    qualitys = [iu.get_data_keyword(f, 'QUALITY') for f in cutouts]

    cut_hduls = [fits.open(f) for f in cutouts]
    cutout_wcss = [wcs.WCS(cuthdul[2].header) for cuthdul in cut_hduls]

    # get the location to put down the apertures
    try:
        xs, ys = [
            c_wcs.all_world2pix(c_obj.icrs.ra, c_obj.icrs.dec, 0)
            for c_wcs in cutout_wcss
        ]
    except Exception as e:
        print('ERR! wcs all_world2pix got {}'.format(repr(e)))
        return

    #
    # plop down a 3 pixel circular aperture at the locations. then make the
    # light curves by doing the sum!
    #
    positions = [(x, y) for x,y in zip(xs, ys)]
    circ_apertures = [
        CircularAperture(position, r=3) for position in positions
    ]

    fluxs = []
    # iterate over sectors
    for img, bkgd, aper in zip(img_fluxs, bkgd_fluxs, circ_apertures):

        img_stack = img - bkgd[:,None,None]

        # iterate over cadences in sector
        s_flux = []
        for _img in img_stack:

            phot_table = aperture_photometry(_img, aper)

            s_flux.append(phot_table['aperture_sum'])

        fluxs.append(np.array(s_flux))

    #
    # FIXME: might want to
    # now take only quality==0 cadences. maybe mask orbit edges too.
    # oh and stitch sectors together!
    #

    # time = time[quality == 0]
    # flux = flux[quality == 0]
    # flux_err = flux_err[quality == 0]

    # normalize each sector by its median
    rel_fluxs = [f/np.nanmedian(f) for f in fluxs]

    rel_flux_errs = [np.sqrt(f)/np.nanmedian(f) for f in fluxs]

    #
    # all except last are output as a 1-dimensional array. time, quality, flux,
    # rel_flux, and rel_flux_err are all of same length. xs and ys are length
    # n_sectors; they are the positions used in the aperture.
    #
    out_dict = {
        'time':np.concatenate(times).flatten(),
        'quality':np.concatenate(qualitys).flatten(),
        'flux':np.concatenate(fluxs).flatten(),
        'rel_flux':np.concatenate(rel_fluxs).flatten(),
        'rel_flux_err':np.concatenate(rel_flux_errs).flatten(),
        'x':np.array(xs).flatten(),
        'y':np.array(ys).flatten(),
        'cutout_wcss': cutout_wcss
    }

    with open(outpath, 'wb') as f:
        pickle.dump(out_dict, f)

    return out_dict
