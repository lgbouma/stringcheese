from astropy.io import fits
import astrobase.imageutils as iu
from astropy import wcs
import numpy as np

import pickle, os

from photutils import aperture_photometry, CircularAperture, CircularAnnulus

from astrobase.lcmath import sigclip_magseries_with_extparams

from stringcheese import lcutils as lcu

def get_lc_given_fficutout(workingdir, cutouts, c_obj, return_pkl=False):
    """
    Do simple aperture photometry on FFI cutouts. Uses world's simplest
    background subtraction -- the cutout median. Imposes an aperture radius of
    3 pixels. Invents the error bars as 1/sqrt(n_counts).

    If multi-sector, each sector is normalized by its median. Sectors are
    stitched together, and only quality==0 cadences are taken. Ridiculous
    outliers are sigma clipped out via a [20sigma, 20sigma] symmetric
    clip on the stitched light curve.

    Saves the lightcurve and related data to a pickle file, in workingdir. If
    the pickle is found to already exist, and return_pkl is True, it is loaded
    and returned

    ----------
    args:

        workingdir (str): directory to which the light curve is written

        cutouts (list of length at least 1): paths to fficutouts. assumed to be
        from different sectors.

        c_obj (astropy.skycoord): coordinate of object, used to project wcs to
        image.

    ----------
    returns:

        if fails to get a good light curve, returns None. else, returns
        dictionary with the following keys.

        'time': time array
        'quality': quality flag array
        'flux': sigma-clipped flux (counts)
        'rel_flux': sigma-clipped relative flux
        'rel_flux_err': sigma-clipped relative flux errors
        'x': location of aperture used to extract light curve
        'y': ditto
        'median_imgs': list of median images of the stack used to extract apertures
        'cutout_wcss': WCSs to the above images
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
    bkgd_fluxs = [np.array([np.nanmedian(img_flux[ix, :, :]) for ix in
                            range(len(img_flux))]) for img_flux in img_fluxs]

    img_flux_errs = [iu.get_data_keyword(f, 'FLUX_ERR') for f in cutouts]

    times = [iu.get_data_keyword(f, 'TIME')+2457000 for f in cutouts]
    qualitys = [iu.get_data_keyword(f, 'QUALITY') for f in cutouts]

    cut_hduls = [fits.open(f) for f in cutouts]
    cutout_wcss = [wcs.WCS(cuthdul[2].header) for cuthdul in cut_hduls]

    # get the location to put down the apertures
    try:
        xs, ys = [], []
        for cutout_wcs in cutout_wcss:
            _x, _y = cutout_wcs.all_world2pix(
                c_obj.icrs.ra, c_obj.icrs.dec, 0
            )
            xs.append(_x)
            ys.append(_y)
    except Exception as e:
        print('ERR! wcs all_world2pix got {}'.format(repr(e)))
        import IPython; IPython.embed()

    #
    # plop down a 3 pixel circular aperture at the locations. then make the
    # light curves by doing the sum!
    #
    positions = [(x, y) for x,y in zip(xs, ys)]
    circ_apertures = [
        CircularAperture(position, r=3) for position in positions
    ]

    fluxs = []
    median_imgs = []
    # iterate over sectors
    for img, bkgd, aper in zip(img_fluxs, bkgd_fluxs, circ_apertures):

        img_stack = img - bkgd[:,None,None]

        # iterate over cadences in sector
        s_flux = []
        for _img in img_stack:

            phot_table = aperture_photometry(_img, aper)

            s_flux.append(phot_table['aperture_sum'])

        fluxs.append(np.array(s_flux))

        median_img = np.nanmedian(img_stack, axis=0)
        median_imgs.append(median_img)

    # normalize each sector by its median
    rel_fluxs = [f/np.nanmedian(f) for f in fluxs]
    rel_flux_errs = [np.sqrt(f)/np.nanmedian(f) for f in fluxs]

    #
    # stitch sectors together and take only quality==0 cadences.
    # sigma clip out any ridiculous outliers via [20sigma, 20sigma] symmetric
    # clip.
    #
    time = np.concatenate(times).flatten()
    quality = np.concatenate(qualitys).flatten()
    flux = np.concatenate(fluxs).flatten()
    rel_flux = np.concatenate(rel_fluxs).flatten()
    rel_flux_err = np.concatenate(rel_flux_errs).flatten()

    sel = (quality == 0)

    time = time[sel]
    flux = flux[sel]
    rel_flux = rel_flux[sel]
    rel_flux_err = rel_flux_err[sel]
    quality = quality[sel]

    stime, srel_flux, srel_flux_err, [sflux, squality] = (
        sigclip_magseries_with_extparams(
        time, rel_flux, rel_flux_err, [flux, quality],
        sigclip=[20,20], iterative=False, magsarefluxes=True)
    )

    #
    # require finite fluxes. then mask gap edges. if you get no finite values,
    # save the dud pickle and return None.
    #
    sel = np.isfinite(srel_flux)

    stime = stime[sel]
    sflux = sflux[sel]
    srel_flux = srel_flux[sel]
    srel_flux_err = srel_flux_err[sel]
    squality = squality[sel]

    if len(stime) == len(sflux) == 0:

        out_dict = {
            'time':stime,
            'quality':squality,
            'flux':sflux,
            'rel_flux':srel_flux,
            'rel_flux_err':srel_flux_err,
            'x':np.array(xs).flatten(),
            'y':np.array(ys).flatten(),
            'median_imgs': median_imgs,
            'cutout_wcss': cutout_wcss
        }

        with open(outpath, 'wb') as f:
            pickle.dump(out_dict, f)

        return None

    if not np.any(median_imgs[0]) and np.any(sflux):
        print('somehow getting nan image but finite flux')
        import IPython; IPython.embed()

    stime, srel_flux, [srel_flux_err, sflux, squality] = (
        lcu.mask_timegap_edges(stime, srel_flux,
                               othervectors=[srel_flux_err, sflux, squality],
                               gap=0.5, padding=12/(24))
    )

    #
    # 1-dimensional arrays:
    # time, quality, flux, rel_flux, and rel_flux_err are all of same length.
    #
    # xs and ys are length n_sectors; they are the positions used in the
    # aperture.
    #
    # median_imgs is list of length n_sectors, for which each entry is the
    # median image in that sector.
    #
    # cutout_wcss is list of length n_sectors, each entry is the WCS
    # corresponding to median_image
    #
    out_dict = {
        'time':stime,
        'quality':squality,
        'flux':sflux,
        'rel_flux':srel_flux,
        'rel_flux_err':srel_flux_err,
        'x':np.array(xs).flatten(),
        'y':np.array(ys).flatten(),
        'median_imgs': median_imgs,
        'cutout_wcss': cutout_wcss
    }

    with open(outpath, 'wb') as f:
        pickle.dump(out_dict, f)

    if len(stime) == len(sflux) == 0:
        return None

    else:
        return out_dict
