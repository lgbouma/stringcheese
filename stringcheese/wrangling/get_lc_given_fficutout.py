from astropy.io import fits
import astrobase.imageutils as iu
from astropy import wcs
import numpy as np
from numpy import array as nparr

import pickle, os

from photutils import aperture_photometry, CircularAperture, CircularAnnulus

from astrobase.lcmath import (
    sigclip_magseries_with_extparams, find_lc_timegroups
)
from numpy.polynomial.legendre import Legendre

from stringcheese import lcutils as lcu

from scipy.stats import iqr

def get_lc_given_fficutout(workingdir, cutouts, c_obj, return_pkl=False):
    """
    Do simple aperture photometry on FFI cutouts. Uses world's simplest
    background subtraction -- the cutout median. Imposes an aperture radius of
    3 pixels. Invents the error bars as 1/sqrt(n_counts).

    To clean up the light curve, the following steps did a decent job:

    If multi-sector, each sector is normalized by its median. Sectors are
    then stitched together, and only quality==0 cadences are taken.
    If any "timegroup" (usually sectors, but not strictly -- here I define it
    by 0.5 day gaps) has much worse interquartile range than the other (>5x the
    median IQR), drop that timegroup. This usually means that the star was on
    bad pixels, or the edge of the detector, or something similar for one
    sector.

    Then, sigma clip out any ridiculous outliers via [7sigma, 7sigma] symmetric
    clip.

    Then, required all fluxes and errors were finite, and for each timegroup
    masked out 0.5 days at the beginning, and 0.5 days at the end. This makes
    the gaps bigger, but mostly throws out ramp systematic-infested data that
    othewise throws off the period measurement.

    Finally, an apparently relatively common long-term systematic in these LCs
    looks just like a linear slope over the course of an orbit. (Maybe b/c
    stars are moving, but aperture center is not? Unclear.) Therefore, to
    detrend, I fit out a LINE in time from each time group, if the group has at
    least two days worth of data. (Not if shorter, to avoid overfitting).

    Then, save the lightcurve and related data to a pickle file, in workingdir.
    If the pickle is found to already exist, and return_pkl is True, it is
    loaded and returned

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
        'rel_flux': sigma-clipped relative flux, fully detrended
        'rel_flux_err': sigma-clipped relative flux errors
        'predetrending_time':
        'predetrending_rel_flux': before fitting out the line, rel flux values
        'predetrending_rel_flux_err':
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
        raise AssertionError('something wrong in tesscut! fix this!')

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
    # concatenate sectors together and take only quality==0 cadences.
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

    #
    # sort everything into time order
    #
    sind = np.argsort(time)

    time = time[sind]
    quality = quality[sind]
    flux = flux[sind]
    rel_flux = rel_flux[sind]
    rel_flux_err = rel_flux_err[sind]

    #
    # if any "timegroup" (usually sectors, but not strictly -- here I define it
    # by 0.5 day gaps) has much worse interquartile range, drop it. This
    # usually means that the star was on bad pixels, or the edge of the
    # detector, or something similar for one sector.
    #
    ngroups, groups = find_lc_timegroups(time, mingap=0.5)

    rel_flux_iqrs = [iqr(rel_flux[group], rng=(25,75)) for group in groups]

    if ngroups >= 3:

        median_iqr = np.nanmedian(rel_flux_iqrs)

        bad_groups = (rel_flux_iqrs > 5*median_iqr)

        if len(bad_groups[bad_groups]) > 0:
            print('WRN! got {} bad time-groups. dropping them.'.
                  format(len(bad_groups[bad_groups])))

            gd_inds = nparr(groups)[~bad_groups]

            time = np.concatenate([time[gd] for gd in gd_inds]).flatten()
            quality = np.concatenate([quality[gd] for gd in gd_inds]).flatten()
            flux = np.concatenate([flux[gd] for gd in gd_inds]).flatten()
            rel_flux = np.concatenate([rel_flux[gd] for gd in gd_inds]).flatten()
            rel_flux_err = np.concatenate([rel_flux_err[gd] for gd in gd_inds]).flatten()

        else:
            # did not find any bad groups
            pass

    else:
        # trickier to detect outlying sectors with fewer groups
        pass

    #
    # sigma clip out any ridiculous outliers via [7sigma, 7sigma] symmetric
    # clip.
    #
    stime, srel_flux, srel_flux_err, [sflux, squality] = (
        sigclip_magseries_with_extparams(
        time, rel_flux, rel_flux_err, [flux, quality],
        sigclip=[7,7], iterative=False, magsarefluxes=True)
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
            'predetrending_time':stime,
            'predetrending_rel_flux':srel_flux,
            'predetrending_rel_flux_err':srel_flux_err,
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
    # Fit out a LINE in time from each time group, if the group has at least
    # two days worth of data. I added this because an apparently relatively
    # common long-term systematic in these LCs looks just like a linear slope
    # over the course of an orbit. (Maybe b/c stars are moving, but aperture
    # center is not? Unclear.)
    #
    predetrending_time = stime
    predetrending_rel_flux = srel_flux
    predetrending_rel_flux_err = srel_flux_err

    ngroups, groups = find_lc_timegroups(stime, mingap=0.5)

    _time, _rflux, _rflux_err = [], [], []

    for group in groups:

        tg_time = stime[group]
        tg_rel_flux = srel_flux[group]
        tg_rel_flux_err = srel_flux_err[group]

        if tg_time.max() - tg_time.min() < 2:

            # don't try fitting out trends in small time groups (risks
            # overfitting).

            _time.append(tg_time)
            _rflux.append(tg_rel_flux)
            _rflux_err.append(tg_rel_flux_err)

            continue

        try:

            p = Legendre.fit(tg_time, tg_rel_flux, 1)
            coeffs = p.coef

            tg_fit_rel_flux = p(tg_time)

            # divide out the linear fit
            tg_dtr_rel_flux = tg_rel_flux/tg_fit_rel_flux

            _time.append(tg_time)
            _rflux.append(tg_dtr_rel_flux)
            _rflux_err.append(tg_rel_flux_err)

        except np.linalg.LinAlgError:
            print('WRN! Legendre.fit failed, b/c bad data for this group. '
                  'Continue.')
            continue

    stime = np.concatenate(_time).flatten()
    srel_flux = np.concatenate(_rflux).flatten()
    srel_flux_err = np.concatenate(_rflux_err).flatten()

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
        'predetrending_time':predetrending_time,
        'predetrending_rel_flux':predetrending_rel_flux,
        'predetrending_rel_flux_err':predetrending_rel_flux_err,
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
