import os
from glob import glob
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from numpy import array as nparr
import matplotlib as mpl

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import LombScargle
from astroquery.mast import Catalogs

import astropy.visualization as vis

from astrobase.lcmath import phase_magseries

def generate_verification_page(lcd, ls, freq, power, cutoutpaths, c_obj,
                               outvppath, outd):
    """
    Make the verification page, which consists of:

    top row: entire light curve (with horiz bar showing rotation period)

    bottom row:
        lomb scargle periodogram  |  phased light curve  |  image w/ aperture

    ----------
    args:

        lcd (dict): has the light curve, aperture positions, some lomb
        scargle results.

        ls: LombScargle instance with everything passed.

        cutoutpaths (list): FFI cutout FITS paths.

        c_obj (SkyCoord): astropy sky coordinate of the target

        outvppath (str): path to save verification page to
    """
    cutout_wcs = lcd['cutout_wcss'][0]

    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'

    plt.close('all')

    fig = plt.figure(figsize=(12,8))

    ax0 = plt.subplot2grid((2, 3), (0, 0), colspan=3)
    ax1 = plt.subplot2grid((2, 3), (1, 0))
    ax2 = plt.subplot2grid((2, 3), (1, 1))
    ax3 = plt.subplot2grid((2, 3), (1, 2), projection=cutout_wcs)

    #  row 0: entire light curve (with horiz bar showing rotation period)
    ax0.scatter(lcd['time'], lcd['rel_flux'], c='k', alpha=1.0, zorder=2, s=10,
                rasterized=True, linewidths=0)

    epoch = np.nanmin(lcd['time']) + lcd['ls_period']
    yval = np.max(lcd['rel_flux']) + 0.5*np.std(lcd['rel_flux'])
    ax0.plot([epoch, epoch+lcd['ls_period']], [yval, yval], color='red', lw=2)

    ax0.set_xlabel('Time [BJD$_{\mathrm{TDB}}$]')
    ax0.set_ylabel('Relative flux')

    # row 1, col 0: lomb scargle periodogram
    ax1.plot(1/freq, power, c='k')
    ax1.set_xscale('log')
    ax1.text(0.03, 0.97, 'FAP={:.1e}\nP={:.1f}d'.format(
        lcd['ls_fap'], lcd['ls_period']), ha='left', va='top',
        fontsize='large', zorder=2, transform=ax1.transAxes
    )
    ax1.set_xlabel('Period [day]')
    ax1.set_ylabel('LS power')

    # row 1, col 1: phased light curve 
    phzd = phase_magseries(lcd['time'], lcd['rel_flux'], lcd['ls_period'],
                           lcd['time'][np.argmin(lcd['rel_flux'])], wrap=False,
                           sort=True)

    ax2.scatter(phzd['phase'], phzd['mags'], c='k', rasterized=True, s=10,
                linewidths=0)

    ax2.set_xlabel('Phase')
    ax2.set_ylabel('Flux')

    #
    # row1, col2: image w/ aperture. put on the nbhr stars as dots too, to
    # ensure the wcs isn't wonky!
    #

    # acquire neighbor stars.
    radius = 2.0*u.arcminute

    nbhr_stars = Catalogs.query_region(
        "{} {}".format(float(c_obj.ra.value), float(c_obj.dec.value)),
        catalog="TIC",
        radius=radius
    )

    try:
        Tmag_cutoff = 15
        px,py = cutout_wcs.all_world2pix(
            nbhr_stars[nbhr_stars['Tmag'] < Tmag_cutoff]['ra'],
            nbhr_stars[nbhr_stars['Tmag'] < Tmag_cutoff]['dec'],
            0
        )
    except Exception as e:
        print('ERR! wcs all_world2pix got {}'.format(repr(e)))
        return

    tmags = nbhr_stars[nbhr_stars['Tmag'] < Tmag_cutoff]['Tmag']

    sel = (px > 0) & (px < 19) & (py > 0) & (py < 19)
    px,py = px[sel], py[sel]
    tmags = tmags[sel]

    ra, dec = float(c_obj.ra.value), float(c_obj.dec.value)
    target_x, target_y = cutout_wcs.all_world2pix(ra,dec,0)

    #
    # finally make it
    #

    img = lcd['median_imgs'][0]

    interval = vis.PercentileInterval(99.9)
    try:
        vmin,vmax = interval.get_limits(img)
    except Exception as e:
        print(e)
        import IPython; IPython.embed()
    norm = vis.ImageNormalize(
        vmin=vmin, vmax=vmax, stretch=vis.LogStretch(1000))

    cset = ax3.imshow(img, cmap='YlGnBu_r', origin='lower', zorder=1,
                      norm=norm)

    ax3.scatter(px, py, marker='x', c='r', s=5, rasterized=True, zorder=2,
                linewidths=1)
    ax3.plot(target_x, target_y, mew=0.5, zorder=5, markerfacecolor='yellow',
             markersize=7, marker='*', color='k', lw=0)

    cb0 = fig.colorbar(cset, ax=ax3, extend='both', fraction=0.046, pad=0.04)

    # overplot aperture
    radius_px = 3
    circle = plt.Circle((target_x, target_y), radius_px,
                         color='C1', fill=False, zorder=5)
    ax3.add_artist(circle)

    name = outd['name']
    group_id = outd['group_id']
    if name=='nan':
        nstr = 'Group {}'.format(group_id)
    else:
        nstr = '{}'.format(name)


    if not np.isfinite(outd['teff']):
        outd['teff'] = 0

    ax3.set_title(
        'Teff={:d}K. {}'.format(int(outd['teff']), nstr)
    )

    #
    # save
    #
    fig.tight_layout(w_pad=1.0, h_pad=1.0)
    fig.savefig(outvppath, dpi=300, bbox_inches='tight')
    print('made {}'.format(outvppath))
