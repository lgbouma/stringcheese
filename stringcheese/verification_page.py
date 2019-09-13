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

from astrobase.lcmath import (
    phase_magseries, find_lc_timegroups, phase_bin_magseries
)

def generate_verification_page(lcd, ls, freq, power, cutoutpaths, c_obj,
                               outvppath, outd, show_binned=True):
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

    fig = plt.figure(figsize=(12,12))

    #ax0 = plt.subplot2grid((3, 3), (0, 0), colspan=3)
    #ax1 = plt.subplot2grid((3, 3), (1, 0), colspan=3)
    #ax2 = plt.subplot2grid((3, 3), (2, 0))
    #ax3 = plt.subplot2grid((3, 3), (2, 1))
    #ax4 = plt.subplot2grid((3, 3), (2, 2), projection=cutout_wcs)

    ax0 = plt.subplot2grid((3, 3), (1, 0), colspan=3)
    ax1 = plt.subplot2grid((3, 3), (2, 0), colspan=3)
    ax2 = plt.subplot2grid((3, 3), (0, 0))
    ax3 = plt.subplot2grid((3, 3), (0, 1))
    ax4 = plt.subplot2grid((3, 3), (0, 2), projection=cutout_wcs)

    #
    # row 0: entire light curve, pre-detrending (with horiz bar showing
    # rotation period). plot model LC too.
    #
    ax0.scatter(lcd['predetrending_time'], lcd['predetrending_rel_flux'],
                c='k', alpha=1.0, zorder=3, s=10, rasterized=True,
                linewidths=0)

    try:
        model_flux = nparr(lcd['predetrending_rel_flux']/lcd['rel_flux'])
    except ValueError:
        model_flux = 0

    if isinstance(model_flux, np.ndarray):
        ngroups, groups = find_lc_timegroups(lcd['predetrending_time'], mingap=0.5)
        for group in groups:
            ax0.plot(lcd['predetrending_time'][group], model_flux[group], c='C0',
                     alpha=1.0, zorder=2, rasterized=True, lw=2)

    # add the bar showing the derived period
    ymax = np.percentile(lcd['predetrending_rel_flux'], 95)
    ymin = np.percentile(lcd['predetrending_rel_flux'], 5)
    ydiff = 1.15*(ymax-ymin)

    epoch = np.nanmin(lcd['predetrending_time']) + lcd['ls_period']
    ax0.plot([epoch, epoch+lcd['ls_period']], [ymax, ymax], color='red', lw=2,
             zorder=4)

    ax0.set_ylim((ymin-ydiff,ymax+ydiff))

    #ax0.set_xlabel('Time [BJD$_{\mathrm{TDB}}$]')
    ax0.set_xticklabels('')
    ax0.set_ylabel('Raw flux')

    name = outd['name']
    group_id = outd['group_id']
    if name=='nan':
        nstr = 'Group {}'.format(group_id)
    else:
        nstr = '{}'.format(name)


    if not np.isfinite(outd['teff']):
        outd['teff'] = 0

    ax0.text(0.98, 0.97,
        'Teff={:d}K. {}'.format(int(outd['teff']), nstr),
             ha='right', va='top', fontsize='large', zorder=2,
             transform=ax0.transAxes
    )

    #
    # row 1: entire light curve (with horiz bar showing rotation period)
    #
    ax1.scatter(lcd['time'], lcd['rel_flux'], c='k', alpha=1.0, zorder=2, s=10,
                rasterized=True, linewidths=0)

    # add the bar showing the derived period
    ymax = np.percentile(lcd['rel_flux'], 95)
    ymin = np.percentile(lcd['rel_flux'], 5)
    ydiff = 1.15*(ymax-ymin)

    epoch = np.nanmin(lcd['time']) + lcd['ls_period']
    ax1.plot([epoch, epoch+lcd['ls_period']], [ymax, ymax], color='red', lw=2)

    ax1.set_ylim((ymin-ydiff,ymax+ydiff))

    ax1.set_xlabel('Time [BJD$_{\mathrm{TDB}}$]')
    ax1.set_ylabel('Detrended flux')

    #
    # row 2, col 0: lomb scargle periodogram
    #
    ax2.plot(1/freq, power, c='k')
    ax2.set_xscale('log')
    ax2.text(0.03, 0.97, 'FAP={:.1e}\nP={:.1f}d'.format(
        lcd['ls_fap'], lcd['ls_period']), ha='left', va='top',
        fontsize='large', zorder=2, transform=ax2.transAxes
    )
    ax2.set_xlabel('Period [day]', labelpad=-1)
    ax2.set_ylabel('LS power')

    #
    # row 2, col 1: phased light curve 
    #
    phzd = phase_magseries(lcd['time'], lcd['rel_flux'], lcd['ls_period'],
                           lcd['time'][np.argmin(lcd['rel_flux'])], wrap=False,
                           sort=True)

    ax3.scatter(phzd['phase'], phzd['mags'], c='k', rasterized=True, s=10,
                linewidths=0, zorder=1)

    if show_binned:
        try:
            binphasedlc = phase_bin_magseries(phzd['phase'], phzd['mags'],
                                              binsize=1e-2, minbinelems=5)
            binplotphase = binphasedlc['binnedphases']
            binplotmags = binphasedlc['binnedmags']

            ax3.scatter(binplotphase, binplotmags, s=10, c='darkorange',
                        linewidths=0, zorder=3, rasterized=True)
        except TypeError as e:
            print(e)
            pass

    xlim = ax3.get_xlim()
    ax3.hlines(1.0, xlim[0], xlim[1], colors='gray', linestyles='dotted',
               zorder=2)
    ax3.set_xlim(xlim)

    ymax = np.percentile(lcd['rel_flux'], 95)
    ymin = np.percentile(lcd['rel_flux'], 5)
    ydiff = 1.15*(ymax-ymin)
    ax3.set_ylim((ymin-ydiff,ymax+ydiff))

    ax3.set_xlabel('Phase', labelpad=-1)
    ax3.set_ylabel('Flux', labelpad=-0.5)

    #
    # row2, col2: image w/ aperture. put on the nbhr stars as dots too, to
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

    # some images come out as nans.
    if np.all(np.isnan(img)):
        img = np.ones_like(img)

    interval = vis.PercentileInterval(99.9)
    vmin,vmax = interval.get_limits(img)
    norm = vis.ImageNormalize(
        vmin=vmin, vmax=vmax, stretch=vis.LogStretch(1000))

    cset = ax4.imshow(img, cmap='YlGnBu_r', origin='lower', zorder=1,
                      norm=norm)

    ax4.scatter(px, py, marker='x', c='r', s=5, rasterized=True, zorder=2,
                linewidths=1)
    ax4.plot(target_x, target_y, mew=0.5, zorder=5, markerfacecolor='yellow',
             markersize=7, marker='*', color='k', lw=0)

    #ax4.coords.grid(True, color='white', ls='dotted', lw=1)
    lon = ax4.coords['ra']
    lat = ax4.coords['dec']

    lon.set_ticks(spacing=1*u.arcminute)
    lat.set_ticks(spacing=1*u.arcminute)

    lon.set_ticklabel(exclude_overlapping=True)
    lat.set_ticklabel(exclude_overlapping=True)

    ax4.coords.grid(True, color='white', alpha=0.3, lw=0.3, ls='dotted')

    #cb0 = fig.colorbar(cset, ax=ax4, extend='neither', fraction=0.046, pad=0.04)

    # overplot aperture
    radius_px = 3
    circle = plt.Circle((target_x, target_y), radius_px,
                         color='C1', fill=False, zorder=5)
    ax4.add_artist(circle)

    #
    # cleanup
    # 
    for ax in [ax0,ax1,ax2,ax3,ax4]:
        ax.get_yaxis().set_tick_params(which='both', direction='in',
                                       labelsize='small', top=True, right=True)
        ax.get_xaxis().set_tick_params(which='both', direction='in',
                                       labelsize='small', top=True, right=True)

    fig.tight_layout(w_pad=0.5, h_pad=0)

    #
    # save
    #
    fig.savefig(outvppath, dpi=300, bbox_inches='tight')
    print('made {}'.format(outvppath))
