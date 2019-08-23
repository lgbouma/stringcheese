import os
from glob import glob
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from numpy import array as nparr
import matplotlib as mpl

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import LombScargle

from astrobase.lcmath import phase_magseries

def generate_verification_page(lcd, ls, freq, power, cutoutpaths):
    """
    Make the verification page, which consists of:

    top row: entire light curve (with horiz bar showing rotation period)

    bottom row:
        lomb scargle periodogram  |  phased light curve  |  image w/ aperture

    args:

        lcd (dict): has the light curve, aperture positions, some lomb
        scargle results.

        ls: LombScargle instance with everything passed.

        cutoutpaths (list): FFI cutout FITS paths.
    """
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'

    plt.close('all')
    fig = plt.figure(figsize=(4,3))

    #  top row: entire light curve (with horiz bar showing rotation period)
    ax0 = plt.subplot2grid((2, 3), (0, 0), rowspan=3)

    #TODO: maybe split by sector? meh.
    ax0.scatter(lcd['time'], lcd['rel_flux'], c='k', alpha=1.0, zorder=2, s=10,
                rasterized=True, linewidths=0)

    # lomb scargle periodogram
    ax1 = plt.subplot2grid((2, 3), (1, 0))

    ax1.plot(1/freq, power, c='k')
    ax1.set_xlabel('period [day]')
    ax1.set_ylabel('LS power')

    # phased light curve 
    ax2 = plt.subplot2grid((2, 3), (1, 1))

    phzd = phase_magseries(lcd['time'], lcd['rel_flux'], lcd['ls_period'],
                           lcd['time'][np.argmin(lcd['rel_flux'])], wrap=False,
                           sort=True)

    ax2.scatter(phzd['phase'], phzd['mags'], c='k', rasterized=True, s=10,
                linewidths=0)


    # image w/ aperture. put on the stars as well!
    ax3 = plt.subplot2grid((2, 3), (1, 2), projection = )
    ...

#FIXME
#FIXME
#FIXME
#FIXME

    ax0 = plt.subplot2grid((3, 3), (0, 0), projection=cutout_wcs)
    #
    # wcs information parsing
    # follow Clara Brasseur's https://github.com/ceb8/tessworkshop_wcs_hack
    #
    radius = 2.0*u.arcminute

    nbhr_stars = Catalogs.query_region(
        "{} {}".format(float(c_obj.ra.value), float(c_obj.dec.value)),
        catalog="TIC",
        radius=radius
    )

    cutout_wcs = cd['cutout_wcs']
    try:
        px,py = cutout_wcs.all_world2pix(
            nbhr_stars[nbhr_stars['Tmag'] < Tmag_cutoff]['ra'],
            nbhr_stars[nbhr_stars['Tmag'] < Tmag_cutoff]['dec'],
            0
        )
    except Exception as e:
        print('ERR! wcs all_world2pix got {}'.format(repr(e)))
        return

    ticids = nbhr_stars[nbhr_stars['Tmag'] < Tmag_cutoff]['ID']
    tmags = nbhr_stars[nbhr_stars['Tmag'] < Tmag_cutoff]['Tmag']

    sel = (px > 0) & (px < 9) & (py > 0) & (py < 9)
    px,py = px[sel], py[sel]
    ticids, tmags = ticids[sel], tmags[sel]

    ra, dec = float(c_obj.ra.value), float(c_obj.dec.value)
    target_x, target_y = cutout_wcs.all_world2pix(ra,dec,0)

    # geometry: there are TWO coordinate axes. (x,y) and (ra,dec). To get their
    # relative orientations, the WCS and ignoring curvature will usually work.
    shiftra_x, shiftra_y = cutout_wcs.all_world2pix(ra+1e-4,dec,0)
    shiftdec_x, shiftdec_y = cutout_wcs.all_world2pix(ra,dec+1e-4,0)


    #
    # ax0: OOT
    #
    vmin = np.min([np.min(cd['m_oot_flux']), np.min(cd['m_intra_flux'])])
    vmax = np.max([np.max(cd['m_oot_flux']), np.max(cd['m_intra_flux'])])

    cset0 = ax0.imshow(cd['m_oot_flux'], cmap='YlGnBu_r', origin='lower',
                       zorder=1, vmin=vmin, vmax=vmax)

    ax0.scatter(cd['ctds_oot'][:,0], cd['ctds_oot'][:,1], marker='o',
                linewidths=0, rasterized=True, c='fuchsia', alpha=0.9,
                zorder=3, s=60)

    ax0.scatter(px, py, marker='x', c='r', s=15, rasterized=True, zorder=2,
                linewidths=1)
    ax0.plot(target_x, target_y, mew=0.5, zorder=5, markerfacecolor='yellow',
             markersize=25, marker='*', color='k', lw=0)

    ax0.set_title('OOT (cyan o: centroid for each OOT window)', fontsize='xx-large')

    cb0 = fig.colorbar(cset0, ax=ax0, extend='neither', fraction=0.046, pad=0.04)

    # DSS is ~1 arcsecond per pixel. overplot apertures on axes 6,7
    for ix, radius_px in enumerate([21,21*1.5,21*2.25]):
        circle = plt.Circle((sizepix/2, sizepix/2), radius_px,
                            color='C{}'.format(ix), fill=False, zorder=5+ix)
        ax6.add_artist(circle)


    ##########################################
    fig.tight_layout(pad=3.0, h_pad=1.7)


