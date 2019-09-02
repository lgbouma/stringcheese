from astrobase import lcmath
import numpy as np, pandas as pd
import os, textwrap
from glob import glob
from datetime import datetime

def mask_timegap_edges(time, flux, othervectors=None, gap=0.5,
                       padding=12/(24), verbose=True):
    """
    Drop the times near the edges of any kind of time gap; this is the obvious
    way of treating the "ramp" systematic.

    args:

        time: times with FINITE fluxes.

        flux: FINITE fluxes corresponding to times

        othervectors: None, or list of vectors with same length as flux.

        gap (float): minimum amount of gap space to be considered a gap.

        padding (float): amount of time at the edge of each time group to drop.
        12 hours seems to be a typical duration for the ramp systematic.

    returns:
        time, flux, othervectors: with `padding` days trimmed out
    """

    if not np.all(np.isfinite(flux)):
        raise AssertionError(
            'mask_timegap_edges requires all fluxes to be finite'
        )

    n_groups, groups = lcmath.find_lc_timegroups(time, mingap=gap)

    if verbose:
        print('mask_timegap_edges: got {} timegroups'.format(n_groups))

    sel = np.zeros_like(time).astype(bool)

    for group in groups:

        tg_time = time[group]

        try:
            start_mask = (np.min(tg_time), np.min(tg_time) + padding)
        except Exception as e:
            print(e)
            import IPython; IPython.embed()
        end_mask = (np.max(tg_time) - padding, np.max(tg_time))

        sel |= (
            (time > max(start_mask)) & (time < min(end_mask))
        )

    return_time = time[sel]
    return_flux = flux[sel]

    if isinstance(othervectors, list):
        return_vectors = []
        for othervector in othervectors:
            return_vectors.append(
                othervector[sel]
            )
    else:
        return_vectors = None

    return return_time, return_flux, return_vectors


def remove_janky_lc_segments(time, flux, othervectors=None, gap=0.5, verbose=True):
    """
    Sometimes a star's light curve in one sector is trash compared to the other
    sectors (e.g., if it turned out to be on the chip edge).

    Given the time, flux, and any ancillary time-series, this function filters
    out such light curve segments.
    """

    if not np.all(np.isfinite(flux)):
        raise AssertionError(
            'mask_timegap_edges requires all fluxes to be finite'
        )

    n_groups, groups = lcmath.find_lc_timegroups(time, mingap=gap)

    print('remove_janky_lc_segments: got {} timegroups'.format(n_groups))

    sel = np.zeros_like(time).astype(bool)

    for group in groups:

        tg_time = time[group]
        tg_flux = flux[group]

        tg_std = np.std(tg_flux)

    #FIXME: WIP finish implementing!
    # for group in groups:

    #     tg_time = time[group]

    #     start_mask = (np.min(tg_time), np.min(tg_time) + padding)
    #     end_mask = (np.max(tg_time) - padding, np.max(tg_time))

    #     sel |= (
    #         (time > max(start_mask)) & (time < min(end_mask))
    #     )

    # return_time = time[sel]
    # return_flux = flux[sel]

    # if isinstance(othervectors, list):
    #     return_vectors = []
    #     for othervector in othervectors:
    #         return_vectors.append(
    #             othervector[sel]
    #         )
    # else:
    #     return_vectors = None

    return return_time, return_flux, return_vectors



def mask_orbit_start_and_end(time, flux, orbitgap=1, expected_norbits=2,
                             orbitpadding=6/(24)):
    """
    Ignore the times near the edges of orbits. (Works if you know the gap
    structure within the orbits).

    args:
        time, flux
    returns:
        time, flux: with `orbitpadding` days trimmed out
    """
    norbits, groups = lcmath.find_lc_timegroups(time, mingap=orbitgap)

    if norbits != expected_norbits:
        errmsg = 'got {} orbits, expected {}. groups are {}'.format(
            norbits, expected_norbits, repr(groups))
        raise AssertionError

    sel = np.zeros_like(time).astype(bool)
    for group in groups:
        tg_time = time[group]
        start_mask = (np.min(tg_time), np.min(tg_time) + orbitpadding)
        end_mask = (np.max(tg_time) - orbitpadding, np.max(tg_time))
        sel |= (
            (time > max(start_mask)) & (time < min(end_mask))
        )

    return_time = time[sel]
    return_flux = flux[sel]

    return return_time, return_flux
