###########
# imports #
###########

import os, subprocess, shlex
import multiprocessing as mp
import numpy as np, pandas as pd
from glob import glob
from itertools import product

from astropy.coordinates import SkyCoord
import astropy.units as u, astropy.constants as const

from numpy import array as nparr

##########
# config #
##########

TICCOLS = ['ID', 'version', 'HIP', 'TYC', 'UCAC', 'TWOMASS', 'SDSS', 'ALLWISE',
           'GAIA', 'APASS', 'KIC', 'objType', 'typeSrc', 'ra', 'dec',
           'POSflag', 'pmRA', 'e_pmRA', 'pmDEC', 'e_pmDEC', 'PMflag', 'plx',
           'e_plx', 'PARflag', 'gallong', 'gallat', 'eclong', 'eclat', 'Bmag',
           'e_Bmag', 'Vmag', 'e_Vmag', 'umag', 'e_umag', 'gmag', 'e_gmag',
           'rmag', 'e_rmag', 'imag', 'e_imag', 'zmag', 'e_zmag', 'Jmag',
           'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag', 'TWOMflag', 'prox',
           'w1mag', 'e_w1mag', 'w2mag', 'e_w2mag', 'w3mag', 'e_w3mag', 'w4mag',
           'e_w4mag', 'GAIAmag', 'e_GAIAmag', 'Tmag', 'e_Tmag', 'TESSflag',
           'SPFlag', 'Teff', 'e_Teff', 'logg', 'e_logg', 'MH', 'e_MH', 'rad',
           'e_rad', 'mass', 'e_mass', 'rho', 'e_rho', 'lumclass', 'lum',
           'e_lum', 'd', 'e_d', 'ebv', 'e_ebv', 'numcont', 'contratio',
           'disposition', 'duplicate_id', 'priority', 'eneg_EBV', 'epos_EBV',
           'EBVflag', 'eneg_Mass', 'epos_Mass', 'eneg_Rad', 'epos_Rad',
           'eneg_rho', 'epos_rho', 'eneg_logg', 'epos_logg', 'eneg_lum',
           'epos_lum', 'eneg_dist', 'epos_dist', 'distflag', 'eneg_Teff',
           'epos_Teff', 'TeffFlag', 'gaiabp', 'e_gaiabp', 'gaiarp', 'e_gaiarp',
           'gaiaqflag', 'starchareFlag', 'VmagFlag', 'BmagFlag', 'splists',
           'e_RA', 'e_Dec', 'RA_orig', 'Dec_orig', 'e_RA_orig', 'e_Dec_orig',
           'raddflag', 'wdflag', 'objID']

TICDTYPE = { 'ID': str, 'version': str, 'HIP': str, 'TYC': str, 'UCAC': str,
            'TWOMASS': str, 'SDSS': str, 'ALLWISE': str, 'GAIA': str, 'APASS':
            str, 'KIC': str, 'objType': str, 'typeSrc': str, 'ra': float,
            'dec': float, 'POSflag': str, 'pmRA': float, 'e_pmRA': float,
            'pmDEC': float, 'e_pmDEC': float, 'PMflag': str, 'plx': float,
            'e_plx': float, 'PARflag': str, 'gallong': float, 'gallat': float,
            'eclong': float, 'eclat': float, 'Bmag': float, 'e_Bmag': float,
            'Vmag': float, 'e_Vmag': float, 'umag': float, 'e_umag': float,
            'gmag': float, 'e_gmag': float, 'rmag': float, 'e_rmag': float,
            'imag': float, 'e_imag': float, 'zmag': float, 'e_zmag': float,
            'Jmag': float, 'e_Jmag': float, 'Hmag': float, 'e_Hmag': float,
            'Kmag': float, 'e_Kmag': float, 'TWOMflag': str, 'prox': float,
            'w1mag': float, 'e_w1mag': float, 'w2mag': float, 'e_w2mag': float,
            'w3mag': float, 'e_w3mag': float, 'w4mag': float, 'e_w4mag': float,
            'GAIAmag': float, 'e_GAIAmag': float, 'Tmag': float, 'e_Tmag':
            float, 'TESSflag': str, 'SPFlag': str, 'Teff': float, 'e_Teff':
            float, 'logg': float, 'e_logg': float, 'MH': float, 'e_MH': float,
            'rad': float, 'e_rad': float, 'mass': float, 'e_mass': float,
            'rho': float, 'e_rho': float, 'lumclass': str, 'lum': float,
            'e_lum': float, 'd': float, 'e_d': float, 'ebv': float, 'e_ebv':
            float, 'numcont': str, 'contratio': float, 'disposition': str,
            'duplicate_id': str, 'priority': float, 'eneg_EBV': float,
            'epos_EBV': float, 'EBVflag': str, 'eneg_Mass': float, 'epos_Mass':
            float, 'eneg_Rad': float, 'epos_Rad': float, 'eneg_rho': float,
            'epos_rho': float, 'eneg_logg': float, 'epos_logg': float,
            'eneg_lum': float, 'epos_lum': float, 'eneg_dist': float,
            'epos_dist': float, 'distflag': str, 'eneg_Teff': float,
            'epos_Teff': float, 'TeffFlag': str, 'gaiabp': float, 'e_gaiabp':
            float, 'gaiarp': float, 'e_gaiarp': float, 'gaiaqflag': str,
            'starchareFlag': str, 'VmagFlag': str, 'BmagFlag': str, 'splists':
            str, 'e_RA': float, 'e_Dec': float, 'RA_orig': float, 'Dec_orig':
            float, 'e_RA_orig': float, 'e_Dec_orig': float, 'raddflag': str,
            'wdflag': str, 'objID': str, }

#######
# run #
#######

def xmatch_sourceids_to_tic8(df, tmag_cut=16):

    df.source_id = df.source_id.astype(str)

    ticdir = '/nfs/phtess2/ar0/TESS/CAT/TIC8'
    ticpaths = np.sort(glob(os.path.join(ticdir, 'tic_*.csv.gz')))[::-1]

    #selcols = ['ID', 'GAIA', 'Vmag', 'Tmag', 'Teff', 'logg', 'rad', 'e_rad',
    #           'mass', 'ebv', 'e_ebv', 'source_id']

    for ix, ticchunk in enumerate(ticpaths):
        print('{}/{}'.format(ix, len(ticpaths)))

        outstr = os.path.basename(ticchunk).replace('.csv.gz','_xmatch.csv')
        outpath = '../data/kc19_ticv8_xmatch/{}'.format(outstr)

        if not os.path.exists(outpath):
            # low_memory=False prevents chunking. force correct data types.
            tic_df = pd.read_csv(
                ticchunk, names=TICCOLS, low_memory=False, dtype=TICDTYPE
            )

            tic_df = tic_df[tic_df.Tmag < tmag_cut]

            mdf = df.merge(tic_df, how='inner', left_on='source_id',
                           right_on='GAIA')

            # outdf = mdf[selcols]
            outdf = mdf

            outdf.to_csv(outpath, index=False)
            print('made {}'.format(outpath))

        else:
            print('found {}'.format(outpath))


if __name__ == "__main__":
    df = pd.read_csv('../data/string_table1.csv')
    xmatch_sourceids_to_tic8(df, tmag_cut=16)
