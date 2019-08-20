import requests
from astroquery.mast import Tesscut

def get_fficutout(c_obj, cutoutdir=None, sector=None):
    # c_obj (SkyCoord): location of target star

    print('beginning download tesscut for {}'.format(repr(c_obj)))
    try:
        tab = Tesscut.download_cutouts(c_obj, size=20, sector=sector,
                                       path=cutoutdir)
    except (requests.exceptions.HTTPError,
            requests.exceptions.ConnectionError) as e:
        print('got {}, try again'.format(repr(e)))
        tab = Tesscut.download_cutouts(c_obj, size=20, sector=sector,
                                       path=cutoutdir)
