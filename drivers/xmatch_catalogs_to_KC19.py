"""
usage:
    python -u xmatch_catalogs_to_KC19.py &> whatever.log &

* ITEM #1: MATCH ROSAT AGAINST KC2019.
* ITEM #1.5: SAME WITH LAMOST (trickier)
* ITEM #2: WISE, SAME IDEA.
* ITEM #3: do a CPM hack.
"""

from astropy import units as u

from stringcheese.crossmatching import (
    match_kc19_to_2rxs,
    plot_kc19_to_2rxs_age_vs_parallax,
    get_kc19_to_2rxs_match_frac,
    match_kc19_to_cs16
)

if __name__ == "__main__":

    # ITEM #1: MATCH ROSAT AGAINST KC2019.
    plot_kc19_to_2rxs_age_vs_parallax(max_sep=12*u.arcsec)
    plot_kc19_to_2rxs_age_vs_parallax(max_sep=20*u.arcsec)
    get_kc19_to_2rxs_match_frac(max_sep=12*u.arcsec)
    get_kc19_to_2rxs_match_frac(max_sep=20*u.arcsec)

    # ITEM #2: MATCH COTTEN & SONG AGAINST KC2019
    match_kc19_to_cs16(max_sep=3*u.arcsec)

