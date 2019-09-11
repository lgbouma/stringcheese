"""
(first step)
parse the Kounkel & Covey (2019) catalogs to get metadata needed to make
"interesting" cuts good for rotation period measurements.
"""

import numpy as np, pandas as pd
from astropy.io.votable import from_table, writeto, parse

getfile = '../data/gaia_archive_kc19_string_table1-result.vot.gz'
vot = parse(getfile)
tab = vot.get_first_table().to_table()
# NOTE: "to pandas" with gaia source IDs is a BAD BAD idea, because it
# literally converts a large fraction of them wrong. (probably int32 vs int64
# type errors? super unclear. i'm also having difficulty making a simple
# working reproducible example for this bug.... regardless, reassigning the
# sourceid seems to work)
gaia_df = tab.to_pandas()
gaia_df['source_id'] = tab['source_id']

# source_id,ra,dec,group_id
df1 = pd.read_csv('../data/string_table1.csv')
df1['source_id'] = np.int64(df1['source_id'])

# group_id,name,age,av,string,width,height,length,l,b,parallax,vl,vb,vr,x,y,z,u,v,w
df2 = pd.read_csv('../data/string_table2.csv')

selcols = 'group_id,name,age,av,string'.split(',')
df2 = df2[selcols]

#
# merge against these
#
# mdf = df1.merge(gaia_df, on='source_id', how='left') #FIXME
mdf = gaia_df.merge(df1, on='source_id', how='left')
print('{} entries from Kounkel table 1'.format(len(df1)))
print('{} entries from gaia match of Kounkel table 1'.format(len(gaia_df)))
print('\n{} entries in merge of the two...'.format(len(mdf)))

#
# get T mag using Stassun+2019 Eq1
#
Tmag_pred = (
    mdf['phot_g_mean_mag']
	- 0.00522555 * (mdf['phot_bp_mean_mag'] - mdf['phot_rp_mean_mag'])**3
	+ 0.0891337 * (mdf['phot_bp_mean_mag'] - mdf['phot_rp_mean_mag'])**2
	- 0.633923 * (mdf['phot_bp_mean_mag'] - mdf['phot_rp_mean_mag'])
	+ 0.0324473
)

mdf['Tmag_pred'] = Tmag_pred

#
# merge on group_id to get ages, and whether it is in a "string" or not
#
amdf = mdf.merge(df2, on='group_id', how='left')

#
# get "absolute G" analog (ignoring extinction). note that it's a bit noisy,
# but it's better than ignoring the parallax entirely. note also that i'm
# ignoring the extinction correction, b/c it's basically not known. (though
# Green+19 might differ).
#
amdf['M_G'] = (
    amdf['phot_g_mean_mag'] +
    + 5*np.log10(amdf['parallax']/1e3)
)


print('\n{} entries after merging outside to get ages...'.format(len(amdf)))

outpath = '../data/kounkel_table1_sourceinfo.csv'
amdf.to_csv(outpath, index=False, header=True)
print('made {}'.format(outpath))
