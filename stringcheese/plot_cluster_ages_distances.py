"""
age vs distance.
"""

import numpy as np, pandas as pd, matplotlib.pyplot as plt

c19_df = pd.read_csv('../data/Curtis19_notable_clusters.csv')

df = pd.read_csv('../data/string_table2.csv')

# manually verified that these show strong rotation signals!
verif_df = pd.read_csv('../data/20190903_verification_pages.csv')

f,ax = plt.subplots(figsize=(1.5*4,1.5*3))

for key, row in c19_df.iterrows():

    ax.scatter(row['distance_pc'], row['age_gyr'], s=10, alpha=.9, c='k')
    ax.annotate(row['name'], xy=(row['distance_pc'], row['age_gyr']),
                fontsize='small')

for key, row in df.iterrows():

    plx_mas = row['parallax']
    dist_pc = float(1/(plx_mas * 1e-3))

    is_str = str(row['string'])

    if dist_pc > 1000:
        continue

    age = 10**(float(row['age']))
    age_gyr = age/(1e9)

    if int(row['group_id']) in np.array(verif_df['group_id']).astype(int):
        color='C1'
        outstr = (
            'verified for: groupid{}: {}. dist: {:.1f} pc, age: {:.3f} Gyr'.
            format(row['group_id'], row['name'], dist_pc, age_gyr)
        )
        print(outstr)
    else:
        color='C0'

    if is_str=='y':
        ax.scatter(dist_pc, age_gyr, s=8, alpha=.7, c=color, lw=0, marker='s')

    else:
        ax.scatter(dist_pc, age_gyr, s=8, alpha=.7, c=color, lw=0, marker='o')


ax.set_xlabel('Distance [pc]')
ax.set_ylabel('Age [Gyr]')

ax.set_xlim([10,2500])

ax.set_title('black: famous clusters/assocs. blue: CK19 table 2 (square=is_str).'
             '\norange: CK19 table 2, and verified to have strong rotators.',
             fontsize='x-small')

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.get_yaxis().set_tick_params(which='both', direction='in')
ax.get_xaxis().set_tick_params(which='both', direction='in')
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize('small')
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize('small')
ax.set_xscale('log')
ax.set_yscale('log')

f.tight_layout()
f.savefig('../results/cluster_ages_distances.png', dpi=300, bbox_inches='tight')
