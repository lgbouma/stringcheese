import numpy as np, pandas as pd, matplotlib.pyplot as plt

df = pd.read_csv('../data/kounkel_table1_sourceinfo.csv')

print('{} before finite Tmag cut'.format(len(df)))
sel = ~pd.isnull(df['Tmag_pred'])

df = df[sel]
print('{} after finite Tmag cut'.format(len(df)))

f,ax = plt.subplots(figsize=(4,3))

magstr = 'Tmag_pred'
bins = np.arange(np.floor(np.min(df[magstr])),
                 np.ceil(np.max(df[magstr]))+1,
                 1)
ax.hist(df[magstr], bins=bins, cumulative=True, color='black', fill=False,
        linewidth=0.5)

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.get_yaxis().set_tick_params(which='both', direction='in')
ax.get_xaxis().set_tick_params(which='both', direction='in')
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize('small')
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize('small')
ax.set_xlabel('predicted TESS magnitude')
ax.set_ylabel('Cumulative number')
ax.set_yscale('log')

f.tight_layout(pad=0.2)
f.savefig('../results/cdf_Tmag_table1_sourceinfo.png', dpi=300)
