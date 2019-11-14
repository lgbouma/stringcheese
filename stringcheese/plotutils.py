from datetime import datetime

def savefig(fig, figpath):
    fig.savefig(figpath, dpi=450, bbox_inches='tight')
    print('{}: made {}'.format(datetime.utcnow().isoformat(), figpath))

    pdffigpath = figpath.replace('.png','.pdf')
    fig.savefig(pdffigpath, bbox_inches='tight', rasterized=True, dpi=450)
    print('{}: made {}'.format(datetime.utcnow().isoformat(), pdffigpath))

def format_ax(ax):
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize('small')
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize('small')
