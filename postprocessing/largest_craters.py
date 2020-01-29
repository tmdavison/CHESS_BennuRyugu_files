'''
A script to produce Figure 3 of Davison et al., "Constraining the age and
strength of Bennu and Ryugu using collisional history modelling" from the
output of the CHESS Monte Carlo model.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import matplotlib.gridspec as gs

# Set the fonts
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \usepackage{amsmath}'
plt.rcParams['font.size'] = 8

# Parameters of interest
strengths = ['1e1', '1e2', '1e3', '1e4', '1e5']
lts = [10, 30, 300, 1000]

colors = [cm.magma(i) for i in [0.8, 0.6, 1, 0.4, 0.2]]
colors2 = [cm.magma(i) for i in [0.8, 0.6, 0.4, 0.2]][::-1]

rbennu = 262.5
sarea = 4. * np.pi * rbennu**2

gridspec = gs.GridSpec(ncols=1, nrows=26)

datasets = []
dlabels = []

def process_model(asteroid, Y, c, lw, lt, ls, medmark):

    # list of surviving (non-disrupted) parent bodies and their largest craters
    slist0, lc0 = np.genfromtxt(
        '{}_{}Pa_{:04d}Myr/survive.txt'.format(asteroid, Y, lt),
        unpack=True, usecols=(1, 3)#, dtype=(int, float)
        )  #[:200]
    slist0 = slist0.astype(int)

    # Number of asteroids
    N = len(slist0)

    # index and final mass fraction of parent bodies
    mlist, fmfrac = np.genfromtxt(
        './{}_{}Pa_{:04d}Myr/final_mass.txt'.format(asteroid, Y, lt),
        unpack=True, usecols=(0, 2))
    mlist = mlist.astype(int)

    print("Processing model: {} {} {}".format(asteroid, Y, lt))
    print("Starting with {} survivors".format(N))
    print("Number of survivors with mass frac > 3 = {}".format(
        (fmfrac > 3).sum()))

    # Filter out any parent bodies which grew too big
    # (i.e. were accreted onto a larger body)
    growlist = mlist[fmfrac > 3]
    slist = np.array([i for i in slist0 if i not in growlist], dtype=int)
    # Scale by 1.3 for final rim-ro-rim diameter
    lc = np.array([i*1.3 for i in lc0 if i not in growlist], dtype=int)
    N = len(slist)

    print("After filtering, processing {} survivors".format(N))

    label = '{:3g} kPa; {:3g} Myr'.format(float(Y) / 1e3, lt)

    if asteroid == 'Bennu':
        threshold = 2e3
    else:
        threshold = 4e3
    datasets.append(lc[lc < threshold])
    dlabels.append(label)

# Process models for the 4 subplots
for Y, c, lw in zip(strengths, colors, (1.5, 1.5, 2.5, 1.5, 1.5)):

    process_model('Bennu', Y, c, lw, 100, '-', 'o')

for lt, c, lw in zip(lts, colors2, (1.5, 1.5, 1.5, 1.5)):

    process_model('Bennu', '1e3', c, lw, lt, '--', 's')


for Y, c, lw in zip(strengths, colors, (1.5, 1.5, 2.5, 1.5, 1.5)):

    process_model('Ryugu', Y, c, lw, 100, '-', 'o')

for lt, c, lw in zip(lts, colors2, (1.5, 1.5, 1.5, 1.5)):

    process_model('Ryugu', '1e3', c, lw, lt, '--', 's')

# Create the figure
fig = plt.figure(figsize=(3.74, 4))

ax2 = fig.add_subplot(gridspec[:6])
ax3 = fig.add_subplot(gridspec[6:12])  #212)
ax4 = fig.add_subplot(gridspec[14:20])  #212)
ax5 = fig.add_subplot(gridspec[20:26])  #212)

fig.subplots_adjust(hspace=0, bottom=0.12, left=0.20, top=0.96, right=0.96)

def plot_ax(ax, order, which, labeltext, xlabels, asteroid, subfig=None):
    ''' Function to plot violins on axis '''

    lb = [dlabels[i].split(';')[which] for i in order]

    ds = [datasets[i] for i in order]

    violin = ax.violinplot(ds, vert=False, showextrema=False, showmedians=True, widths=0.8)

    for i, col in zip(violin['bodies'], colors):
        i.set_facecolor(col)
    violin['cmedians'].set_edgecolors(colors)

    for p in violin['cmedians'].get_paths():
        p.vertices[0, 1] -= 0.13
        p.vertices[1, 1] += 0.13

    ax.set_xscale('log')

    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:4g}'))

    ax.text(0.02, 0.08, labeltext, ha='left', transform=ax.transAxes, fontsize=8)

    ax.set_yticks((1, 2, 3, 4, 5))
    ax.set_yticklabels(lb)

    ax.set_xlim(25, 2500)

    if 'top' in xlabels:
        ax.xaxis.tick_top()
        ax.set_xticklabels([])
    if xlabels == 'bottom':
        ax.set_xlabel('Rim-to-rim crater diameter, $D$ [m]')

    if asteroid is not None:

        ax.text(-0.22, 0, asteroid, rotation=90,
                va='center', ha='center', transform=ax.transAxes)

        ax.text(-0.22, 1, '({})'.format(subfig),
                va='center', ha='center', transform=ax.transAxes)


# Plot the 4 subplots
plot_ax(ax2, [0, 1, 2, 3, 4], 0, 'Fixed 100 Myr age', 'top', 'Bennu', 'a')
plot_ax(ax3, [8, 7, 2, 6, 5], 1, 'Fixed 1 kPa cohesion', 'hello', None)

plot_ax(ax4, [9, 10, 11, 12, 13], 0, 'Fixed 100 Myr age', 'top', 'Ryugu', 'b')
plot_ax(ax5, [17, 16, 11, 15, 14], 1, 'Fixed 1 kPa cohesion', 'bottom', None)

# Plot the largest observed distinct craters on each body
bl_distinct = np.genfromtxt('Walsh2019Distinct.txt', usecols=(0,), unpack=True).max()
bl_all = np.genfromtxt('Walsh2019All.txt', usecols=(0,), unpack=True).max()

rl_12 = np.genfromtxt('Sugita2019_CL-1-2.dat', usecols=(0,), unpack=True).max()
rl_13 = np.genfromtxt('Sugita2019_CL-1-3.dat', usecols=(0,), unpack=True).max()

ax2.axvline(bl_distinct, ls='--', c='#707070', lw=0.75)
ax3.axvline(bl_distinct, ls='--', c='#707070', lw=0.75)

ax4.axvline(rl_12, ls='--', c='#707070', lw=0.75)
ax5.axvline(rl_12, ls='--', c='#707070', lw=0.75)

# Save figure
fig.savefig('largest_crater.pdf')
