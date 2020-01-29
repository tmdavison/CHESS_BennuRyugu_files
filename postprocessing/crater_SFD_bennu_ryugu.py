'''
A script to produce Figure 2 of Davison et al., "Constraining the age and
strength of Bennu and Ryugu using collisional history modelling" from the
output of the CHESS Monte Carlo model.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

import process_model as pm

# Set the fonts
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \usepackage{amsmath}'

# Which models to plot
strengths = ['1e1', '1e2', '1e3', '1e4', '1e5']
lts = [10, 30, 300, 1000]

# Set up the figure
gridspec = gs.GridSpec(ncols=1, nrows=2)
fig = plt.figure(figsize=(5, 5.85))
fig.subplots_adjust(left=0.14, right=0.86, bottom=0.08, top=0.91, hspace=0.10)
ax1 = fig.add_subplot(gridspec[0])
ax2 = fig.add_subplot(gridspec[1])

# Initialise
datasets = []
dlabels = []

# Process and plot the data
datasets, dlabels = pm.plot_loop(
        ['Bennu', 'Ryugu'], [ax1, ax2], [3e3, 3e3],
        strengths, lts, datasets, dlabels, rescale=1.3
        )

# Add the legend
mainlegend = pm.fix_legend(ax1, shift1=69, shift2=61, bbox=(0.98, 1.20), loc=1)

# Plot the observations
pm.plot_bennu(ax1, mleg=mainlegend, y0=0.7, x0=0.98)
pm.plot_ryugu(ax2, mleg=mainlegend, y0=0.8, x0=0.98)

# For scaling
diams = [0.525, 1.]
sareas = [4. * np.pi * (d/2.)**2 for d in diams]

# Customise the axes
for sa, ax in zip(sareas, [ax1, ax2]):

    ax.set_rasterization_zorder(-5)

    axt = ax.twinx()
    axt.set_yscale('log')
    ylim = ax.get_ylim()
    ylim_scaled = [yl/sa for yl in ylim]
    axt.set_ylim(*ylim_scaled)
    axt.set_zorder(-2)

    axt.set_ylabel(r'Cumulative number of\ncraters $\geq D$ per km$^2$')
    axt.yaxis.set_label_coords(1.09, 0.5)

# Save the figure
plt.savefig('crater_SFD_bennu_ryugu_rim-to-rim.pdf')
