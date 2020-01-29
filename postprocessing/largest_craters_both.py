'''
A script to produce Figure 4 of Davison et al., "Constraining the age and
strength of Bennu and Ryugu using collisional history modelling" from the
output of the CHESS Monte Carlo model.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.interpolate as interp

# Set the fonts
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 8
plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \usepackage{amsmath}'

magma = cm.magma

# Models of interest
strengths = ['4e0', '1e1', '1e2', '1e3', '1e4', '1e5']
ages = [0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000]
markers = ['*', 'o', 's', 'D', '^', 'v']

colours = [cm.magma(i) for i in [0.9, 0.75, 0.6, 0.45, 0.3, 0.15]]

# Set up the figure
fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(3.74, 4.6))
fig.subplots_adjust(left=0.17, top=0.92, hspace=0.1, bottom=0.11, right=0.97)

ax1.set_xscale('log')

b_crater = 160.
r_crater = 300.

with open('largest_crater.txt', 'w') as f:

    for asteroid, ax, crater, in zip(
            ['Bennu', 'Ryugu'],
            (ax1, ax2),
            (b_crater, r_crater)
        ):

        for Y, marker, colour in zip(strengths, markers, colours):

            lcmeds = []
            p_ages = []
            lcerr = []
            lc25, lc75 = [], []

            for age in ages:

                try:

                    if age < 1.:
                        this_age = '{:04d}kyr'.format(int(age*1e3))
                    else:
                        this_age = '{:04d}Myr'.format(int(age))

                    slist0, lc0 = np.genfromtxt(
                        '{}_{}Pa_{}/survive.txt'.format(asteroid, Y, this_age),
                        usecols=(1, 3), unpack=True)

                    mlist, fmfrac = np.genfromtxt(
                        './{}_{}Pa_{}/final_mass.txt'.format(
                            asteroid, Y, this_age),
                        unpack=True, usecols=(0, 2))
                    mlist = mlist.astype(int)
                    N = len(slist0)

                    print("Processing model: {} {} {}".format(asteroid, Y, age))
                    print("Starting with {} survivors".format(N))
                    print("Number of survivors with mass frac > 3 = {}".format(
                        (fmfrac > 3).sum()))

                    slist0 = slist0.astype(int)
                    growlist = mlist[fmfrac > 3]
                    slist = np.array([i for i in slist0 if i not in growlist], dtype=int)
                    lc = np.array([i*1.3 for i in lc0 if i not in growlist], dtype=int)
                    N = len(slist)

                    print("After filtering, processing {} survivors".format(N))

                    p_ages.append(age)
                    lcmeds.append(np.median(lc))
                    lcerr.append(np.subtract(*np.percentile(
                        lc, [75, 25])) / 2.)
                    lc25.append(np.percentile(lc, 25))
                    lc75.append(np.percentile(lc, 75))

                except:
                    pass

            if len(lcmeds) > 0:
                ax.semilogx(
                    p_ages, lcmeds, ls='-', marker=marker, c=colour, mfc='None',
                    label='Y = {} kPa'.format(int(Y[0]) * 10**(int(Y[-1])-3))
                    )

                p_ages = np.array(p_ages)
                lcmeds = np.array(lcmeds)
                lcerr = np.array(lcerr)
                lc25 = np.array(lc25)
                lc75 = np.array(lc75)
                age1 = interp.interp1d(lcmeds, p_ages, fill_value='extrapolate')(crater)
                age2 = interp.interp1d(lc25, p_ages, fill_value='extrapolate')(crater)
                age3 = interp.interp1d(lc75, p_ages, fill_value='extrapolate')(crater)

                line = '{} {} : {:8.3f} {:8.3f} {:8.3f} | {:8.3f} + {:8.3f} - {:8.3f}'.format(
                    Y, asteroid, age1, age2, age3, age1, age2-age1, age1-age3)
                print(line)
                f.write(line + '\n')


        ax.set_ylabel('Median largest crater\ndiameter on {} [m]'.format(asteroid))

        if asteroid == 'Bennu':
            ax.set_ylim(0, 600)
            ax.set_xticklabels([])
        else:
            ax.set_ylim(0, 1000)
            ax.set_xlabel('Surface age [Myr]')

        ax.set_xlim(0.20, 1500)

        if asteroid == 'Bennu':
            ax.legend(
                bbox_to_anchor=[0, 1.02, 1., 0.102],
                loc=3, borderaxespad=0, framealpha=1,
                fontsize=7, ncol=3, mode='expand'
                )

        ax.grid(axis='y', alpha=0.7)

        ax.axhline(crater, ls='--', c='k')

ax1.text(-0.15, 1.03, '(a)', ha='right', va='top', transform=ax1.transAxes)
ax2.text(-0.15, 1.03, '(b)', ha='right', va='top', transform=ax2.transAxes)
fig.align_ylabels()
fig.savefig('largest_craters_both.pdf')
