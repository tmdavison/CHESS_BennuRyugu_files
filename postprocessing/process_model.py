import numpy as np
import matplotlib.ticker as ticker
import matplotlib.cm as cm

bennu_all = np.genfromtxt('Walsh2019All.txt')
bennu_distinct = np.genfromtxt('Walsh2019Distinct.txt')

ryugu_all = np.genfromtxt('Sugita2019_CL-1-3.dat')
ryugu_distinct = np.genfromtxt('Sugita2019_CL-1-2.dat')

def process_model(
        asteroid, Y, c, lw, lt, ls, medmark, ax, datasets,
        dlabels, rescale=1., stopafterwrite=False
        ):

    d = np.array([])

    slist0, lc = np.genfromtxt(
        './{}_{}Pa_{:04d}Myr/survive.txt'.format(asteroid, Y, lt),
        unpack=True, usecols=(1, 3)
        )
    slist0 = slist0.astype(int)
    N = len(slist0)

    mlist, fmfrac = np.genfromtxt(
        './{}_{}Pa_{:04d}Myr/final_mass.txt'.format(asteroid, Y, lt),
        unpack=True, usecols=(0, 2)
        )
    mlist = mlist.astype(int)

    if stopafterwrite:
        if len(mlist) < 10000:
            return 'Not'

    print("Processing model: {} {} {}".format(asteroid, Y, lt))
    print("Starting with {} survivors".format(N))
    print("Number of survivors with mass frac > 3 = {}".format((fmfrac > 3).sum()))

    growlist = mlist[fmfrac > 3]
    slist = np.array([i for i in slist0 if i not in growlist], dtype=int)
    N = len(slist)

    print("After filtering, processing {} survivors".format(N))

    try:
        d = np.load('{}_{}Pa_{:04d}Myr_craters.npy'.format(asteroid, Y, lt))
        N = d[0]
        d = d[1:]

        if stopafterwrite:
            return 'Already'

    except IOError:

        print('{}_{}Pa_{:04d}Myr_craters.npy not found. Creating now...'.format(asteroid, Y, lt))

        for i in slist:
            di = np.genfromtxt(
                './{}_{}Pa_{:04d}Myr/craters/craters-{:06d}.txt.gz'.format(
                    asteroid, Y, lt, i), skip_header=1,
                usecols=(3,), unpack=True
                )

            d = np.concatenate((d, di))

            if i%1000 == 0:
                print(Y, i)

        d.sort()
        d = d[::-1]  # reverse sorted

        d *= 1e3  # in m

        if len(mlist) == 10000:
            # simulation has finished running, so we can save the craters file
            np.save('{}_{}Pa_{:04d}Myr_craters.npy'.format(asteroid, Y, lt),
                    np.insert(d, 0, N)  # insert the number of parent bodies at the start
                    )

            if stopafterwrite:
                return ''

    # Convert to rim-to-rim if required
    d *= rescale

    cN = np.arange(1, len(d)+1, 1.) / (float(N))

    print("Ended with    {} survivors".format(int(N)))
    print("Maximum crater size = {:.2f} m".format(d.max()))

    label = '{:3g} kPa; {:4g} Myr'.format(float(Y)/1e3, lt)

    ax.loglog(d, cN, color=c, lw=lw, ls=ls,
              label=label,
              zorder=-10)

    datasets.append(lc[lc < 1e3])
    print(label, (lc > 4e2).sum(), N)
    dlabels.append(label)

    print("\n================================\n")

    return datasets, dlabels

def plot_loop(asteroids, axes, ymaxes, strengths, lts, ds, dl, rescale=1.):

    colors = [cm.magma(i) for i in [0.2, 0.4, 1, 0.6, 0.8]][::-1]
    colors2 = [cm.magma(i) for i in [0.2, 0.4, 0.6, 0.8]][::1]

    for asteroid, ax, ymax in zip(asteroids, axes, ymaxes):

        for Y, c, lw in zip(strengths, colors, (1.5, 1.5, 2.5, 1.5, 1.5)):

            ds, dl = process_model(asteroid, Y, c, lw, 100, '-', 'o', ax, ds, dl, rescale=rescale)

        for lt, c, lw in zip(lts, colors2, (1.5, 1.5, 1.5, 1.5)):

            ds, dl = process_model(
                asteroid, '1e3', c, lw, lt, '--', 's', ax, ds, dl, rescale=rescale)

        ax.axhline(1, color='#AAAAAA', zorder=-20, ls='--', lw=1)

        ax.set_ylabel('Cumulative nunmber of\ncraters $\geq D$ on {}'.format(asteroid))

        ax.set_ylim(3e-1, ymax)
        ax.set_xlim(2e1, 1e3)

        ax.loglog(10, 1e-1, 'w-', label=' ')
        ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:4g}'))

    if len(axes) > 1:
        axes[0].set_xticklabels([])

    axes[-1].set_xlabel('Rim-to-rim crater diameter, $D$ [m]')

    return ds, dl

def fix_legend(ax, shift1=67, shift2=58, bbox=(1.04, 1.07), handlelength=1.5, loc=1):

    handles, labels = ax.get_legend_handles_labels()

    myorder = [0, 1, 2, 3, 4, 8, 7, 9, 6, 5]

    handles = [handles[i] for i in myorder]
    labels = [labels[i] for i in myorder]

    legend1 = ax.legend(
        handles=handles, labels=labels, fontsize=8.5,
        loc=loc, ncol=2, framealpha=1,
        bbox_to_anchor=bbox, borderaxespad=0,
        columnspacing=0.6, handletextpad=0.4,
        handlelength=handlelength)

    shift = (shift1, shift2)
    for tt, t in enumerate(legend1.get_texts()):
        t.set_ha('right')
        if tt < 5:
            t.set_position((shift[0], 0))
        else:
            t.set_position((shift[1], 0))

    ax.add_artist(legend1)

    return legend1

def plot_bennu(ax, mleg=None, y0=0.5, x0=0.95):

    bennu_surface_area = 4. * np.pi * 0.246**2
    b1, = ax.plot(
        bennu_all[:, 0], bennu_all[:, 1] * bennu_surface_area,
        's', mfc='#404040', mec='None', label='Walsh+19, all'
        )
    b2, = ax.plot(
        bennu_distinct[:, 0], bennu_distinct[:, 1] * bennu_surface_area,
        'o', mfc='#808080', mec='None', label='Walsh+19, distinct'
        )

    if mleg is not None:
        bbox = mleg._bbox_to_anchor.inverse_transformed(mleg.axes.transAxes)
        x0 = bbox.x0
    else:
        pass

    ax.legend(handles=[b1, b2], loc=5, bbox_to_anchor=(x0, y0),
              fontsize=8.5, handletextpad=0.4, framealpha=1, borderaxespad=0)

def plot_ryugu(ax, mleg=None, y0=0.5, x0=0.95):

    ryugu_surface_area = 4. * np.pi * 0.448**2
    r1, = ax.plot(
        ryugu_all[:, 0], ryugu_all[:, 1] * ryugu_surface_area,
        's',mfc='#404040', mec='None', label='Sugita+19, CL 1--3'
        )
    r2, = ax.plot(
        
    ryugu_distinct[:, 0], ryugu_distinct[:, 1] * ryugu_surface_area,
        'o', mfc='#808080', mec='None', label='Sugita+19, CL 1--2'
        )

    if mleg is not None:
        bbox = mleg._bbox_to_anchor.inverse_transformed(mleg.axes.transAxes)
        x0 = bbox.x0
    else:
        pass

    ax.legend(handles=[r1, r2], loc=5, bbox_to_anchor=(x0, y0),
              fontsize=8.5, handletextpad=0.4, framealpha=1, borderaxespad=0)
