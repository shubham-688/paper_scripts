import matplotlib.pyplot as plt
import numpy as np

def rf_wiggles_baz_r(str1, str2, tr1, tr2, sta, btyp='baz', tmin=-10., tmax=30,
                   scale=None, save=False, ftitle='Figure_rf_wiggle_baz',
                   wvtype='P', fmt='png'):
    """
    Plots receiver function seismograms sorted by back-azimuth or slowness.
    Args:
        str1 (obspy.stream): Stream 1
        str2 (obspy.stream): Stream 2
        tr1 (obspy.trace):
            Trace 1 (normally obtained from the ``utils.stack_all`` function)
        tr2 (obspy.trace): Trace 2
        sta (str): Station name
        btyp (str, optional): Type of sorting for panel
        tmin (float, optional): Lower bound of time axis (s)
        tmax (float, optional): Upper bound of time axis (s)
        scale (float, optional): Scaling factor
        save (bool, optional): Whether or not to save the figure
        ftitle (str, optional): Title of figure to be saved
        wvtype (str, optional): wave type
    Returns:
        None
    """

    if not (btyp == 'baz' or btyp == 'slow' or btyp == 'dist'):
        raise ValueError('type has to be "baz" or "slow" or "dist"')

    if not fmt in ['png', 'PNG', 'jpg', 'JPG', 'eps', 'EPS', 'pdf', 'PDF']:
        raise ValueError("'fmt' has to be one of 'png', 'jpg', 'eps', 'pdf'")

    print()
    print('Plotting Wiggles by '+btyp)

    # Time axis
    nn = str1[0].stats.npts
    sr = str1[0].stats.sampling_rate
    time = np.arange(-nn/2, nn/2)/sr

    # Initialize figure
    fig = plt.figure(figsize=(3,6))
    plt.clf()

    # Get more control on subplots
    ax1 = fig.add_axes([0.15, 0.825, 0.8, 0.05])
    ax2 = fig.add_axes([0.15, 0.1, 0.8, 0.7])
    # ax3 = fig.add_axes([0.45, 0.825, 0.3, 0.05])
    # ax4 = fig.add_axes([0.45, 0.1, 0.3, 0.7])

    # Plot stack of all traces from str1 on top left
    ax1.fill_between(time, 0., tr1.data, where=tr1.data+1e-6 <= 0.,
                     facecolor='lightsteelblue', linewidth=0)
    ax1.fill_between(time, 0., tr1.data, where=tr1.data+1e-6 >= 0.,
                     facecolor='lightseagreen', linewidth=0)
    ax1.set_ylim(-np.max(np.abs(tr1.data)), np.max(np.abs(tr1.data)))
    ax1.set_yticks(())
    ax1.set_xticks(())
    ax1.set_title('Radial')
    ax1.set_xlim(tmin, tmax)

    # Plot sorted traces from str1 on bottom left panel
    for tr in str1:

        if scale:
            maxval = scale
            # Define y axis
            if btyp == 'baz':
                y = tr.stats.baz
            elif btyp == 'slow':
                y = tr.stats.slow
            elif btyp == 'dist':
                y = tr.stats.slow
        else:
            # Define y axis
            if btyp == 'baz':
                y = tr.stats.baz
                maxval = 180
            elif btyp == 'slow':
                y = tr.stats.slow
                maxval = 0.02
            elif btyp == 'dist':
                y = tr.stats.slow
                maxval = 20

        # Fill positive in red, negative in blue
        ax2.fill_between(time, y, y+tr.data*maxval, where=tr.data+1e-6 <= 0.,
                         facecolor='lightsteelblue', linewidth=0)
        ax2.fill_between(time, y, y+tr.data*maxval, where=tr.data+1e-6 >= 0.,
                         facecolor='mediumaquamarine', linewidth=0)

    ax2.set_xlim(tmin, tmax)

    if btyp == 'baz':
        ax2.set_ylim(-5, 370)
        ax2.set_ylabel('Back-azimuth (deg)')

    elif btyp == 'slow':
        if wvtype == 'P':
            ax2.set_ylim(0.038, 0.082)
        elif wvtype == 'S':
            ax2.set_ylim(0.07, 0.125)
        elif wvtype == 'SKS':
            ax2.set_ylim(0.03, 0.06)
        ax2.set_ylabel('Slowness (s/km)')
    elif btyp == 'dist':
        if wvtype == 'P':
            ax2.set_ylim(28., 92.)
        elif wvtype == 'S':
            ax2.set_ylim(53., 107.)
        elif wvtype == 'SKS':
            ax2.set_ylim(83., 117.)
        ax2.set_ylabel('Distance (deg)')

    ax2.set_xlabel('Time (sec)')
    ax2.grid(ls=':')


    if save:
        plt.savefig(ftitle+'.'+fmt, dpi=300, bbox_inches='tight', format=fmt)
    else:
        plt.show()

    return
