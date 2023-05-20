import warnings

import matplotlib.patheffects as PathEffects
from matplotlib import ticker
from matplotlib.ticker import (MultipleLocator)
from matplotlib.ticker import (AutoMinorLocator, FixedLocator, FixedFormatter,
                               MaxNLocator,MultipleLocator)
import matplotlib.pyplot as plt
import numpy as np


def _label(stream):
    label_fmts = ['{network}.{station}.{location}.{channel}',
                  '{network}.{station}.{location}.{channel:.2}?',
                  '{network}.{station}']
    for label_fmt in label_fmts:
        labelset = {label_fmt.format(**tr.stats) for tr in stream}
        if len(labelset) == 1:
            return labelset.pop()
    return ''

def get_decay_cos(r0,Dt,time):
    m_t=r0**(time/Dt)*np.cos(np.pi*time/Dt)
    return m_t
#######

def plot_rf_edi(stream, dt,Dt,r0,rf_stack_all,rf_pcorr_stack,fname=None, fig_width=7.5, trace_height=0.5,
            stack_height=0.5, dpi=None,
            scale=1, fillcolors=(None, None), trim=None,
            info=(('back_azimuth', u'baz (°)', 'C0'),
                  ('distance', u'dist (°)', 'C3')),
            show_vlines=False):
    """
    Plot receiver functions.

    :param stream: stream to plot
    :param fname: filename to save plot to. Can be None. In this case
        the figure is left open.
    :param fig_width: width of figure in inches
    :param trace_height: height of one trace in inches
    :param stack_height: height of stack axes in inches
    :param dpi: dots per inch for the created figure
    :param scale: scale for individual traces
    :param fillcolors: fill colors for positive and negative wiggles
    :param trim: trim stream relative to onset before plotting using
         `~.rfstream.RFStream.slice2()`
    :param info: Plot one additional axes showing maximal two entries of
        the stats object. Each entry in this list is a list consisting of
        three entries: key, label and color.
        info can be None. In this case no additional axes is plotted.
    :param show_vlines: If True, show vertical alignment grid lines on plot
        at positions of the major x-tick marks.
    """

    if len(stream) == 0:
        return
    if trim:
        stream = stream.slice2(*trim, reftime='onset')
    if info is None:
        info = ()
    N = len(stream)

    N = len(range(0,360,10)) # that is 36
    # calculate axes and figure dimensions
    # big letters: inches, small letters: figure fraction
    H = trace_height
    HS = stack_height
    FB = 0.5
    FT = 0.2
    DW = 0.1
    FH = H * (N + 2) + HS + FB + FT + DW
    h = H / FH
    hs = HS / FH
    fb = FB / FH
    ft = FT / FH
    FL = 0.5
    FR = 0.2
    FW = fig_width
    FW3 = 0.8
    FW2 = FW - FL - FR - (DW + FW3) * bool(info)
    fl = FL / FW
    fr = FR / FW
    fw2 = FW2 / FW
    fw3 = FW3 / FW
    # init figure and axes
    fig = plt.figure(figsize=(FW, FH), dpi=dpi)
    ax1 = fig.add_axes([fl, fb, fw2, h * (N + 2)])
    if info:
        ax3 = fig.add_axes(
            [1 - fr - fw3, fb, fw3, h * (N + 2)], sharey=ax1)
        info = list(info)
        info[0] = [ax3] + list(info[0])
        if len(info) > 1:
            ax4 = ax3.twiny()
            info[1] = [ax4] + list(info[1])
    # plot individual receiver functions

    def _plot(ax, t, d, i):
        c1, c2 = fillcolors
        if c1:
            ax.fill_between(t, d + i, i, where=d >= 0, lw=0., facecolor=c1)
        if c2:
            ax.fill_between(t, d + i, i, where=d < 0, lw=0., facecolor=c2)
        ax.plot(t, d + i, 'k',lw=0.05)
    def _plot_(ax, t, d, i,color):
        c1, c2 = fillcolors
        if c1:
            ax.fill_between(t, d + i, i, where=d >= 0, lw=0, facecolor=c1,alpha=0)
        if c2:
            ax.fill_between(t, d + i, i, where=d < 0, lw=0, facecolor=c2,alpha=0)
        ax.plot(t, d + i,color,lw=1.15,ls='dashed',alpha=0.85)
    xlim = (0, 0)
    max_ = max(np.max(np.abs(tr.data)) for tr in stream)
    for i, bz in enumerate(np.arange(0,360,10)):
        for j, tr in enumerate(stream):
            if tr.stats.back_azimuth == bz+5:
                # print(bz)
                times = tr.times(reftime=tr.stats.onset)
                xlim = (min(xlim[0], times[0]), max(xlim[1], times[-1]))
                _plot(ax1, times, tr.data / max_ * scale, i + 1)
                break
            else:
                continue

    # plot right axes with header information
    for ax, header, label, color in info:
        data = [tr.stats[header] for tr in stream]
        ax.plot(data, 1 + np.arange(len(stream)), '.' + color, mec=color)
        ax.set_xlabel(label, color=color, size='small')
        if header == 'back_azimuth':
            ax.set_xticks(np.arange(5) * 90)
            ax.set_xticklabels(['0', '', '180', '', '360'], size='small')
            #ax.set_yticks(())
            #ax.set_yticklabels('',size='tiny')
            ax.tick_params(axis='y',left=False,labelleft=False,pad=2,length=.1)
        else:
            ax.xaxis.set_major_locator(MultipleLocator(5))
            ax.set_xticklabels(['', '', ''], size='small')
            for l in ax.get_xticklabels():
                l.set_fontsize('small')
            ax.set_yticklabels('')
        ax.xaxis.set_minor_locator(AutoMinorLocator(.5))
    # set x and y limits
    st_label=[' ']
    bz_label=[' ']
    ax1.set_yticks(np.linspace(0,N,N+1))
    ax1.set_xticks([0,5,10,15])
    ax1.xaxis.set_minor_locator(AutoMinorLocator(.5))
    for i, bz in enumerate(np.arange(0,360,10)):
        if (bz+5) %25 == 0:
            bz_label.append(u'{}°'.format(bz+5))
        else:
            bz_label.append(' ')
        for j, tr in enumerate(stream):
            if tr.stats.back_azimuth == bz+5:
                st_label.append(tr.stats.rfs_in_bin)
                break
        if j==len(stream)-1:
            st_label.append(' ')
            # else:
            #     continue
    bz_label[1]=u'5°'
    bz_label[36]=u'355°'

    for i,tr in enumerate(stream):
        st_label.append(tr.stats.rfs_in_bin)
    plt.setp(ax1.get_yticklabels(),fontsize=7)
    ax1.yaxis.set_major_formatter(ticker.FixedFormatter(st_label))
    ax1.set_xlim(*xlim)
    ax1.set_ylim(-0.5, N + 1.5)
    ax1.tick_params(axis='y',left=False,labelleft=True,pad=2,length=.1)
    #ax1.set_yticklabels('')

    ax1.set_xlabel('Time (s)',fontsize=10)###     here
    ax1.set_ylabel('No. of Rf',fontsize=10)####     here


    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax_1T = ax1.twinx()
    ax_1T.set_ylim(-0.5, N + 1.5)
    ax_1T.set_yticks(np.linspace(0,N,N+1))
    ax_1T.tick_params(axis='y',left=False,labelleft=False,labelright=True,pad=2,length=2)
    plt.setp(ax_1T.get_yticklabels(),fontsize=7)
    ax_1T.yaxis.set_major_formatter(ticker.FixedFormatter(bz_label))
    aligner_color = "#a0a0a080"
    # ax_1T.set_ylabel('Back azimuth (°)',fontsize=10)#######     here
    if show_vlines:
        ax1.xaxis.grid(True, color=aligner_color, linestyle=':')
        ###
    # to plot Ps PpPs and PpSs+PsPs 1 deg= 0.0174533 rad
    p=6.4/111.1949
    incl_rad = np.deg2rad(19) # inclination varies between 17 and 27.
    Vp_inv = p/np.sin(incl_rad)
    Vs_inv=Vp_inv*1.7 # k=1.75
    term1 = np.sqrt(Vs_inv ** 2 - p ** 2)
    term2 = np.sqrt(Vp_inv ** 2 - p ** 2)

    t1 = 35 * (term1 - term2)
    t2 = 35 * (term1 + term2)
    t3 = 35 * 2 * term1

    try:
        t1 += dt
        t2 += Dt-dt
        t3 += Dt
    except:
        pass
    #'t1_offset':dt,'t2_offset':Dt - dt,'t3_offset':Dt
    ax1.axvline(t1,0,1, color='navy', ls='dashed',lw=1)
    ax1.axvline(t2,0,1, color='navy', ls='dashed',lw=1)
    ax1.axvline(t3,0,1, color='navy', ls='dashed',lw=1)
    #####
    Vs_inv=Vp_inv*1.7 # k=1.75
    term1 = np.sqrt(Vs_inv ** 2 - p ** 2)
    term2 = np.sqrt(Vp_inv ** 2 - p ** 2)

    t1 = 55 * (term1 - term2)
    t2 = 55 * (term1 + term2)
    t3 = 55 * 2 * term1
    try:
        t1 += dt
        t2 += Dt-dt
        t3 += Dt
    except:
        pass

    ax1.axvline(t1,0,1, color='navy', ls='dashed',lw=1)
    ax1.axvline(t2,0,1, color='navy', ls='dashed',lw=1)
    ax1.axvline(t3,0,1, color='navy', ls='dashed',lw=1)

    # plot stack
    try:
        stack = rf_stack_all
    except ValueError:
        msg = 'Different npts for traces in one RF plot. Do not plot stack.'
        warnings.warn(msg)
    else:
        if len(stack) > 1:
            warnings.warn('Different stations or channels in one RF plot. ' +
                          'Do not plot stack.')
        elif len(stack) == 1:
            ax2 = fig.add_axes([fl, 1 - ft - hs, fw2, hs], sharex=ax1)
            _plot(ax2, times, stack[0].data, 0)
            # ax2 = fig.add_axes([fl, 1 - ft - hs, fw2, hs], sharex=ax1)

            _plot_(ax2, times, rf_pcorr_stack[0].data, 0,'maroon') # plots P-corr Rfs

            # r0 Dt business
            time = np.arange(stack[0].stats.npts) * stack[0].stats.delta
            time=time-2.5 # shifts onset to 0 sec
            if not np.isnan(r0):
                # dt=auto2.stats.tpdelay
                # reso_filter=get_autoC(auto2.stats.r0_mean,auto2.stats.Dt_mean,len(auto2.data),auto2.stats.sampling_rate)
                cos_dec=get_decay_cos(r0,Dt,time+2.5)
                # ax2.plot(time+2.5,cos_dec,color='darkslateblue',lw=1.15,ls='dashed',alpha=0.85)

                _plot_(ax2, time+2.5+dt, cos_dec, 0,'tomato') # plots P-corr Rfs

                roDt = '(r0 $\Delta$t) : (%.2f %.2f)' % (r0, Dt)
                # ax1.annotate(roDt, (fl-.022, 1 - 1.55 * ft),
                #              xycoords='figure fraction', va='top', ha='right',rotation=90,
                #              bbox=bbox, clip_on=False,fontsize=9,color='darkgreen',weight="bold")
                bbox = dict(boxstyle='round', facecolor='white', alpha=0.8, lw=0)
                ax1.annotate(roDt, (fl+.22, 1 - 0.5 * ft),xycoords='figure fraction',\
                            va='top', ha='right',bbox=bbox, clip_on=False,fontsize=8,\
                            color='salmon')
            #######################
            for l in ax2.get_xticklabels():
                l.set_visible(False)
            # ax2.yaxis.set_major_locator(MaxNLocator(1))
            ax2.tick_params(axis='y',left=False,labelleft=False,pad=2,length=.1)
            for l in ax2.get_yticklabels():
                l.set_fontsize('small')
            if show_vlines:
                ax2.xaxis.grid(True, color=aligner_color, linestyle=':')
    # annotate plot with seed id
    bbox = dict(boxstyle='round', facecolor='white', alpha=0.8, lw=0)
    title = '%s Rfs - %s' % (rf_stack_all[0].stats.totalrf, rf_stack_all[0].stats.station)
    ax1.annotate(title, (.9 - .1 * fr, 1 - .5 * ft),
                 xycoords='figure fraction', va='top', ha='right',
                 bbox=bbox, clip_on=False,fontsize=10)
    ax1.annotate('Stack', (fl-.022, 1 - 1.55 * ft),
                 xycoords='figure fraction', va='top', ha='right',rotation=90,
                 bbox=bbox, clip_on=False,fontsize=10,color='darkgreen',weight="bold")
    # save plot
    if fname:
        fig.savefig(fname, dpi=dpi)
        plt.close(fig)
    else:
        return fig
######
def plot_rf_t(stream, dt,Dt,rf_stack_all,fname=None, fig_width=7.5, trace_height=0.5,
            stack_height=0.5, dpi=None,
            scale=1, fillcolors=(None, None), trim=None,
            info=(('back_azimuth', u'baz (°)', 'C0'),
                  ('distance', u'dist (°)', 'C3')),
            show_vlines=False):
    """
    Plot receiver functions.

    :param stream: stream to plot
    :param fname: filename to save plot to. Can be None. In this case
        the figure is left open.
    :param fig_width: width of figure in inches
    :param trace_height: height of one trace in inches
    :param stack_height: height of stack axes in inches
    :param dpi: dots per inch for the created figure
    :param scale: scale for individual traces
    :param fillcolors: fill colors for positive and negative wiggles
    :param trim: trim stream relative to onset before plotting using
         `~.rfstream.RFStream.slice2()`
    :param info: Plot one additional axes showing maximal two entries of
        the stats object. Each entry in this list is a list consisting of
        three entries: key, label and color.
        info can be None. In this case no additional axes is plotted.
    :param show_vlines: If True, show vertical alignment grid lines on plot
        at positions of the major x-tick marks.
    """

    if len(stream) == 0:
        return
    if trim:
        stream = stream.slice2(*trim, reftime='onset')
    if info is None:
        info = ()
    N = len(stream)
    N = len(range(0,360,10)) # that is 36
    # calculate axes and figure dimensions
    # big letters: inches, small letters: figure fraction
    H = trace_height
    HS = stack_height
    FB = 0.5
    FT = 0.2
    DW = 0.1
    FH = H * (N + 2) + HS + FB + FT + DW
    h = H / FH
    hs = HS / FH
    fb = FB / FH
    ft = FT / FH
    FL = 0.5
    FR = 0.2
    FW = fig_width
    FW3 = 0.8
    FW2 = FW - FL - FR - (DW + FW3) * bool(info)
    fl = FL / FW
    fr = FR / FW
    fw2 = FW2 / FW
    fw3 = FW3 / FW
    # init figure and axes
    fig = plt.figure(figsize=(FW, FH), dpi=dpi)
    ax1 = fig.add_axes([fl, fb, fw2, h * (N + 2)])
    if info:
        ax3 = fig.add_axes(
            [1 - fr - fw3, fb, fw3, h * (N + 2)], sharey=ax1)
        info = list(info)
        info[0] = [ax3] + list(info[0])
        if len(info) > 1:
            ax4 = ax3.twiny()
            info[1] = [ax4] + list(info[1])
    # plot individual receiver functions

    def _plot(ax, t, d, i):
        c1, c2 = fillcolors
        if c1:
            ax.fill_between(t, d + i, i, where=d >= 0, lw=0., facecolor=c1)
        if c2:
            ax.fill_between(t, d + i, i, where=d < 0, lw=0., facecolor=c2)
        ax.plot(t, d + i, 'k',lw=.6)
    xlim = (0, 0)
    max_ = max(np.max(np.abs(tr.data)) for tr in stream)
    for i, bz in enumerate(np.arange(0,360,10)):
        for j, tr in enumerate(stream):
            if tr.stats.back_azimuth == bz+5:
                # print(bz)
                times = tr.times(reftime=tr.stats.onset)
                xlim = (min(xlim[0], times[0]), max(xlim[1], times[-1]))
                _plot(ax1, times, tr.data / max_ * scale, i + 1)
                break
            else:
                continue

    # plot right axes with header information
    for ax, header, label, color in info:
        data = [tr.stats[header] for tr in stream]
        ax.plot(data, 1 + np.arange(len(stream)), '.' + color, mec=color)
        ax.set_xlabel(label, color=color, size='small')
        if header == 'back_azimuth':
            ax.set_xticks(np.arange(5) * 90)
            ax.set_xticklabels(['0', '', '180', '', '360'], size='small')
            #ax.set_yticks(())
            #ax.set_yticklabels('',size='tiny')
            ax.tick_params(axis='y',left=False,labelleft=False,pad=2,length=.1)
        else:
            # ax.xaxis.set_major_locator(MaxNLocator(1))
            ax.set_xticklabels(['', '', ''], size='small')
            for l in ax.get_xticklabels():
                l.set_fontsize('small')
            ax.set_yticklabels('')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    # set x and y limits
    st_label=[' ']
    bz_label=[' ']
    ax1.set_yticks(np.linspace(0,N,N+1))
    for i, bz in enumerate(np.arange(0,360,10)):
        if (bz+5) %25 == 0:
            bz_label.append(u'{}°'.format(bz+5))
        else:
            bz_label.append(' ')
        for j, tr in enumerate(stream):
            if tr.stats.back_azimuth == bz+5:
                st_label.append(tr.stats.rfs_in_bin)
                break
        if j==len(stream)-1:
            st_label.append(' ')
            # else:
            #     continue
    bz_label[1]=u'5°'
    bz_label[36]=u'355°'

    for i,tr in enumerate(stream):
        st_label.append(tr.stats.rfs_in_bin)
    plt.setp(ax1.get_yticklabels(),fontsize=5)
    ax1.yaxis.set_major_formatter(ticker.FixedFormatter(st_label))
    ax1.set_xlim(*xlim)
    ax1.set_xticks([0,5,10,15])
    ax1.set_ylim(-0.5, N + 1.5)
    ax1.tick_params(axis='y',left=False,labelleft=True,pad=2,length=.1)
    #ax1.set_yticklabels('')

    ax1.set_xlabel('Time (s)',fontsize=10)###     here
    ax1.set_ylabel('No. of Rf',fontsize=10)####     here

    # ax1.set_xlabel('time (s)')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax_1T = ax1.twinx()
    ax_1T.set_ylim(-0.5, N + 1.5)
    ax_1T.set_yticks(np.linspace(0,N,N+1))
    ax_1T.tick_params(axis='y',left=False,labelleft=False,labelright=True,pad=2,length=2)
    plt.setp(ax_1T.get_yticklabels(),fontsize=5)
    ax_1T.yaxis.set_major_formatter(ticker.FixedFormatter(bz_label))
    aligner_color = "#a0a0a080"
    if show_vlines:
        ax1.xaxis.grid(True, color=aligner_color, linestyle=':')
        ###
    # to plot Ps PpPs and PpSs+PsPs 1 deg= 0.0174533 rad
    p=6.4/111.1949
    incl_rad = np.deg2rad(19) # inclination varies between 17 and 27.
    Vp_inv = p/np.sin(incl_rad)
    Vs_inv=Vp_inv*1.7 # k=1.75
    term1 = np.sqrt(Vs_inv ** 2 - p ** 2)
    term2 = np.sqrt(Vp_inv ** 2 - p ** 2)

    t1 = 35 * (term1 - term2)
    t2 = 35 * (term1 + term2)
    t3 = 35 * 2 * term1

    try:
        t1 += dt
        t2 += Dt-dt
        t3 += Dt
    except:
        pass
    #'t1_offset':dt,'t2_offset':Dt - dt,'t3_offset':Dt
    # ax1.axvline(t1,0,1, color='navy', ls='dashed',lw=1)
    # ax1.axvline(t2,0,1, color='navy', ls='dashed',lw=1)
    # ax1.axvline(t3,0,1, color='navy', ls='dashed',lw=1)
    #####
    Vs_inv=Vp_inv*1.7 # k=1.75
    term1 = np.sqrt(Vs_inv ** 2 - p ** 2)
    term2 = np.sqrt(Vp_inv ** 2 - p ** 2)

    t1 = 55 * (term1 - term2)
    t2 = 55 * (term1 + term2)
    t3 = 55 * 2 * term1
    try:
        t1 += dt
        t2 += Dt-dt
        t3 += Dt
    except:
        pass

    # ax1.axvline(t1,0,1, color='navy', ls='dashed',lw=1)
    # ax1.axvline(t2,0,1, color='navy', ls='dashed',lw=1)
    # ax1.axvline(t3,0,1, color='navy', ls='dashed',lw=1)

    # plot stack
    try:
        stack = rf_stack_all.trim2(-2.5,15.5,'onset')
    except ValueError:
        msg = 'Different npts for traces in one RF plot. Do not plot stack.'
        warnings.warn(msg)
    else:
        if len(stack) > 1:
            warnings.warn('Different stations or channels in one RF plot. ' +
                          'Do not plot stack.')
        elif len(stack) == 1:
            ax2 = fig.add_axes([fl, 1 - ft - hs, fw2, hs], sharex=ax1)
            _plot(ax2, times, stack[0].data, 0)
            for l in ax2.get_xticklabels():
                l.set_visible(False)
            # ax2.yaxis.set_major_locator(MaxNLocator(1))
            ax2.tick_params(axis='y',left=False,labelleft=False,pad=2,length=.1)
            for l in ax2.get_yticklabels():
                l.set_fontsize('small')
            if show_vlines:
                ax2.xaxis.grid(True, color=aligner_color, linestyle=':')
    # annotate plot with seed id
    bbox = dict(boxstyle='round', facecolor='white', alpha=0.8, lw=0)
    title = '%s traces  %s' % (rf_stack_all[0].stats.totalrf, _label(stream))
    ax1.annotate(title, (.9 - .1 * fr, 1 - .5 * ft),
                 xycoords='figure fraction', va='top', ha='right',
                 bbox=bbox, clip_on=False,fontsize=10)
    # save plot
    if fname:
        fig.savefig(fname, dpi=dpi)
        plt.close(fig)
    else:
        return fig
