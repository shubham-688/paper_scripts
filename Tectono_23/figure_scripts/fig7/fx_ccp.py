import warnings

import matplotlib.patheffects as PathEffects
from matplotlib.ticker import (AutoMinorLocator, FixedLocator, FixedFormatter,
                               MaxNLocator,MultipleLocator)
import matplotlib.pyplot as plt
import numpy as np
###################

def _get_geoaxes(crs=None, latlons=None):
    """Return cartopy geoaxis"""
    if crs is None:
        from cartopy.crs import AzimuthalEquidistant
        latlon0 = np.median(latlons, axis=0)
        crs = AzimuthalEquidistant(*latlon0[::-1])
    return plt.axes(projection=crs)


def __pc():
    from cartopy.crs import PlateCarree as PC
    return PC()

def plot_stations(inventory, label_stations=True, ax=None, crs=None, **kwargs):
    """
    Plot stations.

    :param inventory: station inventory
    :param label_stations: weather to label stations
    :param ax: geoaxes (default None: new ax will be created)
    :param crs: coordinate reference system for new geoaxis, (default: None,
        then AzimuthalEquidistant projection with appropriate center is used.)
    :param \*\*kwargs: other kwargs are passed to ax.scatter() call
    """
    try:
        # assume inventory to be stream
        data = [((tr.stats.station_latitude, tr.stats.station_longitude),
                 tr.stats.station) for tr in inventory]
        latlons, names = zip(*list(set(data)))
    except AttributeError:
        # inventory is Inventory
        latlons, names = zip(*[((sta.latitude, sta.longitude), sta.code)
                               for net in inventory for sta in net])
    if ax is None:
        ax = _get_geoaxes(crs=crs, latlons=latlons)
    kw = dict(s=70, marker='v', c='midnightblue',alpha=0.95, linewidth=0.25, zorder=3)
    kw.update(kwargs)
    ax.scatter(*list(zip(*latlons))[::-1], transform=__pc(), **kw)
    if label_stations:
        path_effect = PathEffects.withStroke(linewidth=2, foreground="white")
        kw = {'xycoords': __pc()._as_mpl_transform(ax),'size':8,
              'xytext': (10, 0), 'textcoords': 'offset points', 'zorder': 4,
              'path_effects': [path_effect]}
        for latlon, name in zip(latlons, names):
            ax.annotate(name, latlon[::-1], **kw)
    return ax

def plot_ppoints(ppoints, inventory=None, label_stations=True, ax=None,
                 crs=None, **kwargs):
    """
    Plot piercing points with stations.

    :param ppoints: list of (lat, lon) tuples of piercing points
    :param inventory, label_stations: plot stations, see `plot_stations`
    :param ax: geoaxes (default None: new ax will be created)
    :param crs: coordinate reference system for new geoaxis, (default: None,
        then AzimuthalEquidistant projection with appropriate center is used.)
    :param \*\*kwargs: other kwargs are passed to ax.scatter() call
    """
    if ax is None:
        ax = _get_geoaxes(crs=crs, latlons=ppoints)
    if inventory is not None:
        plot_stations(inventory, label_stations=label_stations, ax=ax)
    kw = dict(s=45, marker='x', linewidth=1,color='brown', alpha=0.15, zorder=2)
    kw.update(kwargs)
    ax.scatter(*list(zip(*ppoints))[::-1], transform=__pc(), **kw)
    return ax

def plot_profile_map(boxes, inventory=None, label_stations=True, ppoints=None,
                     ax=None, crs=None, **kwargs):
    """
    Plot profile map with stations and piercing points.

    :param boxes: boxes created with `~.profile.get_profile_boxes()`
    :param inventory, label_stations: plot stations, see `plot_stations`
    :param ppoints: list of (lat, lon) tuples of piercing points,
        see `plot_ppoints`
    :param ax: geoaxes (default None: new ax will be created)
    :param crs: coordinate reference system for new geoaxis, (default: None,
        then AzimuthalEquidistant projection with appropriate center is used.)
    :param \*\*kwargs: other kwargs are passed to ax.add_geometries() call
    """
    if ax is None:
        latlons = [boxes[len(boxes)//2]['latlon']]
        ax = _get_geoaxes(crs=crs, latlons=latlons)
        ax.set_adjustable("datalim")
        ax.set_aspect(1)

    if inventory is not None:
        plot_stations(inventory, label_stations=label_stations, ax=ax)
    if ppoints is not None:
        plot_ppoints(ppoints, ax=ax)
    kw = dict(facecolor='none', edgecolor='0.8', zorder=1)
    kw.update(kwargs)
    for box in boxes:
        ax.add_geometries([box['poly']], crs=__pc(), **kw)
    return ax
#######
def plot_profile(profile, fname=None, figsize=None, dpi=None,
                 scale=1, fillcolors=('C3', 'C0'),
                 trim=None, top=None, moveout_model='iasp91'):
    """
    Plot receiver function profile.

    :param profile: stream holding the profile
    :param fname: filename to save plot to. Can be None. In this case
        the figure is left open.
    :param figsize: figsize of the created figure
    :param dpi: dots per inch for the created figure
    :param scale: scale for individual traces
    :param fillcolors: fill colors for positive and negative wiggles
    :param trim: trim stream relative to onset before plotting using
         `~.rfstream.RFStream.slice2()`
    :param top: show second axes on top of profile with additional information.
        Valid values: 'hist' - Plot histogram showing the number of receiver
        functions stacked in the corresponding bin
    :param moveout_model: string with model filename. Will be loaded into a
        `~.simple_model.SimpleModel` object to calculate depths for
        tick labels.
    """
    if len(profile) == 0:
        return
    if trim:
        profile = profile.slice2(*trim, reftime='onset')
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.7])
    widths = [tr.stats.box_length for tr in profile]
    pad = max(1, scale) * min(widths)
    xlim = (min(tr.stats.box_pos for tr in profile) - pad,
            max(tr.stats.box_pos for tr in profile) + pad)
    max_ = max(np.max(np.abs(tr.data)) for tr in profile)
    for tr in profile:
        x = tr.stats.box_pos + scale * tr.data / max_ * min(widths)
        y = tr.times() - (tr.stats.onset - tr.stats.starttime)
        ax.plot(x, y, 'k',lw=0.15)
        c1, c2 = fillcolors
        if c1:
            ax.fill_betweenx(y, x, tr.stats.box_pos,
                             where=x >= tr.stats.box_pos, facecolor=c1)
        if c2:
            ax.fill_betweenx(y, x, tr.stats.box_pos,
                             where=x < tr.stats.box_pos, facecolor=c2)
    ax.set_xlabel('Distance (km)',fontsize=14,labelpad=6)
    ax.set_ylim(max(y), min(y))
    ax.set_ylabel('Time (s)',fontsize=14)
    if moveout_model:
        from rf.simple_model import load_model
        model = load_model(moveout_model)
        phase = profile[0].stats.moveout
        slowness = profile[0].stats.slowness
        pd = model.calculate_delay_times(phase=phase, slowness=slowness)
        ax2 = ax.twinx()
        ax.get_shared_y_axes().join(ax, ax2)
        dkm = 50
        if profile[0].stats.endtime - profile[0].stats.onset > 50:
            dkm = 200
        d1 = np.arange(20) * dkm
        d2 = np.arange(100) * dkm / 5
        t1 = np.interp(d1, model.z, pd)
        t2 = np.interp(d2, model.z, pd)
        myLocator = FixedLocator(t1)
        myMinorLocator = FixedLocator(t2)
        myFormatter = FixedFormatter([str(i) for i in d1])
        # ax2.yaxis.set_major_locator(myLocator)
        # ax2.yaxis.set_minor_locator(myMinorLocator)
        # ax2.yaxis.set_major_formatter(myFormatter)
        ax2.set_ylabel('Time (s)',fontsize=14)
        ax2.set_ylim(ax.get_ylim())
    if top is not None:
        ax3 = fig.add_axes([0.1, 0.85, 0.8, 0.1], sharex=ax)
    if top == 'hist':
        left = [tr.stats.box_pos - tr.stats.box_length / 2 for tr in profile]
        height = [tr.stats.num for tr in profile]
        ax3.bar(left, height, widths, color='silver')
        plt.setp(ax3.get_xticklabels(), visible=False)
        ax3.spines['top'].set_color('k')
        ax3.spines['right'].set_color('none')
        ax3.spines['left'].set_color('none')
        ax3.xaxis.set_ticks_position('bottom')
        ax3.yaxis.set_ticks_position('left')
        ax3.set_yticks((0,212))
    elif top is not None:
        raise NotImplementedError("'%s' not supported for top parameter" % top)
    ax.set_xlim(*xlim)
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=10)
    if fname:
        fig.savefig(fname, dpi=dpi)
        plt.close(fig)
    else:
        return fig
