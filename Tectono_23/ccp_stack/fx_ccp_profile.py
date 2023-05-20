import warnings

import matplotlib.patheffects as PathEffects
from matplotlib.ticker import (AutoMinorLocator, FixedLocator, FixedFormatter,
                               MaxNLocator)
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
