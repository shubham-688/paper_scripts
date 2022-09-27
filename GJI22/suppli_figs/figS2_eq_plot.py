import matplotlib.pyplot as plt
import numpy as np
import obspy
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin
from obspy.clients.fdsn import Client
import os
import cartopy.crs as ccrs

import circle as cir_robin
from importlib import reload
reload(cir_robin)
import requests

boundaries = requests.get("https://raw.githubusercontent.com/fraxen/tectonicplates/master/GeoJSON/PB2002_boundaries.json").json()

####
station='AEB15'
sta_lat= -29.3061
sta_long=136.726
client_5g=Client('http://auspass.edu.au:80',user='5g',password='grape71')
inventory = client_5g.get_stations(network="5G",station=station,level='response')
catalog = read_events('AEB15.xml')

# eq_map = Basemap(projection='robin', resolution = 'i', area_thresh = 1000.0,
#               lat_0=sta_lat, lon_0=sta_long)
# eq_map.drawcoastlines()
# eq_map.drawcountries()
# eq_map.fillcontinents(color = 'LightGoldenrodYellow')
# eq_map.drawmapboundary()
# eq_map.drawmeridians(np.arange(0, 360, 60))
# eq_map.drawparallels(np.arange(-90, 90, 30))
# plt.show()

###########


### get eq data
lats, longs = [], []
mags = []
azi=[]

for event in catalog:
    lats.append(event.origins[0].latitude)
    longs.append(event.origins[0].longitude)
    mags.append(event.magnitudes[0].mag)
#

fig, ax=plt.subplots(figsize=(10,6))
# ax = plt.axes(projection=ccrs.Mollweide(central_longitude=sta_long))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=sta_long))


# ax = plt.axes(projection=ccrs.AzimuthalEquidistant(central_longitude=sta_long,central_latitude=sta_lat))

ax.stock_img()
ax.coastlines(color='black', linewidth=.75)
ax.plot(sta_long, sta_lat, color='indigo', marker='^', markersize=7, transform=ccrs.Geodetic())

min_marker_size = 1
for i in range(len(lats)):
    # x,y = eq_map(lon, lat)
    msize = mags[i] * min_marker_size
    # marker_string = get_marker_color(mag)
    ax.plot(longs[i], lats[i],color='black',marker='o',markersize=msize,alpha=.1,transform=ccrs.Geodetic())

X,Y=cir_robin.equi(sta_long, sta_lat, 3330)
X1,Y1=cir_robin.equi(sta_long, sta_lat, 10576)

plt.plot(X,Y,transform=ccrs.Geodetic(),lw=1.25,linestyle='--',c='maroon')
plt.plot(X1,Y1,transform=ccrs.Geodetic(),lw=1.25,linestyle='--',c='maroon')

# Plot boundaries.
for f in boundaries["features"]:
    c = np.array(f["geometry"]["coordinates"])
    lng, lat = c[:, 0], c[:, 1]
    x, y = lng, lat
    mask = np.bitwise_or(np.abs(x) > 1e15, np.abs(y) > 1e15)
    x = np.ma.array(x)
    y = np.ma.array(y)
    x.mask = mask
    y.mask = mask
    plt.plot(x, y, color="Navy", lw=.5,transform=ccrs.Geodetic())

# circle1=plt.Circle((sta_long, sta_lat), 30, fill=False,lw=1, linestyle='--',edgecolor='maroon',transform=ccrs.Geodetic())
# circle2=plt.Circle((sta_long, sta_lat), 95, fill=False,lw=1.5,linestyle='-.', edgecolor='maroon',transform=ccrs.Geodetic())
#
# ax.add_patch(circle1)
# ax.add_patch(circle2)
ax.text(sta_long-10, sta_lat+4, 'AEB15',fontsize=8,fontfamily='serif', color='indigo',transform=ccrs.Geodetic())
ax.text(98, -60, '30°',fontsize=10,fontfamily='serif', color='maroon',transform=ccrs.Geodetic())
ax.text(-125, 0, '95°',fontsize=10,fontfamily='serif', color='maroon',transform=ccrs.Geodetic())

plt.show()
plt.savefig('eq_aeb15_robin.pdf',bbox_inches='tight', pad_inches=0.1)
