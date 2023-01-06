from pylab import *
import os
from scipy.interpolate import interp2d
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf
from matplotlib.pylab import rcParams
rcParams['figure.figsize'] = 12, 6
import xarray as xr
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
from cartopy.util import add_cyclic_point
import cartopy.feature as cf
from mpl_toolkits.basemap import Basemap

#Read the file
#dir1 = "/media/michael/Maxtor/ICE/"
dir1 = "/media/michael/Maxtor/HeatFlux/"
#file_name = ("trendB.220.heatflux.SON.nc")
file_name = ("trendB.220.heatflux.SON.nc")
ftmp= xr.open_dataset(dir1+file_name)


lon1 = ftmp['lon'][:] 
lat1 = ftmp['lat'][:] 
#stemp = ftmp['v'][:]
stemp = ftmp['covvt'][:]

num_year=62
#num_year=116
stemp_mean = np.mean(stemp, axis=0)

stemp_new=stemp_mean*num_year
stemp_10 = stemp_new[5,:,:]    #select a partcular level

#print(stemp_10)
#quit()

#print(np.min(stemp_new))
#print(np.max(stemp_new))
#quit()

fig= plt.figure(figsize=(12,6))

[lons,lats] = np.meshgrid(lon1, lat1)
#x,y = m(lons,lats)
#m = Basemap(projection='npaeqd',boundinglat=10, lon_0=270,resolution='l')
#m = Basemap(projection='robin',lon_0=0,resolution='l')
m = Basemap(projection='mill',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,0], fontsize=15)
m.drawmeridians(np.arange(-180.,181.,60.), labels=[0,0,0,1], fontsize=15)
m.drawmapboundary()
m.fillcontinents()
m.drawcoastlines()
#m.drawcountries()

levels = np.arange(-15,15,0.5)
im = m.contourf(lons,lats, stemp_10,levels,cmap=cm.seismic, latlon=True)
#im = m.pcolormesh(lons,lats, ftmp_10,shading='nearest',latlon=True,cmap="bwr")
cf = m.contour(lons,lats, stemp_10,levels, cmap=cm.binary)
#cbar = m.colorbar(im,pad='10%',format="%.1f",location='bottom',extend = "both")
cbar = m.colorbar(im,fraction=0.05,pad='4%',location='right',extend = "both")
#cbar.set_label(label="[hpa]",size=16)    #controls the size of the title of kthe colorbar
cbar.ax.tick_params(labelsize=16)   #controls the size of the numbers on the legend/colorbar
#fig.tight_layout()
#cf = ax.contour(lons,lats, ftmp_10,  levels , cmap=cm.binary,latlon=True)
# set the limits of the plot to the limits of the data
#axis([lons.min(),lons.max(), lats.max(), lats.min()])


#ax.set_title('SON Temperature Trend(with CO2 forcing): 1900-2015',
#             loc='center', fontsize=16,fontname='Times New Roman', fontweight='bold')
plt.title('SON Heat transpot at 200hpa - SPEEDY: 1950-2021',
             loc='center', fontsize=18, fontweight='bold')
#plt.title('SON Heat transpot at 200hpa - ERA5: 1950-2021',
#             loc='center', fontsize=18, fontweight='bold')
#plt.title('Annual',loc='right', fontsize=18, fontweight='bold')
plt.xticks(fontsize=14, weight='bold')
#ax.set_ylabel('K',fontsize=16)
plt.yticks(fontsize=14, weight='bold')

cbar.set_label(label="[degK m/s]",size=16)    #controls the size of the title of kthe colorbar
cbar.ax.tick_params(labelsize=16)   #controls the size of the numbers on the legend/colorbar
fig.tight_layout()

#pad controls how the color moves close/away from the map
#my_fig = 'Annual'
my_fig ='SON@200hpa_220'
plt.savefig('Output_heat/{}'.format(my_fig))


plt.show()
