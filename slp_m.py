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
#dir1 = "/media/michael/Maxtor/Mixed_Layer_Run/"
dir1 = "/Users/michaeloyelakin/Desktop/Mixed_Layer_Run/"


file_name = ("attm302.nc") 
ftmp= xr.open_dataset(dir1+file_name)


lon1 = ftmp['lon'][:] 
lat1 = ftmp['lat'][:] 
slp1 = ftmp['olr'][:]

#temp1 = ftmp['covvt'][:]

num_year=12
#slp_mean = np.mean(slp1, axis=0)

slp_new=slp1*num_year
ftmp_10 = slp_new[0,0,:,:]

#print(np.min(slp_new))
#print(np.max(slp_new))
#quit()

fig = plt.figure(figsize=(12,6))
m = Basemap(projection='mill',
            llcrnrlat=-90,
            urcrnrlat=90,
            llcrnrlon=-180,
            urcrnrlon=180,
            resolution='c')
[lons,lats] = np.meshgrid(lon1, lat1)
#x,y = m(lons,lats)

m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,0], fontsize=15)
m.drawmeridians(np.arange(-180.,181.,60.), labels=[0,0,0,1], fontsize=15)
m.drawmapboundary()
m.drawcoastlines()
m.drawcountries()

levels = np.arange(-3,3,0.2)
im = m.contourf(lons,lats, ftmp_10,levels,cmap='bwr', latlon=True)
#im = m.pcolormesh(lons,lats, ftmp_10,shading='nearest',latlon=True,cmap="bwr")
cf = m.contour(lons,lats, ftmp_10,levels, cmap=cm.binary)
#cbar = m.colorbar(im,pad='10%',format="%.1f",location='bottom',extend = "both")
cbar = m.colorbar(im,fraction=0.05,pad='10%',location='right',extend = "both")
cbar.set_label(label="[hpa]",size=16)    #controls the size of the title of kthe colorbar
cbar.ax.tick_params(labelsize=16)   #controls the size of the numbers on the legend/colorbar
#fig.tight_layout()
#cf = ax.contour(lons,lats, ftmp_10,  levels , cmap=cm.binary,latlon=True)
# set the limits of the plot to the limits of the data
#axis([lons.min(),lons.max(), lats.max(), lats.min()])


#ax.set_title('SON Temperature Trend(with CO2 forcing): 1900-2015',
#             loc='center', fontsize=16,fontname='Times New Roman', fontweight='bold')
plt.title('SON Mean sea level pressure [attm204]: 1950-2014',
             loc='center', fontsize=16, fontweight='bold')
#plt.title('Annual',loc='right', fontsize=16, fontweight='bold')

#pad controls how the color moves close/away from the map
#my_fig = 'Annual'
my_fig ='SON_mslp_204test'
plt.savefig('Output_Mixed/{}'.format(my_fig))


plt.show()
