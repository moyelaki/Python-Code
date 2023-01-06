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

#Read the file
dir1 = "/media/michael/Maxtor/HeatFlux/"
#dir1 = "/media/michael/Maxtor/Model/HeatFlux/"
#file_name = ("zonmean.220.heatflux.SON.nc") 
#file_name = ("zonmean_ncep_heatflux_SON.nc")
file_name = ("zonmean_attm220_vt_SON.nc")
ftmp= xr.open_dataset(dir1+file_name)


lon1 = ftmp['lon'][:] 
lat1 = ftmp['lat'][:] 
lev1 = ftmp['lev'][:]    
temp1 = ftmp['v'][:]
#temp1 = ftmp['covvt'][:]

nyear=62

#Extract first 10 levels
lev_10 = lev1[:]
ftmp_10 = temp1[0,:,:,0]

#Interpolating
levn = np.linspace(100,1000,num=41)
latn = np.linspace(0,90, num=181)

func = interp2d(lat1, lev_10, ftmp_10, kind='cubic')

# apply to new level and latitude
tnew = func(latn, levn)

#Multiply by number of years
tnew = tnew*nyear
#tnew = tnew

#Visualize the interpolated zonal mean
min_t = floor(np.min(tnew))
max_t = ceil(np.max(tnew))


[lats,levs] = np.meshgrid(latn, levn)
fig, ax = plt.subplots()


#print(tnew.min())
#print(tnew.max())
#quit()

levels = np.arange(-3,2.5,0.1)
#levels=[-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3]

im = ax.contourf( lats,levs, tnew, levels, cmap=cm.seismic)
cf = ax.contour( lats,levs, tnew,  levels , cmap=cm.binary)


ax.set_title('Eddy Heatflux trend [attm220] - SPEEDY: 1960-2021',
             loc='left', fontsize=18, fontweight='bold')
#ax.set_title('Heat flux trend attm209: 1900-2015',
#             loc='center', fontsize=16,fontname='Times New Roman', fontweight='bold')
ax.set_title('Annual',loc='right',fontweight='bold', fontsize=18)
ax.set_xlabel('Latitude [$^o$]', fontsize=16)
plt.xticks(fontsize=14, weight='bold')
ax.set_ylabel('Pressure Level [hPa]',fontsize=16)
plt.yticks(fontsize=14, weight='bold')

# set the limits of the plot to the limits of the data
axis([lats.min(),lats.max(), levs.max(), levs.min()])

#fig.colorbar(im,fraction=0.05, pad=0.05,size=15,label="Temp[K]")

#pad controls how the color moves close/away from the map
cbar = plt.colorbar(im,fraction=0.05, pad=0.04,orientation="vertical",
                    format="%.1f",extend = "both")     #'format' controls the number of decimals you have in your parameter.
cbar.set_label(label="[degK m/s]",size=16)    #controls the size of the title of kthe colorbar
cbar.ax.tick_params(labelsize=16)   #controls the size of the numbers on the legend/colorbar
fig.tight_layout()

#my_fig = 'Annual'
my_fig = 'attm220_vt_ANN'
plt.savefig('Output_heat/{}'.format(my_fig))


plt.show()
