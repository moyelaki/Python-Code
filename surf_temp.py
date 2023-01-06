from pylab import *
import os
#from scipy.interpolate import interp2d
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


num_year_attm220=72     #Expname attm220 ---sea ice+CO2


'''
    Era5 represent the Surface temperature
    ncep represent the near surface air temp
'''

#Read the file
#dir1 = "/media/michael/Maxtor/Mixed_Layer_Run/"
dir1 = "/Users/michaeloyelakin/Desktop/Mixed_Layer_Run/"

#__Annual_______________________________________________________________________________________________
file1a = ("zonmean.mixed.near.surf.air.temp.302.NH.annual.nc")
ncep_ann= xr.open_dataset(dir1+file1a)
lat1 = ncep_ann['lat'][:] 
tempa = ncep_ann['temp0'][:]
tempa_new=tempa*num_year_attm220
ftmpa_10 = tempa_new[0,:,0]

#++++++ERA5
file2a = ("zonmean.mixed.surf.temp.302.NH.annual.nc")
era5_ann= xr.open_dataset(dir1+file2a)
lat2 = era5_ann['lat'][:] 
t2m_ann = era5_ann['st'][:]
t2m_ann_new=t2m_ann*num_year_attm220
t2m_ann_10 = t2m_ann_new[0,:,0]

#__DJF___________________________________________________________________________________________________
file1d = ("zonmean.mixed.near.surf.air.temp.302.NH.DJF.nc")
ncep_djf= xr.open_dataset(dir1+file1d)
tempd = ncep_djf['temp0'][:]
tempd_new=tempd*num_year_attm220
ftmpd_10 = tempd_new[0,:,0]

#++++++ERA5 
file2d = ("zonmean.mixed.surf.temp.302.NH.DJF.nc")
era5_djf= xr.open_dataset(dir1+file2d)
t2m_djf = era5_djf['st'][:]
t2m_djf_new=t2m_djf*num_year_attm220
t2m_djf_10 = t2m_ann_new[0,:,0]

#_MAM___________________________________________________________________________________________________
file1m = ("zonmean.mixed.near.surf.air.temp.302.NH.MAM.nc")
ncep_mam= xr.open_dataset(dir1+file1m)
tempm = ncep_mam['temp0'][:]
tempm_new=tempm*num_year_attm220
ftmpm_10 = tempm_new[0,:,0]

#++++++ERA5 
file2m = ("zonmean.mixed.surf.temp.302.NH.MAM.nc")
era5_mam= xr.open_dataset(dir1+file2m)
t2m_mam = era5_mam['st'][:]
t2m_mam_new=t2m_mam*num_year_attm220
t2m_mam_10 = t2m_ann_new[0,:,0]


#_JJA___________________________________________________________________________________________________
file1j = ("zonmean.mixed.near.surf.air.temp.302.NH.JJA.nc")
ncep_jja= xr.open_dataset(dir1+file1j)
tempj = ncep_jja['temp0'][:]
tempj_new=tempj*num_year_attm220
ftmpj_10 = tempj_new[0,:,0]


#++++++ERA5 
file2j = ("zonmean.mixed.surf.temp.302.NH.JJA.nc")
era5_jja= xr.open_dataset(dir1+file2j)
t2m_jja = era5_jja['st'][:]
t2m_jja_new=t2m_jja*num_year_attm220
t2m_jja_10 = t2m_jja_new[0,:,0]


#SON____________________________________________________________________________________________________
file1s = ("zonmean.mixed.near.surf.air.temp.302.NH.SON.nc")
ncep_son= xr.open_dataset(dir1+file1s)
temps = ncep_son['temp0'][:]
temps_new=temps*num_year_attm220
ftmps_10 = temps_new[0,:,0]

#++++++ERA5 
file2s = ("zonmean.mixed.surf.temp.302.NH.SON.nc")
era5_son= xr.open_dataset(dir1+file2s)
t2m_son = era5_son['st'][:]
t2m_son_new=t2m_son*num_year_attm220
t2m_son_10 = t2m_son_new[0,:,0]
#_____________________________________________________________________________________________________

#fig, (ax1,ax2) = plt.subplots(1,2,sharex=True,sharey=True)    #makes the subplot share thesame x and y axis
fig, (ax1,ax2) = plt.subplots(1,2)
plt.suptitle('Zonal mean Air-temp Change - SPEEDY(Mixed layer): 1950-2021',
             fontsize=18, fontweight='bold')
#plt.title('Annual',loc='right', fontsize=16, fontweight='bold')

ax1.plot(lat1,ftmpd_10, label='DJF')
ax1.plot(lat1,ftmpm_10, label='MAM')
ax1.plot(lat1,ftmpj_10, label='JJA')
ax1.plot(lat1,ftmps_10, label='SON')
ax1.plot(lat1,ftmpa_10, label='ANN')
ax1.set_title('Near-surface air temp- trend: 1950-2021', fontsize=16,fontweight='bold')
ax1.set_xlabel('Latitude [$^o$]',fontsize=16)
#ax1.set_xticks(fontsize=14, weight='bold')
ax1.set_ylabel('K',fontsize=16)
#ax1.set_yticks(fontsize=14, weight='bold')



ax2.plot(lat2,t2m_djf_10, label='DJF')
ax2.plot(lat2,t2m_mam_10, label='MAM')
ax2.plot(lat2,t2m_jja_10, label='JJA')
ax2.plot(lat2,t2m_son_10, label='SON')
ax2.plot(lat2,t2m_ann_10, label='ANN')
ax2.set_title('Surface temperature trend: 1950-2021', fontsize=16,fontweight='bold')
ax2.set_xlabel('Latitude [$^o$]',fontsize=16)
ax2.set_ylabel('K',fontsize=16)
ax2.legend()
#ax2.set_xticks(fontsize=14, weight='bold')
#ax2.set_yticks(fontsize=14, weight='bold')


#for ax in ax.flat:
#    ax.set(xlabel='Latitude [$^o$]',ylabel='K')
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
#for ax in ax1,ax2.flat:
#    ax.label_outer()
#plt.xticks(fontsize=14, weight='bold')
#ax.set_ylabel('K',fontsize=16)
#
#plt.yticks(fontsize=14, weight='bold')
ax1.legend()
ax2.legend()
plt.tight_layout()




#pad controls how the color moves close/away from the map
#my_fig = 'Annual'
plt.show()
my_fig ='surface_temp_attm302'

plt.savefig('Output_mixed/{}'.format(my_fig))



