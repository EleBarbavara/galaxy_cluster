# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 09:27:39 2023

@author: eleobar
"""
#definition
import numpy as np
from astropy.io import fits
from astropy.table import Table
import tabulate
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm
#from tqdm import tqdm
import statistics
from astroplan import Observer, FixedTarget
from astropy.time import Time
import astropy.units as u #u.hour è l'unità in ore
from astropy.coordinates import AltAz, EarthLocation, Angle, get_sun, get_moon, SkyCoord
from astropy.visualization import astropy_mpl_style, quantity_support

#Fot matplotlib, use a nicer set of plot parameters and set up support for plotting/converting quantities
plt.style.use(astropy_mpl_style)
quantity_support()

'''
Suppose you are planning to visit the Sardinia Radio Telescope (Lat. 39°29′34″N - Long. 9°14′42″E - Alt. 600 m). 
You want to plan the observation at 9:00 pm local time, and you want 
to know if the calibrator/cluster/cluster pair will be up.

Minimum observational altitude = 30°
Maximum observational altitude = 80°
'''

# SINGLE CLUSTER
'''
data = Table.read('ACT_DR5_R500.txt', format='ascii.csv', delimiter='\t', header_start=0) #header = #Name	RA_deg	DEC_deg	SNR	z	z_type	z_source	M500c	M500cCal	M500cUncorr	R500_Mpc	R500_arcmin	notes	warnings
# ------------------------------------------------------------------------------------------------- 
#                    LIST NAME COLUMNS and FORMAT of CLUSTER in ACT_DR5_R500.txt
# -------------------------------------------------------------------------------------------------
#SNR_cluster = data['SNR']
#z_cluster = data['z']
#z_type_cluster = data['z_type']
#z_source_cluster = data['z_source']
#M500c_cluster = data['M500c']
#M500cCal_cluster = data['M500cCal']
#M500cUncorr_cluster = data['M500cUncorr']
#R500_Mpc_cluster = data['R500_Mpc']
#sep_cluster = data['R500_arcmin']
#
# -------------------------------------------------------------------------------------------------
name_cluster = data['#Name']
RA_cluster = data['RA_deg']
dec_cluster = data['DEC_deg']
sep_cluster = data['R500_arcmin']
mask_c = sep_cluster < 2
sep_cluster_for_mistral = sep_cluster[mask_c]
name_for_mistral = name_cluster[mask_c]
max_sep_cluster = round(max(sep_cluster_for_mistral), 4) #round(x, n) arrotonda x a n cifre significative
min_sep_cluster = round(min(sep_cluster_for_mistral), 4)
mean_sep_cluster = round(statistics.mean(sep_cluster_for_mistral), 4)
print('----------------MISTRAL observable clusters--------------------')
print('Number of visible cluster (FOV < 4 arcmin) = ', len(sep_cluster_for_mistral),'/',len(sep_cluster))
print(f'Mean angular diameter = {2*mean_sep_cluster:.2f} arcmin')
print(f'Maximum angular diameter = {2*max_sep_cluster:.2f} arcmin')
print(f'Minimum angular diameter = {2*min_sep_cluster:.2f} arcmin')
print('---------------------------------------------------------------')
print(name_for_mistral)
'''

# MULTIPLE SYSTEMS
multiple_systems = 'DR5_multiple-systems_v1.0.fits'
hdul1 = fits.open(multiple_systems, memmap=True)
ms = hdul1[1]
ms_data = Table(ms.data)
hdul1.close()
# ------------------------------------------------------------------------------------------------- 
#                    LIST NAME COLUMNS and FORMAT of MS in DR5_multiple-systems_v1.0.fits
# -------------------------------------------------------------------------------------------------
#    name = 'name'; format = '64A'
#    name = 'meanRADeg'; format = 'D'
#    name = 'meanDecDeg'; format = 'D'
#    name = 'meanRedshift'; format = 'D'
#    name = 'meanRedshift'; format = 'D'
#    name = 'numClusters'; format = 'K'
#    name = 'maxSeparationMpc'; format = 'D'
#    name = 'maxSeparationArcmin'; format = 'D'
#    name = 'includesPhotoRedshift'; format = 'L'
#    name = 'tileName'; format = '17A'
#
# -------------------------------------------------------------------------------------------------
name_tot = ms_data['name']
sep_systems = ms_data['maxSeparationArcmin']
num_systems_tot= ms_data['numClusters']
mask_s = sep_systems < 4
sep_systems_for_mistral = sep_systems[mask_s]
name_for_mistral = name_tot[mask_s]
num_systems = num_systems_tot[mask_s]
ra_obj = ms_data['meanRADeg'][mask_s]
dec_obj = ms_data['meanDecDeg'][mask_s]
max_sep_syst = round(max(sep_systems_for_mistral), 4)
min_sep_syst = round(min(sep_systems_for_mistral), 4)
mean_sep_syst = round(statistics.mean(sep_systems_for_mistral), 4)
print('------------MISTRAL observable multiple systems---------------------')
print('Number of visible multiple system (FOV < 4 arcmin) = ', len(sep_systems_for_mistral),'/',len(sep_systems))
print('Number of clusters for each visible multiple system =', [num_systems[i] for i in range(0, len(sep_systems_for_mistral))])
print(f'Mean angular diameter = {mean_sep_syst:.2f} arcmin')
print(f'Minimum angular diameter = {min_sep_syst:.2f} arcmin')
print(f'Maximum angular diameter = {max_sep_syst:.2f} arcmin')
print('---------------------------------------------------------------------')

SRT_lat = Angle('39°29′34″')
SRT_lon = Angle('9°14′42″E')
SRT_alt = 600
SRT = EarthLocation(lat=SRT_lat, lon=SRT_lon, height=SRT_alt*u.m) #EartLocation.of_site('srt')
loc= Observer(location=SRT)
utcoffset = +2*u.hour  #CEST
time = Time('2024-2-1 21:00:00') - utcoffset
delta_midnight = np.linspace(-12, 12, 1000)*u.hour
time_obs_night = time + delta_midnight
time_obs_range = Time(["2024-2-01 12:00", "2024-2-2 12:00"])
time_range_mistral = Time(["2023-7-01 00:00", "2024-2-29 00:00"])
#%%
#istogramma dell'estensione dei sistemi
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
ax1.hist(sep_systems, bins=30, ec='grey', alpha=0.5, align='mid')
ax1.set_xlabel('Separation [arcmin]')
ax1.set_xticks(ticks=range(0, 100, 10))
ax1.set_ylabel('N')
ax1.grid(False)
#istogramma dei quanti cluster sono in un sistema
counts = np.bincount(num_systems_tot)
bar = ax2.bar(range(5), counts, width=1, align='center', ec='grey', alpha=0.5, label=counts)
ax2.bar_label(bar)
ax2.set_xticks(ticks=range(5), minor=False)
ax2.set_xlim([1.1,5])
ax2.set_ylim([0,160])
ax2.set_xlabel('Number of clusters')
ax2.set_ylabel('N')
ax2.grid(False)
plt.show()

#ec = 'grey' colore bordo barre, fc='white' colore riempimento barra


#%%
coord=[SkyCoord(ra=ra_obj[i]*u.deg, dec=dec_obj[i]*u.deg) for i in range(0, len(sep_systems_for_mistral))]
targets = [FixedTarget(coord[i], name=name_for_mistral[i])  for i in range(0, len(sep_systems_for_mistral))]

from astroplan import AltitudeConstraint#, AirmassConstraint, AtNightConstraint
from astroplan import is_observable, is_always_observable, months_observable, observability_table
constraints = [AltitudeConstraint(30*u.deg, 80*u.deg)]#, AirmassConstraint(5), AtNightConstraint.twilight_civil()]
# Are targets *ever* observable in the time range? -> output: boolean array of lenght=number object for whether or not each target is ever observable in the time range given the constraints 
ever_observable = is_observable(constraints, loc, targets, time_range=time_range_mistral) 

# Are targets *always* observable in the time range?-> output: boolean array of lenght=number object for whether or not each target is ever observable in the time range given the constraints 
always_observable = is_always_observable(constraints, loc, targets, time_range=time_range_mistral) 

# During what months are the targets ever observable? -> output: array of best month in which the target is observable in the current year of time_obs
best_months = months_observable(constraints, loc, targets) 

table_m = observability_table(constraints, loc, targets, time_range=time_range_mistral)

observability_table_m = Table()
observability_table_m['Ever Observable'] = ever_observable
observability_table_m['Always Observable'] = always_observable
observability_table_m['Best Months'] = best_months

'''
ra_c = RA_cluster[mask_c]
dec_c = dec_cluster[mask_c]
gal_c = SkyCoord(ra_c[:], dec_c[:], frame='galactic', unit=u.deg) #coord equatoriali!!!
'''
gal_s = SkyCoord(ms_data['meanRADeg'][:], ms_data['meanDecDeg'][:], frame='galactic', unit=u.deg)
gal_s_v = SkyCoord(ra_obj[:], dec_obj[:], frame='galactic', unit=u.deg) #coord equatoriali!!!
plt.subplot(1,1,1, projection='aitoff')
#plt.scatter(gal_c.l.wrap_at('180d').radian, gal_c.b.radian, s=.1, color='orange')#, alpha=0.05)
plt.scatter(gal_s.l.wrap_at('180d').radian, gal_s.b.radian, c='grey', alpha=.5)
plt.scatter(gal_s_v.l.wrap_at('180d').radian, gal_s_v.b.radian, c='blue')#, c='lightseagreen'), s=1, color='black')#, alpha=0.05)
plt.grid(True)
#plt.title('Map of multiple systems', fontsize=12, x=.5, y=1.1)
#plt.savefig('Map_visible_cluster.png')
plt.show()

print('\n')
print(table_m)
print('\n')
print(observability_table_m)
#%%
from astroplan.plots import plot_sky #, style_sheet = 'dark_style_sheet'
[plot_sky(targets[i], loc, time_obs_night) for i in range(0, len(sep_systems_for_mistral))]
#plot_sky(moon, loc, time_obs_night)
#plot_sky(sun, loc, time_obs_night) 
plt.legend(loc='center left', bbox_to_anchor=(1.1, .5)) #, bbox_to_anchor=(1.25, 0.5) per legenda a destra
#plt.title(f'Target positions in the sky with respect to the observer’s location\n({time_obs_night[0]} to {time_obs_night[-1]})', fontsize=12, x=0.5, y= 1.1)
#[plt.scatter(delta_midnight, targets[i].alt, label=cluster, lw=0, s=8,) for i in range(0, len(sep_systems_for_mistral))]
#plt.savefig('Pos_sky_cluster_pair.png')
plt.show()
#%%
#example of observation of a choosen cluster pair at a certain date at a certain time
print('-------------------------------------------------------------------------')
print(f'             Observation of {name_for_mistral[0]}')
print('-------------------------------------------------------------------------')
obj = coord[0]
coords = obj.transform_to(AltAz(obstime=time,location=SRT)) #trasf in AlzAz coordinates 
print('')
print('Date and Time of observation:', time+utcoffset)
print(f"{name_for_mistral[0][7:]}'s Altitude = {coords.alt:.2f}")
print(f"{name_for_mistral[0][7:]}'s Azimuth = {coords.az:.2f}")
print('----------------------------------------------------------')

#Find the airmass at 100 times evenly spaced between 10pm and 7am CEST:
midnight = Time('2024-2-1 00:00:00') - utcoffset
delta_midnight = np.linspace(-2, 10, 100)*u.hour
frame_night = AltAz(obstime=midnight+delta_midnight, location=SRT)
coords_night = obj.transform_to(frame_night)
airmasss_night = coords_night.secz
'''
plt.plot(delta_midnight, airmasss_night)
plt.xlim(-2, 10)
plt.ylim(1, 4)
plt.xlabel('Hours from CEST Midnight')
plt.ylabel('Airmass [Sec(z)]')
plt.show()
'''

#Find the visibility of a calibrator at 100 times evenly spaced between 12am of a day and 12am of the following day CEST:
delta_midnight = np.linspace(-12, 12, 1000)*u.hour
time_obs_night = midnight + delta_midnight
frame_obs_night = AltAz(obstime=time_obs_night, location=SRT)
coord_obs_night = obj.transform_to(frame_obs_night)

suncoord_obs_night = get_sun(time_obs_night).transform_to(frame_obs_night)
mooncoord_obs_night = get_moon(time_obs_night).transform_to(frame_obs_night)

plt.plot(delta_midnight, suncoord_obs_night.alt, color='gold', label='Sun')
plt.scatter(delta_midnight, mooncoord_obs_night.alt, c=mooncoord_obs_night.az, lw=0, s=8, label='Moon', cmap='gist_gray')# color=[0.75]*3,, ls='--'
plt.colorbar().set_label('Azimuth [deg]')
plt.axhline(30)
plt.axhline(80)
plt.axhspan(0, 30, facecolor='grey', alpha=0.2)
plt.axhspan(80, 90, facecolor='grey', alpha=0.2)
plt.scatter(delta_midnight, coord_obs_night.alt, c=coord_obs_night.az, label=name_for_mistral[0][7:], lw=0, s=8, cmap='viridis')
plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg, suncoord_obs_night.alt < -0*u.deg, color='0.5', zorder=0)
plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg, suncoord_obs_night.alt < -18*u.deg, color='k', zorder=0)
plt.colorbar().set_label('Azimuth [deg]')
plt.legend(loc='upper left', bbox_to_anchor=(0.1, 1.25))
plt.xlim(-12*u.hour, 12*u.hour)
plt.xticks((np.arange(13)*2-12)*u.hour)
plt.ylim(0*u.deg, 90*u.deg)
plt.xlabel('Hours from CEST Midnight')
plt.ylabel('Altitude [deg]')
#plt.title(f'Visibily of {name_for_mistral[0]} (2023/7/1 12:00 to 2023/7/2 12:00)', fontsize=12, x=0.7, y= 1.25)
#plt.savefig('Visibility_cluster.png')
plt.show()

print('-------------------------------------------------------------------------')
print(f'             Observation of {name_for_mistral[1]}')
print('-------------------------------------------------------------------------')
obj = coord[1]
coords = obj.transform_to(AltAz(obstime=time,location=SRT)) #trasf in AlzAz coordinates 
print('')
print('Date and Time of observation:', time+utcoffset)
print(f"{name_for_mistral[0][7:]}'s Altitude = {coords.alt:.2f}")
print(f"{name_for_mistral[0][7:]}'s Azimuth = {coords.az:.2f}")
print('----------------------------------------------------------')

#Find the airmass at 100 times evenly spaced between 10pm and 7am CEST:
midnight = Time('2024-2-1 00:00:00') - utcoffset
delta_midnight = np.linspace(-2, 10, 100)*u.hour
frame_night = AltAz(obstime=midnight+delta_midnight, location=SRT)
coords_night = obj.transform_to(frame_night)
airmasss_night = coords_night.secz
'''
plt.plot(delta_midnight, airmasss_night)
plt.xlim(-2, 10)
plt.ylim(1, 4)
plt.xlabel('Hours from CEST Midnight')
plt.ylabel('Airmass [Sec(z)]')
plt.show()
'''

#Find the visibility of a calibrator at 100 times evenly spaced between 12am of a day and 12am of the following day CEST:
delta_midnight = np.linspace(-12, 12, 1000)*u.hour
time_obs_night = midnight + delta_midnight
frame_obs_night = AltAz(obstime=time_obs_night, location=SRT)
coord_obs_night = obj.transform_to(frame_obs_night)

suncoord_obs_night = get_sun(time_obs_night).transform_to(frame_obs_night)
mooncoord_obs_night = get_moon(time_obs_night).transform_to(frame_obs_night)

plt.plot(delta_midnight, suncoord_obs_night.alt, color='gold', label='Sun')
plt.scatter(delta_midnight, mooncoord_obs_night.alt, c=mooncoord_obs_night.az, lw=0, s=8, label='Moon', cmap='gist_gray')# color=[0.75]*3,, ls='--'
plt.colorbar().set_label('Azimuth [deg]')
plt.axhline(30)
plt.axhline(80)
plt.axhspan(0, 30, facecolor='grey', alpha=0.2)
plt.axhspan(80, 90, facecolor='grey', alpha=0.2)
plt.scatter(delta_midnight, coord_obs_night.alt, c=coord_obs_night.az, label=name_for_mistral[1][7:], lw=0, s=8, cmap='viridis')
plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg, suncoord_obs_night.alt < -0*u.deg, color='0.5', zorder=0)
plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg, suncoord_obs_night.alt < -18*u.deg, color='k', zorder=0)
plt.colorbar().set_label('Azimuth [deg]')
plt.legend(loc='upper left', bbox_to_anchor=(0.1, 1.25))
plt.xlim(-12*u.hour, 12*u.hour)
plt.xticks((np.arange(13)*2-12)*u.hour)
plt.ylim(0*u.deg, 90*u.deg)
plt.xlabel('Hours from CEST Midnight')
plt.ylabel('Altitude [deg]')
#plt.title(f'Visibily of {name_for_mistral[0]} (2023/7/1 12:00 to 2023/7/2 12:00)', fontsize=12, x=0.7, y= 1.25)
#plt.savefig('Visibility_cluster.png')
plt.show()

# prova osservazione triangolo australe
'''
tri = SkyCoord(ra= 16*u.hour, dec=-65*u.deg)
delta_midnight = np.linspace(-12, 12, 1000)*u.hour
time_obs_night = midnight + delta_midnight
frame_obs_night = AltAz(obstime=time_obs_night, location=SRT)
coord_obs_night = tri.transform_to(frame_obs_night)

suncoord_obs_night = get_sun(time_obs_night).transform_to(frame_obs_night)
mooncoord_obs_night = get_moon(time_obs_night).transform_to(frame_obs_night)

plt.plot(delta_midnight, suncoord_obs_night.alt, color='gold', label='Sun')
plt.scatter(delta_midnight, mooncoord_obs_night.alt, c=mooncoord_obs_night.az, lw=0, s=8, label='Moon', cmap='gist_gray')# color=[0.75]*3,, ls='--'
plt.colorbar().set_label('Azimuth [deg]')
plt.axhline(30)
plt.axhline(80)
plt.axhspan(0, 30, facecolor='grey', alpha=0.2)
plt.axhspan(80, 90, facecolor='grey', alpha=0.2)
plt.scatter(delta_midnight, coord_obs_night.alt, c=coord_obs_night.az, label=name_for_mistral[0][7:], lw=0, s=8, cmap='viridis')
plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg, suncoord_obs_night.alt < -0*u.deg, color='0.5', zorder=0)
plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg, suncoord_obs_night.alt < -18*u.deg, color='k', zorder=0)
plt.colorbar().set_label('Azimuth [deg]')
plt.legend(loc='upper left', bbox_to_anchor=(0.1, 1.25))
plt.xlim(-12*u.hour, 12*u.hour)
plt.xticks((np.arange(13)*2-12)*u.hour)
plt.ylim(0*u.deg, 90*u.deg)
plt.xlabel('Hours from CEST Midnight')
plt.ylabel('Altitude [deg]')
plt.title(f'Visibily of {name_for_mistral[0]} (2023/7/1 12:00 to 2023/7/2 12:00)', fontsize=12, x=0.7, y= 1.25)
plt.show()

gal= SkyCoord(tri.ra, tri.dec, frame='galactic', unit=u.deg)
plt.subplot(1,1,1, projection='aitoff')
plt.scatter(gal.l.wrap_at('180d').radian, gal.b.radian, color='orange')
plt.grid(True)
plt.show()

tri = [tri, coord[0]]
tri_name = ['tri', name_for_mistral[0]]
tri = [FixedTarget(tri[i], name=tri_name[i]) for i in range(0, len(tri))]
ever_observable1 = is_observable(constraints, loc, tri, time_range=time_obs_range)
always_observable1 = is_always_observable(constraints, loc, tri, time_range=time_obs_range)
best_months1 = months_observable(constraints, loc, tri)

table_m1 = observability_table(constraints, loc, tri, time_range=time_obs_range)
print(table_m1)
observability_table_m1 = Table()
observability_table_m1['Best Months'] = best_months1
print('\n')
print(observability_table_m1)
'''
#convert ra in degree to ra in hour
'''
RA_hms_cluster = []
for i in tqdm(range(0, len(RA_cluster))):
	ra_c= RA_cluster[i]
	dec_c= dec_cluster[i]
	coord = SkyCoord(ra = ra_c*u.degree, dec = dec_c*u.degree)
	ra_hms = coord.ra.hms
	cs = coord.to_string('hmsdms')
	stop = cs.find(' ')
	ra_hms = cs[0:stop]
	RA_hms_cluster.append(ra_hms)
'''

#2D histogram 
''' plot a 2d histogram in which the colorbar indicates how many cluster are in a certain position (pixel)
NBINS=(60,60) #defines the pixelization of the 2D histogram
img_zero_mpl = plt.hist2d(cluster_data['RADeg'], cluster_data['decDeg'], NBINS, cmap='viridis', norm=LogNorm())  #cmap='viridis' o 'gnuplot' (viole-giallo)
plt.xlabel('RA')
plt.ylabel('Dec')
plt.colorbar()
plt.show()
'''

# ------------------------------------------------------------------------------------------------- 
#                    LIST NAME COLUMNS and FORMAT of CLUSTER in DR5_cluster-catalog_v1.1.fits
#cluster_catalog = 'DR5_cluster-catalog_v1.1.fits'
#hdul0 = fits.open(cluster_catalog, memmap=True)
#cluster = hdul0[1]
#cluster_data = Table(cluster.data)
#hdul0.close()
# -------------------------------------------------------------------------------------------------
#	name = 'name'; format = '19A'                    name = 'M200mUncorr'; format = 'D'
#   name = 'RADeg'; format = 'D'                     name = 'M200mUncorr_errPlus'; format = 'D'
#   name = 'decDeg'; format = 'D'                    name = 'M200mUncorr_errMinus'; format = 'D'
#   name = 'SNR'; format = 'D'                       name = 'M500cUncorr_errMinus'; format = 'D'
#   name = 'y_c'; format = 'D'                       name = 'footprint_DESY3'; format = 'L'
#   name = 'err_y_c'; format = 'D'                   name = 'footprint_HSCs19a'; format = 'L'
#   name = 'fixed_SNR'; format = 'D'                 name = 'footprint_KiDSDR4'; format = 'L'
#   name = 'fixed_y_c'; format = 'D'                 name = 'zCluster_delta'; format = 'D'
#   name = 'fixed_err_y_c'; format = 'D'             name = 'zCluster_errDelta'; format = 'D'
#   name = 'template'; format = '25A'                name = 'zCluster_source'; format = '12A'
#   name = 'tileName'; format = '7A'                 name = 'RM'; format = 'L'
#   name = 'redshift'; format = 'D'                  name = 'RM_LAMBDA'; format = 'D'
#   name = 'redshiftErr'; format = 'D'               name = 'RM_LAMBDA_ERR'; format = 'D'
#   name = 'redshiftType'; format = '1000A'          name = 'RMDESY3'; format = 'L'
#   name = 'redshiftSource'; format = '1000A'        name = 'RMDESY3_LAMBDA_CHISQ'; format = 'D'
#   name = 'M500c'; format = 'D'                     name = 'RMDESY3_LAMBDA_CHISQ_E'; format = 'D'          
#   name = 'M500c_errPlus'; format = 'D'             name = 'CAMIRA'; format = 'L'
#   name = 'M500c_errMinus'; format = 'D'            name = 'CAMIRA_N_mem'; format = 'D'
#   name = 'M500cCal'; format = 'D'                  name = 'opt_RADeg'; format = 'D'
#   name = 'M500cCal_errPlus'; format = 'D'          name = 'opt_decDeg'; format = 'D'
#   name = 'M500cCal_errMinus'; format = 'D'         name = 'opt_positionSource'; format = '10A'  
#   name = 'M200m'; format = 'D'                     name = 'notes'; format = '26A'      
#   name = 'M200m_errPlus'; format = 'D'             name = 'knownLens'; format = '42A'
#   name = 'M200m_errMinus'; format = 'D'            name = 'knownLensRefCode'; format = '14A'
#   name = 'M500cUncorr'; format = 'D'               name = 'warnings'; format = '93A'
#   name = 'M500cUncorr_errPlus'; format = 'D'
#
# -------------------------------------------------------------------------------------------------
