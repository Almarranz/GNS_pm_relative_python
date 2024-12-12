#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 13:24:28 2024

@author: amartinez
"""


import numpy as np
import matplotlib.pyplot as plt
from compare_lists import compare_lists 
from astropy.table import Table
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
import matplotlib.colors as colors_plt
from skimage import data
from skimage import transform
import astroalign as aa
from astroquery.gaia import Gaia
import skimage as ski
import sys
from astropy.stats import sigma_clip
from astropy.io import fits
from astropy.coordinates import SkyCoord
import Polywarp as pw
from astroquery.gaia import Gaia
from astropy import units as u
import cluster_finder
import pandas as pd
import copy
import cluster_finder
from filters import filter_gaia_data
from filters import filter_hosek_data
from filters import filter_gns_data
from filters import filter_vvv_data
# %% 
# %%plotting parametres
from matplotlib import rc
from matplotlib import rcParams
rcParams.update({'xtick.major.pad': '7.0'})
rcParams.update({'xtick.major.size': '7.5'})
rcParams.update({'xtick.major.width': '1.5'})
rcParams.update({'xtick.minor.pad': '7.0'})
rcParams.update({'xtick.minor.size': '3.5'})
rcParams.update({'xtick.minor.width': '1.0'})
rcParams.update({'ytick.major.pad': '7.0'})
rcParams.update({'ytick.major.size': '7.5'})
rcParams.update({'ytick.major.width': '1.5'})
rcParams.update({'ytick.minor.pad': '7.0'})
rcParams.update({'ytick.minor.size': '3.5'})
rcParams.update({'ytick.minor.width': '1.0'})
rcParams.update({'font.size': 20})
rcParams.update({'figure.figsize':(10,5)})
rcParams.update({
    "text.usetex": False,
    "font.family": "sans",
    "font.sans-serif": ["Palatino"]})
plt.rcParams["mathtext.fontset"] = 'dejavuserif'
rc('font',**{'family':'serif','serif':['Palatino']})
plt.rcParams.update({'figure.max_open_warning': 0})# a warniing for matplot lib pop up because so many plots, this turining it of
# %%
pm_folder = '/Users/amartinez/Desktop/PhD/HAWK/pm_gns1_gns2_relative_python/'

field_one, chip_one, field_two, chip_two,t1,t2,max_sig = np.loadtxt('/Users/amartinez/Desktop/PhD/HAWK/GNS_1relative_python/lists/fields_and_chips.txt', 
                                                       unpack=True)
field_one = field_one.astype(int)
chip_one = chip_one.astype(int)
field_two = field_two.astype(int)
chip_two = chip_two.astype(int)

# field_one = 7
# chip_one = 4
# field_two = 7
# chip_two = 1
# max_sig = 0.5#TODO

max_deg = 3
d_m_pm = 2
gns1_pm = Table.read(pm_folder + f'pm_ep1_f{field_one}c{chip_one}_ep2_f{field_two}c{chip_two}deg{max_deg}_dmax{d_m_pm}_sxy%.1f.txt'%(max_sig), format = 'ascii')

Ks_lim = [12,14.5]
# Ks_lim = [0,1000]
Ks_mask = (gns1_pm['Ks1']>Ks_lim[0])&(gns1_pm['Ks1']<Ks_lim[1])
gns1_pm = gns1_pm[Ks_mask]
#Gaia comparison

search_r = 100*u.arcsec

ra_c = np.mean(gns1_pm['ra1'])
dec_c = np.mean(gns1_pm['Dec1'])
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source" # Select early Data Release 3
Gaia.ROW_LIMIT = -1  # Ensure the default row limit.
coord = SkyCoord(ra=ra_c, dec=dec_c, unit='degree', frame='icrs')


j = Gaia.cone_search_async(coord, search_r)

gaia_ = j.get_results()


e_pm = 0.3
gaia_good = filter_gaia_data(
    gaia_table=gaia_,
    astrometric_params_solved=31,
    duplicated_source= False,
    parallax_over_error_min=-10,
    astrometric_excess_noise_sig_max=0.2,
    phot_g_mean_mag_min= 17.5,
    phot_g_mean_mag_max= 13,
    pm_min=0,
    pmra_error_max=e_pm,
    pmdec_error_max=e_pm
    )

gaia_gal = SkyCoord(ra=gaia_good['ra'], dec=gaia_good['dec'],
                        pm_ra_cosdec =gaia_good['pmra'], pm_dec = gaia_good['pmdec'],
                        unit = 'degree',frame = 'icrs', obstime='J2016.0').galactic
gaia_good.add_column(gaia_gal.pm_l_cosb, name = 'pml',index =-1)
gaia_good.add_column(gaia_gal.pm_b,name = 'pmb',index = -1)



# %
fig, ax = plt.subplots(1,1)
ax.scatter(gns1_pm['ra1'],gns1_pm['Dec1'])
ax.scatter(gaia_['ra'],gaia_['dec'])
ax.scatter(gaia_good['ra'],gaia_good['dec'])

# center = (gns1_pm['H1'] - gns1_pm['Ks1'] > 1.3)
# gns1_pm = gns1_pm[center]


gns1_coor = SkyCoord(ra = gns1_pm['ra1'],dec = gns1_pm['Dec1'],unit = 'degree',
                     frame = 'fk5', obstime = 'J2015.4301')
gaia_coord = SkyCoord(ra=gaia_good['ra'], dec=gaia_good['dec'],unit = 'degree',
                     frame = 'icrs', obstime = 'J2016.0')

max_sep = 0.05*u.arcsec#!!!



idx,d2d,d3d = gaia_coord.match_to_catalog_sky(gns1_coor,nthneighbor=1)# ,nthneighbor=1 is for 1-to-1 matchsep_constraint = d2d < max_sep
sep_constraint = d2d < max_sep
gns_match = gns1_pm[idx[sep_constraint]]
# pm_x_match = pm_x[idx[sep_constraint]]
# pm_y_match = pm_y[idx[sep_constraint]]
gaia_match = gaia_good[sep_constraint]

fig, ax = plt.subplots(1,1)
ax.scatter(gns_match['ra1'],gns_match['Dec1'])
ax.scatter(gaia_match['ra'],gaia_match['dec'],s =1)

diff_pmx = gns_match['pmx'] + gaia_match['pml']
diff_pmy = gns_match['pmy'] - gaia_match['pmb']
# diff_pmx = gns_match['pmx'] + gaia_match['pmra']
# diff_pmy = gns_match['pmy'] - gaia_match['pmdec']

mask_pmx, l_lim,h_lim = sigma_clip(diff_pmx, sigma=3, masked = True, return_bounds= True)
mask_pmy, l_lim,h_lim = sigma_clip(diff_pmy, sigma=3, masked = True, return_bounds= True)
diff_pmx_clip = diff_pmx[~mask_pmx.mask]
diff_pmy_clip = diff_pmy[~mask_pmy.mask]

print(np.mean(diff_pmx),np.std(diff_pmx))
print(np.mean(diff_pmy),np.std(diff_pmy))
    # max_sep = 0.08*u.arcsec
    # idx,d2d,d3d = gaia_coord.match_to_catalog_sky(hose_coord,nthneighbor=1)# ,nthneighbor=1 is for 1-to-1 match
    # sep_constraint = d2d < max_sep
    # hose_match = arches[idx[sep_constraint]]
    # gaia_match_hos = gaia_good[sep_constraint]
# %%
fig, (ax,ax1) = plt.subplots(1,2)
ax.set_title('Matching = %s'%(len(gns_match['pmx'])), fontsize = 15)
fig.suptitle('GAIA & GNS1[F%sC%s] GNS2[F%sC%s]'%(field_one, chip_one, field_two, chip_two), fontsize = 15)
ax.hist(diff_pmx_clip, histtype = 'step', label = '$\overline{\Delta x}$ =%.2f\n$\sigma$ = %.2f '%(np.nanmean(diff_pmx_clip),np.nanstd(diff_pmx_clip)))
ax.hist(diff_pmx, histtype = 'step', color ='k', alpha = 0.3)

ax1.hist(diff_pmy_clip, histtype = 'step',color = 'orange', label = '$\overline{\Delta y}$ =%.2f\n$\sigma$ = %.2f '%(np.nanmean(diff_pmy_clip),np.nanstd(diff_pmy_clip)))
ax1.hist(diff_pmy, histtype = 'step',color = 'k', alpha = 0.3, ls = 'dashed')

ax.set_xlabel('$\Delta \mu_{\parallel}$')
ax1.set_xlabel('$\Delta \mu_{\perp}$')
ax.legend()
ax1.legend()
sys.exit(181)
# %%
# Comparison with VVV 
VVV = Table.read('/Users/amartinez/Desktop/PhD/Catalogs/VVV/b333/PMS/b333.dat', format = 'ascii')
# %%
gns1_pm = Table.read(pm_folder + f'pm_ep1_f{field_one}c{chip_one}_ep2_f{field_two}c{chip_two}deg{max_deg}_dmax{d_m_pm}_sxy%.1f.txt'%(max_sig), format = 'ascii')

gns1_coor = SkyCoord(ra = gns1_pm['ra1'],dec = gns1_pm['Dec1'],unit = 'degree',
                     frame = 'fk5', obstime = 'J2015.4301')
rad = search_r.value/3600

delta_ra = (VVV['ra'] - ra_c) * np.cos(np.radians(dec_c))
delta_dec = VVV['dec'] - dec_c

# Calculate the angular distance
angular_distance = np.sqrt(delta_ra**2 + delta_dec**2)

# Select rows where the distance is within the radius
within_radius = angular_distance <= rad

vvv_c = VVV[within_radius]

vvv_c = filter_vvv_data(vvv_c,
                    pmRA = 'good',
                    pmDE = None,
                    epm = 1,
                    ok = 'yes',
                    max_Ks = Ks_lim[0],
                    min_Ks = Ks_lim[1],
                    center = None
                    )



fig, ax = plt.subplots(1,1)
ax.scatter(vvv_c['ra'],vvv_c['dec'])
ax.scatter(gns1_pm['ra1'], gns1_pm['Dec1'], alpha = 0.1)

# Moves coordinates to GNS1 obstime???
tvvv = 2012.29578304
# vvv_c['ra'] = vvv_c['ra'] + (t1-tvvv) * (vvv_c['pmRA']/1000.0)/3600.0 * np.cos(vvv_c['dec']*np.pi/180.)
# vvv_c['dec'] = vvv_c['dec'] + (t1-tvvv) * (vvv_c['pmDEC']/1000.0)/3600.0
vvv_coord = SkyCoord(ra = vvv_c['ra'], dec = vvv_c['dec'], unit = 'degree',
                      frame = 'icrs', obstime = 'J2012.29578304')
# vvv_coord = SkyCoord(ra = vvv_c['ra'], dec = vvv_c['dec'], unit = 'degree',
#                      frame = 'icrs', obstime = 'J2015.4301')




idx,d2d,d3d = vvv_coord.match_to_catalog_sky(gns1_coor,nthneighbor=1)# ,nthneighbor=1 is for 1-to-1 matchsep_constraint = d2d < max_sep
sep_constraint = d2d < max_sep
gns_vvv = gns1_pm[idx[sep_constraint]]
vvv_match = vvv_c[sep_constraint]



sig_cl = 3
dKs = gns_vvv['Ks1']-vvv_match['Ks']
mask_m, l_lim,h_lim = sigma_clip(dKs, sigma=sig_cl, masked = True, return_bounds= True)


fig, (ax, ax1) = plt.subplots(1,2)
ax.scatter(gns_vvv['ra1'],gns_vvv['Dec1'])
ax.scatter(vvv_match['ra'],vvv_match['dec'],s =1)
ax.set_title('Macthes = %s'%(len(vvv_match)))
ax1.hist(dKs,bins = 'auto', label = '$\Delta$Ks = %.2f\n$\sigma$ = %.2f'%(np.mean(dKs),np.std(dKs)))
ax1.axvline(l_lim, ls = 'dashed', color = 'r', label = '$\pm$ %s$\sigma$'%(sig_cl))
ax1.axvline(h_lim, ls = 'dashed', color = 'r')
ax1.legend(loc = 4, fontsize = 12)

vvv_match = vvv_match[~mask_m.mask]
gns_vvv = gns_vvv[~mask_m.mask]
ax1.set_title(f'{sig_cl}$\sigma$ Macthes= %s'%(len(vvv_match)))


vvv_gal = SkyCoord(ra = vvv_match['ra'], dec = vvv_match['dec'], unit = 'degree',
                   pm_ra_cosdec = vvv_match['pmRA']*u.mas/u.yr, pm_dec = vvv_match['pmDEC']*u.mas/u.yr, 
                      frame = 'icrs', obstime = 'J2012.29578304').galactic


vvv_match['l'] = vvv_gal.l
vvv_match['b'] = vvv_gal.b
vvv_match['pml'] = vvv_gal.pm_l_cosb
vvv_match['pmb'] = vvv_gal.pm_b




diff_pmx = gns_vvv['pmx'] + vvv_match['pml']
diff_pmy = gns_vvv['pmy'] - vvv_match['pmb']



mask_pmx, l_lim,h_lim = sigma_clip(diff_pmx, sigma=3, masked = True, return_bounds= True)
mask_pmy, l_lim,h_lim = sigma_clip(diff_pmy, sigma=3, masked = True, return_bounds= True)
diff_pmx_clip = diff_pmx[~mask_pmx.mask]
diff_pmy_clip = diff_pmy[~mask_pmy.mask]

# %
fig, (ax,ax1) = plt.subplots(1,2)
ax.set_title('Matching = %s'%(len(gns_vvv['pmx'])))
ax1.set_title('VVV comparison')
ax.hist(diff_pmx_clip, histtype = 'step', label = '$\overline{\Delta x}$ =%.2f\n$\sigma$ = %.2f '%(np.nanmean(diff_pmx_clip),np.nanstd(diff_pmx_clip)))
ax.hist(diff_pmx, histtype = 'step', color ='k', alpha = 0.3)

ax1.hist(diff_pmy_clip, histtype = 'step',color = 'orange', label = '$\overline{\Delta x}$ =%.2f\n$\sigma$ = %.2f '%(np.nanmean(diff_pmy_clip),np.nanstd(diff_pmy_clip)))
ax1.hist(diff_pmy, histtype = 'step',color = 'k', alpha = 0.3, ls = 'dashed')
ax.legend()
ax1.legend()
ax.set_xlabel('$\Delta \mu_{\parallel}$')
ax1.set_xlabel('$\Delta \mu_{\perp}$')
sys.exit(292)

# %%
 # Comaparison with Hosek
from astropy.io import ascii
catal='/Users/amartinez/Desktop/PhD/Arches_and_Quintuplet_Hosek/'

# choosen_cluster = 'Arches'
choosen_cluster = 'Quintuplet'
center_arc = SkyCoord(ra = '17h45m50.65020s', dec = '-28d49m19.51468s', equinox = 'J2000') if choosen_cluster =='Arches' else SkyCoord('17h46m 14.68579s', '-28d49m38.99169s', frame='icrs',obstime ='J2016.0')#Quintuplet

if choosen_cluster == 'Arches':
    arches = Table.read(catal + 'Arches_from_Article.txt', format = 'ascii')
if choosen_cluster == 'Quintuplet':
    arches = ascii.read(catal + 'Quintuplet_from_Article.txt')
    
RA_DEC = center_arc.spherical_offsets_by(arches['dRA'], arches['dDE'])
RA = RA_DEC.ra
DEC = RA_DEC.dec

arches.add_column(RA,name = 'RA',index = 0)
arches.add_column(DEC,name = 'DEC',index = 1)

e_pm= 0.3

arches = filter_hosek_data(arches,
                      max_e_pos = None,
                      max_e_pm = e_pm,
                       min_mag = 20,
                       max_mag = 10,
                      max_Pclust = 0.00001,
                      # center = 'yes'
                      ) 
# pm_l = (arches['e_pmRA']<e_pm) &(arches['e_pmDE']<e_pm)
# arches = arches[pm_l]
# fig,ax = plt.subplots(1,1)
# ax.hist(arches['t0'],bins ='auto')

gns1_pm = filter_gns_data(gns1_pm, 
                          # center = 'yes',
                           # min_mag = 18,
                           # max_mag = 10
                          )

arches['Jt0'] = [ 'J%s'%(arches['t0'][j]) for j in range(len((arches['t0'])))]

# We add galactic coordinates propermotions to Hosek distributions
arches_gal = SkyCoord(ra=arches['RA'], dec=arches['DEC'],
                    pm_ra_cosdec =arches['pmRA'], pm_dec = arches['pmDE'],
                    unit = 'degree',frame = 'icrs', obstime=arches['Jt0'].value).galactic   

arches['l'] = arches_gal.l
arches['b'] = arches_gal.b
arches['pml'] = arches_gal.pm_l_cosb
arches['pmb'] = arches_gal.pm_b


gns1_coor = SkyCoord(ra = gns1_pm['ra1'],dec = gns1_pm['Dec1'],unit = 'degree',
                     frame = 'fk5', obstime = 'J2015.4301')
# %
fig, ax = plt.subplots(1,1)
ax.scatter(arches['RA'], arches['DEC'])
ax.scatter(gns1_pm['ra1'], gns1_pm['Dec1'], alpha = 0.1)

arches_coord = SkyCoord(ra = arches['RA'], dec = arches['DEC'], unit = 'degree',
                        frame = 'icrs', obstime=arches['Jt0'].value)
max_sep = 0.05*u.arcsec

idx,d2d,d3d = arches_coord.match_to_catalog_sky(gns1_coor,nthneighbor=1)# ,nthneighbor=1 is for 1-to-1 matchsep_constraint = d2d < max_sep
sep_constraint = d2d < max_sep
gns_hos = gns1_pm[idx[sep_constraint]]
arches_match = arches[sep_constraint]

# %
fig, ax = plt.subplots(1,1)
ax.scatter(gns_hos['ra1'],gns_hos['Dec1'])
ax.scatter(arches_match['RA'],arches_match['DEC'],s =1)
# %
diff_pmx = gns_hos['pmx'] + arches_match['pml']
diff_pmy = gns_hos['pmy'] - arches_match['pmb']
# diff_pmx = gns_hos['pmx'] + arches_match['pmRA']
# diff_pmy = gns_hos['pmy']  - arches_match['pmDE']

mask_pmx, l_lim,h_lim = sigma_clip(diff_pmx, sigma=3, masked = True, return_bounds= True)
mask_pmy, l_lim,h_lim = sigma_clip(diff_pmy, sigma=3, masked = True, return_bounds= True)
diff_pmx_clip = diff_pmx[~mask_pmx.mask]
diff_pmy_clip = diff_pmy[~mask_pmy.mask]

# %
fig, (ax,ax1) = plt.subplots(1,2)
ax.set_title('Matching = %s'%(len(gns_hos['pmx'])))
ax1.set_title('HOSEK comparison')
ax.hist(diff_pmx_clip, histtype = 'step', label = '$\overline{\Delta x}$ =%.2f\n$\sigma$ = %.2f '%(np.nanmean(diff_pmx_clip),np.nanstd(diff_pmx_clip)))
ax.hist(diff_pmx, histtype = 'step', color ='k', alpha = 0.3)

ax1.hist(diff_pmy_clip, histtype = 'step',color = 'orange', label = '$\overline{\Delta x}$ =%.2f\n$\sigma$ = %.2f '%(np.nanmean(diff_pmy_clip),np.nanstd(diff_pmy_clip)))
ax1.hist(diff_pmy, histtype = 'step',color = 'k', alpha = 0.3, ls = 'dashed')
ax.legend()
ax1.legend()
ax.set_xlabel('$\Delta \mu_{\parallel}$')
ax1.set_xlabel('$\Delta \mu_{\perp}$')
















