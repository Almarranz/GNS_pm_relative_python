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
max_sig = 0.5#TODO

max_deg = 3
d_m_pm = 2
gns1_pm = Table.read(pm_folder + f'pm_ep1_f{field_one}c{chip_one}_ep2_f{field_two}c{chip_two}deg{max_deg}_dmax{d_m_pm}_sxy%.1f.txt'%(max_sig), format = 'ascii')

# Ks_lim = [12,14.5]
Ks_lim = [0,99]
Ks_mask = (gns1_pm['Ks1'] > Ks_lim[0]) & (gns1_pm['Ks1'] < Ks_lim[1])
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
# %%
fig, (ax,ax1) = plt.subplots(1,2)
ax.scatter(gaia_['phot_g_mean_mag'],gaia_['pmra_error'], label = 'pm_ra')
ax.scatter(gaia_['phot_g_mean_mag'],gaia_['pmdec_error'],label = 'pm_Dec')
ax1.scatter(gaia_['phot_g_mean_mag'],gaia_['ra_error'], label = 'RA')
ax1.scatter(gaia_['phot_g_mean_mag'],gaia_['dec_error'],label = 'Dec')
ax.set_xlabel('G')
ax1.set_xlabel('G')
ax.set_ylabel('$\delta$ pm')
ax1.set_ylabel('$\delta position$')
ax.legend()
ax1.legend()
# %%

e_pm = 0.3
gaia_good = filter_gaia_data(
    gaia_table=gaia_,
    astrometric_params_solved=31,
    duplicated_source= False,
    parallax_over_error_min=-10,
    astrometric_excess_noise_sig_max=2,
    phot_g_mean_mag_min= 17,
    phot_g_mean_mag_max= 13 ,
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

center = (gns1_pm['H1'] - gns1_pm['Ks1'] < 1.3)
gns1_pm = gns1_pm[center]


gns1_coor = SkyCoord(ra = gns1_pm['ra1'],dec = gns1_pm['Dec1'],unit = 'degree',
                     frame = 'fk5', obstime = 'J2015.4301')
gaia_coord = SkyCoord(ra=gaia_good['ra'], dec=gaia_good['dec'],unit = 'degree',
                     frame = 'icrs', obstime = 'J2016.0')

max_sep = 0.05*u.arcsec



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

mask_pmx, l_lim,h_lim = sigma_clip(diff_pmx, sigma=3, masked = True, return_bounds= True)
mask_pmy, l_lim,h_lim = sigma_clip(diff_pmy, sigma=3, masked = True, return_bounds= True)

mask_pmxy = np.logical_and(np.logical_not(mask_pmx.mask), np.logical_not(mask_pmy.mask))

diff_pmx_clip = diff_pmx[mask_pmxy]
diff_pmy_clip = diff_pmy[mask_pmxy]

print(np.mean(diff_pmx),np.std(diff_pmx))
print(np.mean(diff_pmy),np.std(diff_pmy))
    # max_sep = 0.08*u.arcsec
    # idx,d2d,d3d = gaia_coord.match_to_catalog_sky(hose_coord,nthneighbor=1)# ,nthneighbor=1 is for 1-to-1 match
    # sep_constraint = d2d < max_sep
    # hose_match = arches[idx[sep_constraint]]
    # gaia_match_hos = gaia_good[sep_constraint]
# %%
fig, (ax,ax1) = plt.subplots(1,2)
ax.set_title('Matching = %s'%(len(gns_match['pmx'])))
ax1.set_title('GAIA comparison')
ax.hist(diff_pmx_clip, histtype = 'step', label = '$\overline{\Delta x}$ =%.2f\n$\sigma$ = %.2f '%(np.nanmean(diff_pmx_clip),np.nanstd(diff_pmx_clip)))
ax.hist(diff_pmx, histtype = 'step', color ='k', alpha = 0.3)

ax1.hist(diff_pmy_clip, histtype = 'step',color = 'orange', label = '$\overline{\Delta y}$ =%.2f\n$\sigma$ = %.2f '%(np.nanmean(diff_pmy_clip),np.nanstd(diff_pmy_clip)))
ax1.hist(diff_pmy, histtype = 'step',color = 'k', alpha = 0.3, ls = 'dashed')

ax.set_xlabel('$\Delta \mu_{\parallel}$')
ax1.set_xlabel('$\Delta \mu_{\perp}$')
ax.legend()
ax1.legend()

