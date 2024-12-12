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
# field_one = field_one.astype(int)
# chip_one = chip_one.astype(int)
# field_two = field_two.astype(int)
# chip_two = chip_two.astype(int)

field_one = 7
chip_one = 4
field_two = 7
chip_two = 1
max_sig = 0.5#TODO

max_deg = 2
d_m_pm = 2
gns1_pm = Table.read(pm_folder + f'pm_ep1_f{field_one}c{chip_one}_ep2_f{field_two}c{chip_two}deg{max_deg}_dmax{d_m_pm}_sxy%.1f.txt'%(max_sig), format = 'ascii')


# Comaparison with Hosek
from astropy.io import ascii
catal='/Users/amartinez/Desktop/PhD/Arches_and_Quintuplet_Hosek/'

choosen_cluster = 'Arches'
# choosen_cluster = 'Quintuplet'
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














