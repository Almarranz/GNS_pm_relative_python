# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Created on Fri Nov  8 12:42:20 2024

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
import astroalign as aa
import os
import glob
from astropy.modeling.models import Polynomial2D
from astropy.modeling.fitting import LinearLSQFitter
from astropy.modeling import models, fitting
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

field_one, chip_one, field_two, chip_two,t1,t2,max_sig = np.loadtxt('/Users/amartinez/Desktop/PhD/HAWK/GNS_1relative_python/lists/fields_and_chips.txt', 
                                                       unpack=True)
# field_one = field_one.astype(int)
# chip_one = chip_one.astype(int)
# field_two = field_two.astype(int)
# chip_two = chip_two.astype(int)

field_one = 100
chip_one = 2
field_two = 20
chip_two = 2
max_sig = 0.5#TODO




GNS_1relative = '/Users/amartinez/Desktop/PhD/HAWK/GNS_1relative_python/lists/%s/chip%s/'%(field_one, chip_one)
GNS_2relative = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2relative_python/lists/%s/chip%s/'%(field_two, chip_two)
pm_folder = '/Users/amartinez/Desktop/PhD/HAWK/pm_gns1_gns2_relative_python/'



pruebas2 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2relative/pruebas/'
pruebas1 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_1relative/pruebas/'
color = pd.read_csv('/Users/amartinez/Desktop/PhD/python/colors_html.csv')
morralla = '/Users/amartinez/Desktop/morralla/'
strin= color.values.tolist()
indices = np.arange(0,len(strin),1)

# center_only = 'yes'#TODO yes, eliminates foregroud, no, keep them
center_only = 'no'
pix_scale = 0.1064*0.53

max_sep = 0.05 * u.arcsec
sig_cl = 1#!!!
deg = 1#!!!
max_deg = 3
d_m = 0.5#!!! max distance for the fine alignment betwenn GNS1 and 2
d_m_pm = 2#!!! max distance for the proper motions
align_by = 'Polywarp'#!!!
# align_by = '2DPoly'#!!!
# f_mode = 'W'
f_mode = 'WnC'
# f_mode = 'NW'
# f_mode = 'NWnC'

# %%
# 'ra1 0, dec1 1, x1 2, y1 3, f1 4, H1 5, dx1 6, dy1 7, df1 8, dH1 9 ,Ks 10, dKs 11 ID1 12'
# gns1 = np.loadtxt(GNS_1relative + 'stars_calibrated_HK_chip%s_on_gns2_f%sc%s_sxy%s.txt'%(chip_one,field_two,chip_two,max_sig))
gns1 = Table.read(GNS_1relative + 'stars_calibrated_HK_chip%s_on_gns2_f%sc%s_sxy%s.txt'%(chip_one,field_two,chip_two,max_sig), format = 'ascii')
# 'ra1 0, dec1 1, x1 2, y1 3, f1 4, H1 5, dx1 6, dy1 7, df1 8, dH1 9 ID2'
# gns2 = np.loadtxt(GNS_2relative +'stars_calibrated_H_chip%s_on_gns1_f%sc%s_sxy%s.txt'%(chip_two,field_one,chip_one,max_sig))
gns2 = Table.read(GNS_2relative +'stars_calibrated_H_chip%s_on_gns1_f%sc%s_sxy%s.txt'%(chip_two,field_one,chip_one,max_sig), format = 'ascii')


if center_only == 'yes':
    center = np.where(gns1['H1'] - gns1['Ks1'] > 1.3)
elif center_only == 'no':
    center = np.where(gns1['H1'] - gns1['Ks1']  > -1)

# gns1_center = copy.deepcopy(gns1)
gns1 = gns1[center] 

gns1_ra_dec = SkyCoord(ra = gns1['ra1'], dec = gns1['Dec1'], unit ='deg', frame = 'fk5',equinox ='J2000',obstime='J2015.43')
gns2_ra_dec = SkyCoord(ra= gns2['ra2']*u.degree, dec=gns2['Dec2']*u.degree, frame = 'fk5', equinox = 'J2000',obstime='J2022.4')
# 
#I cosider a math if the stars are less than 'max_sep' arcsec away 
# This is for cutting the the overlapping areas of both lists. (Makes astroaling work faster)

idx,d2d,d3d = gns1_ra_dec.match_to_catalog_sky(gns2_ra_dec)# ,nthneighbor=1 is for 1-to-1 match
sep_constraint = d2d < max_sep
gns1_match = gns1[sep_constraint]
gns2_match = gns2[idx[sep_constraint]]

diff_H = gns1_match['H1']-gns2_match['H2']
mask_H, l_lim,h_lim = sigma_clip(diff_H, sigma=sig_cl, masked = True, return_bounds= True)

fig, ax = plt.subplots(1,1)
ax.set_title('Max_dis = %.2f". Matchs = %s'%(max_sep.value, len(gns2_match)))
ax.hist(diff_H, bins = 'auto')
ax.axvline(0.1, color = 'orange', ls = 'dashed')
ax.axvline(-0.1, color = 'orange', ls = 'dashed',label = '$\pm$0.1H')
ax.axvline(l_lim, color = 'red', ls = 'dashed')
ax.axvline(h_lim, color = 'red', ls = 'dashed',label = '$\pm %s\sigma$'%(sig_cl))
ax.legend()
ax.set_xlabel('diff [H]')
loop = 0
comom_ls = []
dic_xy = {}
dic_Kx ={}
dic_xy_final = {}

# xy_1c = np.array((gns1_match['x1'], gns1_match['y1'])).T
# xy_2c = np.array((gns2_match['x2'], gns2_match['y2'])).T
xy_1c = np.array((gns1_match['x1'][~mask_H.mask], gns1_match['y1'][~mask_H.mask])).T
xy_2c = np.array((gns2_match['x2'][~mask_H.mask], gns2_match['y2'][~mask_H.mask])).T

# p = ski.transform.estimate_transform('affine',
#                                              xy_1c, 
#                                              xy_2c)
p = ski.transform.estimate_transform('similarity',
                                              xy_1c, 
                                              xy_2c)
# p = ski.transform.estimate_transform('polynomial',
#                                              xy_1c, 
#                                              xy_2c, order = deg)

gns1_xy = np.array((gns1['x1'],gns1['y1'])).T
gns1_xyt = p(gns1_xy)

s_ls = compare_lists(gns1_xyt, np.array([gns2['x2'],gns2['y2']]).T, d_m)




# %%
fig, ax = plt.subplots(1,1)
ax.scatter(gns2['x2'],gns2['y2'],s=10, color = 'k', alpha = 0.1,label = 'GNS2')
ax.scatter(xy_1c[:,0],xy_1c[:,1], marker = 'x',label = 'GNS1 matched')
# ax.scatter(gns1_xyt[:,0],gns1_xyt[:,1],marker = 'x',color = 'r',label = 'GNS1 transformed')
ax.scatter(xy_2c[:,0],xy_2c[:,1],s=10, label = 'GNS2 matched')

ax.legend(fontsize = 9) 
# ax.set_xlim(2000,2300)
# %%
# sys.exit()
gns1['x1'] = gns1_xyt[:,0]
gns1['y1'] = gns1_xyt[:,1]

gns1.write(pruebas1 + 'gns1_trans.txt', format = 'ascii', overwrite = True)
gns2.write(pruebas2 + 'gns2_trans.txt', format = 'ascii',overwrite = True)

while deg < max_deg:
    loop += 1 
    l1_xy = np.array([gns1['x1'],gns1['y1']]).T
    l2_xy = np.array([gns2['x2'],gns2['y2']]).T
    comp = compare_lists(l1_xy,l2_xy,d_m)
    if len(comom_ls) >1:
        if comom_ls[-1] <= comom_ls[-2]:
            try:
                gns1['x1'] = dic_xy[f'trans_{loop-2}'][:,0]
                gns1['y1'] = dic_xy[f'trans_{loop-2}'][:,1]
                dic_xy_final['xy_deg%s'%(deg)] = np.array([dic_xy[f'trans_{loop-2}'][:,0],dic_xy[f'trans_{loop-2}'][:,1]]).T            
                comom_ls =[]
                dic_xy = {}
                dic_Kx = {}
                deg += 1
                print(f'Number of common star in loop {loop-1} lower tha in loop {loop-2}.\nJupping to degree {deg} ')
                loop = -1
                continue
            except:
                gns1['x1'] = dic_xy_final[f'xy_deg{deg-1}'][:,0]
                gns1['y1'] = dic_xy_final[f'xy_deg{deg-1}'][:,1]
                print(f'Number of common star with polynomial degere {deg} decreases after a single iteration.\nUsing the last iteration of degree {deg -1} ')
                deg = deg
                break
            
    comom_ls.append(len(comp))
    print(f'Common in loop {loop}, degree {deg} = %s'%(len(comp['ind_1'])))
    # if loop == 1:
    #     with open(pruebas + 'sig_and_com.txt', 'a') as file:
    #         file.write('%.1f %.0f\n'%(sig_cl, len(comp['ind_1'])))

    l1_com = gns1[comp['ind_1']]
    l2_com = gns2[comp['ind_2']]
    
    
    diff_mag = l1_com['H1'] - l2_com['H2'] 
    # diff_mag1 = l1_com['IB230_diff'] - l2_com['IB230_diff'] 
    diff_x =  l2_com['x2'] - l1_com['x1'] 
    diff_y =  l2_com['y2'] - l1_com['y1'] 
    diff_xy = (diff_x**2 + diff_y**2)**0.5
    mask_m, l_lim,h_lim = sigma_clip(diff_mag, sigma=sig_cl, masked = True, return_bounds= True)
    
    
    
    # l2_clip = l2_com
    # l1_clip = l1_com
    
    l1_clip = l1_com[~mask_m.mask]
    l2_clip = l2_com[~mask_m.mask]
    
    # fig, (ax,ax1) = plt.subplots(1,2)
    # fig.suptitle(f'Degree = {deg}. Loop = {loop}')
    # ax.set_xlabel('$\Delta$ H')
    # ax.hist(diff_mag, label = 'matches = %s\ndist = %.2f arcsec'%(len(comp['ind_1']), d_m*pix_scale))
    # ax.axvline(l_lim, color = 'red', ls = 'dashed', label = '$\pm$%s$\sigma$'%(sig_cl))
    # ax.axvline(h_lim, color = 'red', ls = 'dashed')
    # ax.legend()
    
    # ax1.hist(diff_x, label = '$\overline{\Delta x} = %.2f\pm%.2f$'%(np.mean(diff_x),np.std(diff_x)), histtype = 'step')
    # ax1.hist(diff_y, label = '$\overline{\Delta y} = %.2f\pm%.2f$'%(np.mean(diff_y),np.std(diff_y)), histtype = 'step')

    # ax1.set_xlabel('$\Delta$ pixel')
    # ax1.legend()
    
    
    
    
    xy_1c = np.array([l1_clip['x1'],l1_clip['y1']]).T
    xy_2c = np.array([l2_clip['x2'],l2_clip['y2']]).T
    
    if align_by == 'Polywarp':
        Kx,Ky=pw.polywarp(xy_2c[:,0],xy_2c[:,1],xy_1c[:,0],xy_1c[:,1],degree=deg)
        
        xi=np.zeros(len(gns1))
        yi=np.zeros(len(gns1))
        
        for k in range(deg+1):
                    for m in range(deg+1):
                        xi=xi+Kx[k,m]*gns1['x1']**k*gns1['y1']**m
                        yi=yi+Ky[k,m]*gns1['x1']**k*gns1['y1']**m
    elif align_by == '2DPoly':
        model_x = Polynomial2D(degree=deg)
        model_y = Polynomial2D(degree=deg)
        
        # Linear least-squares fitter
        fitter = LinearLSQFitter()
        # fitter = fitting.LMLSQFitter()
        
        if f_mode == 'W':
            #==============
            # Fit weighted
            #==============
            fit_xw = fitter(model_x, xy_1c[:,0],xy_1c[:,1], xy_2c[:,0], weights= 1/np.sqrt(l1_clip['dx1']**2 + l1_clip['dy1']**2))  # Fit x-coordinates
            fit_yw = fitter(model_y, xy_1c[:,0],xy_1c[:,1],  xy_2c[:,1],weights= 1/np.sqrt(l1_clip['dx1']**2 + l1_clip['dy1']**2)) 
        
        elif f_mode == 'NW':
            #==============
            # Fit not weighted
            #==============
            fit_xw = fitter(model_x, xy_1c[:,0],xy_1c[:,1], xy_2c[:,0])  # Fit x-coordinates
            fit_yw = fitter(model_y, xy_1c[:,0],xy_1c[:,1],  xy_2c[:,1]) 
        
        elif f_mode == 'WnC':
            #==============
            #Fit weighted and clipped
            #==============
            or_fit = fitting.FittingWithOutlierRemoval(fitter, sigma_clip, niter=3, sigma=3.0)
            
            fit_xw, fm_x = or_fit(model_x, xy_1c[:,0],xy_1c[:,1], xy_2c[:,0], weights= 1/np.sqrt(l1_clip['dx1']**2 + l1_clip['dy1']**2))  # Fit x-coordinates
            fit_yw, fm_y = or_fit(model_y, xy_1c[:,0],xy_1c[:,1],  xy_2c[:,1],weights= 1/np.sqrt(l1_clip['dx1']**2 + l1_clip['dy1']**2)) 
            
            fm_xy = np.logical_not(np.logical_and(np.logical_not(fm_x),np.logical_not(fm_y)))
            
            fig, ax = plt.subplots(1,1)
            plt_per = 20
            ax.set_title(f'Degree = {deg}')
            ax.scatter(xy_1c[:,0][::plt_per],xy_1c[:,1][::plt_per], label = 'Used (%s%%)'%(plt_per))
            ax.scatter(xy_1c[:,0][fm_xy],xy_1c[:,1][fm_xy],s =200, label = 'Clipped')
            ax.legend()
        elif f_mode == 'NWnC':
            or_fit = fitting.FittingWithOutlierRemoval(fitter, sigma_clip, niter=3, sigma=3.0)
            
            fit_xw, fm_x = or_fit(model_x, xy_1c[:,0],xy_1c[:,1], xy_2c[:,0])  # Fit x-coordinates
            fit_yw, fm_y = or_fit(model_y, xy_1c[:,0],xy_1c[:,1],  xy_2c[:,1]) 
            
            fm_xy = np.logical_not(np.logical_and(np.logical_not(fm_x),np.logical_not(fm_y)))
            
            fig, ax = plt.subplots(1,1)
            plt_per = 20
            ax.set_title(f'Degree = {deg}')
            ax.scatter(xy_1c[:,0][::plt_per],xy_1c[:,1][::plt_per], label = 'Used (%s%%)'%(plt_per))
            ax.scatter(xy_1c[:,0][fm_xy],xy_1c[:,1][fm_xy],s =200, label = 'Clipped')
            ax.legend()
            
        # sys.exit(308)
        
        xi = fit_xw(gns1['x1'], gns1['y1'])
        yi = fit_yw(gns1['x1'], gns1['y1'])# Fit y-coordinates
        
    dic_xy[f'trans_{loop+1}'] = np.array([xi,yi]).T
    
    # print(Kx[0][0])
    gns1['x1'] = xi
    gns1['y1'] = yi
# %%
l1_xy = np.array([gns1['x1'],gns1['y1']]).T
l2_xy = np.array([gns2['x2'],gns2['y2']]).T
comp = compare_lists(l1_xy,l2_xy,d_m)
print(len(comp))
# %% 
# Proper motion calculations

l1 = np.array((gns1['x1'],gns1['y1'])).T
l2 = np.array((gns2['x2'],gns2['y2'])).T


l_12 = compare_lists(l1, l2,d_m_pm)

gns1_pm = gns1[l_12['ind_1']]
gns2_pm = gns2[l_12['ind_2']]
delta_t = t2-t1
pm_x = (gns2['x2'][l_12['ind_2']] - gns1['x1'][l_12['ind_1']])*pix_scale*1000/delta_t
pm_y = (gns2['y2'][l_12['ind_2']] - gns1['y1'][l_12['ind_1']])*pix_scale*1000/delta_t
gns1_pm['pmx'] = pm_x
gns1_pm['pmy'] = pm_y

# %%
bins = 20
fig, (ax,ax1)= plt.subplots(1,2)
fig.suptitle(f'GNS1[f{field_one},c{chip_one}] GNS2[f{field_two},c{chip_two}]', ha = 'right')
ax.hist(pm_x, bins = bins, label ='$\overline{\mu_x}$ = %.2f\n$\sigma$ = %.2f'%(np.mean(pm_x),np.std(pm_x)))
ax1.hist(pm_y, bins = bins,label ='$\overline{\mu_y}$ = %.2f\n$\sigma$ = %.2f'%(np.mean(pm_y),np.std(pm_y)))
ax.legend()
ax1.legend()
ax.set_xlabel('$\mu_x$[mas/yr]')
ax1.set_xlabel('$\mu_y$[mas/yr]')

# sys.exit(347)
files_to_remove = glob.glob(os.path.join(pm_folder, f'pm_ep1_f{field_one:.0f}c{chip_one}_ep2_f{field_two}c{chip_two}**.txt'))

# Remove the files
for file in files_to_remove:
    try:
        os.remove(file)
        print(f"Removed: {file}")
    except Exception as e:
        print(f"Error removing {file}: {e}")
gns1_pm.write(pm_folder + f'pm_ep1_f{field_one:.0f}c{chip_one}_ep2_f{field_two}c{chip_two}deg{deg-1}_dmax{d_m_pm}_sxy%.1f.txt'%(max_sig), format = 'ascii', overwrite = True)
print(30*'_' + f'\npm_ep1_f{field_one:.0f}c{chip_one}_ep2_f{field_two}c{chip_two}deg{deg-1}_dmax{d_m_pm}_sxy%.1f.txt\n'%(max_sig)+ 30*'_')
sys.exit(388)


# %%
fig, ax  = plt.subplots(1,1)
ax.scatter(gns1['H1']-gns1['Ks1'],gns1['Ks1'], s = 1)
ax.invert_yaxis()
ax.axvline(1.3)
# %%
# =============================================================================
# knn = 25
# sim_by= 'kernnel'
# lim = 'minimun'
# modes = ['pm_xy_color','pm_xy','pm_color','pm']
# # mode = 'pm_xy'
# # clst = cluster_finder.finder(pm_x, pm_y, l2_c['x'], l2_c['y'], l2_c['RA'], l2_c['Dec'], modes[1], 
# #                              l2_c['H_diff'],l2_c['IB236_diff'],
# #                              knn,sim_by,lim)
# clst = cluster_finder.finder(pm_x, pm_y, l2_c['x2'], l2_c['y2'], l2_c['ra2'], l2_c['Dec2'], modes[0], 
#                               l1_c['H1'],l1_c['Ks1'],
#                               knn,sim_by,lim)
# =============================================================================
# %%
# Checking th alignment with astroalign

# gns1 = Table.read(GNS_1relative + 'stars_calibrated_HK_chip%s_on_gns2_f%sc%s_sxy%s.txt'%(chip_one,field_two,chip_two,max_sig), format = 'ascii')
# gns2 = Table.read(GNS_2relative +'stars_calibrated_H_chip%s_on_gns1_f%sc%s_sxy%s.txt'%(chip_two,field_one,chip_one,max_sig), format = 'ascii')


# if center_only == 'yes':
#     center = np.where(gns1['H1'] - gns1['Ks1'] > 1.3)
# elif center_only == 'no':
#     center = np.where(gns1['H1'] - gns1['Ks1']  > -1)

# gns1_ra_dec = SkyCoord(ra = gns1['ra1'], dec = gns1['Dec1'], unit ='deg', frame = 'fk5',equinox ='J2000',obstime='J2015.43')
# gns2_ra_dec = SkyCoord(ra= gns2['ra2']*u.degree, dec=gns2['Dec2']*u.degree, frame = 'fk5', equinox = 'J2000',obstime='J2022.4')
# # 
# #I cosider a math if the stars are less than 'max_sep' arcsec away 
# # This is for cutting the the overlapping areas of both lists. (Makes astroaling work faster)

# idx,d2d,d3d = gns1_ra_dec.match_to_catalog_sky(gns2_ra_dec)# ,nthneighbor=1 is for 1-to-1 match
# sep_constraint = d2d < max_sep
# gns1_match = gns1[sep_constraint]
# gns2_match = gns2[idx[sep_constraint]]



# xy1 = np.array([gns1_match['x1'],gns1_match['y1']]).T
# xy2 = np.array([gns2_match['x2'],gns2_match['y2']]).T
# m,(c1,c2)= aa.find_transform(xy1,xy2,max_control_points=1000)

# gns2_xy = np.array((gns2['x2'],gns2['y2'])).T
# gns1_xy = np.array((gns1['x1'],gns1['y1'])).T
# gns1_xyt = m(gns1_xy)

# aa_list = compare_lists(gns1_xyt, gns2_xy, d_m)







# %%
# fig, (ax,ax1) = plt.subplots(1,2)
# ax.scatter(gns1['Ks1'], gns1['dx1'],s=1)
# ax1.scatter(gns1['Ks1'], gns1['dy1'],s=1)
# ax.set_ylim(0,0.5)
# ax.set_xlim(10,17)
# ax1.set_ylim(0,0.5)
# ax.grid()
