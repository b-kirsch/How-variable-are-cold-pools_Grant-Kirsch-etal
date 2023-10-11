# -*- coding: utf-8 -*-
"""
@author: Bastian Kirsch (bastian.kirsch@uni-hamburg.de)

Script to perform analyses of FESSTVaL observations and 
ICON-LES simulations for cold pool variogram paper

Dependences on non-standard software:  
- fesstval_routines.py 
- empirical_variogram.py
- variogram_paper_routines.py

Required meta data files:
- cold_pools_fesstval.txt


Last updated: 6 October 2023
"""

import numpy as np
import pandas as pd
import datetime as dt

import fesstval_routines as fst
import empirical_variogram as eva
import variogram_paper_routines as vpr

t_run = dt.datetime.now()

#----------------------------------------------------------------------------
# Paths and meta data files
maindir        = '/Users/bkirs/Desktop/'
cp_file2021    = maindir+'FESSTVaL/cold_pools_fesstval.txt'

datadir_apo_l2 = maindir+'APOLLO_data/level2/'
datadir_wxt_l2 = maindir+'WXT_data/level2/' 
datadir_apo_l3 = maindir+'APOLLO_data/level3/'
datadir_wxt_l3 = maindir+'WXT_data/level3/' 
datadir_icon   = maindir+'ICON_data/'
datadir_io     = maindir+'Variogram_data/'

plotdir        = maindir+'Cold-Pools/Plots/Paper_Variogram/'

#----------------------------------------------------------------------------  
# Basic settings

# Observations
calc_obs   = True
start_time = dt.datetime(2021,6,29,0,0)  #dt.datetime(2021,6,29,12,30)  
end_time   = dt.datetime(2021,6,29,23,59)  #dt.datetime(2021,6,29,16,30)
freq       = 15 # (min)

# ICON data
calc_icon     = True
icon_date     = dt.date(2021,6,29)

# Obs + ICON data
read_grid_data   = True # Obs level3 / ICON re-gridded
write_nc         = False  # Variogram data

# Plots
plot_fig1           = False
plot_figs1_hist     = False
plot_figs2_time     = False


if plot_fig1: read_grid_data = True
#if plot_figs1_hist: freq = 1
#----------------------------------------------------------------------------   
# Advanced settings
# Smoothing of APOLLO data
time_smooth = 30 # 5 (s) Half-window length for smooting of APOLLO T data

# Time resolution of data 
time_res = 10 # (s)

# Distance bins for variogram
dist_min,dist_max,dist_res = 0,15000,500 # (m)

# Duration of cold pool events in observations
dur_event_obs  = 2  # (hours)
dur_event_icon = 3

# Detection limit for CP coverage
# dT_lim_cov = -1  # (K)

# Reference periods for characterization of Jogi (2021-06-29)
cp_onset_obs = dt.datetime(2021,6,29,13,30) # UTC
dt_ref_obs   = slice(cp_onset_obs-dt.timedelta(hours=1),cp_onset_obs-dt.timedelta(minutes=1))
dt_1h_obs    = slice(cp_onset_obs,cp_onset_obs+dt.timedelta(hours=1))
dt_cp_obs    = slice(cp_onset_obs,cp_onset_obs+dt.timedelta(hours=4))

cp_onset_icon = dt.datetime(2021,6,29,12,0) # UTC
dt_ref_icon   = slice(cp_onset_icon-dt.timedelta(hours=1),cp_onset_icon-dt.timedelta(minutes=1))
dt_1h_icon    = slice(cp_onset_icon,cp_onset_icon+dt.timedelta(hours=1))
dt_cp_icon    = slice(cp_onset_icon,cp_onset_icon+dt.timedelta(hours=4))


#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Meta data station network
print('Setting up meta data')
kwargs_l2 = {'datadir_apollo':datadir_apo_l2,'datadir_wxt':datadir_wxt_l2}
kwargs_l3 = {'datadir_apollo':datadir_apo_l3,'datadir_wxt':datadir_wxt_l3}
meta_data_dict = {}
nstat = {}

for yy in [2020,2021]:
    date_meta = (yy,6,29)
    TT,meta_apo = fst.read_fesstval_level2('a','TT',*date_meta,**kwargs_l2,
                                            return_meta=True)
    TT,meta_wxt = fst.read_fesstval_level2('w','TT',*date_meta,**kwargs_l2,
                                            return_meta=True)
    meta_data   = meta_apo.append(meta_wxt).drop(['ALT','HT','NAME','LCZ'],1) 
    
    if yy == 2020:
        lon_ref,lat_ref = 9.99302,53.55074   
    else:    
        lon_ref,lat_ref = meta_data.loc['001Sa','LON'],meta_data.loc['001Sa','LAT']

    x_data,y_data = fst.lonlat_to_xy(meta_data['LON'].values,meta_data['LAT'].values,
                                     lon_ref=lon_ref,lat_ref=lat_ref)
    meta_data['X'] = x_data
    meta_data['Y'] = y_data 
    
    meta_data_dict[yy] = meta_data
    nstat[yy] = meta_data.shape[0]

# Time
days  = pd.date_range(start=start_time,end=end_time,freq='d')
times_obs = pd.date_range(start=start_time,end=end_time,freq=str(freq)+'min')   

# Variograms pairs  
dist_bins = np.arange(dist_min,dist_max+dist_res,dist_res)
dist_mids = dist_bins[:-1]+dist_res/2.

pairs_dict = {}
for yy in [2020,2021]: 
    pairs_dict[yy] = eva.pair_stations(meta_data_dict[yy],bins=dist_bins,
                                       bin_labels=dist_mids)
   
# Cold Pool flag (Obs)    
cp_list = fst.fesstval_cold_pools(cpfile=cp_file2021)
# Change start time of Jogi by 1 h to be consistent with case study analysis
cp_list.rename(index={dt.datetime(2021,6,29,12,30):cp_onset_obs},inplace=True)

cp_flag_obs = pd.Series(index=times_obs,dtype=float)
cp_flag_obs[:] = 0

for t in cp_list.index:    
    cp_time = (times_obs >= t) & (times_obs <= t+dt.timedelta(hours=dur_event_obs))
    cp_flag_obs.loc[cp_time] = 2  
    
    cp_start = times_obs == t
    cp_flag_obs.loc[cp_start] = 1
    

#----------------------------------------------------------------------------
# Reading and analyzing level2 data
variogram_obs_times  = pd.DataFrame(index=times_obs,columns=dist_mids)

valid_pairs_obs = pd.DataFrame(index=pairs_dict[2021].index,columns=times_obs)
    
if calc_obs:
    print('')
    print('Analyzing FESSTVaL observation data')
    for i,d in enumerate(days):
        print(d.strftime('%Y-%m-%d'))
        pairs_year = pairs_dict[d.year]
        
        # Reading level2 data
        date_args = (d.year,d.month,d.day)
        times_day = times_obs[times_obs.date == d.date()]
        TT_apo = fst.read_fesstval_level2('a','TT',*date_args,**kwargs_l2)
        TT_wxt = fst.read_fesstval_level2('w','TT',*date_args,**kwargs_l2)

        # Smooting of raw data  
        TT_apo = TT_apo.rolling(2*time_smooth+1,center=True,\
                                min_periods=time_smooth+1,axis=0).mean()
            
        TT_wxt = TT_wxt.rolling(int(2*(time_smooth/10)+1),center=True,\
                                min_periods=int((time_smooth/10)+1),axis=0).mean()    
        
        TT_obs = TT_apo.loc[start_time:end_time].loc[::time_res].\
                  join(TT_wxt.loc[start_time:end_time]) 
           
        # Outlier station for pre-Jogi period (impacted by another CP)        
        # TT_obs['078Ha'] = np.nan
                  
        # Variogram for single time steps
        for t in times_day:
            variogram_obs_times.loc[t]  = eva.calc_variogram(TT_obs.loc[t],pairs_year)
            valid_pairs_obs[t]          = eva.calc_variogram(TT_obs.loc[t],pairs_year,
                                                             return_valid_pairs=True) 
            
        if d.date() == dt.date(2021,6,29): 
            # Calculation of numbers for overview table (Teams)
            TT_ref_obs  = TT_obs.loc[dt_ref_obs].mean(axis=1).mean()
            dT_mean_obs = (TT_obs.loc[dt_1h_obs] - TT_ref_obs).mean(axis=1).mean()
            dT_min_obs  = (TT_obs.loc[dt_cp_obs] - TT_ref_obs).min().min()
            sT_mean_obs = (TT_obs.loc[dt_1h_obs] - TT_ref_obs).std(axis=1).mean()
            sT_max_obs  = (TT_obs.loc[dt_cp_obs] - TT_ref_obs).std(axis=1).max()
    
            
    if write_nc:
        print('Writing variograms to nc files')
        
        if start_time.date() == end_time.date():
            out_file = 'variogram_obs_'+start_time.strftime('%Y%m%d.nc')
        else:
            out_file = 'variogram_obs_'+start_time.strftime('%Y%m%d-')+end_time.strftime('%Y%m%d.nc')
        vpr.write_netcdf(out_file,times_obs,variogram_obs_times,dist_mids,cp_flag_obs)

        
#----------------------------------------------------------------------------
# Reading and analyzing ICON-LES data
times_icon = pd.date_range(start=icon_date,periods=1440/15,freq='15min')    
domains = range(1,5)

# Accessing specific domain: df.xs(domain,level=1)
index_icon_times   = pd.MultiIndex.from_product([times_icon,domains],
                                                names=['TIME','DOMAIN'])
variogram_icon_times  = pd.DataFrame(index=index_icon_times,columns=dist_mids)

cp_flag_icon = pd.Series(index=times_icon,dtype=float)
cp_flag_icon[:] = 0

cp_flag_icon.loc[cp_onset_icon:cp_onset_icon+dt.timedelta(hours=dur_event_icon)] = 2
cp_flag_icon.loc[cp_onset_icon] = 1
    

valid_pairs_icon = pd.DataFrame(index=pairs_dict[2021].index,columns=domains)
TT_icon_time     = pd.DataFrame(index=times_icon,columns=domains)   

if calc_icon:
    print('')
    print('Analyzing ICON-LES data')
    TT_ref_icon = {}
    
    for i,d in enumerate(days):  
        if d != icon_date: continue
        print(d.strftime('%Y-%m-%d'))
        pairs_year = pairs_dict[d.year]
        date_args = (d.year,d.month,d.day)
        times_day = times_icon[times_icon.date == d.date()]
        
        for dom in domains:
            TT_icon = vpr.read_icon_les(*date_args,dom)
            TT_icon_time[dom] = TT_icon.mean(axis=1)
            
            for t in times_day:
                if TT_icon.loc[t].notnull().sum() == 0: continue
                variogram_icon_times.loc[(t,dom)]  = eva.calc_variogram(TT_icon.loc[t],pairs_year)
     
            # Calculation of numbers for simulation overview table (Teams)
            maxis = 0 if dom <= 2 else 1
                
            TT_ref_icon[dom]  = TT_icon.loc[dt_ref_icon].mean(axis=1).mean()
            dT_mean_icon = (TT_icon.loc[dt_1h_icon] - TT_ref_icon[dom]).mean(axis=maxis).mean()
            dT_min_icon  = (TT_icon.loc[dt_cp_icon] - TT_ref_icon[dom]).min().min()
            sT_mean_icon = (TT_icon.loc[dt_1h_icon] - TT_ref_icon[dom]).std(axis=1).mean()
            sT_max_icon  = (TT_icon.loc[dt_cp_icon] - TT_ref_icon[dom]).std(axis=1).max()
            
            # print('DOM0'+str(dom)+': ','{:1.1f}'.format(TT_ref_icon[dom]),'{:1.1f}'.format(dT_mean_icon),
            #       '{:1.1f}'.format(dT_min_icon),'{:1.2f}'.format(sT_mean_icon),'{:1.2f}'.format(sT_max_icon))

            valid_pairs_icon[dom] = eva.calc_variogram(TT_icon.loc[dt.datetime(2021,6,29,12,0)],
                                                        pairs_year,return_valid_pairs=True)
            
    if write_nc:
        print('Writing variograms to nc files')
        
        dom_ext = {1:'625m',2:'312m',3:'156m',4:'75m'}
        
        for dom in [1,2,3,4]:
            out_file = 'variogram_icon_'+dom_ext[dom]+start_time.strftime('_%Y%m%d.nc')
            vpr.write_netcdf(out_file,times_icon,variogram_icon_times.xs(dom,level=1),
                             dist_mids,cp_flag_icon)
        
#----------------------------------------------------------------------------        
# Reading level3 data (1x1 km2 grid) and calculating CP coverage
if read_grid_data & (len(days) == 1):
    print('Reading level 3 data and re-gridded ICON data')
    date_args = (start_time.year,start_time.month,start_time.day)
    
    # Observations
    TT_obs_l3   = fst.read_fesstval_level3('TT',*date_args,**kwargs_l3) 
    meta_data_l3 = fst.read_fesstval_level3('TT',*date_meta,return_meta=True,**kwargs_l3)    
    
    # # CP Coverage
    # n_cov_max = meta_data_l3['mask_meshgrid'].sum()
    # ii_fin    = np.isfinite(TT_obs_l3).sum(axis=(0,1)) > 0
    
    # n_cov     = ((TT_obs_l3 - TT_ref_obs) <= dT_lim_cov).sum(axis=(0,1))
    # cp_coverage_obs = pd.Series(n_cov/n_cov_max,index=meta_data_l3['time'],dtype=float)
    # cp_coverage_obs[~ii_fin] = np.nan
     
        
    # ICON
    dom_comp  = 3
    TT_icon_grid_dom3 = vpr.read_icon_les(*date_args,dom_comp,type_str='grid')
    meta_icon_grid = vpr.read_icon_les(*date_args,dom_comp,type_str='grid',
                                       return_meta=True) 
    
    # # CP Coverage
    # n_cov_max = np.isfinite(TT_icon_grid_dom3[:,:,:]).sum(axis=(0,1)).max()
    # ii_fin    = np.isfinite(TT_icon_grid_dom3).sum(axis=(0,1)) > 0
    
    # n_cov     = ((TT_icon_grid_dom3 - TT_ref_icon[dom_comp]) <= dT_lim_cov).sum(axis=(0,1))
    # cp_coverage_icon = pd.Series(n_cov/n_cov_max,index=meta_icon_grid['time'],dtype=float)
    # cp_coverage_icon[~ii_fin] = np.nan


#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Plotting    
if plot_fig1 & read_grid_data:
    
    i_obs  = meta_data_l3['time'].get_loc(dt.datetime(2021,6,29,14,15))
    i_icon = meta_icon_grid['time'].get_loc(dt.datetime(2021,6,29,14,0))
    
    obs_data  = TT_obs_l3[:,:,i_obs] - TT_ref_obs
    icon_data = TT_icon_grid_dom3[:,:,i_icon] - TT_ref_icon[3]
    
    rams_data_ideal   = vpr.read_rams_ideal()
    rams_data_oceanic = vpr.read_rams_oceanic(100)
    rams_data_haboob  = vpr.read_rams_haboob('Day')
    rams_data_shear   = vpr.read_rams_shear('upshear')
    
    
    vpr.figure1(obs_data,icon_data,rams_data_ideal,
                rams_data_oceanic,rams_data_haboob,rams_data_shear,
                meta_data,meta_data_l3)


if plot_figs1_hist:
    
    i_obs_min = valid_pairs_obs.sum().iloc[1:].argmin() # Delay by 1 since 0000 UTC misses all WXT data
    col_obs_min = valid_pairs_obs.columns[875+1]
    
    histogram_data = [pairs_dict[start_time.year]['DIST'].values,
                      pairs_dict[start_time.year]['DIST'][valid_pairs_obs[col_obs_min]].values,
                      pairs_dict[start_time.year]['DIST'][valid_pairs_icon[4]].values,
                     ]
    
    vpr.figureS1(histogram_data,dist_bins)
    
    
if plot_figs2_time:
    vpr.figureS2(TT_obs.mean(axis=1),TT_icon_time,cp_onset_obs,cp_onset_icon)    


#----------------------------------------------------------------------------
print(' ')
print('*** Finshed! ***')
fst.print_runtime(t_run)    