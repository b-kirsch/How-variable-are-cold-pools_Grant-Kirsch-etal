# -*- coding: utf-8 -*-
"""
@author: Bastian Kirsch (bastian.kirsch@uni-hamburg.de)

Reading, analysis and plotting routines to be used in variogram_paper_obs_icon.py


Last updated: 5 October 2023
"""

import os
import numpy as np
import datetime as dt
import pandas as pd
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
from cmcrameri import cm as cmc
from netCDF4 import Dataset
import fesstval_routines as fst


#----------------------------------------------------------------------------
# Paths and meta data files
maindir      = '/Users/bkirs/Desktop/'

datadir_icon = maindir+'ICON_data/'
datadir_io   = maindir+'Variogram_data/'
plotdir      = maindir+'Cold-Pools/Plots/Paper_Variogram/'

#----------------------------------------------------------------------------
# Plotting settings

# Fontsize
fs = 12

# Plotting style
mpl.rcParams['font.size'] = fs
mpl.rcParams['font.sans-serif'] = ['Tahoma']
mpl.rcParams['axes.labelpad'] = 8
mpl.rcParams['axes.autolimit_mode'] = 'round_numbers'
mpl.rcParams['axes.grid'] = False
mpl.rcParams['grid.color'] = 'black'
mpl.rcParams['grid.linestyle'] = 'dotted'
mpl.rcParams['grid.alpha'] = 0.25
mpl.rcParams['xtick.labelsize'] = fs
mpl.rcParams['ytick.labelsize'] = fs
mpl.rcParams['legend.handlelength'] = 1.5
mpl.rcParams['legend.handletextpad'] = 0.5 
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top']   = False
    

# Labels
abc = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)']

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Routines

def read_icon_les(yy,mm,dd,dom,var_str='TT',datadir=datadir_icon,type_str='network',
                  return_ds=False,return_meta=False):
    
    var_dict = {'TT':'ta'}
    
    filedate = dt.date(yy,mm,dd)
    filename = datadir+var_dict[var_str]+'_'+type_str+'_icon_dom0'+str(dom)+\
               filedate.strftime('_%Y%m%d.nc')
    
    if os.path.isfile(filename) == False:
        print('File '+filename+' does not exist!')
        return pd.DataFrame()  
    
    ds = xr.open_dataset(filename)

    if return_ds: return ds
    
    val = ds[var_dict[var_str]].values
    if var_str == 'TT': val -= 273.15
    
    if type_str == 'network':
        station_ids = ds['stat_id'].values.astype(str).astype(object)
        data = pd.DataFrame(val,columns=station_ids,index=ds['time'].values)
        meta = {}
    else:
        data = np.transpose(val,(2,1,0))
        time                      = pd.to_datetime(ds['time'].values)
        lon_meshgrid,lat_meshgrid = np.meshgrid(ds['lon'].values,ds['lat'].values)
        x_meshgrid,y_meshgrid     = np.meshgrid(ds['x'].values,ds['y'].values)
        mask                      = np.transpose(ds['mask'].values.astype(bool))
        comment                   = ds['x'].comment
        i_lon                     = comment.find('E')
        lon_ref                   = float(comment[i_lon-8:i_lon])
        i_lat                     = comment.find('N')
        lat_ref                   = float(comment[i_lat-8:i_lat])
        
        meta = {'time'         : time,
                'lon_meshgrid' : lon_meshgrid,
                'lat_meshgrid' : lat_meshgrid,
                'x_meshgrid'   : x_meshgrid,
                'y_meshgrid'   : y_meshgrid,
                'mask_meshgrid': mask,
                'lon_ref'      : lon_ref,
                'lat_ref'      : lat_ref,
                }
        
    ds.close()
    
    if not return_meta:
        return data
    else:
        return meta
    
    
def read_rams_variogram(resstr,timestr,prestr='CP',datadir=datadir_io):
    # RAMS Variogram data processes by Leah
    filename = datadir+'rams/'+resstr+'/'+prestr+'_variogram_gradiogram_'+timestr+'.txt'
    
    if os.path.isfile(filename) == False:
        print('File '+filename+' does not exist!')
        return pd.DataFrame()  
    
    data = pd.read_csv(filename,header=0,delim_whitespace=True,
                       names=['DIST_MIDS','VARIOGRAM','GRADIOGRAM'])
    data['DIST_MIDS'] *= 1000
    data.set_index('DIST_MIDS',inplace=True)

    return data  

def read_rams_ideal(fname='IDEAL-DryBL-50_fig-1-data.nc',datadir=datadir_io):
    filename = datadir+'fig1/'+fname
    
    if os.path.isfile(filename) == False:
        print('File '+filename+' does not exist!')
        return pd.DataFrame()  
    
    ds = xr.open_dataset(filename)
    
    return ds

def read_rams_oceanic(res,fname='CS-TropOce-100_fig-1-data_v04.nc',datadir=datadir_io):
    if res == 100:
        filename = datadir+'fig1/'+fname
    else:   
        filename = datadir+'fig1/'+fname.replace('-100_','-'+str(res)+'_')
    
    if os.path.isfile(filename) == False:
        print('File '+filename+' does not exist!')
        return pd.DataFrame()  
    
    ds = xr.open_dataset(filename)
    
    # Cut domain of network    
    ii_lon = np.where((ds.longitude.values[0,:] >= 118.8) & \
                      (ds.longitude.values[0,:] <= 119.2))[0]
    ii_lat = np.where((ds.latitude.values[:,0] >= 17.8) & \
                      (ds.latitude.values[:,0] <= 18.2))[0]
        
    n_lon,n_lat = len(ii_lon),len(ii_lat)   
    lon_ref,lat_ref = 119,18
    
    slice_lon = slice(ii_lon[0],ii_lon[-1]+1)
    slice_lat = slice(ii_lat[0],ii_lat[-1]+1)
    
    longitude_domain = ds.longitude.values[slice_lat,slice_lon].flatten()
    latitude_domain = ds.latitude.values[slice_lat,slice_lon].flatten()
    
    # Convert to x/y coordinates
    x_domain,y_domain = fst.lonlat_to_xy(longitude_domain,latitude_domain,
                                         lon_ref=lon_ref,lat_ref=lat_ref,
                                         haversine_method=True)
    
    station_x,station_y = fst.lonlat_to_xy(ds.station_lon,ds.station_lat,
                                           lon_ref=lon_ref,lat_ref=lat_ref,
                                           haversine_method=True)
    
    data = {'x':x_domain.reshape(n_lat,n_lon),
            'y':y_domain.reshape(n_lat,n_lon),
            'dT':ds.temperature.values[slice_lat,slice_lon]-ds.T_ref.values[0],
            'station_x':station_x,
            'station_y':station_y,
            }

    return data

def read_rams_haboob(type_str,fname='IDEAL-Haboob-Day_150_fig-1-data_V3.nc',datadir=datadir_io):
    if type_str == 'Day':
        filename = datadir+'fig1/'+fname
    else:   
        filename = datadir+'fig1/'+fname.replace('-Day_','-'+type_str+'_')
    
    if os.path.isfile(filename) == False:
        print('File '+filename+' does not exist!')
        return pd.DataFrame()  
    
    ds = xr.open_dataset(filename)
   
    # Convert to x/y coordinates
    n_lat,n_lon = ds.latitude.values.shape
    lon_ref,lat_ref = ds.station_lon.values[0],ds.station_lat.values[0]
    longitude_domain = ds.longitude.values.flatten()
    latitude_domain = ds.latitude.values.flatten()
    
    x_domain,y_domain = fst.lonlat_to_xy(longitude_domain,latitude_domain,
                                          lon_ref=lon_ref,lat_ref=lat_ref,
                                          haversine_method=True)
    
    station_x,station_y = fst.lonlat_to_xy(ds.station_lon,ds.station_lat,
                                            lon_ref=lon_ref,lat_ref=lat_ref,
                                            haversine_method=True)
    
    data = {'x':x_domain.reshape(n_lat,n_lon),
            'y':y_domain.reshape(n_lat,n_lon),
            'dT':ds.temperature.values-ds.T_ref.values[0],
            'station_x':station_x,
            'station_y':station_y,
            }
    
    return data


def read_rams_shear(type_str,fname='IDEAL-Upshear-100m-fig-1-data_v2.nc',datadir=datadir_io):
    if type_str == 'Upshear':
        filename = datadir+'fig1/'+fname
    else:   
        filename = datadir+'fig1/'+fname.replace('-Upshear-','-'+type_str+'-')
    
    if os.path.isfile(filename) == False:
        print('File '+filename+' does not exist!')
        return pd.DataFrame()  
    
    ds = xr.open_dataset(filename)
    
    # Convert to x/y coordinates
    n_lat,n_lon = ds.latitude.values.shape
    lon_ref,lat_ref = ds.station_lon.values[0],ds.station_lat.values[0]
    longitude_domain = ds.longitude.values.flatten()
    latitude_domain = ds.latitude.values.flatten()
    
    x_domain,y_domain = fst.lonlat_to_xy(longitude_domain,latitude_domain,
                                          lon_ref=lon_ref,lat_ref=lat_ref,
                                          haversine_method=True)
    
    station_x,station_y = fst.lonlat_to_xy(ds.station_lon,ds.station_lat,
                                            lon_ref=lon_ref,lat_ref=lat_ref,
                                            haversine_method=True)
    
    data = {'x':x_domain.reshape(n_lat,n_lon),
            'y':y_domain.reshape(n_lat,n_lon),
            'dT':ds.temperature.values-ds.T_ref.values[0],
            'station_x':station_x,
            'station_y':station_y,
            }
    
    return data     


# Truncate colormap at given relative ranges
def truncate_colormap(cmap,minval=0.0,maxval=1.0,n=100):
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
               'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
               cmap(np.linspace(minval, maxval, n)))
    return new_cmap  
    

def figure1(obs_data,icon_data,rams_data_ideal,
            rams_data_oceanic,rams_data_haboob,rams_data_shear,
            meta_network,grid_data,pdir=plotdir):
    
    labels = ['OBS-Jogi',
              'CS-Jogi-156m',
              'CS-TropOce-100m',
              'IDEAL-DryBL-50m',
              'IDEAL-UpShear-100m',
              'IDEAL-Haboob-20KDay-150m',
              ]
    
    t0_times = [45,60,15,15,15,20]
    
    col_network = 'black'
    ms_network  = 8

    val_min,val_max,val_res = -12,2,0.5 # K
    val_ref = 0
    
    xy_lim = 18 # km
    
    levels = np.arange(val_min,val_max+val_res,val_res)
    
    maxval_truncate = ((val_ref-val_min)+val_max-val_ref)/((val_ref-val_min)*2)
    cmap = truncate_colormap(cmc.vik,maxval=maxval_truncate)
    
    print('Plotting Figure 1')
    
    fig,ax = plt.subplots(2,3,figsize=(10.2,8),dpi=400,constrained_layout=True) #11.5,9
    ax = ax.flatten()
    
    for p in range(2):
        if p == 0: pdata = obs_data 
        if p == 1: pdata = icon_data
        co = ax[p].contourf(grid_data['x_meshgrid']/1000.,grid_data['y_meshgrid']/1000.,
                            pdata,levels,cmap=cmap,extend='both')
        ax[p].scatter(meta_network['X']/1000.,meta_network['Y']/1000.,
                      c=col_network,s=ms_network)
        
    p = 2
    ax[p].contourf(rams_data_oceanic['x']/1000.,rams_data_oceanic['y']/1000.,
                   rams_data_oceanic['dT'],levels,cmap=cmap,extend='both')
    ax[p].scatter(rams_data_oceanic['station_x']/1000.,rams_data_oceanic['station_y']/1000.,
                  c=col_network,s=ms_network)
        
    p = 3
    x_shift = -5
    ax[p].contourf(rams_data_ideal.x.values+x_shift,rams_data_ideal.y.values,
                   rams_data_ideal.temperature.values-rams_data_ideal.T_ref.values,
                   levels,cmap=cmap,extend='both')
    # ax[p].scatter(rams_data_ideal.station_x.values+x_shift,rams_data_ideal.station_y.values,
    #               c=col_network,s=ms_network)
    ax[p].scatter(meta_network['X']/1000.,meta_network['Y']/1000.,
                      c=col_network,s=ms_network)
    
    p = 4
    ax[p].contourf(rams_data_shear['x']/1000.,rams_data_shear['y']/1000.,
                   rams_data_shear['dT'],levels,cmap=cmap,extend='both')
    ax[p].scatter(rams_data_shear['station_x']/1000.,rams_data_shear['station_y']/1000.,
                  c=col_network,s=ms_network)
    
    p = 5    
    ax[p].contourf(rams_data_haboob['x']/1000.,rams_data_haboob['y']/1000.,
                   rams_data_haboob['dT'],levels,cmap=cmap,extend='both')
    ax[p].scatter(rams_data_haboob['station_x']/1000.,rams_data_haboob['station_y']/1000.,
                  c=col_network,s=ms_network)
    
    # Annotations
    p = 0
    stat1,stat2 = '068Ga','086Hw'
    ax[p].annotate('',(meta_network.loc[stat1,'X']/1000.,meta_network.loc[stat1,'Y']/1000.),
                   (meta_network.loc[stat2,'X']/1000.,meta_network.loc[stat2,'Y']/1000.),
                   arrowprops={'arrowstyle':'<|-|>','lw':1.5,'color':'black'})
    ax[p].annotate(r'$d$', xy=(-10.8,-9.5), ha='center', va='center',fontsize=fs)
    ax[p].annotate(r'$T_\mathrm{j}$', xy=(-6.5,-8.0), ha='center', va='center',fontsize=fs)
    ax[p].annotate(r'$T_\mathrm{i}$', xy=(-10.8,-14.5), ha='center', va='center',fontsize=fs)
    
    for p in range(6):
        ax[p].set_title(abc[p]+' '+labels[p],fontsize=fs,loc='center',**{'weight': 'bold'})
        ax[p].set_xlabel('x (km)',fontsize=fs)
        ax[p].set_ylabel('y (km)',fontsize=fs)
        ax[p].set_xticks(np.arange(-20,25,5))
        ax[p].set_yticks(np.arange(-20,25,5))
        ax[p].set_xlim(-xy_lim,xy_lim)
        ax[p].set_ylim(-xy_lim,xy_lim)
        
        if p < 3: 
            ax[p].get_xaxis().set_visible(False) 
            ax[p].spines['bottom'].set_visible(False)
        if p not in [0,3]: 
            ax[p].get_yaxis().set_visible(False) 
            ax[p].spines['left'].set_visible(False)
            
        ax[p].annotate(r'$t_0$+'+str(t0_times[p]), xy=(xy_lim-0.5,xy_lim-0.5), 
                       ha='right', va='top',fontsize=fs)

    cbar = fig.colorbar(co,ax=ax[3:],orientation='horizontal',shrink=0.4,aspect=5)
    cbar.ax.set_xlabel('Temperature perturbation (K)')
    cbar.set_ticks(np.arange(val_min,val_max+val_res,2))
        
    plotname = pdir+'figure1_grl.png'
    fig.savefig(plotname,bbox_inches='tight')
    plt.close()      
    
    print('Plot done!')
    return

def figureS1(dist_data,plot_bins,pdir=plotdir):

    x_res = 3000
    colors = ['darkgrey','royalblue','orange']
    labels = ['Full network','OBS-Jogi (minimum)','CS-Jogi-75m']
    lw = 0.5
      
    print('Plotting distance histogram of station pairs')
    
    fig,ax = plt.subplots(1,1,figsize=(5,4),dpi=400)
    
    i = 0
    ax.hist(dist_data[i],bins=plot_bins,color=colors[i],edgecolor='black',
            lw=0.5,label=labels[i],clip_on=False)
    
    for i in range(1,3):
        ax.hist(dist_data[i],bins=plot_bins,histtype='step',color=colors[i],
                lw=lw+1,label=labels[i],clip_on=False)
        
    ax.set_xlim([0,plot_bins[-1]])
    ax.set_xlabel('Distance (km)',fontsize=fs)
    ax.set_xticks(np.arange(0,plot_bins[-1]+x_res,x_res))
    ax.set_xticklabels((ax.get_xticks()/1000).astype(int).astype(str))
    ax.set_ylim([0,250])
    ax.set_ylabel('Number of station pairs',fontsize=fs)
    ax.grid(visible=False,axis='both')
    ax.legend(fontsize=fs,loc='upper left')
    
    plt.tight_layout()
    plotname = pdir+'network_hist_grl.png'
    fig.savefig(plotname,bbox_inches='tight')
    plt.close()    
    print('Plot done!') 
    return

def figureS2(obs_data,icon_data,cp_times_obs,cp_times_icon,pdir=plotdir):

    print('Plotting time series of observational and simulation data')
    
    lw = 2
    labels = ['OBS-Jogi','CS-Jogi-625m','CS-Jogi-312m','CS-Jogi-156m','CS-Jogi-75m']
    colors = ['black','khaki','goldenrod','darkgoldenrod','saddlebrown']
    
    ymin,ymax = 16,32

    fig,ax = plt.subplots(1,1,figsize=(6,4),dpi=400)

    height = ymax-ymin
    start1 = mpl.dates.date2num(dt.datetime(2021,6,29,2,0))
    end1   = mpl.dates.date2num(dt.datetime(2021,6,29,5,0))
    width1 = end1 - start1
    
    ax.add_patch(mpl.patches.Rectangle((start1,ymin),width1,height-8,color='lightgrey')) 
    # ax.annotate('Nighttime', xy=(start1+width1/2.,height+4), 
    #             ha='center', va='bottom',rotation=90,fontsize=fs+4)
    
    start2 = mpl.dates.date2num(dt.datetime(2021,6,29,22,0))
    end2   = mpl.dates.date2num(dt.datetime(2021,6,30,0,0))
    width2 = end2 - start2
    
    ax.add_patch(mpl.patches.Rectangle((start2,ymin),width2,height,color='lightgrey')) 
    
    i = 0
    ax.plot(obs_data.index+dt.timedelta(hours=2),obs_data.values,
            lw=lw,color=colors[i],label=labels[i])
    obs_start = cp_times_obs+dt.timedelta(hours=2)
    ax.plot([obs_start,obs_start],[ymin,ymax],color=colors[i],linestyle='dashed')
    
    for i in range(1,5):
        pdata = icon_data[i].dropna()
        ax.plot(pdata.index+dt.timedelta(hours=2),pdata,
                lw=lw,color=colors[i],label=labels[i])
    icon_start = cp_times_icon+dt.timedelta(hours=2)
    ax.plot([icon_start,icon_start],[ymin,ymax],color=colors[2],linestyle='dashed')
    
    ax.set_xlim(obs_data.index[0]+dt.timedelta(hours=2),
                obs_data.index[0]+dt.timedelta(hours=24))
    ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
    ax.xaxis.set_major_locator(mpl.dates.HourLocator(byhour=range(0,24,4)))
    
    ax.set_yticks(np.arange(ymin,ymax+2,2))
    ax.set_ylim(ymin,ymax)
    
    ax.set_xlabel(r'Local time',fontsize=fs)
    ax.set_ylabel(r'Station-mean temperature (°C)',fontsize=fs)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    ax.legend(fontsize=fs,loc='upper left')

    plt.tight_layout()
    plotname = plotdir+'timeseries_grl.png'
    fig.savefig(plotname,bbox_inches='tight')
    plt.close()    
    print('Plot done!')   
    return


def write_netcdf(fname,dtime,variogram_data,bins,flag,odir=datadir_io):
    
    filename = odir+fname      
    if os.path.isfile(filename): os.remove(filename)
    
    nbin = len(bins)
    
    ncfile = Dataset(filename,'w',format='NETCDF4')
    
    # Dimensions
    ncfile.createDimension('time',None)
    ncfile.createDimension('distance',nbin)
    
    #Dimension Variables
    utime               = np.array(dtime.astype(np.int64)/1e9,dtype=int)
    time                = ncfile.createVariable('time', 'i4', ('time',))
    time[:]             = utime
    time.long_name      = 'time'
    time.units          = 'seconds since 1970-01-01 00:00:00 UTC'
    time.comment        = 'LT = UTC+2'
    
    distance            = ncfile.createVariable('distance', 'f4', ('distance',)) 
    distance[:]         = bins
    distance.long_name  = 'mid points of distance bins'
    distance.units      = 'm'
    
    variogram           = ncfile.createVariable('variogram','f4',('time','distance',),
                                                 fill_value='nan',zlib=True,
                                                 complevel=1)
    variogram[:]        = variogram_data
    variogram.long_name = 'empirical variogram of near-surface air temperature'
    variogram.units     = 'K2'
    
    cp_flag             = ncfile.createVariable('cp_flag', 'f4', ('time',)) 
    cp_flag[:]          = flag
    cp_flag.long_name   = 'cold pool flag'
    cp_flag.units       = '1'
    cp_flag.comment     = '0 = no cp, 1 = cp start, 2 = cp'
        
    # Global attributes
    ncfile.Processing_date   = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')    
    ncfile.Author            = 'Bastian Kirsch (bastian.kirsch@uni-hamburg.de)'     

    ncfile.close()
    return  