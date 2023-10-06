# -*- coding: utf-8 -*-
"""
@author: Bastian Kirsch (bastian.kirsch@uni-hamburg.de)

Functions used for processing and analyzing FESST@HH 2020 and FESSTVaL 2021 data

Last updated: 13 June 2023
"""

import os
import numpy as np
import datetime as dt
import pandas as pd
import xarray as xr
import geopy as gp
from geopy.distance import distance
import haversine as hav
import scipy.interpolate as sci

maindir = '.'

#----------------------------------------------------------------------------
# Basic functions

def print_runtime(dt_start,return_str=False):
    t_run  = np.round((dt.datetime.now()-dt_start).total_seconds())
    dt_ref = dt.datetime(2000,1,1,0,0,0) + dt.timedelta(seconds=t_run)
    if t_run < 3600: 
        print_str = dt_ref.strftime('Run time: %M:%S min')
    if (t_run >= 3600) & (t_run < 86400): 
        print_str = dt_ref.strftime('Run time: %H:%M h')
    if t_run >= 86400: 
        print_str = dt_ref.strftime('Run time: %d day(s) %H:%M h')
    print(print_str)    
    if return_str:
        return print_str
    else:
        return    

T0 = 273.15

def TK_to_TC(TTK):
    return TTK - T0

def TC_to_TK(TTC):
    return TTC + T0

def datetime_to_unixtime(dtime):
    return np.array(dtime.astype(np.int64)/1e9,dtype=int)

def unixtime_to_datetime(utime):
    return pd.to_datetime(utime*1e9)  

# Coordinates of geographical end point of line for given start coordinates, 
# distance, and direction
def geo_line(lat_start,lon_start,dist_km,dir_angle):
    # dir_angle = 0  --> to North
    # dir_angle = 90 --> to East
    start_point = gp.Point(lat_start,lon_start)
    dist = gp.distance.distance(kilometers=dist_km)
    end_point = dist.destination(point=start_point,bearing=dir_angle)
    lat_end = end_point.latitude
    lon_end = end_point.longitude
    return lat_end,lon_end

# Convert lonlat coordinates to Cartesian coordinates (input and output in m) 
def lonlat_to_xy(lon_data,lat_data,lon_ref=np.nan,lat_ref=np.nan,
                 haversine_method=False):
    
    if type(lon_data) in [int,float,np.float64]: lon_data = np.array([lon_data])
    if type(lat_data) in [int,float,np.float64]: lat_data = np.array([lat_data])
    
    if len(lon_data) != len(lat_data):
        print('LON and LAT coordinates do not have the same length!')
        return np.array([]),np.array([])
    
    npoints = len(lon_data)
    x_data = np.zeros(npoints)*np.nan
    y_data = np.zeros(npoints)*np.nan
    x_sign = np.ones(npoints)
    x_sign[lon_data < lon_ref] = -1
    y_sign = np.ones(npoints)
    y_sign[lat_data < lat_ref] = -1
    
    if np.isnan(lon_ref): lon_ref = np.nanmin(lon_data)
    if np.isnan(lat_ref): lat_ref = np.nanmin(lat_data)
    
    if haversine_method:
        ref     = list(zip(np.array([lat_ref]*npoints),np.array([lon_ref]*npoints)))
        reflon  = list(zip(np.array([lat_ref]*npoints),lon_data))
        reflat  = list(zip(lat_data,np.array([lon_ref]*npoints)))
        x_data = x_sign * hav.haversine_vector(ref,reflon,'km') * 1000
        y_data = y_sign * hav.haversine_vector(ref,reflat,'km') * 1000
        
    else:    
        for i in range(npoints):
            if np.isfinite(lon_data[i]) & np.isfinite(lat_data[i]): 
                x_data[i] = x_sign[i] * distance((lat_data[i],lon_data[i]),
                                                 (lat_data[i],lon_ref)).km * 1000
                y_data[i] = y_sign[i] * distance((lat_data[i],lon_data[i]),
                                                 (lat_ref,lon_data[i])).km * 1000
    return x_data,y_data #(m)


# Convert Cartesian coordinates to lonlat coordinates (input and output in m) 
def xy_to_lonlat(x_data,y_data,lon_ref,lat_ref,meshgrid=False,
                 haversine_method=False):
    # WARNING: No perfect inversion of lonlat_to_xy !!!!
    # x and y (m)
    if type(x_data) in [int,float,np.float64]: x_data = np.array([x_data])
    if type(y_data) in [int,float,np.float64]: y_data = np.array([y_data])
    
    if x_data.shape != y_data.shape:
        print('x and y coordinates do not have the same shape!')
        return np.array([]),np.array([])
    
    if meshgrid == False:
        npoints = len(x_data)
        lon_data = np.zeros(npoints)*np.nan
        lat_data = np.zeros(npoints)*np.nan
        
        if haversine_method:
            for i in range(npoints):
                if np.isnan(x_data[i]) or np.isnan(y_data[i]): continue
                # Dir in radians (=deg * pi/180)
                x_dir = 0.5 * np.pi if x_data[i] >= 0 else 1.5 * np.pi
                y_dir = 0 if y_data[i] >= 0 else np.pi
                lon_data[i] = hav.inverse_haversine((lat_ref,lon_ref),np.abs(x_data[i]/1000.),x_dir,'km')[1]
                lat_data[i] = hav.inverse_haversine((lat_ref,lon_ref),np.abs(y_data[i]/1000.),y_dir,'km')[0] 
        else:     
            for i in range(npoints):
                if np.isnan(x_data[i]) or np.isnan(y_data[i]): continue
                x_dir = 90 if x_data[i] >= 0 else 270
                y_dir = 0 if y_data[i] >= 0 else 180
                lon_data[i] = geo_line(lat_ref,lon_ref,np.abs(x_data[i]/1000.),x_dir)[1]
                lat_data[i] = geo_line(lat_ref,lon_ref,np.abs(y_data[i]/1000.),y_dir)[0]
            
    else:
        ydim, xdim = x_data.shape    
        lon_data = np.zeros([ydim,xdim])*np.nan
        lat_data = np.zeros([ydim,xdim])*np.nan
        
        if haversine_method:
            for j in range(ydim):
                for i in range(xdim):
                    if np.isnan(x_data[j,i]) or np.isnan(y_data[j,i]): continue
                    # Dir in radians (=deg * pi/180)
                    x_dir = 0.5 * np.pi if x_data[j,i] >= 0 else 1.5 * np.pi
                    y_dir = 0 if y_data[j,i] >= 0 else np.pi
                    lon_data[j,i] = hav.inverse_haversine((lat_ref,lon_ref),np.abs(x_data[j,i]/1000.),x_dir,'km')[1]
                    lat_data[j,i] = hav.inverse_haversine((lat_ref,lon_ref),np.abs(y_data[j,i]/1000.),y_dir,'km')[0]   
        else:
            for j in range(ydim):
                for i in range(xdim):
                    if np.isnan(x_data[j,i]) or np.isnan(y_data[j,i]): continue
                    x_dir = 90 if x_data[j,i] >= 0 else 270
                    y_dir = 0 if y_data[j,i] >= 0 else 180
                    lon_data[j,i] = geo_line(lat_ref,lon_ref,np.abs(x_data[j,i]/1000.),x_dir)[1]
                    lat_data[j,i] = geo_line(lat_ref,lon_ref,np.abs(y_data[j,i]/1000.),y_dir)[0]
            
    return lon_data,lat_data

    
#----------------------------------------------------------------------------
# Read functions    

def read_fesstval_level2(stat_type,var_str,yy,mm,dd,ds_version=0,
                         mute=False,return_ds=False,return_meta=False,
                         datadir_apollo=maindir+'APOLLO_data/level2/',
                         datadir_wxt=maindir+'WXT_data/level2/'):
    
    # Read FESST@HH level2 APOLLO and WXT data
    stat_type = stat_type.upper()
    if stat_type == 'A': stat_type = 'APOLLO'
    if stat_type == 'W': stat_type = 'WXT'
    if stat_type not in ['APOLLO','WXT']: stat_type = 'APOLLO'
    campaign = 'fval' 
    if yy == 2019: campaign = 'prefval'
    if yy == 2020: campaign = 'fessthh'
    
    var_dict = {'TT':'ta',
                'PP':'pa',
                'RH':'hur',
                'FF':'wspeed',
                'FB':'wspeed_max',
                'DD':'wdir',
                'RR':'precip',
                'HA':'hail',
                }
    
    if stat_type == 'APOLLO':
        datadir  = datadir_apollo
        var_list = ['TT','PP']

    if stat_type == 'WXT':
        datadir  = datadir_wxt
        var_list = list(var_dict.keys())
        
    if var_str not in var_list:
        print('Variable '+var_str+' is not available! Pick one of '+str(var_list)) 
        return pd.DataFrame()
        
    filedate = dt.date(yy,mm,dd)  # ('%Y/%m/%d/')
    filename = datadir+filedate.strftime('%Y/%m/')+campaign+'_uhh_'+stat_type.lower()+\
               '_l2_'+var_dict[var_str]+'_v'+str(ds_version).zfill(2)+'_'+\
               filedate.strftime('%Y%m%d000000.nc')
                    
    if os.path.isfile(filename) == False:
        if mute == False: print('File '+filename+' does not exist!')
        return pd.DataFrame()  
    
    try:
        ds = xr.open_dataset(filename)
    except ValueError:
        if mute == False: print('File '+filename+' could not be correctly read!')
        return pd.DataFrame() 
    
    if return_ds: return ds
        
    val = ds[var_dict[var_str]].values
    station_ids = ds['station_id'].values.astype(str).astype(object)
    
    if var_str == 'TT': val  = TK_to_TC(val)
    if var_str == 'PP': val *= 0.01
    if var_str == 'RH': val *= 100
    
    data = pd.DataFrame(val,columns=station_ids,index=ds['time'].values)

    if not return_meta:
        return data#.round(2)
    else:
        meta = pd.DataFrame(index=station_ids)
        meta['LAT']   = ds['lat'].values
        meta['LON']   = ds['lon'].values
        meta['ALT']   = ds['zsl'].values
        meta['HT']    = ds['zag'].values
        meta['NAME']  = ds['station_name'].values.astype(str).astype(object)
        #meta['INSTR'] = ds['instrument_id'].values
        meta['LCZ']   = ds['lcz'].values.astype(str).astype(object)
        
        return data,meta #.round(2)
    
    
def read_fesstval_level3(var_str,yy,mm,dd,ds_version=1,check_file=False,
                         mute=False,return_ds=False,return_meta=False,
                         datadir_apollo=maindir+'APOLLO_data/level3/',
                         datadir_wxt=maindir+'WXT_data/level3/'):
    
    # Read FESST@HH level3 APOLLO and WXT data
    campaign = 'fessthh' if yy == 2020 else 'fval'
    
    var_dict = {'DT':'ta',
                'TT':'ta',
                'TT_REF':'ta',
                'DP':'pa',
                }
    var_list = list(var_dict.keys())
    
    if var_str in ['DT','TT','TT_REF','DP']:
        datadir,stat_type = datadir_apollo,'APOLLOWXT' 
    else:
        datadir,stat_type = datadir_wxt,'WXT'
        
    if var_str not in var_list:
        print('Variable '+var_str+' is not available! Pick one of '+str(var_list)) 
        return pd.DataFrame()
        
    filedate = dt.date(yy,mm,dd)
    filename = datadir+filedate.strftime('%Y/%m/')+campaign+'_uhh_'+stat_type.lower()+\
               '_l3_'+var_dict[var_str]+'_v'+str(ds_version).zfill(2)+'_'+\
               filedate.strftime('%Y%m%d.nc')
    
    if check_file:
        return os.path.isfile(filename)           
                    
    if os.path.isfile(filename) == False:
        if mute == False: print('File '+filename+' does not exist!')
        return pd.DataFrame()  
    
    try:
        ds = xr.open_dataset(filename)
    except ValueError:
        if mute == False: print('File '+filename+' could not be correctly read!')
        return pd.DataFrame() 
    
    if return_ds: return ds
    
    if not return_meta:
        if var_str == 'DT':
            return np.transpose(ds['ta_perturbation'].values,(2,1,0))
        if var_str == 'TT_REF':
            return ds['ta_reference'].values-T0
        if var_str == 'TT':
            return np.add(ds['ta_reference'].values,np.transpose(ds['ta_perturbation'].values,(2,1,0)))-T0
        if var_str == 'DP':
            return np.transpose(ds['pa_perturbation'].values,(2,1,0)) * 0.01 #Pa->hPa

    else:
        time                      = pd.to_datetime(ds['time'].values)
        lon_meshgrid,lat_meshgrid = np.meshgrid(ds['lon'].values,ds['lat'].values)
        x_meshgrid,y_meshgrid     = np.meshgrid(ds['x'].values,ds['y'].values)
        mask                      = np.transpose(ds['mask'].values.astype(bool))
        comment                   = ds['x'].comment
        i_lon                     = comment.find('E')
        lon_ref                   = float(comment[i_lon-8:i_lon])
        i_lat                     = comment.find('N')
        lat_ref                   = float(comment[i_lat-8:i_lat])
        event_id                  = ds['event_id'].values
        
        meta = {'time'         : time,
                'lon_meshgrid' : lon_meshgrid,
                'lat_meshgrid' : lat_meshgrid,
                'x_meshgrid'   : x_meshgrid,
                'y_meshgrid'   : y_meshgrid,
                'mask_meshgrid': mask,
                'lon_ref'      : lon_ref,
                'lat_ref'      : lat_ref,
                'event_id'     : event_id,
                }
        
        return meta


# Meta data of FESST@HH 2020 stations
def fessthh_stations(nn,find='NAME',metafile=maindir+'FESSTVaL/FESSTHH/stations_fessthh.txt',
                     serialfile=maindir+'FESSTVaL/APOLLO/APOLLO_Seriennummern.txt',
                     include_ht=False,include_time=False,include_serial=False,
                     include_lcz=False):

    if find.upper() == 'A': find = 'APOLLO'
    if find.upper() == 'W': find = 'WXT'
    if find.upper() == 'N': find = 'NAME'
    if find.upper() == 'S': find = 'SERIAL'
    if find.upper() == 'L': find = 'LCZ'
    
    if find == 'SERIAL': include_serial = True
    if find == 'LCZ': include_lcz = True
    
    if find not in ['STATION','NAME','APOLLO','WXT','SERIAL','LCZ']: find = 'NAME'
    if (find == 'NAME' or find == 'SERIAL' or find == 'LCZ') and (type(nn) != str):
        try:
            nn = str(nn)
        except ValueError:     
            print(find+' search item has to be a string!')
            return pd.DataFrame()
    if (find in ['STATION','APOLLO','WXT']) and (type(nn) != int):
        try:
            nn = int(nn)
        except ValueError:    
            print(find+' search item has to be an integer!')
            return pd.DataFrame()
    
    names = ['STATION','NAME','LAT','LON','ALT','APOLLO','WXT',
             'HT_TT','HT_PP','HT_WXT','START','END','LCZ']
    dtypes = [str,str,float,float,float,str,str,float,float,float,str,str,str]
    data = pd.read_csv(metafile,sep=';',header=0,names=names,dtype=dict(zip(names,dtypes)))

    data['APOLLO'] = data['APOLLO'].apply(lambda x: int(x) if type(x) == str else 0)
    data['WXT']    = data['WXT'].apply(lambda x: int(x) if type(x) == str else 0)
    data['START']  = data['START'].apply(lambda x: pd.to_datetime(x,format='%Y%m%d %H%M'))
    data['END']    = data['END'].apply(lambda x: pd.to_datetime(x,format='%Y%m%d %H%M'))
    data['END'].fillna(dt.datetime(2099,12,31,23,59,59),inplace=True)
    if include_serial:
        data['SERIAL'] = data['APOLLO'].apply(lambda x: apollo_serial(x,serialfile=serialfile) \
                                              if x > 0 else '')
    data['LCZ']    = data['LCZ'].apply(lambda x: x if x != np.nan else '')
        
    if find == 'STATION': 
        ii_find = data[find].str.startswith(str(nn).zfill(3)).values 
    if find == 'NAME': 
        if nn == '':
            # all measurement locations without duplicates due to instrument replacements
            # (last status of network)
            ii_find = (data['STATION'].str[:3].to_frame().duplicated(keep='last') == False)
        if nn == 'l2':
            # all measurement location without duplicates due to instrument replacements
            # that have existed once (level 2 meta data)
            ii_find = (data['STATION'].duplicated() == False)
        elif nn == 'all':
            # all measurement locations with duplicates due to instrument replacements
            ii_find = data['STATION'].notnull()
        else:    
            ii_find = data[find].str.contains(nn,case=False).values
    if find in ['APOLLO','WXT']:
        ii_find = (data[find] == nn).values
    if find == 'SERIAL':
        ii_find = data[find].str.startswith(nn).values
    if find == 'LCZ':
        ii_find = (data[find] == nn).values
    
    if not include_ht:
        data.drop('HT_TT',1,inplace=True)
        data.drop('HT_PP',1,inplace=True)
        data.drop('HT_WXT',1,inplace=True)
    if not include_time:
        data.drop('START',1,inplace=True)
        data.drop('END',1,inplace=True)    
    if not include_lcz:
        data.drop('LCZ',1,inplace=True)    
    if find in ['APOLLO','SERIAL']:    
        data.drop('WXT',1,inplace=True) 
    if find == 'WXT':    
        data.drop('APOLLO',1,inplace=True)     

    return data[ii_find]

# Meta data of FESSTVaL 2021 stations
def fesstval_stations(nn,find='STATION',metafile=maindir+'FESSTVaL/stations_fesstval.txt',
                      serialfile=maindir+'FESSTVaL/APOLLO/APOLLO_Seriennummern.txt',
                      include_time=False,include_serial=False,include_lcz=False):
    
    
    if find.upper() == 'A': find = 'APOLLO'
    if find.upper() == 'W': find = 'WXT'
    if find.upper() == 'S': find = 'SERIAL'
    if find.upper() == 'L': find = 'LCZ'
    if find.upper() == 'N': find = 'NAME'
    
    if find == 'SERIAL': include_serial = True
    if find == 'LCZ': include_lcz = True
    
    if find not in ['STATION','NAME','APOLLO','WXT','SERIAL','LCZ']: find = 'STATION'
    
    names_dtypes = {'STATION':str,
                    'NAME'   :str,
                    'LAT'    :float,
                    'LON'    :float,
                    'ALT'    :float,
                    'APOLLO' :str,
                    'WXT'    :str,
                    'START'  :str,
                    'END'    :str,
                    'LCZ'    :str,
                    'MAINT'  :str,
                   }
    
    names = list(names_dtypes.keys())

    data = pd.read_csv(metafile,sep=';',header=0,names=names,dtype=names_dtypes)
    
    
    data['APOLLO'] = data['APOLLO'].apply(lambda x: int(x) if type(x) == str else 0)
    data['WXT']    = data['WXT'].apply(lambda x: int(x) if type(x) == str else 0)
    data['START']  = data['START'].apply(lambda x: pd.to_datetime(x,format='%Y%m%d %H%M'))
    data['END']    = data['END'].apply(lambda x: pd.to_datetime(x,format='%Y%m%d %H%M'))
    if include_serial:
        data['SERIAL'] = data['APOLLO'].apply(lambda x: apollo_serial(x,serialfile=serialfile) \
                                              if x > 0 else '')   
    
    if find == 'STATION': 
        #ii_find = data[find].str.startswith(str(nn).zfill(3)).values 
        ii_find = data[find].str.contains(nn,case=False).values
        
    if find == 'NAME':
         ii_find = data[find].str.contains(nn,case=False).values
         
    if find in ['APOLLO','WXT']:
        ii_find = (data[find] == nn).values 
        
    if find == 'SERIAL':
        ii_find = data[find].str.startswith(nn).values
        
    if nn == '': 
        ii_find = (data['STATION'].to_frame().duplicated(keep='last') == False)
        
    if nn == 'all':
        ii_find = data['STATION'].notnull() 
        
    if not include_time:
        data.drop('START',1,inplace=True)
        data.drop('END',1,inplace=True)    
    if not include_lcz :
        data.drop('LCZ',1,inplace=True)    
    if find in ['APOLLO','SERIAL']:    
        data.drop('WXT',1,inplace=True) 
    if find == 'WXT':    
        data.drop('APOLLO',1,inplace=True) 
        
    data.drop('MAINT',1,inplace=True)  
    
    return data[ii_find]


# List of cold pool events during FESST@HH 2020
def fessthh_cold_pools(cpfile=maindir+'FESSTVaL/FESSTHH/cold_pools_fessthh.txt'):
    data_columns_type = {'#DATE':str,
                         'TIME(UTC)':str,
                         'NAME':str,                     
                         'DIR_OF_MOVEMENT':str,
                         'DESCRIPTION'  :str,
                         }
    cps   = pd.read_csv(cpfile,sep=';',dtype=data_columns_type)
    dtime = pd.to_datetime(cps['#DATE']+cps['TIME(UTC)'],format='%Y%m%d%H%M')
    cps.set_index(dtime,inplace=True)
    cps.drop(['#DATE','TIME(UTC)'],axis=1,inplace=True)
    return cps

# List of cold pool events during FESSTVaL 2021
def fesstval_cold_pools(cpfile=maindir+'FESSTVaL/cold_pools_fesstval.txt'):
    data_columns_type = {'#DATE':str,
                         'TIME(UTC)':str,
                         'NAME':str,                     
                         'DIR_OF_MOVEMENT':str,
                         'DESCRIPTION'  :str,
                         }
    cps   = pd.read_csv(cpfile,sep=';',dtype=data_columns_type)
    dtime = pd.to_datetime(cps['#DATE']+cps['TIME(UTC)'],format='%Y%m%d%H%M')
    cps.set_index(dtime,inplace=True)
    cps.drop(['#DATE','TIME(UTC)'],axis=1,inplace=True)
    cps['INDEX'] = np.arange(cps.shape[0])
    return cps 


class read_xband_radar:
    
    def __init__(self,yy,mm,dd,hh,datadir=maindir+'Radar_data/level2/',
                 x_res=None,x_min=None,x_max=None,
                 y_res=None,y_min=None,y_max=None,
                 lon_data_ref=None,lat_data_ref=None,
                 hav=True,mute=False):
        
        filetime = dt.datetime(yy,mm,dd,hh)
        loc = 'uhh_wrxHHG' if yy == 2020 else 'fval_uhh_wrxFLK'
        filename = datadir+\
                   filetime.strftime('%Y/%m/'+loc+'_l2_rr_v00_%Y%m%d%H.nc')
        dtimes = pd.date_range(filetime,freq='30s',periods=120)           
        
        ds, valid_ds, valid_grid = xr.Dataset(), False, False
        if not os.path.isfile(filename):
            if not mute: print('File '+filename+' does not exist!')
        else:    
            try:
                ds = xr.open_dataset(filename)
                valid_ds = True
            except ValueError:
                if not mute: print('File '+filename+' could not be correctly read!')
             
        if valid_ds:
            # Convert Lon/Lat data into xy-coordinates
            lon_ref = ds['lon_center'].values*1 if not lon_data_ref else lon_data_ref
            lat_ref = ds['lat_center'].values*1 if not lat_data_ref else lat_data_ref
            
            lon_data = np.ravel(ds['lon'].values)
            lat_data = np.ravel(ds['lat'].values)
            
            x_data,y_data = lonlat_to_xy(lon_data,lat_data,
                                         lon_ref=lon_ref,lat_ref=lat_ref,
                                         haversine_method=hav)
            
            # Set up xy grid
            if not x_res or not y_res:
                x_grid = np.array([])
                y_grid = np.array([])
            else:
                x_res_dim = len(str(int(x_res)))   
                y_res_dim = len(str(int(y_res)))
                
                if not x_min: x_min = np.round(x_data.min(),-x_res_dim+1)-x_res
                if not x_max: x_max = np.round(x_data.max(),-x_res_dim+1)+x_res
                if not y_min: y_min = np.round(y_data.min(),-y_res_dim+1)-y_res
                if not y_max: y_max = np.round(y_data.max(),-y_res_dim+1)+y_res
                
                x_grid = np.arange(x_min,x_max+x_res,x_res,dtype=float)
                y_grid = np.arange(y_min,y_max+y_res,y_res,dtype=float)
                x_meshgrid,y_meshgrid = np.meshgrid(x_grid,y_grid)
                valid_grid = True
                
                self.lon_ref    = lon_ref
                self.lat_ref    = lat_ref
                self.x_data     = x_data
                self.y_data     = y_data
                self.x_meshgrid = x_meshgrid
                self.y_meshgrid = y_meshgrid
                self.hav        = hav
            
        self.ds         = ds
        self.dtimes     = dtimes
        self.valid_grid = valid_grid
        self.mute       = mute
        
        
    def dataset(self):
        return self.ds
    
    def lonlat_polar(self):
        lat = self.ds['lat'].to_pandas().transpose()     
        lon = self.ds['lon'].to_pandas()#.transpose()
        return lon, lat
    
    def xy_reggrid(self):
        if self.valid_grid:
            return self.x_meshgrid,self.y_meshgrid
        else:
            if not self.mute: print('No valid grid defined!')
            return np.array([]),np.array([])
    
    def lonlat_reggrid(self):
        if self.valid_grid:
            lon_meshgrid,lat_meshgrid = xy_to_lonlat(self.x_meshgrid,self.y_meshgrid,
                                                     self.lon_ref,self.lat_ref,
                                                     meshgrid=True,
                                                     haversine_method=self.hav)
            return lon_meshgrid,lat_meshgrid
        else:
            if not self.mute: print('No valid grid defined!')
            return np.array([]),np.array([])    
        
    def mask_meshgrid(self):
        dist_mask_meshgrid = np.zeros_like(self.x_meshgrid,dtype=bool)
        lon_center = self.ds['lon_center'].values*1 
        lat_center = self.ds['lat_center'].values*1 
        
        x_center,y_center = lonlat_to_xy(lon_center,lat_center,
                                         lon_ref=self.lon_ref,lat_ref=self.lat_ref,
                                         haversine_method=self.hav)
        
        dist = np.sqrt(np.add(np.power(np.subtract(x_center,self.x_meshgrid),2),
                              np.power(np.subtract(y_center,self.y_meshgrid),2)))
        range_max = self.ds.range.values[-1] + 30 # half range gate
        dist_mask_meshgrid[dist <= range_max] = True
        return dist_mask_meshgrid    
    
    def itime(self,m,s):
        it = np.where((self.dtimes.minute == m) & (self.dtimes.second == s))[0]
        if len(it) == 0:
            if not self.mute: print('Selected time step is not available!')
            return -1
        else:
            return it[0]
    
    def rr_polar(self,minute,sec):
        it = self.itime(minute,sec)
        if it == -1: return pd.DataFrame()   
        data = self.ds['rr'].isel(time=it).drop_vars(['time']).to_pandas()#.transpose()         
        return data
    
    def rr_reggrid(self,minute,sec):
        if self.valid_grid:
            it = self.itime(minute,sec)
            if it == -1: return pd.DataFrame()  
            rr = np.ravel(self.ds['rr'].values[it,:,:])
            NN_interpol = sci.NearestNDInterpolator(list(zip(self.x_data,self.y_data)),rr)
            rr_interpol = NN_interpol(self.x_meshgrid,self.y_meshgrid)
            rr_mask = self.mask_meshgrid() 
            rr_interpol[~rr_mask] = np.nan
            return rr_interpol
        else:
            if not self.mute: print('No valid grid defined!')
            return np.array([])
        
class read_xband_radar_reggrid:
    
    def __init__(self,yy,mm,dd,hh,datadir=maindir+'Radar_data/level2/',mute=False):
        
        filetime = dt.datetime(yy,mm,dd,hh)
        loc = 'HHG' if yy == 2020 else 'FLK'
        filename = datadir+filetime.strftime('%Y/%m/wrx'+loc+'_rr_reggrid_%Y%m%d%H.nc')
                  
        ds, valid_ds = xr.Dataset(), False
        if not os.path.isfile(filename):
            if not mute: print('File '+filename+' does not exist!')
        else:    
            try:
                ds = xr.open_dataset(filename)
                valid_ds = True
            except ValueError:
                if not mute: print('File '+filename+' could not be correctly read!')
        
        if valid_ds:
            dtimes = pd.DatetimeIndex(ds.time.values)
            x_meshgrid,y_meshgrid = np.meshgrid(ds.x.values.astype(float),
                                                ds.y.values.astype(float))
            lon_ref = float(ds.x.comment.split(' ')[3][:-1])
            lat_ref = float(ds.x.comment.split(' ')[2][:-2])
            
            self.dtimes     = dtimes
            self.lon_ref    = lon_ref
            self.lat_ref    = lat_ref
            self.x_meshgrid = x_meshgrid
            self.y_meshgrid = y_meshgrid
        
        self.valid_ds = valid_ds
        self.ds       = ds
        self.mute     = mute
        
    def dataset(self):
        return self.ds    
    
    def xy_reggrid(self):
        if self.valid_ds:
            return self.x_meshgrid,self.y_meshgrid
        else:
            if not self.mute: print('No valid grid defined!')
            return np.array([]),np.array([])
            
    def itime(self,m,s):
        it = np.where((self.dtimes.minute == m) & (self.dtimes.second == s))[0]
        if len(it) == 0:
            if not self.mute: print('Selected time step is not available!')
            return -1
        else:
            return it[0]    
        
    def rr_reggrid(self,minute,sec):
        it = self.itime(minute,sec)
        if it == -1: return np.array([])   
        data = self.ds['rr'].isel(time=it).values      
        return data
    
    
def read_radar_level3(yy,mm,dd,ds_version=0,mute=False,
                      return_ds=False,return_meta=False,
                      datadir=maindir+'Radar_data/level3/'):
    
    campaign = 'fessthh' if yy == 2020 else 'fval'
        
    filedate = dt.date(yy,mm,dd)
    filename = datadir+filedate.strftime('%Y/%m/')+campaign+'_uhh_wrx_l3_rr_'+\
               'v'+str(ds_version).zfill(2)+'_'+filedate.strftime('%Y%m%d.nc')
                    
    if os.path.isfile(filename) == False:
        if mute == False: print('File '+filename+' does not exist!')
        return pd.DataFrame()  
    
    try:
        ds = xr.open_dataset(filename)
    except ValueError:
        if mute == False: print('File '+filename+' could not be correctly read!')
        return pd.DataFrame() 
    
    if return_ds: return ds
    
    if not return_meta:
        return np.transpose(ds['rr'].values,(2,1,0))
    else:
        time                      = pd.to_datetime(ds['time'].values)
        lon_meshgrid,lat_meshgrid = np.meshgrid(ds['lon'].values,ds['lat'].values)
        x_meshgrid,y_meshgrid     = np.meshgrid(ds['x'].values,ds['y'].values)
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
                'lon_ref'      : lon_ref,
                'lat_ref'      : lat_ref,
                }
        
        return meta    
    