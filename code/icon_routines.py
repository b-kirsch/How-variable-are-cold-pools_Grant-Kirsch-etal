import xarray as xr
import pandas as pd
import datetime as dt
import numpy as np
import os
import scipy.interpolate as sci
from netCDF4 import Dataset
from netCDF4 import stringtochar
import matplotlib as mpl
import matplotlib.pyplot as plt

# import icon_routines as icr
# import importlib
# importlib.reload(icr) # every run

# adding new libraries:
# conda info --envs
# source activate env_fval
# conda install xxx

icondir = '/pool/data/fesstval/SIM/ICON_LES/'
maindir = '/work/um0203/u300510/fesstval/'

grid_dict = {1:'GRIDS/grid_FESSTVaL_0625m.nc',
             2:'GRIDS/grid_FESSTVaL_0312m.nc',
             3:'GRIDS/grid_FESSTVaL_0156m.nc',
             4:'GRIDS/grid_FESSTVaL_0075m.nc'}

var_dict_read = {'TT':'t_2m',
                 'TD':'td_2m',
                 'RH':'rh_2m',
                 'PP':'pres_sfc',
                 'PP_MSL':'pres_msl',
                 'UU':'u_10m',
                 'VV':'v_10m',
                 'SHF':'shfl_s',
                 'LHF':'lhfl_s',
                 'TKE':'tke',
                 'RR':'tot_prec',
                 'CLC':'clct',
                }

var_dict_write = {'TT':'ta',
                 }

def read_les(yy,mm,dd,hh,mi,var_str,domain,datadir=icondir,
             var_dict=var_dict_read,return_ds=False,mute=False):
    
    var_list = list(var_dict.keys())
        
    if var_str not in var_list:
        if not mute: print('Variable '+var_str+' is not available! Pick one of '+str(var_list)) 
        return pd.DataFrame()
    
    filedate = dt.datetime(yy,mm,dd,hh,mi)
    filename = datadir+filedate.strftime('FESSTVaLD2_Lind_corine_%Y%m%d/FESSTVaLD2_00%H%M00_DOM0')+str(domain)+'.nc'
    
    if os.path.isfile(filename) == False:
        if not mute: print('File '+filename+' does not exist!')
        return pd.DataFrame()  
    
    try:
        ds = xr.open_dataset(filename)
    except ValueError:
        if not mute: print('File '+filename+' could not be correctly read!')
        return pd.DataFrame()
    
    if return_ds: return ds
    
    data = ds[var_dict[var_str]].values[0,0,:]
    
    if var_str == 'TT': data -= 273.15
    
    return pd.DataFrame(data)

def check_sim_day(day,datadir=icondir):
    # check if given day was simulated
    daydir = datadir+day.strftime('FESSTVaLD2_Lind_corine_%Y%m%d/')
    return os.path.isdir(daydir)


def read_grid(domain,gridfiles=grid_dict,datadir=icondir,return_ds=False):
    
    if domain not in list(gridfiles.keys()):
        print('Domain '+str(domain)+' not available!')
        return pd.DataFrame()
    filename = datadir+gridfiles[domain]
    
    if os.path.isfile(filename) == False:
        print('File '+filename+' does not exist!')
        return pd.DataFrame()  
    
    try:
        ds = xr.open_dataset(filename)
    except ValueError:
        print('File '+filename+' could not be correctly read!')
        return pd.DataFrame()
    
    if return_ds: return ds

    data = pd.DataFrame()
    for c in ['lon_cell_centre','lat_cell_centre']:
        data[c] = ds[c].values * (180/np.pi)
    
    return data


def read_network(filename='icon/fesstval_network.txt'):
    meta = pd.read_csv(filename,sep=';',header=0)
    meta.set_index('NAME',inplace=True)
    return meta


def interpolate_data(x_data,y_data,data,x_ip,y_ip,ip_type='lin'):
    xy_data = list(zip(x_data,y_data))
    if ip_type == 'lin':
        interp = sci.LinearNDInterpolator(xy_data, data)
    else:    
        interp = sci.NearestNDInterpolator(xy_data, data)
    data_ip = interp(x_ip,y_ip)
    return data_ip


def datetime_to_unixtime(dtime):
    return np.array(dtime.astype(np.int64)/1e9,dtype=int)


def write_nc_network(data,lon_data,lat_data,day,domain,var_str='TT',
                     var_dict=var_dict_write,writedir=maindir+'icon/network_interpolated/'):
    
    print('Writing nc file (network)')
    
    filename = writedir+var_dict[var_str]+'_network_icon_dom0'+str(domain)+day.strftime('_%Y%m%d.nc')
    if os.path.isfile(filename): os.remove(filename)

    ntime,nstat = data.shape
    nchar = 5
    ncfile = Dataset(filename,'w',format='NETCDF4')
    
    # Dimensions
    ncfile.createDimension('time',ntime)
    ncfile.createDimension('stat',nstat)
    ncfile.createDimension('char',nchar)
    
    #Dimension Variables
    utime               = datetime_to_unixtime(data.index)
    time                = ncfile.createVariable('time', 'i4', ('time',))
    time[:]             = utime
    time.standard_name  = 'time'
    time.units          = 'seconds since 1970-01-01 00:00:00 UTC' 
    time.calendar       = 'standard'
    
    lon                 = ncfile.createVariable('lon', 'f4', ('stat',)) 
    lon[:]              = lon_data
    lon.standard_name   = 'longitude'
    lon.units           = 'degrees_east'

    lat                 = ncfile.createVariable('lat', 'f4', ('stat',)) 
    lat[:]              = lat_data
    lat.standard_name   = 'latitude'
    lat.units           = 'degrees_north'
    
    station_id              = ncfile.createVariable('stat_id', 'S1', ('stat','char',))
    station_id[:]           = stringtochar(data.columns.to_numpy(dtype='S'+str(nchar)))
    station_id.long_name    = 'station identifier code'
    
    # Variables
    if var_str == 'TT':
        ta               = ncfile.createVariable('ta','f4',('time','stat',),
                                                 fill_value='nan',zlib=True,
                                                 complevel=1,least_significant_digit=3)
        ta[:,:]          = data + 273.15
        ta.standard_name = 'air_temperature'
        ta.long_name     = 'air temperature'
        ta.units         = 'K'
        ta.comment       = '2-m air temperature'

    # Global attributes
    ncfile.Title             = 'ICON LES FESSTVaLD2 simulation interpolated to station network'
    ncfile.Institution       = 'Meteorological Institute, University of Hamburg (UHH), Germany'
    ncfile.Contact_person    = 'Bastian Kirsch (bastian.kirsch@uni-hamburg.de)'
    ncfile.Source            = 'Simulation data located on Levante at /pool/data/fesstval/SIM/ICON_LES/'
    ncfile.History           = 'None'
    ncfile.Conventions       = 'CF-1.7 where applicable'
    ncfile.Processing_date   = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')    
    ncfile.Author            = 'Bastian Kirsch (bastian.kirsch@uni-hamburg.de)'     
    ncfile.Comments          = 'None'
    ncfile.Licence           = 'This data is licensed under a '+\
                               'Creative Commons Attribution 4.0 '+\
                               'International License (CC BY 4.0).'  
    ncfile.close()
    print('Writing successful')
    return


def write_nc_grid(data,time,grid,day,domain,var_str='TT',
                  var_dict=var_dict_write,writedir=maindir+'icon/grid_interpolated/'):
    
    print('Writing nc file (grid)')
    
    filename = writedir+var_dict[var_str]+'_grid_icon_dom0'+str(domain)+day.strftime('_%Y%m%d.nc')
    if os.path.isfile(filename): os.remove(filename)
    
    lon_data  = grid['lon'][0,:]   
    lat_data  = grid['lat'][:,0]  
    x_data    = grid['x'][0,:]   
    y_data    = grid['y'][:,0]
    mask_data = grid['mask']
    ref_str   = '{:8.5f}E, '.format(grid['lon_ref'])+'{:8.5f}N'.format(grid['lat_ref'])
               
    ny,nx,ntime = data.shape

    ncfile = Dataset(filename,'w',format='NETCDF4')
    
    # Dimensions
    ncfile.createDimension('time',ntime)
    ncfile.createDimension('lon',nx)
    ncfile.createDimension('lat',ny)
    
    #Dimension Variables
    utime               = datetime_to_unixtime(time)
    time                = ncfile.createVariable('time', 'i4', ('time',))
    time[:]             = utime
    time.standard_name  = 'time'
    time.units          = 'seconds since 1970-01-01 00:00:00 UTC' 
    time.calendar       = 'standard'
    
    lon                 = ncfile.createVariable('lon', 'f4', ('lon',)) 
    lon[:]              = lon_data
    lon.standard_name   = 'longitude'
    lon.units           = 'degrees_east'

    lat                 = ncfile.createVariable('lat', 'f4', ('lat',)) 
    lat[:]              = lat_data
    lat.standard_name   = 'latitude'
    lat.units           = 'degrees_north'
    
    x                   = ncfile.createVariable('x', 'f4', ('lon',)) 
    x[:]                = x_data
    x.long_name         = 'x coordinates'
    x.units             = 'm'
    x.comment           = 'origin at '+ref_str

    y                   = ncfile.createVariable('y', 'f4', ('lat',)) 
    y[:]                = y_data
    y.long_name         = 'y coordinates'
    y.units             = 'm'
    y.comment           = 'origin at '+ref_str
    
    mask                = ncfile.createVariable('mask', 'i2', ('lon','lat',)) 
    mask[:,:]           = np.transpose(mask_data,(1,0))
    mask.long_name      = 'interpolation mask'
    mask.comment        = 'binary mask of valid interpolation grid points'
    
    # Variables
    if var_str == 'TT':
        ta               = ncfile.createVariable('ta','f4',('time','lon','lat',),
                                                 fill_value='nan',zlib=True,
                                                 complevel=1,least_significant_digit=3)
        ta[:,:,:]        = np.transpose(data,(2,1,0)) + 273.15
        ta.standard_name = 'air_temperature'
        ta.long_name     = 'air temperature'
        ta.units         = 'K'
        ta.comment       = '2-m air temperature'

    # Global attributes
    ncfile.Title             = 'ICON LES FESSTVaLD2 simulation interpolated to Cartesian grid'
    ncfile.Institution       = 'Meteorological Institute, University of Hamburg (UHH), Germany'
    ncfile.Contact_person    = 'Bastian Kirsch (bastian.kirsch@uni-hamburg.de)'
    ncfile.Source            = 'Simulation data located on Levante at /pool/data/fesstval/SIM/ICON_LES/'
    ncfile.History           = 'None'
    ncfile.Conventions       = 'CF-1.7 where applicable'
    ncfile.Processing_date   = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')    
    ncfile.Author            = 'Bastian Kirsch (bastian.kirsch@uni-hamburg.de)'     
    ncfile.Comments          = 'None'
    ncfile.Licence           = 'This data is licensed under a '+\
                               'Creative Commons Attribution 4.0 '+\
                               'International License (CC BY 4.0).'  
    ncfile.close()
    print('Writing successful')
    return


def plot_les(lon_data,lat_data,data,dtime,var_str,domain,network_stat,
             vmin=None,vmax=None,plotdir=maindir+'icon/plots/'):
    
    print('Plotting '+dtime.strftime('%H:%M'))
    
    var_names = {'TT':'2-m temperature (°C)',
                } 
    titles = {1:'ICON-LES DOM01 (625 m)',
              2:'ICON-LES DOM02 (312 m)',
              3:'ICON-LES DOM03 (156 m)',
              4:'ICON-LES DOM04 (75 m)',
              }
    
    if domain == 1: lon_res,lat_res = 0.005,0.005
    if domain == 2: lon_res,lat_res = 0.002,0.002
    if domain == 3: lon_res,lat_res = 0.001,0.001
    if domain == 4: lon_res,lat_res = 0.001,0.001
    
    # Plot area and interpolation grid
    lon_min,lon_max = lon_data.min(),lon_data.max()
    lat_min,lat_max = lat_data.min(),lat_data.max()
        
    lon_range = np.arange(lon_min,lon_max+lon_res,lon_res,dtype=float)  
    lat_range = np.arange(lat_min,lat_max+lat_res,lat_res,dtype=float) 
    lon_mg,lat_mg = np.meshgrid(lon_range,lat_range)
    
    data_regridded = interpolate_data(lon_data,lat_data,data,lon_mg,lat_mg)
    
    if not vmin: vmin = np.nanmin(data_regridded)
    if not vmax: vmax = np.nanmax(data_regridded)
    
    cmap   = mpl.colormaps['plasma']
    #norm   = mpl.colors.Normalize(vmin,vmax,clip=True)
    levels = np.linspace(vmin,vmax,16)
    norm   = mpl.colors.BoundaryNorm(levels, cmap.N, clip=True)
    
    fig,ax = plt.subplots(1,1,figsize=(5.6,4.5),dpi=200)
    
    # Interpolated data
    pax = ax.pcolormesh(lon_mg,lat_mg,data_regridded,norm=norm,cmap=cmap)
    
    # Station points
    ax.scatter(network_stat['LON'],network_stat['LAT'],c='black',marker='o',s=5) 
    
    ax.set_xlim([lon_min,lon_max])    
    ax.set_ylim([lat_min,lat_max])
    ax.set_xlabel('Lon (°E)')
    ax.set_ylabel('Lat (°N)') 

    # Title and colorbar
    ax.set_title(titles[domain],loc='left',fontsize=10)
    ax.set_title(dtime.strftime('%Y-%m-%d %H:%M UTC'),loc='right',fontsize=10)
    plt.colorbar(pax,label=var_names[var_str])
    
    plotname = plotdir+'dom0'+str(domain)+'/'+var_str+'_icon_dom0'+str(domain)+dtime.strftime('_%Y%m%d_%H%M.png')
    fig.savefig(plotname,bbox_inches='tight')
    plt.close()    
    print('Plot done!')  
    return