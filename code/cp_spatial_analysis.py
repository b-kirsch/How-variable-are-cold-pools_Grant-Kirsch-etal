# -*- coding: utf-8 -*-
"""
@author: Bastian Kirsch (bastian.kirsch@uni-hamburg.de)

Spatial analysis of cold pools based on scattered point observations,
including kriging interpolation and object clustering

Dependences on non-standard software:
- fesstval_routines.py

Last updated: 16 May 2023
"""

import numpy as np
import pandas as pd
import ccl
from scipy import ndimage as ndi
import skimage.segmentation as sks
import skimage.feature as skf
import pykrige as pk
import fesstval_routines as fst


class network_analysis:
    
    def __init__(self,lon_data,lat_data,
                 lon_data_ref=None,lat_data_ref=None,
                 x_res=None,x_min=None,x_max=None,
                 y_res=None,y_min=None,y_max=None,         
                 hav=False):
        
        if lon_data.shape != lat_data.shape:
            print('Lon data and Lat data do not have the same shape!')
        
        if isinstance(lon_data,pd.Series): lon_data = lon_data.to_numpy(float)
        if isinstance(lat_data,pd.Series): lat_data = lat_data.to_numpy(float)
        
        # Convert Lon/Lat data into xy-coordinates
        lon_ref = np.nanmean(lon_data) if not lon_data_ref else lon_data_ref
        lat_ref = np.nanmean(lat_data) if not lat_data_ref else lat_data_ref
        n_points = len(lon_data)
        
        x_data,y_data = fst.lonlat_to_xy(lon_data,lat_data,
                                        lon_ref=lon_ref,lat_ref=lat_ref,
                                        haversine_method=hav)
        
        # Set up xy grid
        if not x_res or not y_res:
            x_grid = np.array([])
            y_grid = np.array([])
            valid_grid = False
        else:
            x_res_dim = len(str(int(x_res)))   
            y_res_dim = len(str(int(y_res)))
            
            if not x_min: x_min = np.round(x_data.min(),-x_res_dim+1)-x_res
            if not x_max: x_max = np.round(x_data.max(),-x_res_dim+1)+x_res
            if not y_min: y_min = np.round(y_data.min(),-y_res_dim+1)-y_res
            if not y_max: y_max = np.round(y_data.max(),-y_res_dim+1)+y_res
            
            x_grid = np.arange(x_min,x_max+x_res,x_res,dtype=float)
            y_grid = np.arange(y_min,y_max+y_res,y_res,dtype=float)
            valid_grid = True
            
        self.lon_ref    = lon_ref
        self.lat_ref    = lat_ref
        self.n_points   = n_points
        self.x_res      = x_res
        self.y_res      = y_res
        self.x_data     = x_data
        self.y_data     = y_data
        self.x_grid     = x_grid
        self.y_grid     = y_grid
        self.valid_grid = valid_grid
        self.hav        = hav
        
        
    def lonlat_ref(self):
        return self.lon_ref,self.lat_ref
    
    def xy_data(self):
        return self.x_data,self.y_data
    
    def xy_gridpoints(self):
        return self.x_grid,self.y_grid
    
    def xy_meshgrid(self):
        return np.meshgrid(self.x_grid,self.y_grid)
    
    def lonlat_meshgrid(self,shift=False):
        x_meshgrid,y_meshgrid = self.xy_meshgrid()
        # Shift grid by half grid spacing for correct plotting
        if shift:
            x_meshgrid -= 0.5*self.x_res
            y_meshgrid -= 0.5*self.y_res
        if len(x_meshgrid) != 0:
            lon_meshgrid,lat_meshgrid = fst.xy_to_lonlat(x_meshgrid,y_meshgrid,
                                                         self.lon_ref,self.lat_ref,
                                                         meshgrid=True,
                                                         haversine_method=self.hav)
            return lon_meshgrid,lat_meshgrid
        else:
            return x_meshgrid,y_meshgrid
    
    # Mask all grid point more than given distance away from any station    
    def mask_meshgrid(self,mask_dist,ii_filter=[]):
        x_meshgrid,y_meshgrid = self.xy_meshgrid()
        dist_mask_meshgrid    = np.zeros_like(x_meshgrid,dtype=bool)
        
        # Define which stations should be considered as valid
        if len(ii_filter) == 0:
            ii_valid = range(self.n_points)
        else:    
            ii_valid = np.arange(self.n_points)[ii_filter].astype(int) 
            
        for i in ii_valid:
            dist = np.sqrt(np.add(np.power(np.subtract(self.x_data[i],x_meshgrid),2),
                                  np.power(np.subtract(self.y_data[i],y_meshgrid),2)))
            dist_mask_meshgrid[dist <= mask_dist] = True
            
        return dist_mask_meshgrid

    def interpolation(self,input_data,x_input=[],y_input=[],x_eval=[],y_eval=[],
                      krig_type='ordinary',sv_model='power',sv_lags=10,sv_plot=False,
                      var_lim=None,weight=False,drift=False,return_var=False,
                      return_sv_params=False):
        
        if not self.valid_grid:
            if (len(x_eval) == 0) or (len(y_eval) == 0):
                print('No defined grid or evaluation points available!')
                return self.xy_meshgrid()[0]
        
        if isinstance(input_data,pd.DataFrame) and input_data.shape[0] == 1:
            input_data = input_data.squeeze()
         
        # Line added 12 Aug 2021 --> still working for 2020?    
        if isinstance(input_data,pd.Series): input_data = input_data.to_numpy(float)
            
        x_input = self.x_data if len(x_input) == 0 else np.array(x_input,dtype=np.float64)
        y_input = self.y_data if len(y_input) == 0 else np.array(y_input,dtype=np.float64)
        
        x_eval = self.x_grid if len(x_eval) == 0 else np.array(x_eval,dtype=np.float64)
        y_eval = self.y_grid if len(y_eval) == 0 else np.array(y_eval,dtype=np.float64)
        
        if (input_data.shape != x_input.shape) or (input_data.shape != y_input.shape):
            print('Input data does not have the same shape as coordinate data!')
            return self.xy_meshgrid()[0]*np.nan
        
        if (np.sum(np.isfinite(input_data)) == 0) or (np.std(input_data) == 0):
            print('No valid input data or variance of input data is zero, no kriging possible!')
            return self.xy_meshgrid()[0]*np.nan
        
        if isinstance(input_data,pd.Series): 
            input_data = input_data.to_numpy(float)
            
        fin = np.isfinite(input_data) 

        krig = pk.ok.OrdinaryKriging(x_input[fin],y_input[fin],input_data[fin],
                                     variogram_model=sv_model,verbose=False,
                                     enable_plotting=sv_plot,nlags=sv_lags,
                                     weight=weight)
        # drift_terms = ['regional_linear'] if drift else []
        # krig = pk.uk.UniversalKriging(x_input[fin],y_input[fin],input_data[fin],
        #                               variogram_model=sv_model,verbose=False,
        #                               enable_plotting=sv_plot,nlags=sv_lags,
        #                               weight=weight,drift_terms=drift_terms)
        krig_exec = krig.execute('grid',x_eval,y_eval)
        krig_data = np.ma.MaskedArray.filled(krig_exec[0])
        krig_var = np.ma.MaskedArray.filled(krig_exec[1])
        
        if var_lim != None:
            krig_data[krig_var > var_lim] = np.nan        
        if not return_var:  
            if not return_sv_params:
                return krig_data
            else:
                return krig_data,krig.variogram_model_parameters
        else:
            return krig_data,krig_var


class cluster_analysis:
    
    def __init__(self,x_meshgrid,y_meshgrid,grid_data,limit_value,
                 limit_bound='upper',method='ccl',connectivity=4,
                 peak_threshold=-2.5,peak_min_dist=5000,peak_number=3): 
        
        if x_meshgrid.shape != y_meshgrid.shape:
            print('x data and y data do not have the same shape!')
        
        if grid_data.shape != x_meshgrid.shape:
            print('Gridded data does not have the same shape as coordinate data!')
        
        if method not in ['ccl','watershed']:
            print('Clustering method '+str(method)+' not available, changed to ccl!')
            method = 'ccl'    
            
        if connectivity not in [4,8]:
            print('Connectivity '+str(connectivity)+' not available, changed to 4!')
            connectivity = 4
            
        if limit_bound not in ['lower','upper']:
            print('Limit bound '+limit_bound+' not available, changed to upper!')    
            limit_bound = 'upper'
        
        with np.errstate(invalid='ignore'):
            if limit_bound == 'upper': binary_mask = grid_data <= limit_value
            if limit_bound == 'lower': binary_mask = grid_data >= limit_value
            
        x_res = x_meshgrid[0,1] - x_meshgrid[0,0]
        y_res = y_meshgrid[1,0] - y_meshgrid[0,0]    
            
        if method == 'ccl':
            cluster_data = ccl.connected_component_labelling(binary_mask,connectivity)
            
        if method == 'watershed':
            min_dist   = int(np.ceil(peak_min_dist/x_res))
            peak_coord = skf.peak_local_max(np.abs(grid_data),
                                            threshold_abs=np.abs(peak_threshold),
                                            min_distance=min_dist,num_peaks=peak_number,
                                            num_peaks_per_label=1)
            peak_mask = np.zeros(grid_data.shape,dtype=bool)
            peak_mask[tuple(peak_coord.T)] = True
            markers, _ = ndi.label(peak_mask)
            cluster_data = sks.watershed(grid_data,markers,mask=binary_mask)    
        
        self.x_meshgrid   = x_meshgrid 
        self.y_meshgrid   = y_meshgrid
        self.x_res        = x_res
        self.y_res        = y_res
        self.grid_data    = grid_data
        self.cluster_data = cluster_data
        
    def data(self):
        return self.cluster_data

    def min_val(self,cluster_val):
        return np.nanmin(self.grid_data[self.cluster_data == cluster_val])
    
    def max_val(self,cluster_val):
        return np.nanmax(self.grid_data[self.cluster_data == cluster_val])
    
    def mean_val(self,cluster_val):
        return np.nanmean(self.grid_data[self.cluster_data == cluster_val])
    
    def std_val(self,cluster_val):
        return np.nanstd(self.grid_data[self.cluster_data == cluster_val])
    
    def perc_val(self,cluster_val,perc):
        # Perc between 0 and 100
        return np.nanpercentile(self.grid_data[self.cluster_data == cluster_val],perc)
    
    def area(self,cluster_val):
        pixel_area   = (self.x_res * self.y_res) * 1e-6 # (km2)
        n_pixel      = (self.cluster_data == cluster_val).sum()
        cluster_area = pixel_area * n_pixel
        return cluster_area # (km2)
    
    def center(self,cluster_val,return_lonlat=False,
               lon_ref=np.nan,lat_ref=np.nan,hav=False):
        # Center of mass
        temp_data = np.ones_like(self.cluster_data) * np.nan
        ii_cl = self.cluster_data == cluster_val
        temp_data[ii_cl] = np.abs(self.grid_data[ii_cl])
        
        x_center = np.nansum((self.x_meshgrid*temp_data/np.nansum(temp_data)))
        y_center = np.nansum((self.y_meshgrid*temp_data/np.nansum(temp_data)))
        
        if not return_lonlat:
            return x_center,y_center
        else:
            if np.isnan(lon_ref) or np.isnan(lat_ref):
                print('Reference point of grid is not defined!')
            lon_center,lat_center = fst.xy_to_lonlat(x_center,y_center,
                                                     lon_ref,lat_ref,
                                                     haversine_method=hav)
            return lon_center[0],lat_center[0]
    
    # Define boundary for given cluster        
    def boundary(self,cluster_val,cl_data,return_perimeter=False):  
        if cl_data.shape != self.x_meshgrid.shape:
            print('Input data does not have the same shape as grid!')
            return np.array([])
        
        bounds = np.zeros_like(self.x_meshgrid,dtype=bool)
        bounds[:,:] = False
        ny,nx = bounds.shape
        perimeter = 0 
        
        for x in range(nx):
            xmin = x-1 if x > 0 else 0
            xmax = x+2 if x < nx else nx+1
            for y in range(ny):
                if cl_data[y,x] != cluster_val: continue
                ymin = y-1 if y > 0 else 0
                ymax = y+2 if y < ny else ny+1
                
                if (np.sum(cl_data[ymin:ymax,xmin:xmax] != cluster_val) > 0) or\
                    (x in [0,nx-1]) or (y in [0,ny-1]):
                    bounds[y,x] = True    
                    perimeter += np.sum(cl_data[ymin:ymax,x] != cluster_val) +\
                                 np.sum(cl_data[y,xmin:xmax] != cluster_val)
                    if x in [0,nx-1]: perimeter += 1   
                    if y in [0,ny-1]: perimeter += 1
        if return_perimeter:            
            return bounds,perimeter 
        else:
            return bounds   
    