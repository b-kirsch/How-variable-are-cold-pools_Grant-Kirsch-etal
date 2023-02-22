# -*- coding: utf-8 -*-
"""
Calculation of empirical variogram for FESSTVaL station network

@author: Bastian Kirsch (bastian.kirsch@uni-hamburg.de)

Last updated: 22 February 2023
"""

import numpy as np
import pandas as pd
from itertools import combinations

# Distance bins for variogram
dist_min,dist_max,dist_res = 0,15000,500 # (m)

dist_bins = np.arange(dist_min,dist_max+dist_res,dist_res)
dist_bin_mids = dist_bins[:-1]+dist_res/2.

# Number of direction bins for variogram anisotropy
number_dir = 2 # 2, 4, or 8


def read_coordinates(filename='network_coordinates.txt'):
    '''
    Parameters
    ----------
    filename : str, optional
        Filepath of meta data file containing station coordinates.

    Returns
    -------
    meta : pandas.DataFrame 
        x and y coordinates (in m) of station locations with columns 
        named 'X' and 'Y' and station names as index
    '''
    meta = pd.read_csv(filename,sep=';',header=0)
    meta.set_index('NAME',inplace=True)
    return meta


def pair_stations(meta,bins=dist_bins,bin_labels=dist_bin_mids,n_dir=number_dir):
    '''
    Parameters
    ----------
    meta : pandas.DataFrame
        x and y coordinates (in m) of station locations with columns 
        named 'X' and 'Y' and station names as index
    bins : numpy.ndarray, optional
        Edges of distance bins (in m)
    bin_labels : numpy.ndarray, optional
        Labels of distance bins, usually mid points of bins (in m)
    n_dir : int, optional
        Number of equally wide direction bins used for anisotropy calculation

    Returns
    -------
    pairs : pandas.DataFrame
        Information on all station pairs (indices, names, distances, 
        distance bins)
    '''
    
    # Number of stations
    n = meta.shape[0]
    # Combinations running station index as pairs with columns named 'I0' and 'I1'
    pairs = pd.DataFrame(list(combinations(range(n),r=2)),columns=['I0','I1'])
    # Setting station names 'S0' and 'S1'
    pairs['S0'] = meta.index[pairs['I0']]
    pairs['S1'] = meta.index[pairs['I1']]
    # Reading x and y coordinates from meta
    x0,x1 = meta.loc[pairs['S0'],'X'],meta.loc[pairs['S1'],'X']
    y0,y1 = meta.loc[pairs['S0'],'Y'],meta.loc[pairs['S1'],'Y']
    # Calculating distance and direction for each station pair
    with np.errstate(invalid='ignore'):
        pairs['DIST'] = np.sqrt((x0.values-x1.values)**2+(y0.values-y1.values)**2)
        pairs['DIR']  = np.arctan((y0.values-y1.values)/(x0.values-x1.values)) * (180./np.pi)
    # Assigning distance bin to each station pair (if provided)
    if len(bins) > 0:
        if len(bin_labels) == 0: bin_labels = np.arange(bins.shape[0]-1)
        pairs['DIST_BIN'] = pd.cut(pairs['DIST'],bins,labels=bin_labels)
    # Assigning direction bin to each station pair (if provided)
    if n_dir in [2,4,8]:
        dir_min,dir_max = -90,90
        dir_binwidth = (dir_max-dir_min)/n_dir
        dir_bins = np.zeros(n_dir+2) * np.nan
        dir_bins[0],dir_bins[-1] = dir_min,dir_max
        dir_bins[1:-1] = np.arange(dir_min+dir_binwidth/2.,dir_max+dir_binwidth/2.,dir_binwidth)
        dir_bin_labels = np.arange(dir_min+dir_binwidth,dir_max+dir_binwidth,dir_binwidth)
        i_labels = np.arange(n_dir+1)
        dict_labels = dict(zip(i_labels,np.append(dir_bin_labels[-1],dir_bin_labels)))
        pairs['DIR_BIN'] = pd.cut(pairs['DIR'],dir_bins,labels=i_labels).replace(dict_labels)
          
    return pairs


def variogram(df):
    '''
    Parameters
    ----------
    df : pandas.DataFrame
        Data of station pairs (usually temperature) for given distance bin 
        with columns named 'T0' and 'T1'

    Returns
    -------
    float
        Variogram value (in K^2)
    '''
    nn = (df['T0'].notnull() & df['T1'].notnull()).sum()
    if nn > 0:
        # factor 1/(2*nn) for semi-variogram
        return (1/nn) * ((df['T0']-df['T1'])**2).sum()
    else:
        return np.nan


def gradiogram(df):
    '''
    Parameters
    ----------
    df : pandas.DataFrame
        Data of station pairs (usually temperature) for given distance bin 
        with columns named 'T0' and 'T1'

    Returns
    -------
    float
        Mean gradient (a.k.a. gradiogram) value (in K/km)
    '''
    nn = (df['T0'].notnull() & df['T1'].notnull()).sum()
    if nn > 0:
        return ((df['T0']-df['T1']).abs()/df['DIST']).mean() * 1000
    else:
        return np.nan      
      
    
def calc_variogram(data,pairs,func=variogram,anisotropy=False):
    '''
    Calculate variogram for single time step 
    
    Parameters
    ----------
    data : pandas.DataFrame
        Data of station network (usually temperature) with station names as 
        index (identical to index of meta)
    pairs : pandas.DataFrame
        Output of pair_stations
    func : function, optional
        Variogram or gradiogram 
    anisotropy : bool, optional
        Indicates if direction-dependent variogram is calculated to check
        for anisotropy

    Returns
    -------
    pandas.DataFrame
        Variogram or gradiogram with mid points of distance bins as index
    '''

    df_calc = pairs.copy(deep=True)
    df_calc['T0'] = data[pairs['S0']].values
    df_calc['T1'] = data[pairs['S1']].values
    
    if not anisotropy:
        return df_calc.groupby('DIST_BIN').apply(func)
    else:
        if 'DIR_BIN' not in pairs.columns:
            print('No direction bins defined! Set n_dir in pair_stations!')
            return
        dir_bins  = np.sort(pairs['DIR_BIN'].unique())
        dist_bins = np.sort(pairs['DIST_BIN'].unique().dropna())
        df_aniso  = pd.DataFrame(index=dir_bins,columns=dist_bins)
        for db in dir_bins:
            ii = df_calc['DIR_BIN'] == db
            df_aniso.loc[db] = df_calc[ii].groupby('DIST_BIN').apply(func)
        return df_aniso 
            
        
def calc_variogram_period(data,pairs,bin_labels=dist_bin_mids,func=variogram):
    '''
    Calculate variogram over given time period

    Parameters
    ----------
    data : pandas.DataFrame
        Data of station network (usually temperature) with time dimension as
        index and station names as columns (identical to index of meta)
    pairs : pandas.DataFrame
        Output of pair_stations
    bin_labels : numpy.ndarray, optional
        Labels of distance bins, usually mid points of bins (in m)
    func : function, optional
        Variogram or gradiogram 

    Returns
    -------
    df_variogram : pandas.DataFrame      
        Variogram or gradiogram with mid points of distance bins as index
    '''

    ntime = data.shape[0]
    df_variogram = pd.Series(index=bin_labels,dtype=float) 
    
    for d_bin in bin_labels:
        df_group = pd.DataFrame(columns=['T0','T1','DIST'])
        group = pairs['DIST_BIN'] == d_bin
        df_group['T0'] = data[pairs.loc[group,'S0']].melt()['value']
        df_group['T1'] = data[pairs.loc[group,'S1']].melt()['value']
        df_group['DIST'] = np.repeat(pairs.loc[group,'DIST'],ntime).values 
        
        df_variogram[d_bin] = func(df_group)
        
    return df_variogram
