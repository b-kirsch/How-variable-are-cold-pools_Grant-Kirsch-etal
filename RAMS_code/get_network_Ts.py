import h5py
import numpy as np
import datetime
import pandas as pd
import os
from scipy.interpolate import griddata


def get_network_ts(ramsfolder,startt,num_sim_times,sim_output_freq,file_format,network_cent_lat,network_cent_lon,main_dir,network_file):

    if not os.path.isdir(main_dir):
        os.mkdir(main_dir)
    network_info_sname = main_dir+'network_info.h5'
    network_Ts_sfolder = main_dir+'network_Ts/'
    if not os.path.isdir(network_Ts_sfolder):
        os.mkdir(network_Ts_sfolder)

    network_dataframe = pd.read_csv(network_file,delimiter=';')
    network_array = np.array(network_dataframe)
    network_delta_x = network_array[:,1].astype(float)
    network_delta_y = network_array[:,2].astype(float)

    network_delta_lat = network_delta_y/110574.0
    network_lat = network_cent_lat+network_delta_lat
    network_delta_lon = network_delta_x/(111320.0*np.cos(np.deg2rad(network_lat)))
    network_lon = network_cent_lon+network_delta_lon

    saved_file = h5py.File(network_info_sname, 'w')
    saved_file.create_dataset('network_lat',data=network_lat)
    saved_file.create_dataset('network_lon',data=network_lon)
    saved_file.create_dataset('network_delta_x',data=network_delta_x)
    saved_file.create_dataset('network_delta_y',data=network_delta_y)
    saved_file.create_dataset('network_cent_lat',data=network_cent_lat)
    saved_file.create_dataset('network_cent_lon',data=network_cent_lon)
    saved_file.close()

    all_times = pd.date_range(startt, periods=num_sim_times, freq=sim_output_freq)
    for currt in all_times:

        ct_str = currt.strftime("%Y-%m-%d-%H%M%S")
        rams_file = h5py.File(ramsfolder+file_format.format(ct_str), 'r')
        lats = np.array(rams_file['GLAT'])
        lons = np.array(rams_file['GLON'])
        theta = np.array(rams_file['THETA'][1,:,:])
        pi = np.array(rams_file['PI'][1,:,:])/1004
        rams_file.close()
        temperature = theta*pi

        network_Ts = griddata((lats.flatten(),lons.flatten()),temperature.flatten(),
                        (network_lat,network_lon),method='linear')
        
        saved_file = h5py.File(network_Ts_sfolder+ct_str+'.h5', 'w')
        saved_file.create_dataset('network_Ts',data=network_Ts)
        saved_file.close()

        print(currt)

'''
ramsfolder = '/sumatra/nfalk/CAMP2EX_forecast-like_simulations/SF16/'
startt = datetime.datetime(2019,9,28,22,0)
num_sim_times = 49
sim_output_freq = '15min'
file_format = "a-A-{0}-g1.h5"
network_cent_lat = 18
network_cent_lon = 119
main_dir = './SF16/'
network_file = './network_coordinates.txt'
'''


ramsfolder = '/sumatra/nfalk/HR/HR_29_600_T12/'
startt = datetime.datetime(2019,9,29,0,0)
num_sim_times = 121
sim_output_freq = '5min'
file_format = "a-L-{0}-g1.h5"
network_cent_lat = 18
network_cent_lon = 119
main_dir = './HRT12/'
network_file = './network_coordinates.txt'


get_network_ts(ramsfolder,startt,num_sim_times,sim_output_freq,file_format,
    network_cent_lat,network_cent_lon,main_dir,network_file)