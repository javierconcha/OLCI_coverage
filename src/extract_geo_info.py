#!/usr/bin/env python3
# coding: utf-8
"""
Created on Wed Jul  8 19:57:47 2020

@author: javier.concha
"""
"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from netCDF4 import Dataset
import pandas as pd
import subprocess
import sys
import os

# user defined functions
def create_list_products(path_source,path_out):
    cmd = f'find {path_source} -name "*OL_2_WFR*trim_MED*"> {path_out}/file_list.txt'
    prog = subprocess.Popen(cmd, shell=True,stderr=subprocess.PIPE)
    out, err = prog.communicate()
    if err:
        print(err)  
    path_to_list = f'{path_out}/file_list.txt'
    return path_to_list

def extract_geo_info(path_to_list,df):
    with open(path_to_list,'r') as file:
        for cnt, line in enumerate(file):   
            path_im = line[:-1]
            coordinates_filename = 'geo_coordinates.nc'
            filepah = os.path.join(path_im,coordinates_filename)
            nc_f0 = Dataset(filepah,'r')
            lat = nc_f0.variables['latitude'][:,:]
            lon = nc_f0.variables['longitude'][:,:]
            
            UL_lat = lat[0,0]
            UL_lon = lon[0,0]
            UR_lat = lat[0,-1]
            UR_lon = lon[0,-1]
            LL_lat = lat[-1,0]
            LL_lon = lon[-1,0]
            LR_lat = lat[-1,-1]
            LR_lon = lon[-1,-1]
            
            #%% create csv
            sensor = path_im.split('/')[-1].split('_')[0]
            datetimestr = path_im.split('/')[-1].split('_')[7]
            date =  datetimestr.split('T')[0]
            time = datetimestr.split('T')[1]
            doy = path_im.split('/')[-2]
            filename = path_im.split('/')[-1]
            
            granule = {
                    'sensor': sensor,
                    'datetimestr': datetimestr,
                    'date': date,
                    'time': time,
                    'doy': doy,
                    'UL_lat': UL_lat,
                    'UL_lon': UL_lon,
                    'UR_lat': UR_lat,
                    'UR_lon': UR_lon,
                    'LL_lat': LL_lat,
                    'LL_lon': LL_lon,
                    'LR_lat': LR_lat,
                    'LR_lon': LR_lon,
                    'filename': filename,
                    'filepah':   path_im
                    }
            
            df = df.append(granule,ignore_index=True) 
    
    return df
#%%
def main():
    if sys.platform == 'linux':
        path_source = '/dst04-data1/OC/OLCI/trimmed_sources'
        path_out = '/home/Javier.Concha/OLCI_coverage/Figures'  
    elif sys.platform == 'darwin':
        path_source = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/Images/OLCI/trimmed_sources'
        path_out = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/Call_ESA_MED/Figures'
    else:
        print('Error: host flag is not either mac or vm')
    
    # create Dataframe    
    cols = ['sensor','datetimestr','date','time','doy','UL_lat','UL_lon','UR_lat','UR_lon','LL_lat','LL_lon','LR_lat','LR_lon','filename','filepah']
    df = pd.DataFrame(columns = cols) 
    
    year_str = '2020'
    sdoy = 1
    edoy = 366
    
    for doy in range(sdoy,edoy+1):
        print('-----------')
        
        doy_str = str(doy)
        if doy<100:
            if doy<10:
                doy_str = '00'+ doy_str
            else:
                doy_str = '0'+ doy_str
        path_to_doy = os.path.join(path_source,year_str,doy_str)
        if os.path.isdir(path_to_doy):
            print(f'Extracting {year_str} {doy_str}')
            path_to_list = create_list_products(path_to_doy,path_out)
            df = extract_geo_info(path_to_list,df)
        else:
            print(f'Warning: Directory {year_str} {doy_str} does not exist.')
            
        
    # print(df[['sensor','datetimestr']])
    
    csv_filename = os.path.join(path_out,'OLCI_geo_info.csv')
    df.to_csv(csv_filename,index=False)
#%%
if __name__ == '__main__':
    main()

