#!/usr/bin/env python3
# coding: utf-8
"""
Created on Mon Jul  6 12:33:11 2020

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
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import os.path
import os
import sys

import matplotlib.patches
import shapely.geometry
import cartopy.io.shapereader as shapereader
import pandas as pd
from matplotlib.lines import Line2D

# import cartopy.crs as ccrs
# from descartes import PolygonPatch

def create_map():
    lat_min = 10
    lat_max = 60
    step_lat = 10
    lon_min = -30
    lon_max = 60
    step_lon = 10
    
    m = Basemap(llcrnrlat=lat_min,urcrnrlat=lat_max,\
    	llcrnrlon=lon_min,urcrnrlon=lon_max, resolution='l')
    m.drawparallels(range(lat_min, lat_max, step_lat),labels=[1,0,0,1],color='grey',linewidth=0.1)
    m.drawmeridians(range(lon_min, lon_max, step_lon),labels=[1,0,0,1],color='grey',linewidth=0.1)
    m.drawcoastlines(linewidth=0.1)
    m.fillcontinents(color='grey')
    return m

def draw_polygon(lats,lons,m,sensor):
    x, y = m( lons, lats )
    xy = [(x[0],y[0]),(x[1],y[1]),(x[2],y[2]),(x[3],y[3])]
    
    if sensor == 'S3A':
        fc = 'red'
    elif sensor == 'S3B':
        fc = 'blue'
    poly = matplotlib.patches.Polygon( xy, facecolor=fc, alpha=0.2,closed=True, ec='k', lw=1)
    plt.gca().add_patch(poly)
    poly_all = shapely.geometry.Polygon(xy)
    return m, poly_all



    
def create_MED_shp(path_to_shp):
    #find the MED polygon. 
    for sea in shapereader.Reader(path_to_shp).records(): 
        if sea.attributes['NAME']=='Mediterranean Sea - Western Basin': 
            MedW = sea.geometry 
        if sea.attributes['NAME']=='Mediterranean Sea - Eastern Basin': 
            MedE = sea.geometry   
        if sea.attributes['NAME']=='Adriatic Sea': 
            Adri = sea.geometry
        if sea.attributes['NAME']=='Aegean Sea': 
            Aege = sea.geometry
        if sea.attributes['NAME']=='Ligurian Sea': 
            Ligu = sea.geometry   
        if sea.attributes['NAME']=='Ionian Sea': 
            Ioni = sea.geometry
        if sea.attributes['NAME']=='Tyrrhenian Sea': 
            Tyrr = sea.geometry 
        if sea.attributes['NAME']=='Balearic Sea': 
            Bale = sea.geometry 
        if sea.attributes['NAME']=='Alboran Sea': 
            Albo = sea.geometry  
      
    Med = MedE.union(MedW).union(Adri).union(Aege).union(Ligu).union(Ioni).union(Tyrr).union(Bale).union(Albo) 
    return Med           
        
def coverage(path_out,df_doy,df_coverage):
    # create empty polygons per each sensor
    S3Apoly = shapely.geometry.Polygon()
    S3Bpoly = shapely.geometry.Polygon()
    m = create_map()
    
    for row in df_doy.itertuples(index=True, name='Pandas'):
        # extract data from dataframe
        UL_lat = row.UL_lat
        UL_lon = row.UL_lon
        UR_lat = row.UR_lat
        UR_lon = row.UR_lon
        LL_lat = row.LL_lat
        LL_lon = row.LL_lon
        LR_lat = row.LR_lat
        LR_lon = row.LR_lon    
        sensor = row.sensor
        date = row.date   

        # draw map
        lats_poly =[UL_lat,UR_lat,LR_lat,LL_lat]
        lons_poly =[UL_lon,UR_lon,LR_lon,LL_lon]
        m, poly = draw_polygon(lats_poly,lons_poly,m,sensor)
        if sensor == 'S3A':
            S3Apoly = S3Apoly.union(poly)
        elif sensor == 'S3B':
            S3Bpoly = S3Bpoly.union(poly)
        plt.gcf()
        custom_lines = [Line2D([0], [0], color='red', lw=4),
            Line2D([0], [0], color='blue', lw=4)]

        plt.legend(custom_lines, ['S3A', 'S3B'],loc='upper left')
            
    # calculate coverage percentage
    AB_union_perc, AB_inter_perc = coverage_calc(date,S3Apoly, S3Bpoly)
    
    # save figure
    plt.gcf()
    ofname = os.path.join(path_out,str(date)+'_coverage.pdf')
    plt.savefig(ofname, dpi=200)
    
    plt.close()
    
    row_percentage = {
        'date': str(date),
        'A union B': AB_union_perc,
        'A interc B': AB_inter_perc
        }
    
    df_coverage = df_coverage.append(row_percentage,ignore_index=True) 
        
    return df_coverage,S3Apoly, S3Bpoly

def coverage_calc(date,S3Apoly, S3Bpoly):
    # #create figure
    # plt.figure(figsize=(10,10))  
    # PLT = plt.axes(projection=ccrs.PlateCarree())
    # PLT.set_extent([-30,60,20,50])
    # PLT.gridlines()
    
    #import and display shapefile
    path_to_shp = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/Call_ESA_MED/src/World_Seas_IHO_v1/World_Seas.shp'
    
    # PLT.add_patch(PolygonPatch(S3Apoly,  fc='pink', ec='#555555', alpha=0.1, zorder=5))
    # PLT.add_patch(PolygonPatch(S3Bpoly,  fc='cyan', ec='#555555', alpha=0.1, zorder=5))
    
    Med = create_MED_shp(path_to_shp)
    
    #calculate coverage 
    AB_union = S3Apoly.union(S3Bpoly).intersection(Med)
    AB_inter = S3Apoly.intersection(S3Bpoly).intersection(Med)
    # PLT.add_patch(PolygonPatch(AB_union, fc='black', alpha=1)) 
    # PLT.add_patch(PolygonPatch(AB_inter, fc='red', alpha=1)) 
    
    AB_union_perc = AB_union.area/Med.area*100
    AB_inter_perc = AB_inter.area/Med.area*100
    
    plt.gcf()
    plt.title(f'{date}; A$\cup$B={AB_union_perc:.2f}%; A$\cap$B={AB_inter_perc:.2f}%')
    
    print(f'Coverage A+B: {AB_union_perc:.2f}%')
    print(f'Coverage A intersection B: {AB_inter_perc:.2f}%')
    
    return AB_union_perc, AB_inter_perc       
#%%
# def main():
if sys.platform == 'linux': 
    path_source = '/dst04-data1/OC/OLCI/trimmed_sources'
    path_out = '/home/Javier.Concha/OLCI_coverage/Figures'  
elif sys.platform == 'darwin':
    path_source = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/Images/OLCI/trimmed_sources'
    path_out = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/Call_ESA_MED/Figures'
else:
    print('Error: host flag is not either mac or vm')   

    # create Dataframe    
cols = ['date','A union B','A interc B']
df_coverage = pd.DataFrame(columns = cols) 
  
csv_filename = os.path.join(path_out,'OLCI_geo_info_ALL.csv')
df = pd.read_csv(csv_filename)

year_str = '2020'
sdoy = 148
edoy = 183

for doy in range(sdoy,edoy+1):
    
    
    if not df.loc[df['doy'] == doy].empty:
        df_doy = df.loc[df['doy'] == doy]
        print('-----------')
        print(df_doy)
        df_coverage, S3Apoly, S3Bpoly = coverage(path_out,df_doy,df_coverage)  

csv_coverage_filaname = os.path.join(path_out,'OLCI_coverage_info.csv')
df_coverage.round({'A union B':1,'A interc B':1}).to_csv(csv_coverage_filaname)
#%%
# if __name__ == '__main__':
#     main()
