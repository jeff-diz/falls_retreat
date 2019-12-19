# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 12:58:00 2019

@author: Jeff Disbrow
"""
import logging
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, LineString, MultiLineString
from shapely.ops import split, linemerge
from scipy.stats import linregress

from lib import create_points, create_plot_size, snap2closest_line, \
                iterative_split_multilinestring, find_drainage_by_loc, \
                y_fmt
from saf_retreat import saf_retreat_date


#### PROJECT DIRECTORIES AND PATHS
## Source directories
PRJ_DIR = r'C:\Users\Jeff Disbrow\Documents\UMN-Drive\Fall19\ESCI8701\project'
GIS_DIR = r'C:\GIS'
## Vector paths
FALLS_EXCEL = 'falls.xlsx'
FALLS_SHP = os.path.join(PRJ_DIR, 'falls.shp')
FALLS_SNAPPED_SHP = 'falls_snapped.shp'
RIVERS_PATH = os.path.join(PRJ_DIR, 'dnr_rivers_metro.shp')
## Drainage basins paths
DRAINAGE_PATH_GRASS = r'C:\GIS\esci8701\r_water_outlet\all_outlet_basins.gpkg' # GRASS delineated basins
DRAINAGE_PATH_SSTAT = r'C:\GIS\esci8701\streamstats\all_outletbasins_ss.gpkg'  # USGS streamstats basins
# DRAINAGE_PATH = DRAINAGE_PATH_GRASS # Choose which basins to use



## PARAMS -- string literals
# Excel column names
NAME = 'name'
LOC = 'loc'
PRES_LAT = 'present_lat'
PRES_LON = 'present_lon'
# Excel field values
LOC_PRES = 'present'
LOC_CONF = 'confluence'
LOC_MISS = 'miss_int'
# Shapefile field names
RIVER_TYPE = 'Strm_type_'
RIVER = 'Centerline (River)'
KITTLE_NAME = 'KITTLE_NAM'
KITTLE_NBR = 'KITTLE_NBR'
# Shapefile field values
MISS_RIVER = 'Mississippi River'
# Created column names
INCIS_YEARS = 'incis_date_yrs_bp'
RET_DIST = 'ret_distance'
RETREAT = 'retreat_rate'
DRAINAGE = 'drainage'
DRAIN_AREA = 'drainage_area_km'
PLOT_SIZE = 'plot_size'
RIVER_GEOM = 'river_geom'
# CRS codes
WGS84 = {'init':'epsg:4326'}


#### LOAD DATA
# Read excel with falls lat, lon, convert to points and geodataframe
falls_p = os.path.join(PRJ_DIR, FALLS_EXCEL)
falls = pd.read_excel(falls_p, header=0)
falls['geometry'] = falls.apply(lambda x: create_points(x, PRES_LAT, PRES_LON), axis=1)
falls = gpd.GeoDataFrame(falls, crs=WGS84)

# Read shapefile with creek and river lines
rivers = gpd.read_file(RIVERS_PATH)
# Dissolve on kittle name to only have 1 line per river or creek -- Better way?
rivers = rivers.dissolve(by=[KITTLE_NBR])
rivers.reset_index(inplace=True)
# Convert falls to UTM (crs of rivers)
falls = falls.to_crs(rivers.crs)
falls.to_file(FALLS_SHP)


#### FIND DISTANCES
## Snap falls to closest river line
falls_snapped = snap2closest_line(falls, rivers, 1000, line_cols='all')
## Create new dataframe with paired points on one row
paired_points = []
for name in falls_snapped[NAME].unique():
    present_loc = falls_snapped[(falls_snapped[NAME]==name) & (falls_snapped[LOC]==LOC_PRES)]['geometry'].values[0]
    conf_loc = falls_snapped[(falls_snapped[NAME]==name) & (falls_snapped[LOC]==LOC_CONF)]['geometry'].values[0]
    miss_int = falls_snapped[(falls_snapped[NAME]==name) & (falls_snapped[LOC]==LOC_MISS)]['geometry'].values[0]
    paired_points.append([name, present_loc, conf_loc, miss_int])
paired_df = pd.DataFrame(paired_points, columns=[NAME, LOC_PRES, LOC_CONF, LOC_MISS])

## Find only lines that points are on, to then merge and calc distance
matching_river_geoms = []
for i, pf in paired_df.iterrows():
    # Point geometries
    p1 = paired_df.iloc[i][LOC_PRES]
    p2 = paired_df.iloc[i][LOC_CONF]
    # Find all river lines that intersect current point
    pf_matches = []
    for row in rivers.itertuples():
        if row.geometry.distance(p1) < 0.001 or row.geometry.distance(p2) < 0.001:
            if row.geometry.geom_type == 'MultiLineString':
                lines = iterative_split_multilinestring(row.geometry)
                pf_matches.extend(lines)
            else:        
                pf_matches.append(row.geometry)
                
    if len(pf_matches) > 1:
        multiline = MultiLineString([x for x in pf_matches])
        matched_line = linemerge(multiline)
    else:
        matched_line = pf_matches[0]
    matching_river_geoms.append(matched_line)
paired_df[RIVER_GEOM] = matching_river_geoms

## Find distances between present and confluence points
distances = []
for i, pf in paired_df.iterrows():
    p1_proj = pf.river_geom.project(pf[LOC_PRES])
    p2_proj = pf.river_geom.project(pf[LOC_CONF])
    dist = p1_proj - p2_proj
    distances.append(abs(dist))
paired_df[RET_DIST] = distances


#### INTERPOLATE DATE SAF PASSED EACH CONFLUENCE (MISS RET)
paired_df[INCIS_YEARS] = paired_df[LOC_MISS].apply(saf_retreat_date)
paired_df[RETREAT] = paired_df[RET_DIST] / paired_df[INCIS_YEARS]


##### Identify watershed areas
## Loop drainage determination methods
for DRAINAGE_PATH in [DRAINAGE_PATH_GRASS, DRAINAGE_PATH_SSTAT]:
    ## Final spreadsheet path
    if DRAINAGE_PATH == DRAINAGE_PATH_GRASS:
        version = 'GRASS'
    elif DRAINAGE_PATH == DRAINAGE_PATH_SSTAT:
        version = 'stream_stat'
    RESULTS = 'falls_retreat_results_{}.xlsx'.format(version)
    
    ## Read in and rename goemetry column
    drainages = gpd.read_file(DRAINAGE_PATH)
    drainages.rename(columns={'geometry': DRAINAGE}, inplace=True)
    drainages.set_geometry(DRAINAGE, inplace=True)
    ## Convert crs
    drainages = drainages.to_crs(rivers.crs)
    
    falls_ct_before = len(paired_df)
    falls_names_before = list(paired_df.name)
    paired_df = paired_df.merge(drainages, how='inner', 
                                left_on='name', right_on='Name')
    # Check all drainages were found, if not print warning
    if len(paired_df) != falls_ct_before:
        print('Drainages were not found for all points')
        print('Missing falls:\n{}'.format([fall for fall in falls_names_before if fall not in list(paired_df.name)]))
    # This finds drainages by location
    # paired_df[DRAINAGE] = paired_df[LOC_PRES].apply(lambda x: find_drainage_by_loc(x, drainages))
    DRAIN_AREA = '{}_{}'.format(DRAIN_AREA, version)
    paired_df[DRAIN_AREA] = paired_df[DRAINAGE].apply(lambda x: x.area*1e-6)


#### Plotting
def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.annotate(s=str(point['val']), xy=(point['x'], point['y']), 
                    xytext=(6,-6),
                    textcoords='offset pixels')


plt.style.use('ggplot')
#### PLOTTING
## Scatter of drainage area vs retreat rate
fig, ax = plt.subplots(1,1, figsize=(8,8))
paired_df.plot.scatter(x=DRAIN_AREA, y=RETREAT, s=35, ax=ax)
lm = linregress(paired_df[DRAIN_AREA], paired_df[RETREAT])
ax.plot(paired_df[DRAIN_AREA], lm.slope*paired_df[DRAIN_AREA] +lm.intercept, color='red')
ax.set_xlabel('Drainage Area ($km^2$)')
ax.set_ylabel('Retreat Rate ($m / 1000 years$)')
# ax.set_title('Impact of Drainage Area on Waterfall Retreat Rate')
ax.xaxis.set_major_formatter(FuncFormatter(y_fmt))
ax.annotate(s="""Slope: {:.2f}\nIntercept: {:.2f}\n$R^2$: {:.2f}\np-val: {:.4f}\nstd. err: {:.2f}\n"""
                .format(lm.slope, 
                        lm.intercept,
                        lm.rvalue**2,
                        lm.pvalue,
                        lm.stderr),
             xy=(0.75, 0.05),
             xycoords='axes fraction',
             fontsize=10
             )
label_point(paired_df[DRAIN_AREA], 
            paired_df[RETREAT],
            paired_df[NAME], ax=ax)
plt.tight_layout()

## Map
# Plot falls on rivers
# map_fig, map_ax = plt.subplots(1,1, figsize=(10,10))
# drainages.plot(alpha=0.5, ax=map_ax, color='y', edgecolor='black')
# rivers[PLOT_SIZE] = rivers.apply(lambda x: create_plot_size(x, RIVER_TYPE, RIVER), axis=1)
# rivers['plot_color'] = np.random.randint(1, 1000, rivers.shape[0])
# rivers[rivers[RIVER_TYPE]!=RIVER].plot(ax=map_ax, color='#61c0ff', linewidth=rivers[PLOT_SIZE]) #column='plot_color',
# rivers[rivers[RIVER_TYPE]==RIVER].plot(ax=map_ax, color='blue', linewidth=rivers[PLOT_SIZE])

# falls.plot(ax=map_ax, color='g', markersize=20)
# falls_snapped.plot(ax=map_ax, color='r', markersize=30)

# plt.tight_layout()

#### Write created features to shapefiles - final dataframe to csv
# falls_snapped.to_file(os.path.join(PRJ_DIR, FALLS_SNAPPED_SHP))
# rivers.to_file(os.path.join(PRJ_DIR, 'rivers_dissolve.shp'))
# falls_final = gpd.GeoDataFrame(paired_df, geometry=LOC_PRES, crs=falls.crs)
# falls_final = falls_final.drop(columns=[LOC_CONF,
#                                         LOC_MISS,
#                                         RIVER_GEOM,
#                                         DRAINAGE])
# falls_final.to_file(FALLS_SHP)
drop_cols = ['present', 'confluence', 'miss_int', 'river_geom', 
             'layer', 'cat', 'Name', 'drainage', 'path']
drop_cols = [dc for dc in drop_cols if dc in list(paired_df)]
paired_df.drop(columns=drop_cols, inplace=True)
paired_df.rename(columns={'name':'Falls', 
                          'ret_distance':'Retreat Distance',
                          'incis_date_yrs_bp':'Time Incising (thousands years BP)',
                          'retreat_rate':'Retreat Rate',
                          'drainage_area_km':'Drainage Area (km^2)'},
                 inplace=True)
paired_df['Falls'] = paired_df['Falls'].apply(lambda x: x.title())
paired_df.round(decimals=2).to_excel(os.path.join(PRJ_DIR, RESULTS))