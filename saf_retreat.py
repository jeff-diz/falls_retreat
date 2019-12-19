# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 13:55:48 2019

@author: Jeff Disbrow
"""

import matplotlib.pyplot as plt
import numpy as np

import geopandas as gpd
import fiona
from shapely.geometry import Point, shape, MultiLineString, LineString
from shapely.ops import split, linemerge

from lib import snap2closest_line, create_splitter

def saf_retreat_date(interp_point, write_river_section=None):
    """
    Returns the date st anthony falls passed interp_point
    """
    utm_15n = {'init':'epsg:26915'}
    wgs84   = {'init':'epsg:4326'}
    
    #### Read in Mississippi River
    mr_path = r'C:\Users\Jeff Disbrow\Documents\UMN-Drive\Fall19\ESCI8701\project\miss_river.shp'
    mr = shape(next(iter(fiona.open(mr_path)))['geometry'])
    mr = linemerge(mr)
    mr = MultiLineString([line for line in mr])
    mr = linemerge(mr)
    
    #### Split MR at SAF, Hidden Falls, St. Paul
    ## Point locations
    saf = Point(-93.2476297, 44.9791802)
    hf  = Point(-93.1958309, 44.9063822)
    stp = Point(-93.07419, 44.94705)
    points = gpd.GeoDataFrame({'loc': ['saf', 'hf', 'stp']}, geometry=[saf, hf, stp], crs=wgs84)
    points = points.to_crs(utm_15n)
    
    mr_df = gpd.GeoDataFrame(geometry=[l for l in mr])
    
    # Snap points to MR for splitting, then get back lines as shapely geometries
    mr_df = gpd.GeoDataFrame(geometry=[mr])
    pts_snapped = snap2closest_line(points, mr_df, tolerance=10000)
    saf = pts_snapped[pts_snapped['loc']=='saf'].geometry.values[0]
    hf  = pts_snapped[pts_snapped['loc']=='hf'].geometry.values[0]
    stp = pts_snapped[pts_snapped['loc']=='stp'].geometry.values[0]
    
    ## Split
    # Select only relevent section of MR (one LineString from MultiLineString)
    mr = [l for l in mr if l.distance(hf)<100][0]
    # Create splitters
    hf_splitter = create_splitter(hf.x, hf.y, frac=0.0005)
    saf_splitter = create_splitter(saf.x, saf.y, frac=0.0005)
    stp_splitter = create_splitter(stp.x, stp.y, frac=0.0005)
    # Split to get just sections: stp->hf->saf
    below_saf, _above_saf = split(mr, saf_splitter)
    _below_stp, above_stp = split(below_saf, stp_splitter)
    below_hf, above_hf = split(above_stp, hf_splitter)
    
    if write_river_section is not None:
        river_section = gpd.GeoDataFrame(geometry=[below_hf, above_hf], crs=utm_15n)
        river_section.to_file(write_river_section)
    
    #### Interpolate position
    ## 12 -> 13.4 below HF -> stp
    ## direction starts at stp
    len_below = below_hf.length
    dist_arr_below = [0, len_below]
    time_arr_below = [13.4, 12]
    
    ## 12 -> 0 above HF -> SAF
    ## direction starts at hf
    len_above = above_hf.length
    dist_arr_above = [0, len_above]
    time_arr_above = [12, 0]
    
    ## Determine if iterp point is above or below hf
    if above_hf.distance(interp_point) < below_hf.distance(interp_point):
        date = np.interp(above_hf.project(interp_point), dist_arr_above, time_arr_above)
    else:
        date = np.interp(below_hf.project(interp_point), dist_arr_below, time_arr_below)
    
    return date

### Debug plotting
#fig, ax = plt.subplots(1,1,figsize=(8,8))
#splitter = gpd.GeoDataFrame(geometry=[hf_splitter, saf_splitter, stp_splitter])
#splits = gpd.GeoDataFrame(geometry=[above_hf])
##tmr = gpd.GeoDataFrame(geometry=[r for r in mr], crs=utm_15n)
#mr_df['rand'] = np.random.randint(1, 6, mr_df.shape[0])
#splits['rand'] = np.random.randint(1, 6, splits.shape[0])
#
##mr_df.plot(column='rand', ax=ax)
#splitter.plot(ax=ax)
#splits.plot(column='rand', ax=ax)
#points.plot(ax=ax)
#test = above_hf.interpolate(100)
#t = gpd.GeoDataFrame(geometry=[test])
#t.plot(ax=ax)
##pts_snapped.plot(ax=ax, color='blue')