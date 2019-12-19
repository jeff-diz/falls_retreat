# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:25:32 2019

@author: Jeff Disbrow
"""

import logging
import os
import matplotlib.pyplot as plt
import numpy as np

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, LineString, MultiLineString
from shapely.ops import split, linemerge
import sympy


def create_points(row, lat, lon):
    """
    Create point from longitude and latitude.
    """
    geom = Point(row[lon], row[lat])
    
    return geom
    

def create_plot_size(row, RIVER_TYPE, RIVER):
    """
    Creates plotting sizes for rivers, with main (river) channels being bigger.
    """
    if row[RIVER_TYPE] == RIVER:
        plot_size = 5
    else:
        plot_size = 2
    
    return plot_size


def snap2closest_line(pts, lines, tolerance, line_cols=None):
    """
    Snaps each point to the closest line.
    pts:      (geodataframe) points to move
    lines:    (geodataframe) lines to snap to
    tolerance (int)          maximum distance 
                             to move a point - in units of projection
    line_cols (list)         list of columns from lines gdf to keep
    Returns
    geodataframe: original points geodataframe with new, snapped geometry
    """
    #### Initial search for close lines
    ## Create bounding box for spatial index
    bbox = pts.bounds + [-tolerance, -tolerance, tolerance, tolerance]
    ## Find river lines that intersect bboxes
    hits = bbox.apply(lambda row: list(lines.sindex.intersection(row)), axis=1)
    ## Collect all intersections (pot. multiple per point, and respective indicies)
    tmp = pd.DataFrame({"pt_idx": np.repeat(hits.index, hits.apply(len)),
                        "line_i": np.concatenate(hits.values)})
    ## Join back to lines, then to corresponding points
    tmp = tmp.join(lines.reset_index(drop=True), on="line_i")
    tmp = tmp.join(pts.geometry.rename("point"), on="pt_idx")
    ## Convert back to geodataframe
    tmp = gpd.GeoDataFrame(tmp, geometry="geometry", crs=pts.crs)
    
    #### Determine closest line
    ## Find distance to bbox intersecting line(s)
    tmp['snap_dist'] = tmp.geometry.distance(gpd.GeoSeries(tmp['point']))
    ## Discard lines further away than tolerance
    tmp = tmp.loc[tmp['snap_dist'] <= tolerance]
    ## Sort, so closest is first, then keep only first match for each point
    tmp = tmp.sort_values(by=["snap_dist"])
    closest = tmp.groupby("pt_idx").first()
    ## Geodataframe of only closest lines
    closest = gpd.GeoDataFrame(closest, geometry='geometry')
    
    #### Find distance along line to nearest point on line
    pos = closest.geometry.project(gpd.GeoSeries(closest['point']))
    ## Find new point location - distance along closest line
    new_pts = []
    for i, row in closest.iterrows():
        new_pts.append(row.geometry.interpolate(pos.iloc[i]))
    
    ### Copy back to original pts dataframe, using snapped locations
    ## Line cols to keep
    if line_cols == None:
        line_cols = []
    elif line_cols == 'all':
        line_cols = list(lines)
    ## Geodataframe of relevant line columns and new snapped geometries
    snapped = gpd.GeoDataFrame(closest[line_cols], geometry=new_pts)
    ## Join to original pts, switching to snapped geometry
    pts_snapped = pts.drop(columns=['geometry']).join(snapped)
    ## Drop points that were farther than tolerance
    pts_snapped = pts_snapped.dropna(subset=['geometry'])
    
    return pts_snapped


def split_multilinestring(MultiLineString_object):
    """
    Split a MultiLineString into it's parts.

    Parameters
    ----------
    MultiLineString_object : shapely.geometry.MultiLineString
        Shapely MultiLineString object.

    Returns
    -------
    lines : list
        of shapely LineStrings and/or MultiLineStrings.

    """
    lines = []
    for line in MultiLineString_object:
        lines.append(line)
    return lines


def iterative_split_multilinestring(MultiLineString_object):
    """
    Iteratively split a MultiLineString object into a list of LineString 
    objects and/or MultiLineString objects. If result includes MultiLineString
    splits it, and returns it LineString parts until only LineStrings remain.
    MultiLineString.    

    Parameters
    ----------
    MultiLineString_object : shapely.geometry.MultiLineString
        shapley MultiLineStringObject.

    Returns
    -------
    lines : list
        of shapely LineString onjects.
    """
    multi_present = True
    while multi_present == True:
        lines = split_multilinestring(MultiLineString_object)
        geom_types = [x.geom_type for x in lines]
        if 'MultiLineString' in geom_types:
            multi_present = True
        else:
            multi_present = False
    
    return lines


def find_drainage_by_loc(pt, drainages):
    """
    Performs point in polygon to identify appropriate drainage
    """
    selected_drainage = [x for x in drainages.geometry if pt.within(x)][0]

    return selected_drainage


def create_splitter(x, y, frac=0.001):
    # create 'splitters' - lines that extend above and below point to ensure intersection with line
    splitter = LineString([(x, y-y*frac), (x, y), (x, y+y*frac)])

    return splitter


def stream_power_m(R, A, K=1.0):
    """
    Solve the stream power equation for coefficient
    setting the power law scaling between retreat
    rate and drainage area.

    Parameters
    ----------
    R : float
        Retreat rate (mm / yr).
    A : float
        Area (km^2).
    K : float
        bedrock erodibility coefficient

    Returns
    -------
    float
        m: coeffiicient setting power law scaling
           between retreat rate and drainage area

    """
    m = sympy.Symbol('m')
    calculated_m = sympy.solve(K*(A**m)-R, m)
    if len(calculated_m) == 1:
        calculated_m = float(calculated_m[0])
    else:
        calculated_m = None
        
    return calculated_m


#### PLOTTING FUNCTIONS
def y_fmt(y, pos):
    '''
    Formatter for y axis of plots. Returns the number with appropriate suffix
    y: value
    pos: *Not needed?
    '''
    decades = [1e9, 1e6, 1e3, 1e0, 1e-3, 1e-6, 1e-9 ]
    suffix  = ["G", "M", "k", "" , "m" , "u", "n"  ]
    if y == 0:
        return str(0)
    for i, d in enumerate(decades):
        if np.abs(y) >=d:
            val = y/float(d)
            signf = len(str(val).split(".")[1])
            if signf == 0:
                return '{val:d} {suffix}'.format(val=int(val), suffix=suffix[i])
            else:
                if signf == 1:
                    if str(val).split(".")[1] == "0":
                       return '{val:d}{suffix}'.format(val=int(round(val)), suffix=suffix[i]) 
                tx = "{"+"val:.{signf}f".format(signf = signf) +"} {suffix}"
                return tx.format(val=val, suffix=suffix[i])
    return y