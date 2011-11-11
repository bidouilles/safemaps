#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
# Copyright (C) 2011  Lionel Bergeret
#
# ----------------------------------------------------------------
# The contents of this file are distributed under the CC0 license.
# See http://creativecommons.org/publicdomain/zero/1.0/
# ----------------------------------------------------------------

import os
import time
import cPickle
from optparse import OptionParser
import numpy as np
from matplotlib.mlab import griddata
from shapely.geometry import Point
from shapely.ops import cascaded_union
import safecastCommon
import csv

try:
   import psyco
   psyco.full()
except ImportError:
   print "Psyco plugin missing, will run slower"
   pass

# Japan limits
lat_min = 32.17023
lon_min = 136.7910
lat_max = 46.06
lon_max = 147.45

# -----------------------------------------------------------------------------
# Pre-compute the griddata interpolation
# -----------------------------------------------------------------------------
def PreCompute(safecastDataset, safecastDatasetCPMThreashold, safecastGridsize):
    # Setup: loading data...
    nx, ny = safecastGridsize, safecastGridsize # grid size
    npts, x, y, z, zraw = LoadSafecastData(safecastDataset, safecastDatasetCPMThreashold, lat_min, lat_max, lon_min, lon_max)

    # Compute area with missing data
    print "Compute area with missing data"
    measures = np.vstack((x,y)).T
    points = [Point(a,b) for a, b in measures]
    spots = [p.buffer(0.04) for p in points] # 0.04 degree ~ 1km radius
    # Perform a cascaded union of the polygon spots, dissolving them into a 
    # collection of polygon patches
    missing = cascaded_union(spots)

    # Create the grid
    print "Create the grid"
    xil = np.linspace(x.min(), x.max(), nx)
    yil = np.linspace(y.min(), y.max(), ny)
    xi, yi = np.meshgrid(xil, yil)

    # Calculate the griddata
    print "Calculate the griddata (%d x %d)" % (nx, ny)
    t1 = time.clock()
    zi = griddata(x,y,z,xi,yi,interp='nn')
    grid = zi.reshape((ny, nx))
    print "done in",time.clock()-t1,'seconds.'

    toSave = [npts, x, y, z, zraw, xil, yil, grid, missing]
    cPickle.dump(toSave,open('safecast.pickle','wb'),-1)
    print "Griddata saved (safecast.pickle)."

# -----------------------------------------------------------------------------
# Load the Safecast csv data file
# -----------------------------------------------------------------------------
def LoadSafecastData(filename, CPMclip, latmin, latmax, lonmin, lonmax):
    # Load data
    data = csv.reader(open(filename))

    # Read the column names from the first line of the file
    fields = data.next()

    x = []
    y = []
    z = []
    zraw = []
    
    # Process data
    for row in data:
      # Zip together the field names and values
      items = zip(fields, row)

      # Add the value to our dictionary
      item = {}
      for (name, value) in items:
         item[name] = value.strip()

      # Ignore if outside limits
      if not ((float(item["lat_avg"])>latmin) and (float(item["lat_avg"])<latmax) and (float(item["lon_avg"])>lonmin) and (float(item["lon_avg"])<lonmax)):
        continue

      cpm = float(item["cpm_avg"])
      zraw.append(cpm)
      
      if (cpm>CPMclip): cpm=CPMclip # clip

      x.append(float(item["lon_avg"]))
      y.append(float(item["lat_avg"]))
      z.append(float(cpm))

    npts = len(x)
    print "%s measurements loaded." % npts
   
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    zraw = np.array(zraw)
    
    return npts, x, y, z, zraw

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    parser = OptionParser("Usage: safecastPrecompute [options] <safecast-csv-file>")
    parser.add_option("-c", "--clip",
                      type="string", dest="clip", default="350",
                      help="specify the clipping value in CPM")
    parser.add_option("-g", "--gridsize",
                      type="string", dest="gridsize", default="1500",
                      help="specify the grid size for interpolation")

    (options, args) = parser.parse_args()
    
    if len(args) != 1:
        parser.error("Wrong number of arguments")

    PreCompute(args[0], int(options.clip), int(options.gridsize))
