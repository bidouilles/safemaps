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
import re
import cPickle
from optparse import OptionParser

# Import shapefile informations
from shapelib import ShapeFile
import dbflib

# Shapely
from shapely.geometry import Polygon
from shapely.ops import cascaded_union

# Numpy and matplotlib
import numpy as np
from matplotlib.nxutils import points_inside_poly

try:
   import psyco
   psyco.full()
except ImportError:
   print "Psyco plugin missing, will run slower"
   pass

def main(shapefile, picklefile):
    if picklefile:
      [npts, x, y, z, zraw, xil, yil, grid, missing] = cPickle.load(open(picklefile,'rb'))
      points = np.vstack((x,y)).T

    # Load administrative area
    shp = ShapeFile(shapefile)
    dbf = dbflib.open(shapefile)

    coastline = []

    # Process every shape from the ShapeFile
    print "Processing shapes ..."
    for npoly in range(shp.info()[0]):
        shp_object = shp.read_object(npoly)
        shp_dict = dbf.read_record(npoly)
        verts = shp_object.vertices()

        if "NAME_1" in shp_dict:
          name = "%s" % (shp_dict["NAME_1"])
        else:
          name = "Unknown"

        print "Processing %s" % (name)
        # Extract city polygon vertices (ring per ring)
        for ring in verts:
          vx = []
          vy = []
          for point in ring:
            vx.append(point[0])
            vy.append(point[1])

          # Only process big enough rings
          if len(vx) > 256: # big enough
            poly_verts = zip(vx,vy)

            if picklefile:
              # Compute intersections with the city
              intersection = points_inside_poly(points, poly_verts)
              npts = sum(1 for x in points_inside_poly(points, poly_verts) if x)
            else:
              npts = 1 # Add this polygon

            # Add the ring to the coastine if measures inside
            if npts > 0:
              polygon = Polygon(poly_verts)
              if not polygon.is_empty and polygon.is_valid:
                print "- Add polygon (%d)" % (len(vx))
                coastline.append(polygon)
            else:
                print "- Skip polygon (%d)" % (len(vx))
    
    print "Union of %d polygons" % len(coastline)
    coast = cascaded_union(coastline)
    cPickle.dump(coast,open('coastline.pickle','wb'),-1)
    print "Done."

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    parser = OptionParser("Usage: safecastCoastline <shapefile>")
    parser.add_option("-s", "--safecast", dest="scfilename",
                      help="provice the safecast.pickle file for intersections.", metavar="FILE")

    (options, args) = parser.parse_args()
    
    if len(args) != 1:
        parser.error("Wrong number of arguments")

    main(args[0], options.scfilename)
