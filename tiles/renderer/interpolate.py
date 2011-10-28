#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
# Copyright (C) 2011  Lionel Bergeret
#
# ----------------------------------------------------------------
# The contents of this file are distributed under the CC0 license.
# See http://creativecommons.org/publicdomain/zero/1.0/
# ----------------------------------------------------------------

# CGI modules
import cgi, cgitb
#cgitb.enable()  # for troubleshooting

# System
import cPickle
from optparse import OptionParser
import os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__))+ "/../..")

# Matplotlib
import matplotlib
matplotlib.use('Agg') # for CGI script (no display)
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

# Scipy, Numpy and pylab
from scipy import interpolate
import numpy as np

# Shapely
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import cascaded_union

# PIL
import Image, ImageMath

# Safecast common
from safecastCommon import GoogleProjection, cmap_discretize, mask_outside_polygons, makeColorTransparent

# Default parameter 
rendererFolder = "./"
tileFolder = "cached"
shapefile = "../../data/JPN_adm1"

uncovered = False # generate uncovered yellow areas only

googleWaterColorHtml = "#a5bfdd"
googleWaterColor = (165, 191, 221, 255)

def DrawTile(precision, lon_min, lon_width, lat_min, lat_height, gx, gy, gzoom, gsize_lon, gsize_lat, x, y, z, xi, yi, grid, missing, uncovered, coastline, waterbodies, name):
    fig = Figure()
    canvas = FigureCanvas(fig)

    # Precompute gsize x gsize tiles in one shot
    lon_max = lon_min + gsize_lon * lon_width
    lat_max = lat_min + gsize_lat * lat_height

    # Create the basemap and load the shapefile
    m = Basemap(projection='merc', llcrnrlon=lon_min ,llcrnrlat=lat_min, urcrnrlon=lon_max ,urcrnrlat=lat_max, resolution='i')

    # Cleanup all axes, title, ... from the canvas
    dpi = 64
    m.ax = fig.add_axes([0, 0, 1, 1], axis_bgcolor=(1.0,1.0,1.0,0.0), frameon=False)
    m.ax.get_frame().set_linewidth(0.0)
    fig.set_size_inches((256*gsize_lon)/dpi, (256*gsize_lat)/dpi) # we need 256x256 pixel images
    fig.set_facecolor((1.0,1.0,1.0,0.0))
    fig.figurePatch.set_alpha(0.0)
    fig.gca().axesPatch.set_alpha(0.0)
    
    #japan_shp_info = m.readshapefile("%s/%s" % (rendererFolder, shapefile),"%s/%s" % (rendererFolder, shapefile), color='b', linewidth = 0.5)

    # Compute Safecast color map
    cmap = cmap_discretize(cm.RdYlBu_r, 16, 0.)
    normCPM = colors.Normalize(vmin=0,vmax=350)

    # Exclude uncovered areas
    if uncovered:
      for patch in missing.geoms:
        assert patch.geom_type in ['Polygon']
        assert patch.is_valid

        if patch.area > 0.0016: # more than (0.04 degree x 0.04 degree) ~ (1km x 1km) area
          # Fill and outline each patch
          x, y = patch.exterior.xy
          x, y = m(x, y)
          m.ax.fill(x, y, color='#FFFF00', aa=True, alpha=1.0, hatch="x") 
          m.plot(x, y, color=googleWaterColorHtml, aa=True, lw=1.0, alpha=0.0) # needed for basemap to scale/crop the area

    # Draw countour interpolation map
    if not uncovered:
      xim, yim = m(*np.meshgrid(xi, yi))
      levels = range(0, 400, 350/16)
      m.contourf(xim,yim, grid, levels, cmap=cmap, norm=normCPM, vmin=0, vmax=350)

    # Draw Safecast data on the map
    if not uncovered:
      lon,lat = m(x,y)
      m.scatter(lon,lat,s=2, c=z, cmap=cmap, norm=normCPM, linewidths=0.2, alpha=1.0)

    # Clip outside coastlines area and water bodies
    if len(coastline) > 0 and len(waterbodies) > 0:
      polygonsToClip = []
      for patch in coastline.geoms:
        if not patch.is_empty and patch.is_valid:
          vx, vy = patch.exterior.xy
          vx.reverse()
          vy.reverse()
          mvx, mvy = m(vx,vy)
          polygonsToClip.append(zip(mvx, mvy))

      for patch in waterbodies.geoms:
        if not patch.is_empty and patch.is_valid:
          vx, vy = patch.exterior.xy
          mvx, mvy = m(vx,vy)
          polygonsToClip.append(zip(mvx, mvy))

      mask_outside_polygons(polygonsToClip, googleWaterColorHtml, ax = m.ax)

    # Save the result
    if not os.path.exists("%s/%s/%s" % (rendererFolder, tileFolder, gzoom)):
      os.makedirs("%s/%s/%s" % (rendererFolder, tileFolder, gzoom))

    #
    # Multiple tiles
    #
    if (gsize_lon > 1):
      tileNameTemp = "%s/%s/%s/%s" % (rendererFolder, tileFolder, gzoom, name)
      canvas.print_figure(tileNameTemp, dpi=dpi)

      # Start the tile cutting process
      image = Image.open(tileNameTemp)
      tile_width = 256
      tile_height = 256

      # Cut the tiles 
      currentx = 0
      currenty = 0
      googlex = gx
      googley = gy
      while currenty < image.size[1]:
        while currentx < image.size[0]:
          tile = image.crop((currentx,currenty,currentx + tile_width,currenty + tile_height))
          tilename = "safecast_griddata_%s_%s_%s.png" %(googlex, googley, gzoom)
          filename = "%s/%s/%s/%s" % (rendererFolder, tileFolder, gzoom, tilename)

          if not os.path.exists(filename):
            makeColorTransparent(tile, googleWaterColor).save(filename)

          currentx += tile_width
          googlex += 1

        currentx = 0 
        currenty += tile_height

        googley += 1
        googlex = gx
      os.remove(tileNameTemp)

    #
    # One tile
    #
    else:
      name = "safecast_griddata_%s_%s_%s.png" %(gx, gy, gzoom)
      tileName = "%s/%s/%s/%s" % (rendererFolder, tileFolder, gzoom, name)
      canvas.print_figure(tileName, dpi=dpi)
      makeColorTransparent(Image.open(tileName), googleWaterColor).save(tileName)

    # Clear the plot (free the memory for the other threads)
    plt.clf()

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
  # Check if running on web
  #if "HTTP_HOST" in os.environ.keys():
  if len(sys.argv) == 1 :
     running_on_web = True
  else:
     running_on_web = False

  if running_on_web:
     # Create instance of FieldStorage 
     form = cgi.FieldStorage() 

     # Get data from fields
     gx = int(form.getvalue('x'))
     gy = int(form.getvalue('y'))
     gzoom = int(form.getvalue('zoom'))
     gsize = int(form.getvalue('size'))
  else:
     parser = OptionParser("Usage: interpolate <x> <y> <zoom> <size>")

     (options, args) = parser.parse_args()
    
     if len(args) != 4:
        parser.error("Wrong number of arguments")

     # Get data from command line arguments
     gx = int(args[0])
     gy = int(args[1])
     gzoom = int(args[2])
     gsize = int(args[3])

  # Create the cached filename
  tilename = "safecast_griddata_%s_%s_%s.png" %(gx, gy, gzoom)
  filename = "%s/%s/%s/%s" % (rendererFolder, tileFolder, gzoom, tilename)

  projection = GoogleProjection()
  tilearea = projection.convert(gx,gy, gzoom)

  try:
    # Open requested tile
    f = open(filename, 'rb')
  except:
    # Load precomputed safecast griddata
    [npts, x, y, z, zraw, xi, yi, grid, missing] = cPickle.load(open('%s/safecast.pickle' % rendererFolder,'rb'))

    # Load precomputed coastline if any
    try:
      coastline = cPickle.load(open("%s/coastline.pickle" % rendererFolder,'rb'))
    except:
      coastline = []
      pass

    # Load precomputed waterbodies if any
    try:
      waterbodies = cPickle.load(open("%s/waterbodies.pickle" % rendererFolder,'rb'))
    except:
      waterbodies = []
      pass

    # All uncovered areas
    difference = coastline.difference(missing)
    missing = difference
    
    # Draw the new tiles
    DrawTile(60, tilearea[0], tilearea[1], tilearea[2], tilearea[3], gx, gy, gzoom, gsize, 1, x, y, z, xi, yi, grid, missing, uncovered, coastline, waterbodies, "safecast_temp_%s_%s_%s.png" %(gx, gy, gzoom))
    # Open requested tile
    f = open(filename, 'rb')
  
  if running_on_web:  
    # Send the png file to the client
    print "Content-type: image/png\n"
    print f.read()

  f.close()
    
