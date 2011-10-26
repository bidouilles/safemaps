#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
# Copyright (C) 2011  Lionel Bergeret
#
# ----------------------------------------------------------------
# The contents of this file are distributed under the CC0 license.
# See http://creativecommons.org/publicdomain/zero/1.0/
# ----------------------------------------------------------------

# system libraries
import os, sys
import re
import cPickle
from optparse import OptionParser

# matplotlib libraries
import matplotlib
matplotlib.use('Agg') # for CGI script (no display)
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import matplotlib.patheffects as PathEffects
from mpl_toolkits.basemap import Basemap as Basemap

# mathematical libraries
import numpy as np
import pylab as pl
import math

# Import shapefile informations
from shapelib import ShapeFile
import dbflib

# Shapely
from shapely.geometry import Polygon
from shapely.ops import cascaded_union

# Multiprocess
from multiprocessing import Pool

# Safecast common
from safecastCommon import mask_outside_polygons, cmap_discretize

try:
   import psyco
   psyco.full()
except ImportError:
   print "Psyco plugin missing, will run slower"
   pass

# Japan limits
lat_min = 25.29
lon_min = 120.91
lat_max = 46.06
lon_max = 147.45

# Default parameters
dataFolder = "data"
renderedMapsFolder = "SafecastMap"

administrativeShapefile = "%s/JPN_adm1" % dataFolder
locationShapefile = "%s/asia_eastern_asia_japan_location" % dataFolder
lakesPickle = "%s/waterbodies.pickle" % dataFolder
coastlinePickle = "%s/coastline.pickle" % dataFolder

renderCPM = True
distanceAngle = 45
safecastDatasetName = "100m"

def main(pickleName, renderCities, sieverts, uncovered):
    try:
      [npts, x, y, z, zraw, xil, yil, grid, missing] = cPickle.load(open(pickleName,'rb'))
      lakes = cPickle.load(open(lakesPickle,'rb'))
      if uncovered:
        coastline = cPickle.load(open(coastlinePickle,'rb'))
      print "Pickle file loaded (%d points)" % npts
    except:
      print "An error occurs loading the pickle datas !!!"
      sys.exit(1)

    if renderCities:
      administrativeShapefile = "%s/JPN_adm2" % dataFolder
      distanceAngle = 5
    else:
      administrativeShapefile = "%s/JPN_adm1" % dataFolder
      distanceAngle = 45

    if sieverts:
      print "Rendering using uSv/h as unit"
      renderCPM = False
    else:
      print "Rendering using CPM as unit"
      renderCPM = True

    # Load locations names and position
    print "Load locations"
    shp = ShapeFile(locationShapefile)
    dbf = dbflib.open(locationShapefile)

    lx = []
    ly = []
    lname = []
    for name in range(shp.info()[0]):
        shp_object = shp.read_object(name)
        shp_dict = dbf.read_record(name)
        cx , cy = shp_object.vertices()[0]
        names = re.split('[()]', shp_dict["NAME"])
        if len(names)>1:
          cityname = names[1]
        else:
          cityname = ""
        lx.append(cx)
        ly.append(cy)
        lname.append(cityname)

    if uncovered:
      # Compute all uncovered areas
      difference = coastline.difference(missing)
      missing = difference

    # Create the thread pool
    pool = Pool()

    # Load Japan administrative area
    shp = ShapeFile(administrativeShapefile)
    dbf = dbflib.open(administrativeShapefile)

    # Process every shape from the ShapeFile
    print "Processing shapes ..."
    for npoly in range(shp.info()[0]):
        shpsegs = []
        shpinfo = []
        vx = []
        vy = []
 
        shp_object = shp.read_object(npoly)
        shp_dict = dbf.read_record(npoly)
        verts = shp_object.vertices()

        # Start building the city map
        if renderCities:
          name = "%s" % (shp_dict["NAME_2"])
          folder = "%s" % (shp_dict["NAME_1"])
        else:
          name = "%s" % (shp_dict["NAME_1"])
          folder = "All"

        #if (name not in ["Koriyama", "Fukushima", "Tokyo"]):
        #if (folder not in ["Fukushima"]):
        #if (shp_dict["NAME_1"] not in ["Fukushima", "Miyagi", "Tokyo", "Chiba", "Ibaraki", "Shizuoka", "Iwate", "Kanagawa"]):
        #   print "Skipping %s" % name
        #   continue

        # Extract city polygon vertices
        # Biggest ring only
        biggestRing = 0
        biggestRingSize = 0
        for ring in verts:
          if len(ring) > biggestRingSize:
             biggestRingSize = len(ring)
             biggestRing = ring

        for point in biggestRing:
            vx.append(point[0])
            vy.append(point[1])

        poly_verts = zip(vx,vy)

        # Compute intersections with the city
        from matplotlib.nxutils import points_inside_poly
        points = np.vstack((x,y)).T
        intersection = points_inside_poly(points, poly_verts)

        # Compute a small statistics for the city
        measures = []
        for i in range(len(intersection)-1):
          if (intersection[i]):
            measures.append(zraw[i])
        measures = np.array(measures)

        npts, minCPM, maxCPM, medianCPM = (0,0,0,0)
        npts = np.size(measures)
        if npts:
          minCPM = measures.min()
          maxCPM = measures.max()
          medianCPM = np.median(measures)

        # Only keep cities information inside the region
        cities_verts = zip(lx,ly)
        cities_inside = points_inside_poly(cities_verts, poly_verts)
        cities = [(lx[i], ly[i], lname[i]) for i in range(len(cities_inside)-1) if (cities_inside[i])]

        if npts>0:
          title = 'Safecast %s - Griddata (%s points) - %s [%s]\n(min, median, max) = (%d, %d, %d) CPM' % (safecastDatasetName, npts, name, folder, minCPM, medianCPM, maxCPM)
          DrawMap(title, 2, vx, vy, min(vx), max(vx), min(vy), max(vy), npts, x, y, z, xil, yil, grid, name, folder, cities, missing, uncovered, lakes)
          #pool.apply_async(DrawMap, (title, 2, vx, vy, min(vx), max(vx), min(vy), max(vy), npts, x, y, z, xil, yil, grid, name, folder, cities, missing, uncovered, lakes))
        else:
           print "Skipping %s" % name

    # Wait the pool to be completed
    pool.close()
    pool.join()

def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.
 
    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")
 
    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)
 
    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):
 
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu
 
    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
 
    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
 
    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi
 
    return (glon2, glat2, baz)
 
def equi(m, centerlon, centerlat, radius, xlim, ylim, *args, **kwargs):
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])
 
    X,Y = m(X,Y)
    plt.plot(X,Y,**kwargs)
    
    # Add label only on visible circles
    angle = 0
    for xpt,ypt in zip(X,Y):
       if (xpt>xlim[0] and xpt<xlim[1]) and (ypt>ylim[0] and ypt<ylim[1]):
         if (angle%distanceAngle == 0):
           label = plt.text(xpt+5,ypt+5,"%dkm" % radius, fontsize=10, zorder = 21)
           plt.setp(label, path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")], zorder = 21)
         angle+=1

def DrawMap(title, landwidth, vx, vy, lon_min, lon_max, lat_min, lat_max, npts, x, y, z, xi, yi, grid, name, folder, cities, missing, uncovered, lakes):
    print "Generating map for %s - %s (%d points) [%f,%f - %f,%f]" % (name, folder, npts, lat_max, lon_min, lat_min, lon_max)

    # Create the basemap
    m = Basemap(projection='merc', llcrnrlon=lon_min ,llcrnrlat=lat_min, urcrnrlon=lon_max ,urcrnrlat=lat_max, resolution='i')

    # Load and draw the shapefile
    japan_shp_info = m.readshapefile(administrativeShapefile,administrativeShapefile, color='b', linewidth = landwidth)
    
    # Compute the scale for the parallels and meridians lines
    scale = min((lon_max-lon_min), (lat_max-lat_min))/4

    # Draw parallels and meridians lines
    m.drawparallels(np.arange(y.min(),y.max(),scale),labels=[1,0,0,0],color='black',dashes=[1,0],labelstyle='+/-',linewidth=0.2) # draw parallels
    m.drawmeridians(np.arange(x.min(),x.max(),scale),labels=[0,0,0,1],color='black',dashes=[1,0],labelstyle='+/-',linewidth=0.2) # draw meridians

    # Compute Safecast color map
    cmap = cmap_discretize(cm.RdYlBu_r, 16, 0.)
    if renderCPM:
      normCPM = colors.Normalize(vmin=0,vmax=350)
    else:
      normCPM = colors.Normalize(vmin=0.0,vmax=1.0)

    # Draw Safecast data on the map
    lon,lat = m(x,y)
    m.scatter(lon,lat,s=2, c=z, cmap=cmap, linewidths=0.5, alpha=0.4, norm=normCPM, zorder = 5)

    # Draw contour color map
    xim, yim = m(*np.meshgrid(xi, yi))
    levels = range(0, 400, 350/16)
    #CS = m.contour(xim,yim, grid, levels, colors='k', linewidths=0.5)
    m.contourf(xim,yim, grid, levels, cmap=cmap, norm=normCPM)

    #plt.hexbin(lon,lat, C=z, cmap=cmap, gridsize=100, norm=normCPM, marginals=False)

    # Exclude uncovered areas
    if uncovered:
      for patch in missing.geoms:
        assert patch.geom_type in ['Polygon']
        assert patch.is_valid

        # Fill and outline each patch
        x, y = patch.exterior.xy
        x, y = m(x, y)
        plt.fill(x, y, color='#FFFF00', aa=True, alpha=1.0, hatch="x") 

    # Clip the water bodies area (white hatch)
    for patch in lakes.geoms:
        assert patch.geom_type in ['Polygon']
        assert patch.is_valid

        # Fill and outline each patch
        x, y = patch.exterior.xy
        x, y = m(x, y)
        plt.fill(x, y, color="#FFFFFF", aa=True, alpha=1.0, hatch="o") 
        m.plot(x, y, color="#FFFFFF", aa=True, lw=1.0, alpha=0.0) # needed for basemap to scale/crop the area

    # Add city names
    for city in cities:
      x,y = m(city[0], city[1])
      m.plot(x, y, 'bo')
      cnametext = plt.text(x+50, y+50, city[2].decode('utf-8'), fontsize=5, zorder = 30)
      plt.setp(cnametext, path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")], zorder = 30)

    # Add prefecture contour
    if len(vx) > 0:
      # Find main polygon
      vertices = zip(vx, vy)
      final = 0
      for i in reversed(range(1,len(vertices)-1,1)):
        x, y = vertices[i-1]
        nx, ny = vertices[len(vertices)-1]
        if (x==nx) and (y==ny):
          final = i
          break

      ax = plt.gca() # get current axes instance 
      #vxx, vyy = m(vx, vy) # project the points
      vxx, vyy = m(vx[final-1:], vy[final-1:]) # project the points
      poly_verts = zip(vxx,vyy) # zip the result
      
      # Reverse polygon vertices
      org = poly_verts
      poly_verts = []
      for a in reversed(org):
         poly_verts.append(a)

      # Mask out
      mask_outside_polygons([poly_verts], "white")

      # Draw contour
      from matplotlib.patches import Polygon
      from matplotlib.patches import Shadow
      poly = Polygon(poly_verts, edgecolor="white", fill=False, label=name, linewidth=4, zorder = 10) # facecolor="b"
      ax.add_patch(poly)
      poly = Polygon(poly_verts, edgecolor="black", fill=False, label=name, linewidth=2, zorder = 10) # facecolor="b"
      ax.add_patch(poly)

    #clbls = plt.clabel(CS, inline=1, fontsize=5, fmt='%1.0f CPM', use_clabeltext=True)
    #plt.setp(clbls, path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")], zorder = 10)

    DefaultSize = plt.gcf().get_size_inches()

    MaxSize = max(DefaultSize[0], DefaultSize[1])
    SizeFactor = 16.5/max(DefaultSize[0], DefaultSize[1])
    #print "Page size [%f, %f]" % (DefaultSize[0] * SizeFactor, DefaultSize[1] * SizeFactor)
    plt.gcf().set_size_inches( (DefaultSize[0] * SizeFactor, DefaultSize[1] * SizeFactor) )

    plt.title(title)

    # Add distance circle from Daiichi power plant
    daiichilat = 37.425252
    daiichilon = 141.033247
    radii = [20, 30, 60, 100, 160, 250, 400]
    ax = plt.gca() # get current axes instance 
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    #for radius in radii:
    for radius in range(20, 400, 10):
      equi(m, daiichilon, daiichilat, radius, xlim, ylim, linewidth=1, zorder = 20)

    plt.colorbar() # Legend

    # Save the result
    if not os.path.exists("%s/%s" % (renderedMapsFolder, folder)):
      os.makedirs("%s/%s" % (renderedMapsFolder, folder))

    # Cleanup the output filename
    name.replace("'","_")
    name.replace(" ","_")

    # Save png file
    plt.savefig("%s/%s/safecast_%s_%s.png" % (renderedMapsFolder, folder, name, folder), dpi = (200))
    plt.clf() # clear the plot (free the memory for the other threads)

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    parser = OptionParser("Usage: safecastGridMaps [options] <safecast-pickle>")
    parser.add_option("-c", "--cities",
                      action="store_true", dest="cities", default=False,
                      help="render the cities maps (default is prefectures)")
    parser.add_option("-s", "--sieverts",
                      action="store_true", dest="sieverts", default=False,
                      help="render the maps using uSv/h data (default in CPM)")
    parser.add_option("-u", "--uncovered",
                      action="store_true", dest="uncovered", default=False,
                      help="add uncovered yellow area to the map (default is none)")

    (options, args) = parser.parse_args()
    
    if len(args) != 1:
        parser.error("Wrong number of arguments")

    main(args[0], options.cities, options.sieverts, options.uncovered)
