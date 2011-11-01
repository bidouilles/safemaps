#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
# Copyright (C) 2011  Lionel Bergeret
#
# ----------------------------------------------------------------
# The contents of this file are distributed under the CC0 license.
# See http://creativecommons.org/publicdomain/zero/1.0/
# ----------------------------------------------------------------

import cgi, cgitb 

# http://localhost:8888/renderer/interpolate.py?x=7275&y=3226&zoom=13

import cPickle
from optparse import OptionParser
import os
import sys
import time
sys.path.append(os.path.dirname(os.path.abspath(__file__))+ "/../..")

# External commands 
import shlex, subprocess

# Download URL
import urllib

# Multiprocess
from multiprocessing import Pool

# Safecast common
from safecastCommon import GoogleProjection

# Build the tiles using the python cgi or not
running_on_web = False

# Japan limits
jp_lat_min = 34.6332
jp_lon_min = 135.0076
jp_lat_max = 41.2778
jp_lon_max = 142.0532

gridSize = 4
concurrency = 4

def download(url, output):
    print "Downloading tile %s" % (url)
    urllib.urlretrieve(url, output)

def buildTile(command):
    print "Building tile %s" % (command)
    args = shlex.split(command)
    retcode = subprocess.call(args) # Success!

# from http://stackoverflow.com/questions/6728236/exception-thrown-in-multiprocessing-pool-not-detected
# Pool debugging

if __name__ == "__main__":
 parser = OptionParser("Usage: safecastPrecompute [options] <safecast-csv-file>")
 parser.add_option("-u", "--uncovered",
                      action="store_true", dest="uncovered", default=False,
                      help="generate uncovered yellow area to the map (default is none)")

 (options, args) = parser.parse_args()

 projection = GoogleProjection()

 for gzoom in range(4, 14):
  print "Zoom level", gzoom

  gx0 , gy0 = projection.fromLLtoPixel((jp_lon_min, jp_lat_max),gzoom) # top right
  gx1 , gy1 = projection.fromLLtoPixel((jp_lon_max, jp_lat_min),gzoom) # bottom left
  
  gx0 = int(gx0/256)
  gy0 = int(gy0/256)

  gx1 = int(gx1/256)
  gy1 = int(gy1/256)
  
  # Optimal size for the slice computation
  gridSize = gx1 - gx0 + 1
  if gridSize == 0: gridSize = 1
  if gridSize > 82: gridSize = 64

  pool = Pool(processes = concurrency)

  # Start looping
  gx = gx0
  gy = gy0
  while (gx <= gx1):
   while (gy <= gy1):
     if running_on_web:
       # Create the cached filename
       url = "http://localhost/renderer/interpolate.py?x=%s&y=%s&zoom=%s&size=%s" % (gx, gy, gzoom, gridSize)

       # Draw the new tiles
       #download(url, "/dev/null")
       pool.apply_async(download, (url, "/dev/null"))
     else:
       # Create the cached filename
       cmd = "./interpolate.py %s %s %s %s %d" % (gx, gy, gzoom, gridSize, options.uncovered)

       # Draw the new tiles
       #buildTile(cmd)
       pool.apply_async(buildTile, (cmd,))

     gy += 1
     
   gy = gy0
   gx += gridSize
  
  # Wait the pool to be completed
  pool.close()
  pool.join()
