#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
# Copyright (C) 2011  Lionel Bergeret
#
# ----------------------------------------------------------------
# The contents of this file are distributed under the CC0 license.
# See http://creativecommons.org/publicdomain/zero/1.0/
# ----------------------------------------------------------------

# mathematical libraries
import numpy as np
from scipy import interpolate

# matplotlib
from matplotlib import cm, colors
import matplotlib.pyplot as plt

# Math
from math import pi,cos,sin,log,exp,atan
DEG_TO_RAD = pi/180
RAD_TO_DEG = 180/pi

# -----------------------------------------------------------------------------
# Discretize a colormap
# -----------------------------------------------------------------------------
def cmap_discretize(cmap, N, crop = 0):
    """Return a discrete colormap from the continuous colormap cmap"""

    cdict = cmap._segmentdata.copy()
    # N colors
    colors_i = np.linspace(crop,1.,N)
    # N+1 indices
    indices = np.linspace(0,1.,N+1)
    for key in ('red','green','blue'):
        # Find the N colors
        D = np.array(cdict[key])
        I = interpolate.interp1d(D[:,0], D[:,1])
        acolors = I(colors_i)
        # Place these colors at the correct indices.
        A = np.zeros((N+1,3), float)
        A[:,0] = indices
        A[1:,1] = acolors
        A[:-1,2] = acolors
        # Create a tuple for the dictionary.
        L = []
        for l in A:
            L.append(tuple(l))
        cdict[key] = tuple(L)

    # Add dark gray at the end
    for key in ('red','green','blue'):
      L = list(cdict[key])
      L[len(L)-1] = (1.0,0.3,0.3) # gray
      cdict[key] = tuple(L)

    return colors.LinearSegmentedColormap('colormap',cdict,1024)

# -----------------------------------------------------------------------------
# Perform a google projection
# -----------------------------------------------------------------------------
# from http://svn.openstreetmap.org/applications/rendering/mapnik/generate_tiles.py
#
class GoogleProjection:
    def __init__(self,levels=18):
        self.Bc = []
        self.Cc = []
        self.zc = []
        self.Ac = []
        c = 256
        for d in range(0,levels):
            e = c/2;
            self.Bc.append(c/360.0)
            self.Cc.append(c/(2 * pi))
            self.zc.append((e,e))
            self.Ac.append(c)
            c *= 2

    def minmax (self,a,b,c):
        a = max(a,b)
        a = min(a,c)
        return a
                
    def fromLLtoPixel(self,ll,zoom):
         d = self.zc[zoom]
         e = round(d[0] + ll[0] * self.Bc[zoom])
         f = self.minmax(sin(DEG_TO_RAD * ll[1]),-0.9999,0.9999)
         g = round(d[1] + 0.5*log((1+f)/(1-f))*-self.Cc[zoom])
         return (e,g)
     
    def fromPixelToLL(self,px,zoom):
         e = self.zc[zoom]
         f = (px[0] - e[0])/self.Bc[zoom]
         g = (px[1] - e[1])/-self.Cc[zoom]
         h = RAD_TO_DEG * ( 2 * atan(exp(g)) - 0.5 * pi)
         return (f,h)

    def convert(self, gx, gy, zoom):
        # Calculate pixel positions of bottom-left & top-right
        p0 = (gx * 256, (gy + 1) * 256)
        p1 = ((gx + 1) * 256, gy * 256)

        # Convert to LatLong (EPSG:4326)
        l0 = self.fromPixelToLL(p0, zoom);
        l1 = self.fromPixelToLL(p1, zoom);

        # Get tile width and height in degrees
        lonWidth = (l1[0]-l0[0])
        latHeight = (l1[1]-l0[1])

        # Convert main tile position to LatLong (EPSG:4326)
        latlon = self.fromPixelToLL(((gx)*256, (gy+1)*256), zoom); # top-left

        return (latlon[0], lonWidth, latlon[1], latHeight)

# -----------------------------------------------------------------------------
# Mask outside polygons
# -----------------------------------------------------------------------------
# Original from http://stackoverflow.com/questions/3320311/fill-outside-of-polygon-mask-array-where-indicies-are-beyond-a-circular-boundar
# Modified to add multiple polygons support (Lionel)
#
def mask_outside_polygons(polygons, pcolor, ax=None):
    """
    Plots a mask on the specified axis ("ax", defaults to plt.gca()) such that
    all areas outside of the polygon specified by "poly_verts" are masked.  

    "poly_verts" must be a list of tuples of the verticies in the polygon in
    counter-clockwise order.

    Returns the matplotlib.patches.PathPatch instance plotted on the figure.
    """
    import matplotlib.patches as mpatches
    import matplotlib.path as mpath

    if ax is None:
        ax = plt.gca()

    # Get current plot limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Verticies of the plot boundaries in clockwise order
    bound_verts = [(xlim[0], ylim[0]), (xlim[0], ylim[1]), 
                   (xlim[1], ylim[1]), (xlim[1], ylim[0]), 
                   (xlim[0], ylim[0])]

    # A series of codes (1 and 2) to tell matplotlib whether to draw a line or 
    # move the "pen" (So that there's no connecting line)
    bound_codes = [mpath.Path.MOVETO] + (len(bound_verts) - 1) * [mpath.Path.LINETO]
    poly_codes = []
    poly_verts = []
    for poly in polygons:
      poly_codes += [mpath.Path.MOVETO] + (len(poly) - 1) * [mpath.Path.LINETO]
      poly_verts += poly

    # Plot the masking patch
    path = mpath.Path(bound_verts + poly_verts, bound_codes + poly_codes)
    patch = mpatches.PathPatch(path, facecolor=pcolor, edgecolor='none')
    patch = ax.add_patch(patch)

    # Reset the plot limits to their original extents
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    return patch

# -----------------------------------------------------------------------------
# Make a color transparent from an image
# -----------------------------------------------------------------------------
def makeColorTransparent(image, color):
    """
    Replace the color from the input image by transparent
    """
    img = image
    img = img.convert("RGBA")
    pixdata = img.load()

    for y in xrange(img.size[1]):
       for x in xrange(img.size[0]):
          if pixdata[x, y] == color:
              pixdata[x, y] = (255, 255, 255, 0)

    return img

