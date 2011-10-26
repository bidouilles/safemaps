#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
# Copyright (C) 2011  Lionel Bergeret
#
# ----------------------------------------------------------------
# The contents of this file are distributed under the CC0 license.
# See http://creativecommons.org/publicdomain/zero/1.0/
# ----------------------------------------------------------------

import CGIHTTPServer
import BaseHTTPServer

class Handler(CGIHTTPServer.CGIHTTPRequestHandler):
    cgi_directories = ["/renderer"]

# -----------------------------------------------------------------------------
# Simple tile server started at localhost port 8888
# -----------------------------------------------------------------------------
if __name__ == '__main__':
  PORT = 8888
  httpd = BaseHTTPServer.HTTPServer(("", PORT), Handler)
  print "The server is ready and can be accessed at http://localhost:%s/index.html" % PORT
  httpd.serve_forever()
