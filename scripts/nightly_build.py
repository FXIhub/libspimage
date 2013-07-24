#!/usr/bin/env python

import os
import socket

os.chdir(os.path.dirname(os.path.realpath(__file__))+"/../build-"+socket.gethostname())
os.system("make Nightly")
