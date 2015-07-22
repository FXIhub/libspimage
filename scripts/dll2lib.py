#!/usr/bin/env python

import os
import sys
import subprocess

if len(sys.argv) < 2:
    print "Usage dll2lib.py <file.dll>"
    sys.exit(0)


dll_name = sys.argv[1]
def_name = dll_name[:-4]+".def"
lib_name = dll_name[:-4]+".lib"
output = subprocess.check_output(['dumpbin','/exports',dll_name])
output = output.split('\r\n')
fdef = open(def_name,'w')
fdef.write('EXPORTS\r\n')
for line in output:
    fields = line.split()
    if len(fields) == 4:
        try:
            int(fields[0])
            int(fields[1],16)
            fdef.write("%s\r\n" % (fields[3]))
        except ValueError:
            pass
fdef.close()

os.system("lib /def:%s /out:%s /machine:x64" % (def_name, lib_name))
