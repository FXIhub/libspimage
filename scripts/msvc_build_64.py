#!/usr/bin/env python

import os
import subprocess
import urllib

def dll2lib(dll_name):
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
    def2lib(def_name)
    # os.system("lib /def:%s /out:%s /machine:x64" % (def_name, lib_name))

def def2lib(def_name):
    lib_name = def_name[:-4]+".lib"
    os.system("lib /def:%s /out:%s /machine:x64" % (def_name, lib_name))

def _reporthook(numblocks, blocksize, filesize, url=None):
#print "reporthook(%s, %s, %s)" % (numblocks, blocksize, filesize)
    base = os.path.basename(url)
#XXX Should handle possible filesize=-1.
    try:
        percent = min((numblocks*blocksize*100)/filesize, 100)
    except:
        percent = 100
    if numblocks != 0:
        sys.stdout.write("\b"*70)
        sys.stdout.write("%-66s%3d%%" % (base, percent))

def geturl(url, dst):
    print "get url '%s' to '%s'" % (url, dst)
    if sys.stdout.isatty():
        urllib.urlretrieve(url, dst,
                           lambda nb, bs, fs, url=url: _reporthook(nb,bs,fs,url))
        sys.stdout.write('\n')
    else:
        urllib.urlretrieve(url, dst)

env_map = os.environ;


# Install cmake
cmake_url = "http://www.cmake.org/files/v3.2/cmake-3.2.3-win32-x86.exe"
cmake_exe = "cmake-3.2.3-win32-x86.exe";
if(not os.path.exists(cmake_exe)):
    geturl(cmake_url,cmake_exe);
    os.system(cmake_exe);

# Install git
git_url = "https://github.com/msysgit/msysgit/releases/download/Git-1.9.5-preview20150319/Git-1.9.5-preview20150319.exe"
git_exe = "Git-1.9.5-preview20150319.exe"
if(not os.path.exists(git_exe)):
    geturl(git_url,git_exe);
    os.system(git_exe);

# Install msys2 (64-bit)
msys2_url = "http://sourceforge.net/projects/msys2/files/Base/x86_64/msys2-x86_64-20150512.exe/download"
msys2_exe = "msys2-x86_64-20150512.exe"
if(not os.path.exists(msys2_exe)):
    geturl(msys2_url,msys2_exe);
    os.system(msys2_exe);

# Install Visual Studio Community 2013
vs_url = "https://go.microsoft.com/fwlink/?LinkId=532606&clcid=0x409"
vs_exe = "vs_community.exe"
if(not os.path.exists(vs_exe)):
    geturl(vs_url,vs_exe);
    os.system(vs_exe);

# Download msys64 packages (libpng, libtiff, zlib, gsl)
print 'Open msys2 terminal and run'
print 'pacman -Ss libpng'
print 'to search for libpng'
print 'and something like pacman -S mingw-w64-x86_64-libpng to install it'
print 'Repeat for libtiff zlib and gsl'

# Convert msys64 .dll to .lib
dll2lib('c:\msys2\mingw64\bin\libpng16-16.dll')
dll2lib('c:\msys2\mingw64\bin\libtiff-5.dll')
dll2lib('c:\msys2\mingw64\bin\zlib1.dll')
dll2lib('c:\msys2\mingw64\bin\libgsl-0.dll')
dll2lib('c:\msys2\mingw64\bin\libgslcblas-0.dll')

# Install latest HDF5
hdf5_url = "http://www.hdfgroup.org/ftp/HDF5/current/bin/windows/hdf5-1.8.15-patch1-win64-vs2013-shared.zip"
hdf5_zip = "hdf5-1.8.15-patch1-win64-vs2013-shared.zip"

if(not os.path.exists(hdf5_zip)):
    geturl(hdf5_url,hdf5_zip);
    file = zipfile.ZipFile(hdf5_zip, "r");
    file.extractall("c:/hdf5");
    os.system("c:/hdf5/HDF5-1.8.15-win64");

# Install fftw3 (not from msys64)

fftw_zip = "fftw-3.3.4-dll64.zip";
fftw_url = "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4-dll64.zip";

if(not os.path.exists(fftw_zip)):
    geturl(fftw_url,fftw_zip);
    file = zipfile.ZipFile(fftw_zip, "r");
    file.extractall("c:/fftw3");

def2lib('c:/fftw3/libfftw3-3.def')
def2lib('c:/fftw3/libfftw3f-3.def')
def2lib('c:/fftw3/libfftw3l-3.def')


# Install CUDA 7
cuda_url = 'http://developer.download.nvidia.com/compute/cuda/7_0/Prod/network_installers/cuda_7.0.28_windows_network.exe'
cuda_exe = 'cuda_7.0.28_windows_network.exe'
if(not os.path.exists(cuda_exe)):
    geturl(cuda_url,cuda_exe);
    os.system(cuda_exe);


# Run cmake

p = subprocess.Popen(['cmake',
                  '-G','Visual Studio 12 Win64',
                  '-DTIFF_LIBRARY:FILEPATH=c:/msys64/mingw64/bin/libtiff-5.lib',
                  '-DTIFF_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DPNG_LIBRARY:FILEPATH=c:/msys64/mingw64/bin/libpng16-16.lib',
                  '-DPNG_PNG_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DZLIB_LIBRARY:FILEPATH=c:/msys64/mingw64/bin/zlib1.lib',
                  '-DZLIB_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DFFTW3_FFTWF_LIBRARY:FILEPATH=c:/fftw3/libfftw3f-3.lib',
                  '-DFFTW3_FFTW_LIBRARY:FILEPATH=c:/fftw3/libfftw3-3.lib',
                  '-DFFTW3_FFTWL_LIBRARY:FILEPATH=c:/fftw3/libfftw3l-3.lib',
                  '-DFFTW3_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DHDF5_LIBRARY:FILEPATH=c:/Program Files/HDF_Group/HDF5/1.8.15/lib/hdf5.lib',
                  '-DHDF5_INCLUDE_DIR:PATH=c:/Program Files/HDF_Group/HDF5/1.8.15/include',
                  '-DGSL_GSL_LIBRARY:FILEPATH=c:/msys64/mingw64/bin/libgsl-0.lib',
                  '-DGSL_GSLCBLAS_LIBRARY:FILEPATH=c:/msys64/mingw64/bin/libgslcblas-0.lib',
                  '-DGSL_INCLUDE_DIR:PATH=c:/msys64/mingw64/include',
                  '-DCMAKE_VERBOSE_MAKEFILE:BOOL=FALSE',
                  '..'], env=env_map)

# Run vcvarsall.bat to get msbuild in the path

# Build things with msbuild spimage.sln

