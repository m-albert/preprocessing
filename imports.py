# this file handles the dependencies and adapts to the platform

import sys,os,inspect,subprocess,shutil,copy,ast,time,itertools,re
from itertools import combinations

importsDict = dict()

curDir = os.path.dirname(os.path.abspath(__file__))

timestamp = time.localtime()
timestamp = '%02d%02d%02d_%02d%02d%02d' %tuple([timestamp[i] for i in range(6)])
importsDict['timestamp'] = timestamp

relPathAppends = []
if sys.platform[:3] == 'lin':
    relPathAppends.append('dependencies_linux/SimpleITKcurrent')
    # if sys.version[:3] == '2.6':
    #     relPathAppends.append('dependencies_linux/SimpleITKcurrent')
    #     #relPathAppends.append('dependencies_linux/SimpleITK-0.9.0.dev235-py2.6-linux-x86_64.egg')
    # elif sys.version[:3] == '2.7':
    #     #relPathAppends.append('dependencies_linux/')
    #     relPathAppends.append('dependencies_linux/SimpleITK-0.6.1-py2.7-linux-x86_64.egg_FILES')
    elastixDir = os.path.join(curDir,'dependencies_linux/elastix_linux64_v4.7/bin/elastix')
    transformixDir = os.path.join(curDir,'dependencies_linux/elastix_linux64_v4.7/bin/transformix')
    tmpFolder = os.path.join('/tmp','fusionTmpFolder_%s' %timestamp)
    #tmpFolder = os.path.join('/dev/shm','fusionTmpFolder_%s' %timestamp)
    os.environ['LD_LIBRARY_PATH'] = os.path.join(curDir,'dependencies_linux/elastix_linux64_v4.7/lib')

if sys.platform[:3] == 'win':
    relPathAppends.append('dependencies_windows\\SimpleITK-0.6.1-py2.7-win-amd64')
    elastixDir = os.path.join(curDir,'dependencies_windows\\elastix_windows64_v4.6\\elastix.exe')
    transformixDir = os.path.join(curDir,'dependencies_windows\\elastix_windows64_v4.6\\transformix.exe')
    tmpFolder = os.path.join('C:\\Windows\\Temp','fusionTmpFolder_%s' %timestamp) #untested!

if sys.platform[:3] == 'dar':
    elastixDir = os.path.join(curDir,'../../software/elastix_macosx64_v4.8/bin/elastix')
    transformixDir = os.path.join(curDir,'../../software/elastix_macosx64_v4.8/bin/transformix')
    os.environ['LD_LIBRARY_PATH'] = os.path.join(curDir,'../../software/elastix_macosx64_v4.8/lib')
    tmpFolder = os.environ['TMPDIR']

for relPathAppend in relPathAppends:
    sys.path.append(os.path.join(curDir,relPathAppend))

print 'importing SimpleITK09 as sitk09!'
import SimpleITK as sitk
#import SimpleITK09 as sitk09
import numpy as n
from scipy import ndimage,optimize,spatial
import pdb
import tifffile
import pickle
import h5py

import extras
import beads
import fusion
import filters
import transformations


from configuration import config

importsDict['elastixDir'] = elastixDir
importsDict['transformixDir'] = transformixDir
if not os.path.exists(tmpFolder): os.mkdir(tmpFolder)
importsDict['tmpDir'] = tmpFolder

# some aliases
sitk.gifa = sitk.GetImageFromArray
sitk.gafi = sitk.GetArrayFromImage

os.environ['SITK_SHOW_EXTENSION']='.tif'
os.environ['SITK_SHOW_COMMAND']='fiji'
