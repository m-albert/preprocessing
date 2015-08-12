__author__ = 'malbert'

from imports import *
import sys

if len(sys.argv) < 5 or sys.argv[-2][-2:] != 'h5':
    print 'usage: python <inFile> <labelIndex> <multiplicationFactor> <outDir>'

# process input arguments
inFile = sys.argv[-4]
labelIndex = sys.argv[-3]
multiplicationFactor = sys.argv[-2]
outDir = sys.argv[-1]

# take only label with indicated index
im = n.array(h5py.File(inFile)['DS1']).squeeze()[:,:,:,labelIndex]

# convert to tif friendly uint16
im = (im*multiplicationFactor).astype(n.uint16)

# write image
tifffile.imsave(os.path.join(outDir,inFile.split('.')[0]+'.tif'),im)

