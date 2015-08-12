__author__ = 'malbert'

from imports import *

#####################
#####################
#####################
# process arguments

inFile1 = sys.argv[-2]
inFile2 = sys.argv[-1]


if inFile1[-2:] == 'h5':
    im1 = n.array(h5py.File(inFile1)['Data'])
    im2 = n.array(h5py.File(inFile2)['Data'])
    im1 = sitk.gifa(im1)
    im2 = sitk.gifa(im2)
elif inFile1[-2:] in ['if','ff']:
    im1 = sitk.ReadImage(inFile1)
    im2 = sitk.ReadImage(inFile2)

im1 = beads.scaleStack([3,3,3],im1)
im2 = beads.scaleStack([3,3,3],im2)

segmentedBeads = beads.segmentBeadsNils([im1,im2],1)
params = list(beads.findParametersNils(segmentedBeads))
params[1] = [1,-0.002,-0.29,0,-1,0,0,0,0.93,280,2037,154.3]
# try:
#     params = list(beads.findParametersSheared(segmentedBeads))
# except:
#     raise(Exception('could not find bead parameters'))

x = -n.mean([params[1][0],params[1][8]])

shearMatrixLeft = n.array([1,0,x,0,1,0,0,0,1])

im21 = beads.transformStack(params[1],im2)

tim1 = beads.transformStack(shearMatrixLeft,im1)
tim2 = beads.transformStack(shearMatrixLeft,im21)

sitk.WriteImage(im1,os.path.join(outputFolder,'stackf1.tif'))
sitk.WriteImage(tim2,os.path.join(outputFolder,'stackf2.tif'))

sitk.WriteImage(sitk.MaximumProjection(tim1,2),os.path.join(outputFolder,'stackf1_proj2.tif'))
sitk.WriteImage(sitk.MaximumProjection(tim2,2),os.path.join(outputFolder,'stackf2_proj2.tif'))

