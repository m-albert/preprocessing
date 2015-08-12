__author__ = 'malbert'

from imports import *

#####################
#####################
#####################
# process arguments

print '\n\n\ninput args:\n%s\n\n\n' %sys.argv

inFilePattern1                                  = sys.argv[-13] #expects %05d...
inFilePattern2                                  = sys.argv[-12]
startTime                                       = int(sys.argv[-11])
endTime                                         = int(sys.argv[-10])
mexicanHatFilterSize                            = float(sys.argv[-9])
inFile1                                         = sys.argv[-8]
inFile2                                         = sys.argv[-7]
direction1                                      = float(sys.argv[-6]) # 1
direction2                                      = float(sys.argv[-5]) # 0
diagdz                                          = float(sys.argv[-4]) # relative diagonal z spacing
isotropicSpacingFactor                          = float(sys.argv[-3]) # multiply xy spacing by this factor for obtaining
                                                                      # final isotropic spacing
outputFolder                                    = sys.argv[-2]
option                                          = int(sys.argv[-1])

dz = diagdz*n.cos(n.pi/4.)
outSpacing = n.array([isotropicSpacingFactor,isotropicSpacingFactor,isotropicSpacingFactor/dz])

if direction1 == 0:
    direction1 = -1.
if direction2 == 0:
    direction2 = -1.

if inFile1[-2:] == 'h5':
    im1 = n.array(h5py.File(inFile1,'r')['Data'])
    im2 = n.array(h5py.File(inFile2,'r')['Data'])
    im1 = sitk.gifa(im1)
    im2 = sitk.gifa(im2)
elif inFile1[-2:] in ['if','ff']:
    im1 = sitk.ReadImage(inFile1)
    im2 = sitk.ReadImage(inFile2)

im1.SetSpacing([1,1,dz])
im2.SetSpacing([1,1,dz])

tmpShape = n.array(im1.GetSize())
outShape = tmpShape/outSpacing

shearParams1 = n.zeros(12,dtype=n.float)
tmpShear1 = transformations.shear_matrix(n.pi/4.,[direction1,0,0],[0,0,0],[0,0,1])
shearParams1[:9] = n.array(tmpShear1[:3,:3].flatten())
shearParams1[9:] = n.zeros(3)

shearParams2 = n.zeros(12,dtype=n.float)
tmpShear2 = transformations.shear_matrix(n.pi/4.,[direction2,0,0],[0,0,0],[0,0,1])

shearParams2[:9] = n.array(tmpShear2[:3,:3].flatten())
shearParams2[9:] = n.zeros(3)

tim1 = beads.transformStack(shearParams1,im1)
tim2 = beads.transformStack(shearParams2,im2)

tim1.SetSpacing([1.,1,1])
tim2.SetSpacing([1.,1,1])
tim1.SetOrigin([0.,0,0])
tim2.SetOrigin([0.,0,0])

if option == 1:

    dz = diagdz*n.cos(n.pi/4.)
    segmentedBeads = beads.segmentBeadsNils([tim1,tim2],mexicanHatFilterSize)
    params = list(beads.findParametersNils(segmentedBeads))

    tmpParams = fusion.invertParams(params[1])
    print tmpParams
    tim1 = beads.transformStack([1.,0,0,0,1,0,0,0,1,0,0,0],tim1,outSpacing=outSpacing,outShape=outShape)
    tim2 = beads.transformStack(tmpParams,tim2,outSpacing=outSpacing,outShape=outShape)

sitk.WriteImage(tim1,os.path.join(outputFolder,'beads_1.tif'))
sitk.WriteImage(tim2,os.path.join(outputFolder,'beads_2.tif'))

sitk.WriteImage(sitk.MaximumProjection(tim1,2),os.path.join(outputFolder,'beads_1_proj2.tif'))
sitk.WriteImage(sitk.MaximumProjection(tim2,2),os.path.join(outputFolder,'beads_2_proj2.tif'))

averageBeads = (tim1+tim2)/2
sitk.WriteImage(averageBeads,os.path.join(outputFolder,'beads_average.tif'))
sitk.WriteImage(sitk.MaximumProjection(averageBeads,2),os.path.join(outputFolder,'beads_proj2_average.tif'))

averageBeads = sitk.gafi(averageBeads) #z,y,x
averageBeads = averageBeads.swapaxes(0,2)
averageBeads = sitk.gifa(averageBeads)

sitk.WriteImage(averageBeads,os.path.join(outputFolder,'beads_resliced_average.tif'))
sitk.WriteImage(sitk.MaximumProjection(averageBeads,2),os.path.join(outputFolder,'beads_resliced_proj2_average.tif'))



for time in range(startTime,endTime+1):

    inFile1 = inFilePattern1 %time
    inFile2 = inFilePattern2 %time

    if inFile1[-2:] == 'h5':
        im1 = n.array(h5py.File(inFile1,'r')['Data'])
        im2 = n.array(h5py.File(inFile2,'r')['Data'])
        im1 = sitk.gifa(im1)
        im2 = sitk.gifa(im2)
    elif inFile1[-2:] in ['if','ff']:
        im1 = sitk.ReadImage(inFile1)
        im2 = sitk.ReadImage(inFile2)

    im1.SetSpacing([1,1,diagdz])
    im2.SetSpacing([1,1,diagdz])

    tim1 = beads.transformStack(shearParams1,im1)
    tim2 = beads.transformStack(shearParams2,im2)

    tim1.SetSpacing([1.,1,1])
    tim2.SetSpacing([1.,1,1])
    tim1.SetOrigin([0.,0,0])
    tim2.SetOrigin([0.,0,0])

    if option == 1:
        tmpParams = fusion.invertParams(params[1])
        print tmpParams

        tim1 = beads.transformStack([1.,0,0,0,1,0,0,0,1,0,0,0],tim1,outSpacing=outSpacing,outShape=outShape)
        tim2 = beads.transformStack(tmpParams,tim2,outSpacing=outSpacing,outShape=outShape)
        # tim2 = beads.transformStack(tmpParams,tim2)

    sitk.WriteImage(tim1,os.path.join(outputFolder,'time%06d_1.tif' %time))
    sitk.WriteImage(tim2,os.path.join(outputFolder,'time%06d_2.tif' %time))

    sitk.WriteImage(sitk.MaximumProjection(tim1,2),os.path.join(outputFolder,'time%06d_proj2_1.tif' %time))
    sitk.WriteImage(sitk.MaximumProjection(tim2,2),os.path.join(outputFolder,'time%06d_proj2_2.tif' %time))

    sitk.WriteImage((tim1+tim2)/2,os.path.join(outputFolder,'time%06d_average.tif' %time))
    sitk.WriteImage(sitk.MaximumProjection((tim1+tim2)/2,2),os.path.join(outputFolder,'time%06d_proj2_average.tif' %time))

# write script call

tmpFile = open(os.path.join(outputFolder,'call.txt'),'w')
tmpFile.write(str(sys.argv))
tmpFile.close()