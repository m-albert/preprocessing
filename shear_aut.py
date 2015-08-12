__author__ = 'malbert'

from imports import *
import re
#####################
#####################
#####################
# process arguments

print '\n\n\ninput args:\n%s\n\n\n' %sys.argv

inFilePattern1                                  = sys.argv[-10] #expects %05d...
inFilePattern2                                  = sys.argv[-9]
startTime                                       = int(sys.argv[-8])
endTime                                         = int(sys.argv[-7])
mexicanHatFilterSize                            = float(sys.argv[-6])
inFile1                                         = sys.argv[-5]
inFile2                                         = sys.argv[-4]
isotropicSpacingFactor                          = float(sys.argv[-3]) # multiply xy spacing by this factor for obtaining
                                                                      # final isotropic spacing
outputFolder                                    = sys.argv[-2]
option                                          = int(sys.argv[-1])



def getMeta(fileName):

    md = h5py.File(fileName,'r').attrs.items()[2][1]
    ss = md[md.find('\"Stack\"'):]

    # pdb.set_trace()
    z1 = float(re.search("\"Z1\":(-?[0-9]+),",ss).groups()[0])
    z2 = float(re.search("\"Z2\":(-?[0-9]+),",ss).groups()[0])
    N = float(re.search("\"N\":(-?[0-9]+)}",ss).groups()[0])

    leftOrRight = re.search("Selected View::(\\D+),",md).groups()[0]


    if leftOrRight == 'Left': leftOrRight = 1
    else: leftOrRight = 0

    if leftOrRight==1:
        if n.sign(z2-z1)>0:
            direction = 1.
        else:
            direction = -1.
    else:
        if n.sign(z2-z1)>0:
            direction = -1.
        else:
            direction = 1.

    diagdz = (z2-z1)/float(N-1)/0.65

    return z1,z2,diagdz,direction

def getShearParams(direction):
    shearParams = n.zeros(12,dtype=n.float)
    tmpShear = transformations.shear_matrix(n.pi/4.,[direction,0,0],[0,0,0],[0,0,1])
    shearParams[:9] = n.array(tmpShear[:3,:3].flatten())
    shearParams[9:] = n.zeros(3)
    return shearParams




metab1 = getMeta(inFile1)
metab2 = getMeta(inFile2)

direction1 = getMeta(inFile1)[3]
direction2 = getMeta(inFile2)[3]
diagdz = getMeta(inFile2)[2]

dz = diagdz*n.cos(n.pi/4.)
outSpacing = n.array([isotropicSpacingFactor,isotropicSpacingFactor,isotropicSpacingFactor/dz])

# stack = re.search("\"Stack\":{\"X\":(\\d+),\"Y\":(\\d+),\"Z1\":(\\d+),\"Z2\":(\\d+),\"R\":(\\d+),\"N\":(\\d+)}",
#                   ss)

print 'extracted from metadat:\ndiagdz=%s\ndz=%s\ndirection1=%s\ndirection2=%s\nassuming: dx/dy=0.65um' %(diagdz,dz,direction1,direction2)


im1 = n.array(h5py.File(inFile1,'r')['Data'])
im2 = n.array(h5py.File(inFile2,'r')['Data'])
im1 = sitk.gifa(im1)
im2 = sitk.gifa(im2)

im1.SetSpacing([1,1,dz])
im2.SetSpacing([1,1,dz])

tmpShape = n.array(im1.GetSize())
outShape = tmpShape/outSpacing

shearParams1 = getShearParams(direction1)
shearParams2 = getShearParams(direction2)

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

    metad1 = getMeta(inFile1)
    metad2 = getMeta(inFile2)

    im1 = n.array(h5py.File(inFile1,'r')['Data'])
    im2 = n.array(h5py.File(inFile2,'r')['Data'])
    im1 = sitk.gifa(im1)
    im2 = sitk.gifa(im2)

    im1.SetSpacing([1,1,diagdz])
    im2.SetSpacing([1,1,diagdz])

    shearParams1 = getShearParams(metad1[3])
    shearParams2 = getShearParams(metad2[3])

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