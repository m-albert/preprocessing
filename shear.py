__author__ = 'malbert'

__author__ = 'malbert'

from imports import *
import re
#####################
#####################
#####################
# process arguments

print '\n\n\ninput args:\n%s\n\n\n' %sys.argv

inFilePattern1                                  = sys.argv[-11] #expects %05d...
inFilePattern2                                  = sys.argv[-10]
startTime                                       = int(sys.argv[-9])
endTime                                         = int(sys.argv[-8])
stepTime                                        = int(sys.argv[-7])
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

    diagdz = (z2-z1)/float(N-1)/6.5*30.
    dz = diagdz*n.cos(n.pi/4.)

    if direction == -1:
        deltaOrigin = n.zeros(3)
    else:
        deltaOrigin = n.array([(N-1)*diagdz*n.sin(n.pi/4.)*(-1)*direction,0,0.])
    deltaShape  = n.abs(n.array([(N-1)*diagdz*n.sin(n.pi/4.)*direction,0,0.])).astype(n.uint16)#/6.5*30

    return z1,z2,diagdz,direction,dz,deltaOrigin,deltaShape

def getShearParams(direction):
    shearParams = n.zeros(12,dtype=n.float)
    tmpShear = transformations.shear_matrix(n.pi/4.,[direction,0,0],[0,0,0],[0,0,1])
    shearParams[:9] = n.array(tmpShear[:3,:3].flatten())
    shearParams[9:] = n.zeros(3)
    return shearParams


metab1 = getMeta(inFile1)
metab2 = getMeta(inFile2)

diagdz = metab1[2]
dz = metab1[4]



print 'extracted from metadat:\ndiagdz=%s\ndz=%s\ndirection1=%s\ndirection2=%s\nassuming: dx/dy=0.65um' %(diagdz,dz,metab1[3],metab2[3])


im1 = n.array(h5py.File(inFile1,'r')['Data'])
im2 = n.array(h5py.File(inFile2,'r')['Data'])
im1 = sitk.gifa(im1)
im2 = sitk.gifa(im2)

im1.SetSpacing([1,1,metab1[4]])
im2.SetSpacing([1,1,metab2[4]])

outSpacing = n.array([isotropicSpacingFactor,isotropicSpacingFactor,isotropicSpacingFactor/dz])
tmpShape = n.array(im1.GetSize())
outShape = tmpShape/outSpacing

shearParams1 = getShearParams(metab1[3])
shearParams2 = getShearParams(metab2[3])



tim1 = beads.transformStack(shearParams1,im1,outShape=n.array(im1.GetSize())+metab1[6],outOrigin=n.array(im1.GetOrigin())+metab1[5])
tim2 = beads.transformStack(shearParams2,im2,outShape=n.array(im2.GetSize())+metab2[6],outOrigin=n.array(im2.GetOrigin())+metab2[5])

if metab1[1]<metab1[0]:
    tim1 = tim1[:,:,::-1]
    tim1.SetDirection([1,0,0,0,1,0,0,0,1])
    tmpz1 = metab1[0]
    metab1[0] = metab1[1]
    metab1[1] = tmpz1

if metab2[1]<metab2[0]:
    tim2 = tim2[:,:,::-1]
    tim2.SetDirection([1,0,0,0,1,0,0,0,1])
    tmpz1 = metab2[0]
    metab2[0] = metab2[1]
    metab2[1] = tmpz1

tim1.SetSpacing([1.,1,1])
tim2.SetSpacing([1.,1,1])
tim1.SetOrigin([0.,0,0])
tim2.SetOrigin([0.,0,0])

if option in [1]:

    segmentedBeads = beads.segmentBeadsNils([tim1,tim2],mexicanHatFilterSize)
    params = list(beads.findParametersNils(segmentedBeads))

    tmpParams = fusion.invertParams(params[1])
    #tmpParams,metab1,metab2 = pickle.load(open('/data/malbert/tmp/nilsfusion/beadResults.pc','r'))
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



for time in range(startTime,endTime+1,stepTime):

    inFile1 = inFilePattern1 %time
    inFile2 = inFilePattern2 %time

    metad1 = getMeta(inFile1)
    metad2 = getMeta(inFile2)

    im1 = n.array(h5py.File(inFile1,'r')['Data'])
    im2 = n.array(h5py.File(inFile2,'r')['Data'])
    im1 = sitk.gifa(im1)
    im2 = sitk.gifa(im2)

    # pdb.set_trace()

    outSpacing = n.array([isotropicSpacingFactor,isotropicSpacingFactor,isotropicSpacingFactor/metad1[4]])
    tmpShape = n.array(im1.GetSize())
    outShape = tmpShape/outSpacing

    im1.SetSpacing([1,1,metad1[4]])
    im2.SetSpacing([1,1,metad2[4]])

    shearParams1 = getShearParams(metad1[3])
    shearParams2 = getShearParams(metad2[3])

    tim1 = beads.transformStack(shearParams1,im1,outShape=n.array(im1.GetSize())+metad1[6],outOrigin=n.array(im1.GetOrigin())+metad1[5])
    tim2 = beads.transformStack(shearParams2,im2,outShape=n.array(im2.GetSize())+metad2[6],outOrigin=n.array(im2.GetOrigin())+metad2[5])

    if metad1[1]<metad1[0]:
        tim1 = tim1[:,:,::-1]
        tim1.SetDirection([1,0,0,0,1,0,0,0,1])
        tmpz1 = metad1[0]
        metad1[0] = metad1[1]
        metad1[1] = tmpz1

    if metad2[1]<metad2[0]:
        tim2 = tim2[:,:,::-1]
        tim2.SetDirection([1,0,0,0,1,0,0,0,1])
        tmpz1 = metad2[0]
        metad2[0] = metad2[1]
        metad2[1] = tmpz1

    tim1.SetSpacing([1.,1,metad1[4]/metab1[4]])
    tim2.SetSpacing([1.,1,metad2[4]/metab2[4]])
    tim1.SetOrigin([0.,0,metad1[0]-metab1[0]])
    tim2.SetOrigin([0.,0,metad2[0]-metab2[0]])

    outOrigin = tim1.GetOrigin()

    if option == 1:
        #tmpParams = fusion.invertParams(params[1])
        print tmpParams

        tim1 = beads.transformStack([1.,0,0,0,1,0,0,0,1,0,0,0],tim1,outSpacing=outSpacing,outShape=outShape,outOrigin=outOrigin)
        tim2 = beads.transformStack(tmpParams,tim2,outSpacing=outSpacing,outShape=outShape,outOrigin=outOrigin)
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
