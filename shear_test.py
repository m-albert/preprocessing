__author__ = 'malbert'

from imports import *
import re
from shearing import *

#####################
#####################
#####################
# process arguments

lastArg = sys.argv[-1]
if lastArg[-2:] == 'py':
    configurationFile = lastArg
else:
    print '\nWrong arguments:\n\n Usage: python shear.py <configurationFile.py>\n '
    sys.exit()

#configurationFile = sys.argv[-1] # use file passed as argument
print '\n\nstarting preprocessing pipeline using configuration file:\n\n    %s\n' %configurationFile

#####################
#####################
#####################
# read and handle configuration file
configFile = open(configurationFile)
configCode = configFile.read()
configCode = configCode.replace('\r','')
configFile.close()
if configCode[:39] == '# dual view SPIM preprocessing configuration file':
    exec(configCode)
else:
    print '\n\n Error: specified configuration file \n\n %s \n\n is not a proper configuration file (first line must state: \'# SPIM preprocessing configuration file\')\n\n' %configurationFile
    sys.exit()

#####################
#####################
#####################
# process arguments

inFilePattern1                                  = config['inFilePattern1']
inFilePattern2                                  = config['inFilePattern2']
startTime                                       = config['startTime']
endTime                                         = config['endTime']
stepTime                                        = config['stepTime']
mexicanHatFilterSize                            = config['mexicanHatFilterSize']
inFile1                                         = config['inFile1']
inFile2                                         = config['inFile2']
isotropicSpacingFactor                          = config['isotropicSpacingFactor'] # multiply xy spacing by this factor for obtaining
                                                                      # final isotropic spacing
outputFolder                                    = config['outputFolder']
option                                          = config['option']
tenXObjective                                   = config['tenXObjective']
beadParams                                      = config['beadParams']


metab1 = getMeta(inFile1,tenXObjective)
metab2 = getMeta(inFile2,tenXObjective)

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

if option in [1,2]:

    if config['beadParams'] is None:
        segmentedBeads = beads.segmentBeadsNils([tim1,tim2],mexicanHatFilterSize)
        params = list(beads.findParametersNils(segmentedBeads,dz=dz))
    else:
        params = config['beadParams']

    tmpParams = fusion.invertParams(params[1])
    #tmpParams,metab1,metab2 = pickle.load(open('/data/malbert/tmp/nilsfusion/beadResults.pc','r'))
    print tmpParams
    tim1 = beads.transformStack([1.,0,0,0,1,0,0,0,1,0,0,0],tim1,outSpacing=outSpacing,outShape=outShape)
    tim2 = beads.transformStack(tmpParams,tim2,outSpacing=outSpacing,outShape=outShape)

tim1.SetOrigin((0,0,0))
tim2.SetOrigin((0,0,0))
tim1.SetSpacing((1,1,1))
tim2.SetSpacing((1,1,1))


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

    metad1 = getMeta(inFile1,tenXObjective)
    metad2 = getMeta(inFile2,tenXObjective)

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

    if option == 1:

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

        #tmpParams = fusion.invertParams(params[1])
        print tmpParams

        tim1 = beads.transformStack([1.,0,0,0,1,0,0,0,1,0,0,0],tim1,outSpacing=outSpacing,outShape=outShape,outOrigin=outOrigin)
        tim2 = beads.transformStack(tmpParams,tim2,outSpacing=outSpacing,outShape=outShape,outOrigin=outOrigin)
        # tim2 = beads.transformStack(tmpParams,tim2)

    # tim1.SetOrigin((0,0,0))
    # tim2.SetOrigin((0,0,0))
    # tim1.SetSpacing((1,1,1))
    # tim2.SetSpacing((1,1,1))

    sitk.WriteImage(tim1,os.path.join(outputFolder,'time%06d_1.tif' %time))
    sitk.WriteImage(tim2,os.path.join(outputFolder,'time%06d_2.tif' %time))

    sitk.WriteImage(sitk.MaximumProjection(tim1,2),os.path.join(outputFolder,'time%06d_proj2_1.tif' %time))
    sitk.WriteImage(sitk.MaximumProjection(tim2,2),os.path.join(outputFolder,'time%06d_proj2_2.tif' %time))

    # sitk.WriteImage((tim1+tim2)/2,os.path.join(outputFolder,'time%06d_average.tif' %time))
    # sitk.WriteImage(sitk.MaximumProjection((tim1+tim2)/2,2),os.path.join(outputFolder,'time%06d_proj2_average.tif' %time))

# write script call

tmpFile = open(os.path.join(outputFolder,'call.txt'),'w')
tmpFile.write(str(sys.argv))
tmpFile.close()


if option == 2:
    parameters = [[1.,0,0,0,1,0,0,0,1,0,0,0],tmpParams]
    produceXML([tim1,tim2],parameters)