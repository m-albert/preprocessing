from imports import *

import beads
import transformations
import filters

def calculateOpposingWeights(size,sigmoidHalfWidthRelativeToStackSize=0.05):
    print 'calculating weights for opposing stacks fusion'

    size = (size[0],size[1],size[2]+2)
    halfY = 0.9
    nZPlanes = size[2]
    sigmoidHalfWidthInPixels = nZPlanes*sigmoidHalfWidthRelativeToStackSize
    halfY = 0.9 # half width is defined as the x value corresponding to y=halfY
    sigmoidA = n.log(1./halfY-1.)/sigmoidHalfWidthInPixels

    zs = n.linspace(0,size[2]-1,size[2])
    zsigmoid1 = n.abs(1/(1+n.exp(sigmoidA*(zs-size[2]/2.))))
    #zsigmoid2 = n.abs(1/(1+n.exp(sigmoidA*(-zs+size[2]/2.))))
    zyplane = n.ones(size[1:])
    zyplane1 = zyplane*zsigmoid1
    #zyplane2 = zyplane*zsigmoid2
    zyplane1 = zyplane1.astype(n.float32)
    #zyplane2 = zyplane2.astype(n.float32)
    w1 = (n.ones(size,dtype=n.uint16)*zyplane1).swapaxes(0,1).swapaxes(1,2).swapaxes(0,1)
    #w2 = (n.ones(size,dtype=n.uint16)*zyplane2).swapaxes(0,1).swapaxes(1,2).swapaxes(0,1)
    w1 = sitk.GetImageFromArray(w1)
    w2 = beads.transformStack([1.,0,0,0,1,0,0,0,-1,0,0,size[2]-2],w1,outShape=(size[0],size[1],size[2]-2))
    w1 = beads.transformStack([1.,0,0,0,1,0,0,0,1,0,0,2],w1,outShape=(size[0],size[1],size[2]-2))
    return [w1,w2]


def calculateOpposingWeightsOld(size,sigmoidHalfWidthRelativeToStackSize=0.05):
    print 'calculating weights for opposing stacks fusion'

    halfY = 0.9
    nZPlanes = size[2]
    sigmoidHalfWidthInPixels = nZPlanes*sigmoidHalfWidthRelativeToStackSize
    halfY = 0.9 # half width is defined as the x value corresponding to y=halfY
    sigmoidA = n.log(1./halfY-1.)/sigmoidHalfWidthInPixels

    zs = n.linspace(0,size[2]-1,size[2])
    zsigmoid1 = n.abs(1/(1+n.exp(sigmoidA*(zs-size[2]/2.))))
    zsigmoid2 = n.abs(1/(1+n.exp(sigmoidA*(-zs+size[2]/2.))))
    zyplane = n.ones(size[1:])
    zyplane1 = zyplane*zsigmoid1
    zyplane2 = zyplane*zsigmoid2
    zyplane1 = zyplane1.astype(n.float32)
    zyplane2 = zyplane2.astype(n.float32)
    w1 = (n.ones(size,dtype=n.uint16)*zyplane1).swapaxes(0,1).swapaxes(1,2).swapaxes(0,1)
    w1 = sitk.GetImageFromArray(w1)
    w2 = (n.ones(size,dtype=n.uint16)*zyplane2).swapaxes(0,1).swapaxes(1,2).swapaxes(0,1)
    w2 = sitk.GetImageFromArray(w2)
    return [w1,w2]

def fuseOpposingStacks(p1,p2,s1,s2,axis=0):
    s1 = beads.transformStack(invertParams(p1),s1)
    s2 = beads.transformStack(invertParams(p2),s2)
    size = s1.GetSize()
    s1 = sitk.GetArrayFromImage(s1)
    s2 = sitk.GetArrayFromImage(s2)

    # halfY = 0.9
    nZPlanes = size[2]
    sigmoidHalfWidthInPixels = nZPlanes*0.05
    halfY = 0.9 # half width is defined as the x value corresponding to y=halfY
    sigmoidA = n.log(1./halfY-1.)/sigmoidHalfWidthInPixels

    length = s1.shape[axis]
    s1 = s1.swapaxes(0,axis)
    s2 = s2.swapaxes(0,axis)
    for z in range(length):
        w1 = n.abs(1/(1+n.exp(sigmoidA*(z-size[2]/2.))))
        w2 = n.abs(1/(1+n.exp(sigmoidA*(-z+size[2]/2.))))
        s1[z,:,:] = (w1*s1[z,:,:]+w2*s2[z,:,:])/(w1+w2)
    return sitk.GetImageFromArray(s1.swapaxes(0,axis))

def fuseOpposingStacksGus(s1,s2,axis=0):
    # s1 = beads.transformStack(invertParams(p1),s1)
    # s2 = beads.transformStack(invertParams(p2),s2)
    size = s1.GetSize()
    s1 = sitk.GetArrayFromImage(s1)
    s2 = sitk.GetArrayFromImage(s2)

    # halfY = 0.9
    nZPlanes = size[2]
    sigmoidHalfWidthInPixels = nZPlanes*0.05
    halfY = 0.9 # half width is defined as the x value corresponding to y=halfY
    sigmoidA = n.log(1./halfY-1.)/sigmoidHalfWidthInPixels

    length = s1.shape[axis]
    s1 = s1.swapaxes(0,axis)
    s2 = s2.swapaxes(0,axis)
    for z in range(length):
        w1 = n.abs(1/(1+n.exp(sigmoidA*(z-size[2]/2.))))
        w2 = n.abs(1/(1+n.exp(sigmoidA*(-z+size[2]/2.))))
        s1[z,:,:] = (w1*s1[z,:,:]+w2*s2[z,:,:])/(w1+w2)
        s1[z,:,:] = (w1*s1[z,:,:])
        s2[z,:,:] = (w2*s2[z,:,:])

    return sitk.GetImageFromArray(s1.swapaxes(0,axis)),sitk.GetImageFromArray(s2.swapaxes(0,axis))

def fuseOpposingStacksOld(p1,p2,s1,s2,normalizeIntensities=True,axis=0):
    print 'fusing (normalize intensities %s)' %(['OFF','ON'][int(normalizeIntensities)])
    s1 = beads.transformStack(invertParams(p1),s1)
    s2 = beads.transformStack(invertParams(p2),s2)
    s1 = sitk.GetArrayFromImage(s1)
    s2 = sitk.GetArrayFromImage(s2)
    if normalizeIntensities:
        n1,n2 = n.mean(s1),n.mean(s2)
    else:
        n1,n2 = 1.,1.
    length = s1.shape[axis]
    s1 = s1.swapaxes(0,axis)
    s2 = s2.swapaxes(0,axis)
    for z in range(length):
        w1 = 1/(1+n.exp(-40./length*(z-length/2.)))
        w2 = 1/(1+n.exp(40./length*(z-length/2.)))
        s1[z,:,:] = (w1*s1[z,:,:]+w2*s2[z,:,:])/(w1+w2)
    return sitk.GetImageFromArray(s1.swapaxes(0,axis))


def fuseStacks(s,weights,p=None):
    if p != None:
        for iis in range(len(s)):
            s[iis] = beads.transformStack(p[iis],s[iis],outShape=s[0].GetSize())
    # pdb.set_trace()
    if weights is None:
        print 'fusing stacks by averaging'
        N = len(s)
        tmp = sitk.Cast(s[0],6)/N
        for i in range(N-1):
            tmp += sitk.Cast(s[i+1],6)/N
        return sitk.Cast(tmp,3)
    print 'fusing stacks by using weight matrices'
    for iweight in range(len(weights)):
            weights[iweight].SetSpacing(s[iweight].GetSpacing())
            weights[iweight].SetOrigin(s[iweight].GetOrigin())
    #for iweight in range(len(weights)):    
        #weights[iweight] = weights[iweight]*sitk.Cast(s[iweight],6)
    result = weights[0]*sitk.Cast(s[0],6)
    for iweight in range(1,len(weights)):
        result += weights[iweight]*sitk.Cast(s[iweight],6)
    result = sitk.Cast(result,3)
    return result

def calculateFinalWeightMatricesOld(params,size,spacings,isotropicSpacing,outShape=None,outSpacing=None,sigmoidHalfWidthRelativeToImageSize=0.05):

    size = n.ceil(size).astype(n.uint16)

    def calcAngles(Y,Z):
        norm = n.sqrt(n.power(Y,2)+n.power(Z,2))
        norm[norm==0] = 0.001
        angles = []
        angles.append(n.arccos(Y/norm))
        angles.append(2*n.pi-n.arccos(Y/norm))
        angles.append(n.pi-n.arccos(-Y/norm))
        angles.append(n.pi+n.arccos(-Y/norm))
        finalAngles = n.zeros(Y.shape,dtype=n.float64)
        criteria = [(Y>=0)*(Z>=0),(Y>=0)*(Z<0),(Y<0)*(Z>=0),(Y<0)*(Z<0)]
        for icrit,crit in enumerate(criteria):
            finalAngles[crit] = angles[icrit][crit]
        return finalAngles

    center = n.array(size)[1:][::-1]/2.
    bestAngles = []
    for iparam,param in enumerate(params):
        point = n.array([[0,size[1]/2.,size[0]/2.]])#,size[1]/2.,0.]])
        point[0] = point[0]*spacings[iparam]
        refPoint = beads.backTransform(param,point)[0][1:]
        print refPoint
        refPoint = refPoint-center
        refPoint = refPoint/n.linalg.norm(refPoint)
        angle = calcAngles(n.array([refPoint[0]]),n.array([refPoint[1]]))
        bestAngles.append(angle)

    Y,Z = n.mgrid[0:size[2],0:size[1]]
    Y = Y-center[0]
    Z = Z-center[1]

    norm = n.sqrt(n.power(Y,2)+n.power(Z,2))

    sigmoidHalfWidthInPixels = n.mean(size[1:])*sigmoidHalfWidthRelativeToImageSize
    halfY = 0.9 # half width is defined as the x value corresponding to y=halfY
    norm[norm==0] = 0.0001
    sigmoidA = n.log(1./halfY-1.)/n.arctan(sigmoidHalfWidthInPixels/norm)
    
    finalAngles = calcAngles(Y,Z)

    # reference angles
    viewAngles = []
    for i in range(len(bestAngles)):
        viewAngle = []
        viewAngle.append(bestAngles[i])
        viewAngle.append(bestAngles[i]+2*n.pi)
        viewAngle.append(bestAngles[i]-2*n.pi)
        viewAngles.append(viewAngle)

    # calculate difference to ref angles
    weights = []
    for viewAngle in viewAngles:
        tmp = n.min([n.abs(finalAngles-i)%(n.pi) for i in viewAngle],0)
        tmpCrit = tmp>=n.pi/2.
        tmp[tmpCrit] = tmp[tmpCrit] - n.pi/2.
        weights.append(tmp)

    weightsDiff = []
    for iview,view in enumerate(viewAngles):
        tmpDiff = n.min([-weights[iview]+weights[i] for i in n.delete(range(len(viewAngles)),iview,0)],0)
        weightsDiff.append(1/(1+n.exp(-sigmoidA*(tmpDiff))))

    weights = []
    weightsDiffSum = n.sum(weightsDiff,0)
    for iview,view in enumerate(viewAngles):
        weights.append(weightsDiff[iview]/weightsDiffSum)

    finalWeights = []
    w = n.ones(size[::-1],dtype=n.uint16).swapaxes(0,2).swapaxes(1,2)#*w1).swapaxes(0,2).swapaxes(0,1)
    for i in range(len(weights)):
        tmp = (w*weights[i]).swapaxes(0,2).swapaxes(0,1)
        tmp = sitk.GetImageFromArray(tmp)
        finalWeights.append(sitk.Cast(tmp,6))

    return finalWeights


def calculateFinalWeightMatrices(params,
                                 stacks,
                                 axisOfRotation=1,
                                 sigmoidHalfWidthRelativeToImageSize=0.05,
                                 adaptiveCenterFile=None,
                                 sampleIntensityThreshold=15
                                 ):
    # old function arguments:
    # params,size,spacings,isotropicSpacing,outShape=None,outSpacing=None,sigmoidHalfWidthRelativeToImageSize=0.05)
    """

    :rtype : object
    - expects params and size, axisOfRotation, adaptiveCenterFile(sitk) relative to axes definition: (x,y,z)
    - rotation around x
    """

    #origins = infoDict['origins']
    #positions = infoDict['positions']
    #spacing = infoDict['spacing']
    #centerOfRotation = infoDict['centerOfRotation']
    #axisOfRotation = infoDict['axisOfRotation']
    #sizes = infoDict['sizes']

    params = n.array(params)
    finalSize = n.ceil(stacks[0].GetSize()).astype(n.int64)
    # finalSpacing = n.array(stacks[0].GetSpacing()).astype(n.float64)

    # tmpMatrix = n.zeros((4,4)).astype(n.float)
    # tmpMatrix[:3,:3] = params[1][:9].reshape((3,3))
    # direction = transformations.rotation_from_matrix(tmpMatrix)[1]
    axisOfRotation = n.argmax(n.abs(params[1][[0,4,8]]))
    print 'determining axis of rotation from parameters... %s' %axisOfRotation


    axes = range(3)
    axes.remove(axisOfRotation)
    firstAxis = axes[0]
    lastAxis = axes[1]

    def calcAngles(Y,Z):
        norm = n.sqrt(n.power(Y,2)+n.power(Z,2))
        norm[norm==0] = 0.001
        angles = []
        angles.append(n.arccos(Y/norm))
        angles.append(2*n.pi-n.arccos(Y/norm))
        #angles.append(n.pi-n.arccos(-Y/norm))
        #angles.append(n.pi+n.arccos(-Y/norm))
        finalAngles = n.zeros(Y.shape,dtype=n.float64)
        #criteria = [(Y>=0)*(Z>=0),(Y>=0)*(Z<0),(Y<0)*(Z>=0),(Y<0)*(Z<0)]
        criteria = [(Z>=0),(Z<0)]
        for icrit,crit in enumerate(criteria):
            finalAngles[crit] = angles[icrit][crit]
        return finalAngles

    if adaptiveCenterFile == None:
        centers = n.array([[finalSize[lastAxis],finalSize[firstAxis]]])/2.
    else:
        print 'calculating stack centers for fusion weights (assuming sample intensity of >%s)' %sampleIntensityThreshold
        #adaptiveCenterFileSpacing = n.array(adaptiveCenterFile.GetSpacing())
        centers = []
        Z, Y = n.mgrid[0:finalSize[lastAxis], 0:finalSize[firstAxis]]
        Z, Y = sitk.GetImageFromArray(Y.astype(n.uint16)), sitk.GetImageFromArray(Z.astype(n.uint16))
        for x in range(adaptiveCenterFile.GetSize()[axisOfRotation]):
            if axisOfRotation == 0:
                xImage = adaptiveCenterFile[x]
            elif axisOfRotation == 1:
                xImage = adaptiveCenterFile[:,x]
            elif axisOfRotation == 2:
                xImage = adaptiveCenterFile[:,:,x]
            #pdb.set_trace()
            Z.SetOrigin(xImage.GetOrigin())
            Z.SetSpacing(xImage.GetSpacing())
            Y.SetOrigin(xImage.GetOrigin())
            Y.SetSpacing(xImage.GetSpacing())
            xImage = sitk.Cast(xImage>sampleIntensityThreshold,7)*1000+0.01 #careful, could create overflow
            totalIntensity = sitk.SumProjection(sitk.SumProjection(xImage,0),1).GetPixel(0,0)
            z = sitk.SumProjection(sitk.SumProjection(sitk.Cast(Z,7)*xImage,0),1).GetPixel(0,0)/totalIntensity
            y = sitk.SumProjection(sitk.SumProjection(sitk.Cast(Y,7)*xImage,0),1).GetPixel(0,0)/totalIntensity
            # totalIntensity = sitk.SumProjection(sitk.SumProjection(sitk.GetImageFromArray(n.ones(xImage.GetSize())),0),1).GetPixel(0,0)
            #pdb.set_trace()
            # z = sitk.SumProjection(sitk.SumProjection(Z,0),1).GetPixel(0,0)/totalIntensity
            # y = sitk.SumProjection(sitk.SumProjection(Y,0),1).GetPixel(0,0)/totalIntensity
            # tmp = sitk.GetArrayFromImage(xImage)
            # tmp[int(y),int(z)] = 2000
            # plt.imshow(tmp)
            # plt.colorbar()
            # plt.show()
            center = [y,z]
            centers.append(center)

    # smoothen center line
    centers = n.array(centers)
    #centers = n.array([ndimage.gaussian_filter(centers[:,i],finalSize[[firstAxis,lastAxis][i]]/30.) for i in range(2)])
    #centers = centers.swapaxes(0,1)

    # get angles
    refCenter = n.array(stacks[0].TransformIndexToPhysicalPoint([int(i) for i in n.array(stacks[0].GetSize())/2.]))
    #origCenter = n.array([refCenter[i] for i in [lastAxis,firstAxis]])

    bestAngles = []
    refPoints = []
    print 'refCenter %s' %refCenter
    refPoints.append(refCenter)
    for iparam,param in enumerate(params):
        # point2 is refCenter in rotated view
        point2 = beads.transform(param,[refCenter])[0]
        # shift point2 in z (positive here since later compared to positive again)
        point2[lastAxis] += 100
        # transform back to reference view
        refPoint = beads.backTransform(param,[point2])[0]
        refPoints.append(refPoint)
        # calc direction vector
        refVector = refPoint-refCenter
        refVector = n.array([refVector[i] for i in [lastAxis,firstAxis]])
        refVector = refVector/n.linalg.norm(refVector)
        # calc angle with respect to [1,1]
        angle = calcAngles(n.array([refVector[0]]),n.array([refVector[1]]))
        bestAngles.append(angle)

    #pdb.set_trace()
    refPoints = n.array(refPoints)
    # plt.plot(refPoints[:,0],refPoints[:,2],'o')
    # plt.show()

    zCenter = finalSize[[lastAxis,firstAxis]]/2.
    lower = n.array([0,0])
    upper = finalSize[[lastAxis,firstAxis]]
    templateLower = lower - n.max(centers,0)+zCenter
    templateUpper = upper - n.min(centers,0)+zCenter
    print templateLower,templateUpper
    deltaplus = n.ceil(templateUpper-zCenter)
    deltaminus = n.ceil(zCenter-templateLower)
    templateSize = deltaminus+deltaplus

    Y,Z = n.mgrid[0:templateSize[0],0:templateSize[1]]
    Y = Y-deltaminus[0]
    Z = Z-deltaminus[1]
    norm = n.sqrt(n.power(Y,2)+n.power(Z,2))

    sigmoidHalfWidthInPixels = n.mean(finalSize[[firstAxis,lastAxis]])*sigmoidHalfWidthRelativeToImageSize
    halfY = 0.9 # half width is defined as the x value corresponding to y=halfY
    norm[norm==0] = 0.0001
    sigmoidA = n.log(1./halfY-1.)/n.arctan(sigmoidHalfWidthInPixels/norm)

    finalAngles = calcAngles(Y,Z)
    # plt.imshow(finalAngles)
    # plt.title('finalAngles')
    # plt.colorbar()
    # plt.show()

    # add cameras here
    viewAngles = []
    for i in range(len(bestAngles)):
        viewAngle = []
        viewAngle.append(bestAngles[i])
        viewAngle.append((bestAngles[i]+n.pi)%(2*n.pi))
        viewAngles.append(viewAngle)

    print viewAngles
    # calculate difference to ref angles
    weights = []
    for iview,viewAngle in enumerate(viewAngles):
        tmp = n.min([n.min([n.abs(finalAngles-i),n.abs(2*n.pi-n.abs(finalAngles-i))],0) for i in viewAngle],0)
        weights.append(tmp.astype(n.float32))
        # plt.imshow(weights[iview])
        # plt.title('weights1 %s' %iview)
        # plt.colorbar()
        # plt.show()
        # #tifffile.imsave('/home/malbert/delme/test%s.tif' %(iview),weights[iview])

    weightsDiff = []
    for iview,view in enumerate(viewAngles):
        weightsDiff.append(n.min([-weights[iview]+weights[i] for i in n.delete(range(len(viewAngles)),iview,0)],0))
        #tifffile.imsave('/home/malbert/delme/test%s.tif' %(iview),weightsDiff[iview])
        #sitk.WriteImage(sitk.Cast(sitk.GetImageFromArray(weightsDiff[iview]),3),'/home/malbert/delme/test%s.tif' %(iview))
        # plt.imshow(weightsDiff[iview])
        # plt.title('weightsDiff1 %s' %iview)
        # plt.colorbar()
        # plt.show()
    for iview,view in enumerate(viewAngles):
        weightsDiff[iview] = (1/(1+n.exp(sigmoidA*(weightsDiff[iview]))))
        # plt.imshow(weightsDiff[iview])
        # plt.title('weightsDiff2 %s' %iview)
        # plt.colorbar()
        # plt.show()
        #pdb.set_trace()
        #sitk.WriteImage(sitk.Cast(sitk.GetImageFromArray(weightsDiff[iview]*100),3),'/home/malbert/delme/test%s.tif' %(iview))

    weights = []
    weightsDiffSum = n.sum(weightsDiff,0)
    for iview,view in enumerate(viewAngles):
        weights.append(weightsDiff[iview]/weightsDiffSum)
        #sitk.WriteImage(sitk.Cast(sitk.GetImageFromArray(weights[iview]*1000),3),'/home/malbert/delme/test%s.tif' %(iview))

    #sitk.WriteImage(sitk.Cast(finalWeights[iview]*1000,3),'/home/malbert/delme/test0.tif')

    xWeights = []
    for icenter, center in enumerate(centers):
        xWeight = []
        for iview,view in enumerate(viewAngles):
            # pdb.set_trace()
            # plt.imshow(weights[iview])
            # plt.title('test')
            # plt.colorbar()
            # plt.show()
            xW = sitk.GetImageFromArray(weights[iview])
            #xW.SetOrigin((templateLower[0], templateLower[1]))
            xW.SetOrigin([0,0])
            tmpTransform = sitk.Transform(2,1)
            tmpSize = finalSize[[lastAxis,firstAxis]].astype(n.int64)
            tmpSize = [int(i) for i in tmpSize[::-1]]
            tmpOrigin = n.array([-center[i]+deltaminus[i] for i in [1,0]])
            #pdb.set_trace()
            tmpTransform.SetParameters(tmpOrigin)
            #pdb.set_trace()
            xW = sitk.Resample(xW,
                               tmpSize,
                               tmpTransform,
                               sitk.sitkLinear,
                               #tmpOrigin
                               )
            #pdb.set_trace()
            xW = sitk.GetArrayFromImage(xW)
            xWeight.append(xW)
            # plt.imshow(xW)
            # plt.title('weights %s %s' %(iview,center))
            # plt.colorbar()
            # plt.show()

        xWeights.append(xWeight)


    finalWeights = []
    for i in range(len(params)):
        #finalWeights.append(sitk.GetImageFromArray(xWeights[0][i]))
        if adaptiveCenterFile is None:
            tmp = n.array([xWeights[0][i]]*finalSize[axisOfRotation])
        else:
            tmp = n.array([xWeights[j][i] for j in range(finalSize[axisOfRotation])])
        if axisOfRotation == 0: numpyAxisOfRotation = 2
        elif axisOfRotation == 1: numpyAxisOfRotation = 1
        elif axisOfRotation == 2: numpyAxisOfRotation = 0
        tmp = tmp.swapaxes(0,numpyAxisOfRotation)
        if numpyAxisOfRotation == 2:
            otherAxes = range(3)
            otherAxes.remove(numpyAxisOfRotation)
            tmp = tmp.swapaxes(otherAxes[0],otherAxes[1])
        #tmp = tmp.swapaxes(0, 1).swapaxes(1, 2).astype(n.float32)
        tmp = tmp.astype(n.float32)
        finalWeights.append(sitk.GetImageFromArray(tmp))
    #sitk.WriteImage(sitk.Cast(finalWeights[iview]*1000,3),'/home/malbert/delme/test.tif')

    return finalWeights



def invertParams(params,scaling=''):
    params = n.array(params)
    newParams = n.array(tuple(params[:9].reshape((3,3))[::-1,::-1].flatten())+tuple(params[9:][::-1]))
    if scaling !='': newParams[9:] = newParams[9:]/scaling
    return newParams


def eliminateAffineScaling(params):
    # adopt parameters to image spacing
    params = n.array(params)
    oldMatrix = n.diag((1.,1.,1.,1.))
    oldMatrix[:3,:3] = params[:9].reshape((3,3))
    #print oldMatrix
    decomp = transformations.decompose_matrix(oldMatrix)
    #print decomp[0]
    newMatrix = transformations.compose_matrix(scale=n.ones(3),shear=decomp[1],angles=decomp[2])
    newParams = newMatrix[:3,:3].flatten()
    newTrans = params[9:]*n.array([decomp[0][0],1.,1.])
    newParams = n.concatenate([newParams,newTrans],0)
    return newParams


def adaptParametersToSpacing(params,spacingRef,spacingFinal):
    # adopt parameters to new image spacing
    params = n.array(params)
    spacingRef = n.diag(spacingRef)
    spacingFinal = n.diag(spacingFinal)
    #print spacingRef,spacingFinal
    newParams = n.dot(n.dot(spacingFinal,params[:9].reshape((3,3))),n.linalg.inv(spacingRef)).flatten()
    newTrans = n.dot(spacingFinal,params[9:])
    newParams = n.concatenate([newParams,newTrans],0)
    return newParams
    

def composeAffineTransformations(paramsList):
    """
    compose parameters from subsequent affine transformations given in list paramsList
    """
    currentParams = n.array([1,0,0,0,1,0,0,0,1,0,0,0])
    for params in paramsList:
        params = n.array(params)
        tmpMatrix = n.dot(currentParams[:9].reshape((3,3)),params[:9].reshape((3,3)))
        tmpTrans = n.dot(currentParams[:9].reshape((3,3)),params[9:])+currentParams[9:]
        currentParams = n.concatenate([tmpMatrix.flatten(),tmpTrans],0)
    return currentParams

def composeAffineTransformationsNew(paramsList):
    """
    compose parameters from subsequent affine transformations given in list paramsList
    """
    currentParams = n.array([1,0,0,0,1,0,0,0,1,0,0,0])
    for params in paramsList:
        params = n.array(params)
        tmpMatrix = n.dot(params[:9].reshape((3,3)),currentParams[:9].reshape((3,3)))
        tmpTrans = n.dot(params[:9].reshape((3,3)),currentParams[9:])+params[9:]
        currentParams = n.concatenate([tmpMatrix.flatten(),tmpTrans],0)
    return currentParams

def registerWithElastix(fi,mi,config,initialParams=None,tmpFolder=None):

    if tmpFolder is None:
        tmpFolder = '/tmp'
    if initialParams == None:
        #params = n.array([1,0,0,0,-1,0,0,0,1,0,fi.GetSize()[1],0])
        # params = n.array([-1,0,0,0,1,0,0,0,1,fi.GetSize()[0],0,0])
        #params = invertParams(n.array([0.99728742898194889, -0.00038156465789470365, -0.00045136981651977868, -0.00078167254385604035, 1.0043850721125063, -0.014534116681577268, -0.0073061921565782829, -0.015172334919639853, -1.003839698888308, 0.90114927532248734, 137.0586870286873, 2173.3420588855224]))
        params = n.array([1,0,0,0,1,0,0,0,1,0,0,0])
    else:
        params = initialParams

    #fi = fi[fi.GetSize()[0]

    fiPath = os.path.join(tmpFolder,'tmpFixedImage.mhd')
    miPath = os.path.join(tmpFolder,'tmpMovingImage.mhd')
    sitk.WriteImage(filters.constantBackgroundSubtraction(fi,100),fiPath)
    sitk.WriteImage(filters.constantBackgroundSubtraction(mi,100),miPath)

    paramDict = dict()
    paramDict['el'] = config['elastixDir']
    parameterPath = os.path.join(tmpFolder,'elastixParameters.txt')
    initialTransformPath = os.path.join(tmpFolder,'elastixInitialTransform.txt')
    paramDict['out'] = tmpFolder
    extras.createInitialTransformFile(n.array([1,1,1.]),params,config['initialTransformTemplateString'],initialTransformPath)
    extras.createParameterFile(n.array([1,1,1.]),initialTransformPath,config['parameterTemplateStringAffine'],parameterPath)
    paramDict['params'] = parameterPath
    paramDict['f1'] = fiPath
    paramDict['f2'] = miPath
    paramDict['initialTransform'] = initialTransformPath

    # run elastix
    cmd = ('%(el)s -f %(f1)s -m %(f2)s -p %(params)s -t0 %(initialTransform)s -out %(out)s' %paramDict).split(' ')
    print "\n\nCalling elastix for image based registration with arguments:\n\n %s\n\n" %cmd
    subprocess.Popen(cmd).wait()

    # read output parameters from elastix output file
    rawOutParams = open(os.path.join(paramDict['out'],'TransformParameters.0.txt')).read()
#    outParams = rawOutParams.split('\n')[2][:-1].split(' ')[1:]
#    outParams = n.array([float(i) for i in outParams])
#    outCenterOfRotation = rawOutParams.split('\n')[19][:-1].split(' ')[1:]
#    outCenterOfRotation = n.array([float(i) for i in outCenterOfRotation])

#    # alternatively: use rigid transform instead of affine
#    affineOutParams = transformations.euler_matrix(outParams[0],outParams[1],outParams[2])

#    # correct for different center of rotation used in elastix
#    #affineOutParams = outParams
#    totalTranslate = outParams[-3:] - n.dot(affineOutParams.flatten()[:9].reshape((3,3)),outCenterOfRotation) + outCenterOfRotation
#    affineOutParams = n.concatenate([affineOutParams.flatten()[:9],totalTranslate],0)
#    #pdb.set_trace()
#    finalParams = composeAffineTransformations([affineOutParams,params])


    outParams = rawOutParams.split('\n')[2][:-1].split(' ')[1:]
    outParams = n.array([float(i) for i in outParams])
    outCenterOfRotation = rawOutParams.split('\n')[19][:-1].split(' ')[1:]
    outCenterOfRotation = n.array([float(i) for i in outCenterOfRotation])

    # alternatively: use rigid transform instead of affine
    #affineOutParams = transformations.euler_matrix(outParams[0],outParams[1],outParams[2])

    # correct for different center of rotation used in elastix
    # affineOutParams = outParams
    totalTranslate = outParams[-3:] - n.dot(outParams[:9].reshape((3,3)),outCenterOfRotation) + outCenterOfRotation
    affineOutParams = n.concatenate([outParams[:9],totalTranslate],0)
    finalParams = composeAffineTransformations([affineOutParams,params])

    return finalParams
