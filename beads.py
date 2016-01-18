#import filing
from imports import *
import affine_fit
import ransac

def loadFiles(files,config=None):
    print 'loading files...'
    ims = []
    for ifi,fi in enumerate(files):
        if fi[-3:] == 'tif' or fi[-3:] == 'iff':
            tmpImage = sitk.ReadImage(fi)
        elif fi[-3:] == '.h5':
            tmpFile = h5py.File(fi,mode='r')
            #pdb.set_trace()
            tmpImage = sitk.GetImageFromArray(tmpFile[config['h5path']].value)
            #timeStamp = tmpFile.attrs['TimeStamp']
            tmpFile.close()
            tmpFile = 0
            #tmpFile = open(os.path.join(config['fusedDir'],'timestamps.txt'),'a')
            #tmpFile.write(fi+': \t'+timeStamp+'\n')
            #tmpFile.close()
        tmpSize = n.array(tmpImage.GetSize())
        print tmpSize
        # make sure largest axes go first, f.e. resulting size should be (2000,800,200)
        if not(tmpSize[2]<=tmpSize[1] and tmpSize[1]<=tmpSize[0]):
            axes = n.argsort(tmpSize)[::-1]
            transformationParams = n.zeros((3,3))
            for iax,ax in enumerate(axes):
                transformationParams[iax][ax] = 1
            transformationParams = list(transformationParams.flatten())+[0,0,0]
            tmpImage = transformStack(transformationParams,tmpImage,outShape=tmpSize[axes])
            #print 'scaled!'
            #tmpImage = transformStack(None,tmpImage,outShape=tmpSize[axes]/2,outSpacing=n.array(tmpImage.GetSpacing())*2)
        if config['cropInputImages']:
            if not ifi%2:
                lowerCropBoundary = n.array(config['cropOffsetLC']).astype(n.int32)
                upperCropBoundary = n.array(tmpImage.GetSize())-n.array(config['cropSizeLC'])-lowerCropBoundary
            else:
                lowerCropBoundary = n.array(config['cropOffsetRC']).astype(n.int32)
                upperCropBoundary = n.array(tmpImage.GetSize())-n.array(config['cropSizeRC'])-lowerCropBoundary
            upperCropBoundary = upperCropBoundary.astype(n.uint32)
            if lowerCropBoundary.min()<0 or upperCropBoundary.min()<0:
                tmpTuple = (config['cropOffsetLC'],config['cropOffsetRC'],config['cropSizeLC'],config['cropSizeRC'],tmpImage.GetSize())
                raise Exception("\nwrong crop settings:\ncropOffsets: LC %s RC %s\ncropSize: LC %s RC %s\nimage size: %s" %tmpTuple)
            tmpImage = sitk.Crop(tmpImage,tuple(int(i) for i in lowerCropBoundary),tuple(int(i) for i in upperCropBoundary))
        ims.append(tmpImage)
    return ims


def segmentBeads(ims,mexicanHatKernelSize):
    #thresholds = [n.mean(i)+2*n.std(i) for i in ims]
    print 'segmenting beads (threshold: mean+std)'
    rs = []
    # pdb.set_trace()
    for imi,im in enumerate(ims):
        #im = sitk.GetImageFromArray(im)
        im = sitk.SmoothingRecursiveGaussian(im,mexicanHatKernelSize)-sitk.SmoothingRecursiveGaussian(im,mexicanHatKernelSize*2)
        stats = sitk.Statistics(im).values()
        threshold = stats[1]+5*stats[3]
        l = sitk.ConnectedComponent(im>threshold)
        l = sitk.GetArrayFromImage(l)
        im = sitk.GetArrayFromImage(im)
        centers = ndimage.center_of_mass(im,l,range(1,n.max(l)+1))
        rs.append(n.array(centers))
        del l
        print 'found %s beads in image %s' %(len(centers),imi)
    return rs

def segmentBeadsFaster(ims,mexicanHatKernelSize):
    #thresholds = [n.mean(i)+2*n.std(i) for i in ims]
    print 'segmenting beads (threshold: mean+std)'
    rs = []
    for imi,im in enumerate(ims):
        #im = sitk.GetImageFromArray(im)
        im = sitk.SmoothingRecursiveGaussian(im,mexicanHatKernelSize)-sitk.SmoothingRecursiveGaussian(im,mexicanHatKernelSize*2)
        stats = sitk.Statistics(im).values()
        threshold = stats[1]+5*stats[3]
        l = sitk.ConnectedComponent(im>threshold)
        l = sitk.GetArrayFromImage(l)
        im = sitk.GetArrayFromImage(im)
        centers = ndimage.center_of_mass(im,l,range(1,n.max(l)+1))
        rs.append(n.array(centers))
        del l
        print 'found %s beads in image %s' %(len(centers),imi)
    return rs

def segmentBeadsNils(ims,mexicanHatKernelSize):
    #thresholds = [n.mean(i)+2*n.std(i) for i in ims]
    print 'segmenting beads (threshold: mean+std)'
    rs = []
    for imi,im in enumerate(ims):
        im = sitk.Cast(im,3)
        zeroRegion = sitk.Cast(im==0,sitk.sitkUInt8)
        for iter in range(5):
            zeroRegion = sitk.BinaryDilate(zeroRegion)
        zeroRegion = sitk.Cast(zeroRegion==0,sitk.sitkUInt8)
        #im = sitk.GetImageFromArray(im)
        im = sitk.Cast(im,sitk.sitkInt32)
        im = sitk.Abs(im-sitk.Cast(sitk.SmoothingRecursiveGaussian(im,mexicanHatKernelSize*5),im.GetPixelID()))
        #pdb.set_trace()
        for iter in range(1):
            im = sitk.SmoothingRecursiveGaussian(im,mexicanHatKernelSize)-sitk.SmoothingRecursiveGaussian(im,mexicanHatKernelSize*2)
        im = im*sitk.Cast(zeroRegion,im.GetPixelID())
        stats = sitk.Statistics(im).values()
        threshold = stats[1]+10*stats[3]
        l = sitk.ConnectedComponent(im>threshold)
        l = sitk.GetArrayFromImage(l)
        im = sitk.GetArrayFromImage(im)
        centers = ndimage.center_of_mass(im,l,range(1,n.max(l)+1))
        rs.append(n.array(centers))
        del l
        print 'found %s beads in image %s' %(len(centers),imi)
    return rs

def segmentBeadsWithSample(ims,mexicanHatKernelSize,sampleSignalLevel=220,ignoreSizeRelativeToPixelSize=1/20.,stdFactor=20.,excludeSize=None):
    #thresholds = [n.mean(i)+2*n.std(i) for i in ims]
    print 'segmenting beads (threshold: mean+std)'
    rs = []
    brightnesses = []
    for imi,im in enumerate(ims):
        sums = sitk.SumProjection(sitk.SumProjection(im,0),1)
        sums = sitk.GetArrayFromImage(sums).squeeze()
        #pdb.set_trace()
        goodPlaneEnd = n.where(sums==0)[0]
        goodPlaneStart = n.where(sums>0)[0]
        goodPlaneEnd = goodPlaneEnd[goodPlaneEnd>n.min(goodPlaneStart)]
        if not len(goodPlaneEnd): goodPlaneEnd = im.GetSize()[2]
        else: goodPlaneEnd = n.min(goodPlaneEnd)
        upperCrop = tuple([int(i) for i in [im.GetSize()[0],im.GetSize()[1],goodPlaneEnd]])
        im = sitk.Slice(im,(0,0,0),upperCrop)
        cc = sitk.ConnectedComponent(im>sampleSignalLevel)
        #statsFilter = sitk.LabelStatisticsImageFilter()
        #statsFilter.Execute(cc,cc)
        cc2 = sitk.GetArrayFromImage(cc).flatten()
        cc2 = cc2[cc2.nonzero()]
        sizes = n.zeros(len(set(cc2)))
        for e in cc2:
            sizes[int(e)-1] += 1
        sizes = n.array(sizes)
        imsize = n.array(im.GetSize())
        sizeThreshold = n.product(imsize)*n.power(ignoreSizeRelativeToPixelSize,3)
        sampleLabels = n.where(sizes>sizeThreshold)[0]
        #pdb.set_trace()
        for sampleLabel in sampleLabels:
            sampleMask = (cc==(sampleLabel+1))
            sampleMask = sitk.BinaryFillhole(sampleMask)
            im = im*sitk.Cast(sitk.Abs(sitk.Cast(sampleMask,sitk.sitkInt8)-1),3)+sampleSignalLevel*sitk.Cast(sampleMask,3)
            del sampleMask
        del cc,cc2
        sitk.WriteImage(im,'/home/malbert/delme/delme1_%s.tif' %imi)
        iterations = 1
        if len(mexicanHatKernelSize) == 1:
            for i in range(iterations):
                im = sitk.SmoothingRecursiveGaussian(im,mexicanHatKernelSize)-sitk.SmoothingRecursiveGaussian(im,mexicanHatKernelSize*2)
        else:
            for i in range(iterations):
                im1 = sitk.Cast(im,im.GetPixelID())
                im2 = sitk.Cast(im,im.GetPixelID())
                for axis in range(len(mexicanHatKernelSize)):
                    if mexicanHatKernelSize[axis] == 0: continue
                    im1 = sitk.RecursiveGaussian(im1,mexicanHatKernelSize[axis],direction=axis)
                    im2 = sitk.RecursiveGaussian(im2,mexicanHatKernelSize[axis]*2,direction=axis)
                im = im1-im2

        # eliminate border effects
        origSize = n.array(im.GetSize())
        im = sitk.Crop(im,(1,1,1),(1,1,1))
        im = sitk.Resample(im,tuple(int(i) for i in origSize))

        sitk.WriteImage(im,'/home/malbert/delme/delme2_%s.tif' %imi)
        stats = sitk.Statistics(im).values()
        threshold = stats[1]+stdFactor*stats[3]
        #pdb.set_trace()
        sitk.WriteImage(sitk.Threshold(im,threshold,64000),'/home/malbert/delme/delme3_%s.tif' %imi)
        #sitk.WriteImage(im,'/home/malbert/delme/delme4_%s.tif' %imi)
        print threshold
        #sitk.OtsuMultipleThresholds(im, 15)
        #l = sitk.ConnectedComponent(sitk.OtsuThreshold(im, 5))
        l = sitk.ConnectedComponent(sitk.Threshold(im,threshold,64000)>0,fullyConnected=True)

        cc2 = sitk.GetArrayFromImage(l).flatten()
        cc2 = cc2[cc2.nonzero()]
        sizes = n.zeros(len(set(cc2)))
        for e in cc2:
            sizes[int(e)-1] += 1
        sizes = n.array(sizes)

        if excludeSize is not None:
            sizeThreshold = n.min([excludeSize,sizeThreshold])
            sizeThreshold = excludeSize
            sampleLabels = n.where(sizes<=sizeThreshold)[0]+1
            l = sitk.GetArrayFromImage(l)
        else:
            l = sitk.GetArrayFromImage(l)
            sampleLabels = range(1,n.max(l)+1)

        im = sitk.GetArrayFromImage(im)
        centers = ndimage.center_of_mass(im,l,sampleLabels)
        #pdb.set_trace()
        rs.append(n.array(centers))
        bright = []
        for ii in range(len(centers)):
            bright.append(ims[imi][int(centers[ii][2]),int(centers[ii][1]),int(centers[ii][0])])
        brightnesses.append(bright)
        del l
        print 'found %s beads in image %s' %(len(centers),imi)
    return rs,[sizes,brightnesses]

def segmentBeadsWithSampleRot(ims,mexicanHatKernelSize,sampleSignalLevel=220,ignoreSizeRelativeToPixelSize=1/20.,stdFactor=20.,excludeSize=None,rotationAxis=1):
    #thresholds = [n.mean(i)+2*n.std(i) for i in ims]
    print 'segmenting beads (threshold: mean+std)'
    rs = []
    brightnesses = []
    for imi,im in enumerate(ims):
        sums = sitk.SumProjection(sitk.SumProjection(im,0),1)
        sums = sitk.GetArrayFromImage(sums).squeeze()
        #pdb.set_trace()
        goodPlaneEnd = n.where(sums==0)[0]
        goodPlaneStart = n.where(sums>0)[0]
        goodPlaneEnd = goodPlaneEnd[goodPlaneEnd>n.min(goodPlaneStart)]
        if not len(goodPlaneEnd): goodPlaneEnd = im.GetSize()[2]
        else: goodPlaneEnd = n.min(goodPlaneEnd)
        upperCrop = tuple([int(i) for i in [im.GetSize()[0],im.GetSize()[1],goodPlaneEnd]])
        im = sitk.Slice(im,(0,0,0),upperCrop)
        cc = sitk.ConnectedComponent(im>sampleSignalLevel)
        #statsFilter = sitk.LabelStatisticsImageFilter()
        #statsFilter.Execute(cc,cc)
        cc2 = sitk.GetArrayFromImage(cc).flatten()
        cc2 = cc2[cc2.nonzero()]
        sizes = n.zeros(len(set(cc2)))
        for e in cc2:
            sizes[int(e)-1] += 1
        sizes = n.array(sizes)
        imsize = n.array(im.GetSize())
        sizeThreshold = n.product(imsize)*n.power(ignoreSizeRelativeToPixelSize,3)
        sampleLabels = n.where(sizes>sizeThreshold)[0]
        #pdb.set_trace()
        for sampleLabel in sampleLabels:
            sampleMask = (cc==(sampleLabel+1))
            sampleMask = sitk.BinaryFillhole(sampleMask)
            im = im*sitk.Cast(sitk.Abs(sitk.Cast(sampleMask,sitk.sitkInt8)-1),3)+sampleSignalLevel*sitk.Cast(sampleMask,3)
            del sampleMask
        del cc,cc2
        sitk.WriteImage(im,'/home/malbert/delme/delme1_%s.tif' %imi)
        iterations = 3
        for i in range(iterations):
            im = sitk.RecursiveGaussian(im,mexicanHatKernelSize,direction=rotationAxis)-sitk.RecursiveGaussian(im,mexicanHatKernelSize*2,direction=rotationAxis)

        # eliminate border effects
        origSize = n.array(im.GetSize())
        im = sitk.Crop(im,(1,1,1),(1,1,1))
        im = sitk.Resample(im,tuple(int(i) for i in origSize))

        sitk.WriteImage(im,'/home/malbert/delme/delme2_%s.tif' %imi)
        stats = sitk.Statistics(im).values()
        threshold = stats[1]+stdFactor*stats[3]
        #pdb.set_trace()
        sitk.WriteImage(sitk.Threshold(im,threshold,64000),'/home/malbert/delme/delme3_%s.tif' %imi)
        #sitk.WriteImage(im,'/home/malbert/delme/delme4_%s.tif' %imi)
        print threshold
        #sitk.OtsuMultipleThresholds(im, 15)
        #l = sitk.ConnectedComponent(sitk.OtsuThreshold(im, 5))
        l = sitk.ConnectedComponent(sitk.Threshold(im,threshold,64000)>0,fullyConnected=True)

        cc2 = sitk.GetArrayFromImage(l).flatten()
        cc2 = cc2[cc2.nonzero()]
        sizes = n.zeros(len(set(cc2)))
        for e in cc2:
            sizes[int(e)-1] += 1
        sizes = n.array(sizes)

        if excludeSize is not None:
            sizeThreshold = n.min([excludeSize,sizeThreshold])
            sizeThreshold = excludeSize
            sampleLabels = n.where(sizes<=sizeThreshold)[0]+1
            l = sitk.GetArrayFromImage(l)
        else:
            l = sitk.GetArrayFromImage(l)
            sampleLabels = range(1,n.max(l)+1)

        im = sitk.GetArrayFromImage(im)
        centers = ndimage.center_of_mass(im,l,sampleLabels)
        #pdb.set_trace()
        rs.append(n.array(centers))
        bright = []
        for ii in range(len(centers)):
            bright.append(ims[imi][int(centers[ii][2]),int(centers[ii][1]),int(centers[ii][0])])
        brightnesses.append(bright)
        del l
        print 'found %s beads in image %s' %(len(centers),imi)
    return rs,[sizes,brightnesses]

def checkBeads(rs,images,beadIntensity=2000):
    newims = [sitk.Cast(i,i.GetPixelID()) for i in images]
    for ii in range(len(rs)):
        for ip,p in enumerate(rs[ii]):
            newims[ii][int(p[0]),int(p[1]),int(p[2])] = beadIntensity
    return newims


def scaleBeads(beads,scaling=None,im=None,verbose=False):
    if scaling == None:
        scaling = estimateScaling(beads,im)
    if verbose: print 'scaling beads with factors %s' %scaling
    beads = n.array(beads).swapaxes(0,1)
    beads = n.array([beads[i]*scaling[i] for i in range(3)]).swapaxes(0,1)
    return beads

def calcFeatures3(rs,N=4):
    print 'calculating invariant features...'
    tree = spatial.cKDTree(rs)
    features = []
    boundaries = n.max(rs,0)
    for ii,i in enumerate(rs):
        ns = tree.query(i,N+1)
        neighInds = ns[1][1:]
        neighDists = ns[0][1:]
        relCoords = n.array([rs[j]-i for j in neighInds])
        scalarProducts = n.array([n.dot(relCoords[ii1],relCoords[ii2])/n.linalg.norm(relCoords[ii1])/n.linalg.norm(relCoords[ii2]) for (ii1,ii2) in combinations(range(N),2)])
        relatives = n.min([i,boundaries-i],0)/boundaries
        feature = n.concatenate([scalarProducts,relatives])
        features.append(feature)
    return n.array(features)

def calcFeatures2(rs,N=4):
    print 'calculating invariant features...'
    tree = spatial.cKDTree(rs)
    features = []
    for ii,i in enumerate(rs):
        ns = tree.query(i,N+1)
        neighInds = ns[1][1:]
        neighDists = ns[0][1:]
        features.append(neighDists)
    return n.array(features)


def calcFeatures(rs,N=4):
    #print 'calculating invariant features...'
    tree = spatial.KDTree(rs)
    features = []
    for ii,i in enumerate(rs):
        ns = tree.query(i,N+1)
        neighInds = ns[1][1:]
        neighDists = ns[0][1:]
        relCoords = n.array([rs[j]-i for j in neighInds])
        scalarProducts = n.array([n.dot(relCoords[ii1],relCoords[ii2])/n.linalg.norm(relCoords[ii1])/n.linalg.norm(relCoords[ii2]) for (ii1,ii2) in combinations(range(N),2)])
        features.append(scalarProducts)
    return n.array(features)

def calcFeaturesFallback(rs,N=4):
    #print 'calculating invariant features...'
    tree = spatial.KDTree(rs)
    features = []
    for ii,i in enumerate(rs):
        ns = tree.query(i,4+1)
        neighInds = ns[1][1:]
        neighDists = ns[0][1:]
        neighCoords = rs[neighInds]
        dataM = n.append(neighCoords,n.ones((4,1)),1).T
        colinear = False
        try:
            barycoords = n.linalg.solve(dataM,n.array([i[0],i[1],i[2],1]))
        except:
            colinear = True
        if colinear:
            feature = n.zeros(4)
        else:
            feature = barycoords
        features.append(feature)
    return n.array(features)

def calcFeaturesOld(rs,N=4):
    print 'calculating invariant features...'
    tree = spatial.cKDTree(rs)
    features = []
    for ii,i in enumerate(rs):
        ns = tree.query(i,N+1)
        neighInds = ns[1]
        neighDists = ns[0]
        relCoords = n.array([rs[neighInds[j1]]-rs[neighInds[j2]] for (j1,j2) in combinations(range(len(neighInds)),2)])
        scalarProducts = n.array([n.dot(relCoords[i],1./relCoords[i+1]) for i in range(len(relCoords)-1)])
        features.append(scalarProducts)
    return n.array(features)

def bestMatchNormal(features1,features2):
    print 'matching beads...'
    # returns index order relative to nearestneigh(f1) and range(len(f2))
    tree = spatial.cKDTree(features1)
    ns = tree.query(features2,1)
    indices = ns[1]
    distances = ns[0]
    order = n.argsort(distances)
    pairs = n.array([indices[order],order])
    return pairs,distances[order]

def bestMatch(features1,features2):
    print 'matching beads...'
    # returns index order relative to nearestneigh(f1) and range(len(f2))
    i1 = n.where(n.sum(n.abs(features1),1))[0]
    i2 = n.where(n.sum(n.abs(features2),1))[0]
    #pdb.set_trace()
    tree = spatial.KDTree(features1[i1])
    ns = tree.query(features2[i2],1)
    indices = i1[ns[1]]
    distances = ns[0]
    order = n.argsort(distances)
    pairs = n.array([indices[order],order])
    return pairs,distances[order]

def transform(p,rs):
    p = n.array(p)
    rs = n.array(rs)
    newrs = n.zeros(rs.shape)
    for ii,i in enumerate(rs):
        tmp = n.dot(p[:9].reshape((3,3)),i)+p[9:]
        newrs[ii] = tmp
    return n.array(newrs)

def backTransform(p,rs):
    rs = n.array(rs)
    newrs = n.zeros(rs.shape)
    for ii,i in enumerate(rs):
        m = n.linalg.inv(p[:9].reshape((3,3)))
        tmp = n.dot(m,(i-p[9:]))
        newrs[ii] = tmp
    return n.array(newrs)

def errfuncFull(p,rs,refrs):
    dists = n.zeros(len(refrs))
    newrs = transform(p,rs)
    for ii,i in enumerate(refrs):
        dists[ii] = n.linalg.norm(i-newrs[ii])
    #print n.sum(dists)
    return dists

def checkim(im,beads):
    outim = im.copy()
    for bead in beads:
        tmp = tuple(bead.astype(n.int64))
        broken = False
        for i in range(3):
            if tmp[i]<0 or tmp[i]>im.shape[i]-1:
                broken = True
                break
        if broken: continue
        outim[tuple(bead.astype(n.int64))] = 20000
    return outim

def selectGoodInds(dists,ngauss,fstd=1,qbins=10,nmax=50,nsamples=1):
    maxind = n.where(n.isinf(distances)==False)[0][-1]
    dists = dists[:maxind+1]
    maxind = n.where(dists<=n.percentile(dists,95))[0][-1]
    dists = dists[:maxind+1]
    steps = len(dists)/qbins
    x = n.linspace(dists.min(),dists.max(),steps)
    hist = ndimage.histogram(dists,dists.min(),dists.max(),steps)
    def fitfunc(p,x,ngauss):
        gausses = n.zeros(len(x))
        for igauss in range(ngauss):
            A = p[igauss*len(p)/ngauss+0]
            mu = p[igauss*len(p)/ngauss+1]
            std = p[igauss*len(p)/ngauss+2]
            gausses += A*n.exp(-n.power((x-mu)/std,2))
        return gausses
    def errfunc(p,x,ngauss,hist):
        return fitfunc(p,x,ngauss)-hist
    if not ngauss-1: x0 = [10,x[0],x.std()/2]
    else:
        x0 = []
        for igauss in range(ngauss):
            x0+=[qbins,x[0]+igauss/(ngauss-1)*(x.max()-x.min())/ngauss,x.std()/2]
    print x0
    fit,succ = optimize.leastsq(errfunc,x0,args=(x,ngauss,hist),maxfev=1000000)
    print fit
    gaussoffset = len(x0)/ngauss*n.argmin([fit[i*len(x0)/ngauss+1] for i in range(ngauss)])
    goodinds = n.where(dists<fit[gaussoffset+1]+fstd*n.abs(fit[gaussoffset+2]))[0]
    print gaussoffset
    print fit[gaussoffset+1]+fstd*fit[gaussoffset+2]
    samples = []
    for sample in range(nsamples):
        samples.append(n.random.permutation(goodinds)[:nmax])
    #pdb.set_trace()
    return samples,goodinds,fit,hist,x,fitfunc(fit,x,ngauss)

def getRotationMatrix(angle,scaling=n.ones(3)):
    matrix = n.diag(scaling)
    matrix[0][0] = n.cos(angle)
    matrix[0][1] = -n.sin(angle)
    matrix[1][0] = n.sin(angle)
    matrix[1][1] = n.cos(angle)
    return matrix

def estimateScalingOld(beads):
    #side = n.min(im.shape)
    beads = n.array(beads)
    side = n.min([n.max(i) for i in beads.swapaxes(0,1)])
    result = []
    for i in range(3):
        inds = range(3)
        inds.remove(i)
        tmp1 = n.where(beads.swapaxes(0,1)[inds[0]]<=side)[0]
        tmp2 = n.where(beads.swapaxes(0,1)[inds[1]]<=side)[0]
        tmp = len(set(tmp1).intersection(set(tmp2)))
        #tmp = tmp/float(im.shape[i])
        #tmp = tmp/float(n.max(beads[:,i]))
        result.append(tmp)
    #result = 1./n.array(result)
    #result = result/result[2]
    return result,tmp


def findParametersOld(rs,N=4,useFeatures=True):
    # problem in bestMatch() solved?
    params =[[1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.]]
    for ibeads in range(1,len(rs)):
        tmpWorked = False
        for featureFunc in [calcFeatures,calcFeaturesFallback]:
            if tmpWorked: break
            for scale in n.linspace(1,10,50):
                scale = n.array([scale,1.,1.])
                #print scale
                tmpbeads = scaleBeads(rs[ibeads],scale,verbose=True)
                refbeads = scaleBeads(rs[0],scale)
                if useFeatures:
                    f0 = featureFunc(refbeads,N=N)
                    # find initial matches and calculate approximate transform
                    tmpf = featureFunc(tmpbeads,N=N)
                else:
                    f0 = refbeads
                    tmpf = tmpbeads
                indices,distances = bestMatch(f0,tmpf)
                number = 20
                number = n.max([len(indices[0])/10,number])
                number = n.min([len(indices[0]),number])
                model = ransac.AffineTransformationModel(rs[0][indices[0]],rs[ibeads][indices[1]])
                data = n.array(range(number))
                print number
                try:
                    #pdb.set_trace()
                    fit,inliers = ransac.ransac(data,model,4,number*10,50,number/2,return_all=1)
                    tmpWorked = True
                    print 'OK: found parameters'
                    break
                except ValueError:
                    print 'try again with different z scale factor...'
                    pass

        if not tmpWorked:
            raise Exception('Could not register bead images.')
        
        # find better transform using redefined matches
        indices2,distances2 = bestMatch(scaleBeads(transform(fit,rs[0]),scale),scaleBeads(rs[ibeads],scale))
        inds = n.where(distances2<5.)[0]
        if len(inds)<number:
            print '#'*10+'\ncheck bead registration since only a few beads are well aligned\n'+'#'*10
            inds = n.array(range(number))
        fit2 = affine_fit.Affine_Fit(rs[0][indices2[0,inds]],rs[ibeads][indices2[1,inds]]).Matrix()
        params.append(fit2)
        indices3,distances3 = bestMatch(transform(fit2,rs[0]),rs[ibeads])
        print 'mean bead distance in pixels after alignment is\n%s' %n.mean(distances3)
    return params

def findParameters(rs,N=4,useFeatures=True):
    # problem in bestMatch() solved?
    params =[[1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.]]
    for ibeads in range(1,len(rs)):
        tmpWorked = False
        # for featureFunc in [calcFeatures,calcFeaturesFallback]:
        for ifeaturefunc,featureFunc in enumerate([calcFeatures,
                                                   calcFeatures,
                                                   calcFeatures,
                                                   calcFeaturesFallback]):
            if tmpWorked: break
            if useFeatures:
                if ifeaturefunc in [0,2]: scales = n.linspace(1,10,50)
                else: scales = n.linspace(1,10,500)
            else: scales = [1.]
            for scale in scales:
                scale = n.array([scale,1.,1.])
                #print scale
                tmpbeads = scaleBeads(rs[ibeads],scale,verbose=True)
                refbeads = scaleBeads(rs[0],scale)
                if useFeatures:
                    f0 = featureFunc(refbeads,N=N)
                    # find initial matches and calculate approximate transform
                    tmpf = featureFunc(tmpbeads,N=N)
                else:
                    f0 = refbeads
                    tmpf = tmpbeads
                if len(f0)>=len(tmpf): indList = [0,1]
                else: indList = [1,0]
                indices,distances = bestMatch(*tuple([[f0,tmpf][i] for i in indList]))
                number = 20
                number = n.max([len(indices[0])/10,number])
                number = n.min([len(indices[0]),number])
                model = ransac.AffineTransformationModel(rs[0][indices[indList[0]]],rs[ibeads][indices[indList[1]]])
                data = n.array(range(number))
                #print indices
                #print number
                if useFeatures:
                    if ifeaturefunc in [0,1]: maxIterations = number*10
                    else: maxIterations = number*500
                else: maxIterations = number*1000
                try:
                    #pdb.set_trace()
                    fit,inliers = ransac.ransac(data,model,4,maxIterations,50,number/2,return_all=1)
                    tmpWorked = True
                    print 'OK: found parameters'
                    break
                except ValueError:
                    print 'try again with different z scale factor...'
                    pass

        if not tmpWorked:
            raise Exception('Could not register bead images.')

        # find better transform using redefined matches
        tmpArgs = [[scaleBeads(transform(fit,rs[0]),scale),scaleBeads(rs[ibeads],scale)][i] for i in indList]
        indices2,distances2 = bestMatch(*tuple(tmpArgs))
        inds = n.where(distances2<5.)[0]
        if len(inds)<number:
            print '#'*10+'\ncheck bead registration since only a few beads are well aligned\n'+'#'*10
            inds = n.array(range(number))
        fit2 = affine_fit.Affine_Fit(rs[0][indices2[indList[0],inds]],rs[ibeads][indices2[indList[1],inds]]).Matrix()
        params.append(fit2)
        indices3,distances3 = bestMatch(transform(fit2,rs[0]),rs[ibeads])
        print 'mean bead distance in pixels after alignment is\n%s' %n.mean(distances3)
    return params

def findParametersSheared(rs,N=4,shearxs=n.linspace(0.2,1,10)):
    # problem in bestMatch() solved?
    rs[0] = rs[0][:,[2,1,0]]
    rs[1] = rs[1][:,[2,1,0]]
    tmpWorked=False
    params =[[1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.]]
    shearMatricesL = [n.array([1.,0,x,0,1,0,0,0,1,0,0,0]) for x in range(2*len(shearxs))]
    shearMatricesR = [n.array([1.,0,x,0,1,0,0,0,1,0,0,0]) for x in range(2*len(shearxs))]
    for ix,x in enumerate(shearxs):
        shearMatricesL[2*ix][2] = x
        shearMatricesL[2*ix+1][2] = -x
        shearMatricesR[2*ix][2] = -x
        shearMatricesR[2*ix+1][2] = +x

    for ix in range(len(shearMatricesL)):
        refbeads = transform(shearMatricesL[ix],rs[0])
        tmpbeads = transform(shearMatricesR[ix],rs[1])
        print shearMatricesL[ix]

        f0 = calcFeatures(refbeads,N=N)
        # find initial matches and calculate approximate transform
        tmpf = calcFeatures(tmpbeads,N=N)

        if len(f0)>=len(tmpf): indList = [0,1]
        else: indList = [1,0]
        indices,distances = bestMatch(*tuple([[f0,tmpf][i] for i in indList]))
        number = 20
        number = n.max([len(indices[0])/10,number])
        number = n.min([len(indices[0]),number])
        model = ransac.AffineTransformationModel(rs[0][indices[indList[0]]],rs[1][indices[indList[1]]])
        data = n.array(range(number))
        #print indices
        #print number
        maxIterations = number*10
        try:
            fit,inliers = ransac.ransac(data,model,4,maxIterations,50,number/2,return_all=1)
            tmpWorked = True
            print 'OK: found parameters'
            break
        except ValueError:
            print 'try again with different shearing parameter'
            pass

    if not tmpWorked:
        raise Exception('Could not register bead images.')

    # find better transform using redefined matches
    tmpArgs = [[scaleBeads(transform(fit,rs[0]),scale),scaleBeads(rs[ibeads],scale)][i] for i in indList]
    indices2,distances2 = bestMatch(*tuple(tmpArgs))
    inds = n.where(distances2<5.)[0]
    if len(inds)<number:
        print '#'*10+'\ncheck bead registration since only a few beads are well aligned\n'+'#'*10
        inds = n.array(range(number))
    fit2 = affine_fit.Affine_Fit(rs[0][indices2[indList[0],inds]],rs[ibeads][indices2[indList[1],inds]]).Matrix()
    params.append(fit2)
    indices3,distances3 = bestMatch(transform(fit2,rs[0]),rs[ibeads])
    print 'mean bead distance in pixels after alignment is\n%s' %n.mean(distances3)
    return params

def findParametersNilsNew(rs,N=4,useFeatures=True,dz=None):
    # problem in bestMatch() solved?
    params =[[1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.]]
    for ibeads in range(1,len(rs)):
        tmpWorked = False
        for featureFunc in [calcFeatures,calcFeaturesFallback]:
            if tmpWorked: break
            if useFeatures and dz is None:
                scales = n.append(n.linspace(0.3,1,10),n.linspace(1,10,30),0)
            elif not dz is None:
                scales = [float(dz)]
            else: scales = [1.]
            for scale in scales:
                scale = n.array([scale,1.,1.])
                #print scale
                tmpbeads = scaleBeads(rs[ibeads],scale,verbose=True)
                refbeads = scaleBeads(rs[0],scale)
                if useFeatures:
                    f0 = featureFunc(refbeads,N=N)
                    # find initial matches and calculate approximate transform
                    tmpf = featureFunc(tmpbeads,N=N)
                else:
                    f0 = refbeads
                    tmpf = tmpbeads
                if len(f0)>=len(tmpf): indList = [0,1]
                else: indList = [1,0]
                indices,distances = bestMatch(*tuple([[f0,tmpf][i] for i in indList]))
                number = 20
                number = n.max([len(indices[0])/10,number])
                number = n.min([len(indices[0]),number])
                model = ransac.AffineTransformationModel(rs[0][indices[indList[0]]],rs[ibeads][indices[indList[1]]])
                data = n.array(range(number))
                #print indices
                #print number
                if useFeatures: maxIterations = number*10
                else: maxIterations = number*1000
                try:
                    fit,inliers = ransac.ransac(data,model,4,maxIterations,50,number/2,return_all=1)
                    tmpWorked = True
                    print 'OK: found parameters'
                    break
                except ValueError:
                    print 'try again with different z scale factor...'
                    pass

        if not tmpWorked:
            raise Exception('Could not register bead images.')

        # find better transform using redefined matches
        tmpArgs = [[scaleBeads(transform(fit,rs[0]),scale),scaleBeads(rs[ibeads],scale)][i] for i in indList]
        indices2,distances2 = bestMatch(*tuple(tmpArgs))
        inds = n.where(distances2<5.)[0]
        if len(inds)<number:
            print '#'*10+'\ncheck bead registration since only a few beads are well aligned\n'+'#'*10
            inds = n.array(range(number))
        fit2 = affine_fit.Affine_Fit(rs[0][indices2[indList[0],inds]],rs[ibeads][indices2[indList[1],inds]]).Matrix()
        params.append(fit2)
        indices3,distances3 = bestMatch(transform(fit2,rs[0]),rs[ibeads])
        print 'mean bead distance in pixels after alignment is\n%s' %n.mean(distances3)
    return params

# def findParametersNils(rs,dz,N=4):
def findParametersNils(rs,dz):
    # Ns = [4,3,5,6]
    Ns = [4]
    useFeatures = True
    # problem in bestMatch() solved?
    params =[[1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.]]
    for ibeads in range(1,len(rs)):
        tmpWorked = False
        for N in Ns:
            if tmpWorked: break
            for featureFunc in [calcFeatures,calcFeaturesFallback]:
                if tmpWorked: break
                if useFeatures:
                    scales = n.linspace(0.75*dz,1.25*dz,30)
                else: scales = [1.]
                for scale in scales:
                    scale = n.array([scale,1.,1.])
                    #print scale
                    tmpbeads = scaleBeads(rs[ibeads],scale,verbose=True)
                    refbeads = scaleBeads(rs[0],scale)
                    if useFeatures:
                        f0 = featureFunc(refbeads,N=N)
                        # find initial matches and calculate approximate transform
                        tmpf = featureFunc(tmpbeads,N=N)
                    else:
                        f0 = refbeads
                        tmpf = tmpbeads
                    if len(f0)>=len(tmpf): indList = [0,1]
                    else: indList = [1,0]
                    indices,distances = bestMatch(*tuple([[f0,tmpf][i] for i in indList]))
                    number = 20
                    number = n.max([len(indices[0])/10,number])
                    number = n.min([len(indices[0]),number])
                    model = ransac.AffineTransformationModel(rs[0][indices[indList[0]]],rs[ibeads][indices[indList[1]]])
                    data = n.array(range(number))
                    #print indices
                    #print number
                    if useFeatures: maxIterations = number*10
                    else: maxIterations = number*1000
                    try:
                        fit,inliers = ransac.ransac(data,model,4,maxIterations,50,number/2,return_all=1)
                        tmpWorked = True
                        print 'OK: found parameters'
                        break
                    except ValueError:
                        print 'try again with different z scale factor...'
                        pass

        if not tmpWorked:
            raise Exception('Could not register bead images.')

        # find better transform using redefined matches
        tmpArgs = [[scaleBeads(transform(fit,rs[0]),scale),scaleBeads(rs[ibeads],scale)][i] for i in indList]
        indices2,distances2 = bestMatch(*tuple(tmpArgs))
        inds = n.where(distances2<5.)[0]
        if len(inds)<number:
            print '#'*10+'\ncheck bead registration since only a few beads are well aligned\n'+'#'*10
            inds = n.array(range(number))
        fit2 = affine_fit.Affine_Fit(rs[0][indices2[indList[0],inds]],rs[ibeads][indices2[indList[1],inds]]).Matrix()
        params.append(fit2)
        indices3,distances3 = bestMatch(transform(fit2,rs[0]),rs[ibeads])
        print 'mean bead distance in pixels after alignment is\n%s' %n.mean(distances3)
    return params

def findParametersNilsWorking(rs,N=4,useFeatures=True):
    # problem in bestMatch() solved?
    params =[[1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.]]
    for ibeads in range(1,len(rs)):
        tmpWorked = False
        for featureFunc in [calcFeatures,calcFeaturesFallback]:
            if tmpWorked: break
            if useFeatures:
                scales = n.append(n.linspace(0.3,1,10),n.linspace(1,10,30),0)
            else: scales = [1.]
            for scale in scales:
                scale = n.array([scale,1.,1.])
                #print scale
                tmpbeads = scaleBeads(rs[ibeads],scale,verbose=True)
                refbeads = scaleBeads(rs[0],scale)
                if useFeatures:
                    f0 = featureFunc(refbeads,N=N)
                    # find initial matches and calculate approximate transform
                    tmpf = featureFunc(tmpbeads,N=N)
                else:
                    f0 = refbeads
                    tmpf = tmpbeads
                if len(f0)>=len(tmpf): indList = [0,1]
                else: indList = [1,0]
                indices,distances = bestMatch(*tuple([[f0,tmpf][i] for i in indList]))
                number = 20
                number = n.max([len(indices[0])/10,number])
                number = n.min([len(indices[0]),number])
                model = ransac.AffineTransformationModel(rs[0][indices[indList[0]]],rs[ibeads][indices[indList[1]]])
                data = n.array(range(number))
                #print indices
                #print number
                if useFeatures: maxIterations = number*10
                else: maxIterations = number*1000
                try:
                    fit,inliers = ransac.ransac(data,model,4,maxIterations,50,number/2,return_all=1)
                    tmpWorked = True
                    print 'OK: found parameters'
                    break
                except ValueError:
                    print 'try again with different z scale factor...'
                    pass

        if not tmpWorked:
            raise Exception('Could not register bead images.')

        # find better transform using redefined matches
        tmpArgs = [[scaleBeads(transform(fit,rs[0]),scale),scaleBeads(rs[ibeads],scale)][i] for i in indList]
        indices2,distances2 = bestMatch(*tuple(tmpArgs))
        inds = n.where(distances2<5.)[0]
        if len(inds)<number:
            print '#'*10+'\ncheck bead registration since only a few beads are well aligned\n'+'#'*10
            inds = n.array(range(number))
        fit2 = affine_fit.Affine_Fit(rs[0][indices2[indList[0],inds]],rs[ibeads][indices2[indList[1],inds]]).Matrix()
        params.append(fit2)
        indices3,distances3 = bestMatch(transform(fit2,rs[0]),rs[ibeads])
        print 'mean bead distance in pixels after alignment is\n%s' %n.mean(distances3)
    return params

def findParametersTest(rs,N=4,useFeatures=True,infoDict=None):

    origins = infoDict['origins']
    positions = infoDict['positions']
    spacing = infoDict['spacing']
    centerOfRotation = infoDict['centerOfRotation']
    axisOfRotation = infoDict['axisOfRotation']
    sizes = infoDict['sizes']

    # problem in bestMatch() solved?
    params =[[1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.]]
    for ibeads in range(1,len(rs)):
        tmpWorked = False
        for featureFunc in [calcFeatures,calcFeaturesFallback]:
            if tmpWorked: break
            for scale in n.linspace(1,10,50):
                scale = n.array([scale,1.,1.])
                #print scale
                tmpbeads = scaleBeads(rs[ibeads],scale,verbose=True)
                refbeads = scaleBeads(rs[0],scale)
                if useFeatures:
                    f0 = featureFunc(refbeads,N=N)
                    # find initial matches and calculate approximate transform
                    tmpf = featureFunc(tmpbeads,N=N)
                else:
                    f0 = refbeads
                    tmpf = tmpbeads
                indices,distances = bestMatch(f0,tmpf)
                number = n.max([len(indices[0])/10,20])
                model = ransac.AffineTransformationModel(rs[0][indices[0]],rs[ibeads][indices[1]])
                data = n.array(range(number))
                try:
                    fit,inliers = ransac.ransac(data,model,4,number*10,50,number/2,return_all=1)
                    tmpWorked = True
                    print 'OK: found parameters'
                    break
                except:
                    print 'try again with different z scale factor...'
                    pass

        if not tmpWorked:
            raise Exception('Could not register bead images.')
        
        # find better transform using redefined matches
        indices2,distances2 = bestMatch(scaleBeads(transform(fit,rs[0]),scale),scaleBeads(rs[ibeads],scale))
        inds = n.where(distances2<5.)[0]
        if len(inds)<5:
            print '#'*10+'\ncheck bead registration since only a few beads are well aligned\n'+'#'*10
            inds = n.array(range(5))
        fit2 = affine_fit.Affine_Fit(rs[0][indices2[0,inds]],rs[ibeads][indices2[1,inds]]).Matrix()
        params.append(fit2)
        indices3,distances3 = bestMatch(transform(fit2,rs[0]),rs[ibeads])
        print 'mean bead distance in pixels after alignment is\n%s' %n.mean(distances3)
    return params


def findParametersLessAccurate(rs,N=4):
    # less accurate
    beads = []
    features = []
    #distancess = []
    params =[[1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.]]
    for ibeads in range(1,len(rs)):
        for scale in n.linspace(1,10,50):
            scale = n.array([scale,1.,1.])
            #print scale
            tmpbeads = scaleBeads(rs[ibeads],scale,verbose=True)
            refbeads = scaleBeads(rs[0],scale)
            f0 = calcFeatures(refbeads,N=N)
            # find initial matches and calculate approximate transform
            tmpf = calcFeatures(tmpbeads,N=N)
            indices,distances = bestMatch(f0,tmpf)
            #distancess.append(distances)
            number = n.max([len(indices[0])/5,40])
            inds = n.array(range(number))
            model = ransac.AffineTransformationModel(rs[0][indices[0]],rs[ibeads][indices[1]])
            data = n.array(range(number))
            try:
                fit,inliers = ransac.ransac(data,model,4,number*10,10,number/2,return_all=1)
                print 'OK: found parameters'
                break
            except:
                print 'try again with different z scale factor...'
                pass
                
        # find better transform using redefined matches
        #indices2,distances2 = bestMatch(scaleBeads(transform(fit,rs[0]),scale),scaleBeads(rs[ibeads],scale))
        #inds = n.where(distances2<5.)[0]
        #fit2 = affine_fit.Affine_Fit(rs[0][indices2[0,inds]],rs[ibeads][indices2[1,inds]]).Matrix()
        params.append(fit)
    return params

def transformStack(p,stack,outShape=None,outSpacing=None,outOrigin=None):
    # can handle composite transformations (len(p)%12)
    # 20140326: added outOrigin option
    numpyarray = False
    if type(stack)==n.ndarray:
        numpyarray = True
        stack = sitk.GetImageFromArray(stack)
    transf = sitk.Transform(3,8)
    if not (p is None):
        for i in range(len(p)/12):
            transf.AddTransform(sitk.Transform(3,6))
        #p = n.array(p)
        #p = n.concatenate([p[:9].reshape((3,3))[::-1,::-1].flatten(),p[9:][::-1]])
        transf.SetParameters(n.array(p,dtype=n.float64))
    if outShape is None: shape = stack.GetSize()
    else:
        shape = n.ceil(n.array(outShape))
        shape = [int(i) for i in shape]
    if outSpacing is None: outSpacing = stack.GetSpacing()
    else: outSpacing = n.array(outSpacing)
    if outOrigin is None: outOrigin = stack.GetOrigin()
    else: outOrigin = n.array(outOrigin)
    print stack.GetSize(),shape
    #pdb.set_trace()
    newim = sitk.Resample(stack,shape,transf,sitk.sitkLinear,outOrigin,outSpacing)
    if numpyarray:
        newim = sitk.GetArrayFromImage(newim)
    return newim

def transformStack2(p,stack,refStack):
    # can handle composite transformations (len(p)%12)
    numpyarray = False
    if type(stack)==n.ndarray:
        numpyarray = True
        stack = sitk.GetImageFromArray(stack)
    transf = sitk.Transform(3,8)
    if p != None:
        for i in range(len(p)/12):
            transf.AddTransform(sitk.Transform(3,6))
    newim = sitk.Resample(stack,refStack.GetSize(),transf,sitk.sitkLinear,refStack.GetOrigin(),refStack.GetSpacing())
    if numpyarray:
        newim = sitk.GetArrayFromImage(newim)
    return newim

def transformStackAndRef(p,stack,refStack):
    # can handle composite transformations (len(p)%12)
    newim = transformStack(p,stack)
    if not (refStack is None):
        newim = sitk.Resample(newim,refStack)
    return newim

def scaleStack(p,stack):
    # expects three values
    p = n.array(p,dtype=n.float64)[::-1]
    numpyarray = False
    if type(stack)==n.ndarray:
        numpyarray = True
        stack = sitk.GetImageFromArray(stack)
    transf = sitk.Transform(3,2)
    transf.SetParameters(p)
    #print tuple((n.array(oldim.GetSize())/p).astype('uint64'))
    newSize = tuple([int(i) for i in n.array(stack.GetSize())/p])
    newim = sitk.Resample(stack,newSize,transf,sitk.sitkLinear)
    #newim = sitk.Resample(oldim,oldim.GetSize(),transf)
    if numpyarray:
        newim = sitk.GetArrayFromImage(newim)
    return newim

if __name__ == "__main__":

    #files = ['Time000001_00_00.tif','Time000001_00_90.tif']
    #files = ['hard/Time000000_00.tif','hard/Time000000_90.tif']
    scaling = n.array([8.,1.,1.])
    files = ['rotation2/Time000000_00.tif','rotation2/Time000000_90.tif']
    angles = [0,n.pi/2]
    ims = loadFiles(files)
    rs = segmentBeads(ims)
    #refbeads = scaleBeads(rs[0],scaling)
    #f0 = calcFeatures(refbeads,N=N)

    beads = []
    features = []
    distancess = []
    params = []
    for ibeads in range(1,len(rs)):
        for scale in n.linspace(1,10,10):
            scale = n.array([scale,1.,1.])
            print scale
            tmpbeads = scaleBeads(rs[ibeads],scale)
            refbeads = scaleBeads(rs[0],scale)
            f0 = calcFeatures(refbeads,N=N)
            # find initial matches and calculate approximate transform
            tmpf = calcFeatures(tmpbeads,N=N)
            indices,distances = bestMatch(f0,tmpf)
            distancess.append(distances)
            number = n.max([len(indices)/10,20])
            #number = len(indices[0])
            inds = n.array(range(number))
            model = ransac.AffineTransformationModel(rs[0][indices[0]],rs[1][indices[1]])
            data = n.array(range(number))
            try:
                fit,inliers = ransac.ransac(data,model,4,number*10,50,number/2,return_all=1)
                break
            except:
                print 'did not work'
                pass
                
        # find better transform using redefined matches
        indices2,distances2 = bestMatch(transform(fit,rs[0]),rs[1])
        number = n.max([len(indices2[0])/2,10])
        inds = n.array(range(number))
        model = ransac.AffineTransformationModel(rs[0][indices2[0]],rs[1][indices2[1]])
        data = n.array(range(number))
        fit2 = ransac.ransac(data,model,number/5,100,2,number/2)
        
        params.append(fit2)

#def scaleImage(im,scaleFactors):
    

#        param = n.concatenate([fit2[:9].reshape((3,3))[::-1,::-1].flatten(),fit2[9:][::-1]])
#        oldim = sitk.GetImageFromArray(ims[1])
#        transf = sitk.Transform(3,6)
#        transf.SetParameters(param)
#        newim = sitk.Resample(oldim,oldim.GetSize(),transf)
#        newim = sitk.GetArrayFromImage(newim)

        #filing.arrayToTif(newim,files[ibeads][:-4]+'_transformed.tif')

    #os.system('./AffineTransform %s 1 1 1 %s %s %s %s %s %s %s %s %s %s %s %s %s 0 0 0 0.000000 1.000000 1.000000 1.000000' %((files[1],'test.tif')+tuple(nparams[:9].reshape((3,3))[::-1,::-1].flatten())+tuple(nparams[9:][::-1])))
