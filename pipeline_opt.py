# Main script for spim preprocessing pipeline
# author: Marvin Albert (marvin.albert@gmail.com)

from imports import *

#####################
#####################
#####################
# process arguments

lastArg = sys.argv[-1]
if lastArg[-2:] == 'py':
    configurationFile = lastArg
    onlyExecutePart = 0
elif int(lastArg[-2:]) in [0,1,2]:
    configurationFile = sys.argv[-2]
    onlyExecutePart = int(lastArg[-2:])
else:
    print '\nWrong arguments:\n\n Usage: python pipeline.py <configurationFile.py> [split flag: 0 (execute all), 1 (only opposing view fusion), 2 (only image based fusion)]\n '
    sys.exit()

#configurationFile = sys.argv[-1] # use file passed as argument
print '\n\nstarting preprocessing pipeline using configuration file:\n\n    %s\n' %configurationFile

#####################
#####################
#####################
# read and handle configuration file
configFile = open(configurationFile)
configCode = configFile.read()
configFile.close()
if configCode[:39] == '# SPIM preprocessing configuration file':
    exec(configCode)
else:
    print '\n\n Error: specified configuration file \n\n %s \n\n is not a proper configuration file (first line must state: \'# SPIM preprocessing configuration file\')\n\n' %configurationFile
    sys.exit()

# create handler for configuration and files
handler = extras.handler(config)

# use modified config file
config = handler.config

# copy config file to pipe folder
configCopyPath = 'configuration_%s.py' %handler.timestamp
shutil.copy2(configurationFile,os.path.join(handler.config['pipeDir'],configCopyPath))

# initialize dict to contain time measurements
timeDict = dict()
timeDict['preparation'] = []
timeDict['elastix'] = []
timeDict['finalization'] = []
tmpTimeI = time.time()

#############################
# Bead Alignment
#############################

# check whether bead registration is already calculated, else calculate it

print '\n\n'+30*'-'+'\nBEAD ALIGNMENT start\n'
beadParamsFile = os.path.join(config['pipeDir'],'affineParams.py')
if os.path.exists(beadParamsFile):
    cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],config['pipeDir'])))
    if cmd_subfolder not in sys.path:
        sys.path.insert(0, cmd_subfolder)
    exec('from %s import alignParams') %os.path.basename(beadParamsFile[:-3])
    print 'found alignment parameters in %s' %beadParamsFile
else:
    print 'bead alignment...'
    beadFiles = [[config['beads_LC'][i],config['beads_RC'][i]] for i in range(len(config['beads_LC']))]
    if config['useFileCaching']:
        origBeadFiles = [[handler.origConfig['beads_LC'][i],handler.origConfig['beads_RC'][i]] for i in range(len(config['beads_LC']))]
        targetFiles = [os.path.join(dataDir,config['filePattern'] %{'t':0}) for dataDir in sum(beadFiles,[])]
        sourceFiles = [os.path.join(dataDir,handler.origConfig['filePattern'] %{'t':0}) for dataDir in sum(origBeadFiles,[])]
        beadsDownload = handler.copyFiles(sourceFiles,targetFiles,isDownload=True)
        handler.finalizeProcesses(beadsDownload)
    segmentedBeads = [beads.segmentBeads(
                      beads.loadFiles([os.path.join(config['beads_LC'][i],config['filePattern'] %{'t':0}),
                                       os.path.join(config['beads_RC'][i],config['filePattern'] %{'t':0})]),
                                       config['mexicanHatKernelSize']) for i in range(len(config['beads_LC']))]
    alignParams = []
    for iview,view in enumerate(beadFiles):
        tmpParams = []
        tmpParams.append(list(beads.findParameters(segmentedBeads[iview])))
        tmpParams.append(list(beads.findParameters([segmentedBeads[0][0],segmentedBeads[iview][0]])))
        tmpParams = [[list(i) for i in j] for j in tmpParams]
        alignParams.append(tmpParams)
    commentString = '# parameter structure:\n# p=[rot0,...,rotN]\n# roti=[[param(LC->LC),param(RC->LC)],[param(LC[0]->param(LC[0]),param(LC->LC[0])]]\n\n\n'
    newParamsFile = open(beadParamsFile,'w')
    newParamsFile.write(commentString+'alignParams = '+alignParams.__str__())
    newParamsFile.close()
alignParams = n.array(alignParams)

print '\nBEAD ALIGNMENT end\n'+30*'-'

tmpTimeF = time.time()
timeDict['beads'] = tmpTimeF-tmpTimeI
tmpTimeI = tmpTimeF

############################
# Fuse rotational views
############################

print 30*'-'+'\nFUSIONING start\n'

# reorganizing
paramDict = dict(el=elastixDir)

dataFiles = list([[config['data_LC'][i],config['data_RC'][i]] for i in range(len(config['data_LC']))])
numberOfViews = len(dataFiles)

if config['processAdditionalTimepoints']:
    addDataFilesList = []
    addDataFiles = list([[config['addData_LC'][i],config['addData_RC'][i]] for i in range(len(config['addData_LC']))])
    numberOfAdditionalChannels = len(addDataFiles)/float(numberOfViews)
    if not numberOfAdditionalChannels == int(numberOfAdditionalChannels):
        print '\n\nadditionalChannels not defined correctly\n\n'
        sys.exit()
    numberOfAdditionalChannels = int(numberOfAdditionalChannels)
    for iadd in range(numberOfAdditionalChannels):
        addDataFilesList.append(addDataFiles[iadd*numberOfViews:(iadd+1)*numberOfViews])
        #extras.checkDir(os.path.join(config['fusedDir'],'additional%s'%iadd),'additional%s'%iadd)

relZScales = []
tmpMatrix = n.diag((1.,1.,1.,1.))
for iviews,iview in enumerate(dataFiles):
    if not iviews:
        if len(dataFiles) == 1: relZScale = 1
        else:
            tmpMatrix[:3,:3] = n.array(alignParams[1][1][1])[:9].reshape((3,3))
            tmpScaleFactors = transformations.decompose_matrix(tmpMatrix)[0]
            relZScale = tmpScaleFactors[0] #axis!
    else:
        tmpMatrix[:3,:3] = n.array(alignParams[iviews][1][1])[:9].reshape((3,3))
        tmpScaleFactors = transformations.decompose_matrix(tmpMatrix)[0]
        relZScale = 1./tmpScaleFactors[1] #axis!
    relZScales.append(relZScale)

stackSpacings = []
for i in relZScales:
    stackSpacings.append(n.array([1.,1.,i]))

# eliminate scaling in affine parameters (assume X and Y spacing as 1)
newAlignParams = []
for iviews,views in enumerate(dataFiles):
    tmpViews = []
    for imod,mod in enumerate(alignParams[iviews]):
        tmpMods = []
        for icam,cam in enumerate(mod):
            tmpParams = fusion.invertParams(cam)
            tmpParams = fusion.adaptParametersToSpacing(tmpParams,spacingRef=stackSpacings[0],spacingFinal=stackSpacings[iviews])
            tmpMods.append(tmpParams)
        tmpViews.append(tmpMods)
    newAlignParams.append(tmpViews)
newAlignParams = n.array(newAlignParams)

# define isotropic scaling value in terms of x and y scaling
isotropicSpacing = n.ones(3)*float(config['isotropicOutputSpacing'])

# weights are calculated during first iteration
opposingWeights = None
finalWeights = None

# initialize dictionary for obtained parameters
obtainedParameters = dict()

# initial alignment parameters are taken from beads
if handler.config['ignoreBeadsForRotatedViews']:
    initialParameters = n.array([1.,0.,0,0,1,0,0,0,1,0,0,0])
else:
    initialParameters = newAlignParams[:,1,1]

# rescale factors dict
rescaleFactors = dict()

# loop over timepoints
for itp,tp in enumerate(config['tps']):
    print 'processing timepoint %s' %tp

    tempFiles = []
    if handler.alreadyFused(tp):
        print '\n\tfound existing file in fusion directory, so skipping timepoint %s\n\n' %tp
        continue

    # initialize vars
    opposingFusions = []
    if config['processAdditionalTimepoints']:
        addOpposingFusionsList = [[] for i in range(int(numberOfAdditionalChannels))]
    viewFileNames,viewFileNamesFiltered = [],[]

    # handle files
    if config['useFileCaching']:
        if not globals().has_key('dataDownloads'):
            if onlyExecutePart in [0,1]:
                dataDownloads = handler.downloadDataTimepoint(tp)
            else:
                dataDownloads = handler.downloadIntermediateDataTimepoint(tp)
        handler.finalizeProcesses(dataDownloads)
        #if len(config['tps']) > itp+1:
        #    if onlyExecutePart in [0,1]:
        #        dataDownloads = handler.downloadDataTimepoint(config['tps'][itp+1])
        #    else:
        #        dataDownloads = handler.downloadIntermediateDataTimepoint(config['tps'][itp+1])
        if len(config['tps']) > itp+1:
            for nextTp in range(config['tps'][itp+1],config['tps'][-1]+1):
                if handler.alreadyFused(nextTp):
                    continue
                else:
                    if onlyExecutePart in [0,1]:
                        dataDownloads = handler.downloadDataTimepoint(nextTp)
                    else:
                        dataDownloads = handler.downloadIntermediateDataTimepoint(nextTp)
                    break

    # bead based alignment of opposing stacks, loop over views
    masks = []
    maskFileNames = []
    for iviews,views in enumerate(dataFiles):

        # arrange filenames
        tmpFileName = os.path.join(config['alignDirs'][iviews],config['filePattern'].split('.')[0] %{'t':tp}+'.mhd')
        tmpFileNameFiltered = os.path.join(config['alignDirs'][iviews],config['filePattern'].split('.')[0] %{'t':tp}+'_filtered.mhd')
        #tmpFileNameMask = os.path.join(config['alignDirs'][iviews],config['filePattern'].split('.')[0] %{'t':tp}+'_mask.mhd')
        tmpFileNameMask = os.path.join(config['pipeDir'],config['filePattern'].split('.')[0] %{'t':tp}+'_mask_rot%s.mhd'%iviews)
        maskFileNames.append(tmpFileNameMask)
        viewFileNames.append(tmpFileName)
        viewFileNamesFiltered.append(tmpFileNameFiltered)

        # fuse opposing stacks
        print 'fusing and filtering opposing stacks from view %s' %iviews
        if onlyExecutePart not in [0,1]:
            tmpims = sitk.ReadImage(os.path.join(config['alignDirs'][iviews],config['filePattern'] %{'t':tp}))
            tmpims.SetSpacing(stackSpacings[iviews])
        else:
            tmpSuccess = False
            while not tmpSuccess:
                try:
                    tmpims = beads.loadFiles([os.path.join(view,config['filePattern']) %{'t':tp} for view in views])
                    tmpSuccess = True
                except:
                    time.sleep(5)
                    pass
            #if opposingWeights == None: opposingWeights = fusion.calculateOpposingWeights(tmpims[0].GetSize())
            tmpims = fusion.fuseOpposingStacks(alignParams[iviews][0][0],alignParams[iviews][0][1],tmpims[0],tmpims[1],axis=0)
            #tmpims = fusion.fuseStacks([tmpims[0],tmpims[1]],opposingWeights,[alignParams[iviews][0][0],alignParams[iviews][0][1]])
            tmpims.SetSpacing(stackSpacings[iviews])
            tmpims = filters.constantBackgroundSubtraction(tmpims,config['backgroundLevel'])
            if onlyExecutePart == 1:
                sitk.WriteImage(tmpims,os.path.join(config['alignDirs'][iviews],config['filePattern'] %{'t':tp}))

        opposingFusions.append(tmpims)
        del tmpims

        # fuse opposing views of additional channels if specified
        if config['processAdditionalTimepoints']:
            for ichannel in range(numberOfAdditionalChannels):
                print 'fusing opposing stacks from additional channel %s' %ichannel
                #tmpAddFiles = [(os.path.join(view,handler.config['filePattern']) %{'t':tp}).split('.')[0]+'.tif' for view in addDataFilesList[ichannel][iviews]]
                tmpAddFile = os.path.join(handler.config['additionalAlignDirs'][ichannel][iviews],config['filePattern']%{'t':tp})
                if onlyExecutePart in [2]:
                    tmpims = sitk.ReadImage(tmpAddFile)
                else:
                    tmpSuccess = False
                    while not tmpSuccess:
                        try:
                            tmpims = beads.loadFiles([os.path.join(addDataFiles[iviews][iview],config['filePattern']) %{'t':tp} for iview in range(len(addDataFiles[iviews]))])
                            tmpSuccess = True
                        except:    
                            time.sleep(5)
                            pass
                    #tmpims = beads.loadFiles([os.path.join(addDataFiles[iviews][iview],config['filePattern']%{'t':tp}) for iview in range(len(addDataFiles[iviews]))])
                    #tmpims = beads.loadFiles([view %{'t':tp} for view in addDataFilesList[ichannel][iviews]])
                    #tmpims = fusion.fuseStacks([tmpims[0],tmpims[1]],opposingWeights,[alignParams[iviews][0][0],alignParams[iviews][0][1]])
                    tmpims = fusion.fuseOpposingStacks(alignParams[iviews][0][0],alignParams[iviews][0][1],tmpims[0],tmpims[1])
                    tmpims = filters.constantBackgroundSubtraction(tmpims,config['backgroundLevel'])
                tmpims.SetSpacing(stackSpacings[iviews])
                if onlyExecutePart in [1]:
                    sitk.WriteImage(tmpims,tmpAddFile)
                addOpposingFusionsList[ichannel].append(tmpims)
                del tmpims

        # rescaleBeforeImageBasedFusion
        if handler.config['rescaleBeforeImageBasedFusion'] and onlyExecutePart not in [1]:
            if not rescaleFactors.has_key('rot%s' %iviews):
                tmpFile = open(os.path.join(handler.config['pipeDir'],handler.config['rescaleFactorsFiles'][iviews]))
                rescaleFactors['rot%s' %iviews] = extras.parseRescaleFactors(tmpFile.read())
                tmpFile.close()
            if rescaleFactors['rot%s' %iviews].has_key(str(tp)):
                tmpFactor = rescaleFactors['rot%s' %iviews][str(tp)]
                print 'multiplying rot%s by rescale factor %s' %(iviews,tmpFactor)
                if tmpFactor != 1.:
                    opposingFusions[iviews] = sitk.Cast(sitk.Cast(opposingFusions[iviews],6)*tmpFactor,3)
            if config['processAdditionalTimepoints']:
                for ichannel in range(numberOfAdditionalChannels):
                    if not rescaleFactors.has_key('rot%s_add%s' %(iviews,ichannel)):
                        tmpIndex = (ichannel+1)*numberOfViews+iviews
                        tmpFile = open(os.path.join(handler.config['pipeDir'],handler.config['rescaleFactorsFiles'][tmpIndex]))
                        rescaleFactors['rot%s_add%s' %(iviews,ichannel)] = extras.parseRescaleFactors(tmpFile.read())
                        tmpFile.close()
                    if rescaleFactors['rot%s_add%s' %(iviews,ichannel)].has_key(str(tp)):
                        tmpFactor = rescaleFactors['rot%s_add%s' %(iviews,ichannel)][str(tp)]
                        print 'multiplying rot%s add%s by rescale factor %s' %(iviews,ichannel,tmpFactor)
                        if tmpFactor != 1.:
                            addOpposingFusionsList[ichannel][iviews] = sitk.Cast(sitk.Cast(addOpposingFusionsList[ichannel][iviews],6)*tmpFactor,3)

        if numberOfViews > 1:
            # filter image which will serve as input to elastix
            tmpimsFiltered = filters.gaussianBackgroundSubtraction(opposingFusions[iviews],config['gaussianFilterKernelSize'])
            sitk.WriteImage(tmpimsFiltered,tmpFileNameFiltered)
            tempFiles.append(tmpFileNameFiltered)
            tempFiles.append(tmpFileNameFiltered[:-3]+'raw')
        
            # create mask for elastix
            if not iviews or handler.config['normalizeRotatedViews']:
                print 'creating masks for elastix input'
                if os.path.exists(tmpFileNameMask):
                    mask = sitk.ReadImage(tmpFileNameMask)
                    mask.SetSpacing(n.round(n.array(mask.GetSpacing()),3))
                    #mask = beads.transformStack(None,mask,
                    #            outSpacing=opposingFusions[iviews].GetSpacing(),
                    #            outShape=n.array(opposingFusions[0].GetSize()))
                else:
                    mask = filters.createMask(tmpimsFiltered,config['maskThreshold'])
                    mask.SetSpacing(n.round(n.array(mask.GetSpacing()),3))
                    sitk.WriteImage(mask,tmpFileNameMask)
                mask = beads.transformStack(None,mask,
                                outSpacing=opposingFusions[iviews].GetSpacing(),
                                outShape=n.array(opposingFusions[0].GetSize()))
                if onlyExecutePart in [0]:
                    tempFiles.append(tmpFileNameMask)
                    tempFiles.append(tmpFileNameMask[:-3]+'raw')

                masks.append(mask)
                del mask
            del tmpimsFiltered

    if onlyExecutePart in [1]:
        if config['useFileCaching']:
            if globals().has_key('dataUploads'):
                handler.finalizeProcesses(dataUploads)
                #handler.eliminateTempFiles(origOrFused=True)
            handler.eliminateDataTimepoint(tp)
            handler.eliminateTempFiles(origOrFused=False)
            dataUploads = handler.uploadIntermediateDataTimepoint(tp)
        else:
            if not config['keepTemporaryImageFiles']:
                for file in tempFiles:
                    try:
                        os.remove(file)
                    except:
                        pass
        continue

    # normalize
    if handler.config['normalizeRotatedViews'] and onlyExecutePart in [1,2]:
        normRatioFile = os.path.join(handler.config['pipeDir'],'normalizationRatios.pc')
        # calculate ratio before
        ratioKnown = False
        if os.path.exists(normRatioFile):
            normalizationRatiosBefore = pickle.load(open(normRatioFile))
            if normalizationRatiosBefore.has_key(str(tp)):
                normalizationRatiosBefore = normalizationRatiosBefore[str(tp)]
                ratioKnown = True
        if not ratioKnown:
            if onlyExecutePart in [2]: print '#'*10+'\nratio not known! continuing without normalization\n'+'#'*10
            normalizationRatiosBefore = []
            refValue = n.sum(sitk.GetArrayFromImage(opposingFusions[0]*masks[0]))
            for iview in range(1,len(opposingFusions)):
                currValue = n.sum(sitk.GetArrayFromImage(opposingFusions[iview]*masks[iview]))
                normalizationRatiosBefore.append(currValue/float(refValue))
            if os.path.exists(normRatioFile): 
                ratioDict = pickle.load(open(normRatioFile))
            else:
                ratioDict = dict()
            ratioDict[str(tp)] = normalizationRatiosBefore
            pickle.dump(ratioDict,open(normRatioFile,'w'))
        # calculate ratio after
        normalizationRatiosAfter = []
        refValue = n.sum(sitk.GetArrayFromImage(opposingFusions[0]*masks[0]))
        for iview in range(1,len(opposingFusions)):
            currValue = n.sum(sitk.GetArrayFromImage(opposingFusions[iview]*masks[iview]))
            normalizationRatiosAfter.append(currValue/float(refValue))

        print '\nnormalization ratio before: %s\nnormalization ratio after: %s\n' %(normalizationRatiosBefore,normalizationRatiosAfter)
        # perform normalization
        for iview in range(1,len(opposingFusions)):
            tmpFactor = normalizationRatiosBefore[iview-1]/normalizationRatiosAfter[iview-1]
            print 'multiplying opposing view %s by factor %s' %(iview,tmpFactor)
            opposingFusions[iview] = sitk.Cast(sitk.Cast(opposingFusions[iview],6)*tmpFactor,3)
    else:
        print '#'*10+'\ncontinuing without normalization\n'+'#'*10

    # image based alignment using elastix
    alignedViews = []
    if numberOfViews > 1:
        paramDict['f1'] = viewFileNamesFiltered[0]
        tmpAlignedView = beads.transformStack(newAlignParams[0][1][1],opposingFusions[0],
                                outSpacing=isotropicSpacing,
                                outShape=opposingFusions[0].GetSpacing()/isotropicSpacing*opposingFusions[0].GetSize())
        alignedViews.append(tmpAlignedView)
        if config['processAdditionalTimepoints']:
            addAlignedViewsList = [[] for i in range(numberOfAdditionalChannels)]
            for ichannel in range(numberOfAdditionalChannels):
                tmpAddAlignedView = beads.transformStack(newAlignParams[0][1][1],addOpposingFusionsList[ichannel][0],
                                        outSpacing=isotropicSpacing,
                                        outShape=addOpposingFusionsList[ichannel][0].GetSpacing()/isotropicSpacing*addOpposingFusionsList[ichannel][0].GetSize())
                addAlignedViewsList[ichannel].append(tmpAddAlignedView)

    # loop over views
    finalParamsList = [[1,0,0,0,1,0,0,0,1,0,0,0]]
    for iof,of in enumerate(dataFiles[1:]):

        # set elastix parameters
        parameterPath = os.path.join(os.path.join(config['elastixDirs'][iof+1],'Time%06d' %tp),'elastixParameters.txt')
        initialTransformPath = os.path.join(os.path.join(config['elastixDirs'][iof+1],'Time%06d' %tp),'elastixInitialTransform.txt')
        params = initialParameters[iof+1]
        paramDict['out'] = handler.checkDir(os.path.join(config['elastixDirs'][iof+1],'Time%06d' %tp),'elastix_time',ask=False)
        extras.createInitialTransformFile(stackSpacings[1],params,config['initialTransformTemplateString'],initialTransformPath)
        extras.createParameterFile(stackSpacings[1],initialTransformPath,config['parameterTemplateString'],parameterPath)
        paramDict['params'] = parameterPath
        paramDict['f2'] = viewFileNamesFiltered[iof+1]
        paramDict['initialTransform'] = initialTransformPath
        paramDict['fMask'] = maskFileNames[0]

        tmpTimeF = time.time()
        timeDict['preparation'].append(tmpTimeF-tmpTimeI)
        tmpTimeI = tmpTimeF

        # run elastix
        cmd = ('%(el)s -f %(f1)s -m %(f2)s -fMask %(fMask)s -p %(params)s -t0 %(initialTransform)s -out %(out)s' %paramDict).split(' ')
        #cmd = ('%(el)s -f %(f1)s -m %(f2)s -p %(params)s -t0 %(initialTransform)s -out %(out)s' %paramDict).split(' ')
        subprocess.Popen(cmd).wait()

        tmpTimeF = time.time()
        timeDict['elastix'].append(tmpTimeF-tmpTimeI)
        tmpTimeI = tmpTimeF

        # read output parameters from elastix output file
        rawOutParams = open(os.path.join(paramDict['out'],'TransformParameters.0.txt')).read()
        outParams = rawOutParams.split('\n')[2][:-1].split(' ')[1:]
        outParams = n.array([float(i) for i in outParams])
        outCenterOfRotation = rawOutParams.split('\n')[19][:-1].split(' ')[1:]
        outCenterOfRotation = n.array([float(i) for i in outCenterOfRotation])

        # alternatively: use rigid transform instead of affine
        #affineOutParams = transformations.euler_matrix(outParams[0],outParams[1],outParams[2])

        # correct for different center of rotation used in elastix
        affineOutParams = outParams
        totalTranslate = outParams[-3:] - n.dot(outParams[:9].reshape((3,3)),outCenterOfRotation) + outCenterOfRotation
        affineOutParams = n.concatenate([outParams[:9],totalTranslate],0)
        finalParams = fusion.composeAffineTransformations([affineOutParams,params])

        # align views using the obtained parameters
        tmpAlignedView = beads.transformStack(finalParams,opposingFusions[iof+1],
                            outSpacing=isotropicSpacing,
                            outShape=opposingFusions[0].GetSpacing()/isotropicSpacing*opposingFusions[0].GetSize())
        alignedViews.append(tmpAlignedView)

        # align additional channels if specified
        if config['processAdditionalTimepoints']:
            for ichannel in range(numberOfAdditionalChannels):
                tmpAddAlignedView = beads.transformStack(finalParams,addOpposingFusionsList[ichannel][iof+1],
                            outSpacing=isotropicSpacing,
                            outShape=addOpposingFusionsList[ichannel][0].GetSpacing()/isotropicSpacing*addOpposingFusionsList[ichannel][0].GetSize())
                addAlignedViewsList[ichannel].append(tmpAddAlignedView)

        finalParamsList.append(finalParams)

    if numberOfViews > 1:
        # set initial params for next timepoint to output params from current
        initialParams = n.array(finalParamsList)

        # add obtained parameters to dict
        obtainedParameters['%s' %tp] = n.array(finalParamsList)

        # calculate weight matrices if not already done
        if finalWeights == None:
            finalWeights = fusion.calculateFinalWeightMatrices([newAlignParams[0][1][1],newAlignParams[1][1][1]],
                    opposingFusions[0].GetSpacing()/isotropicSpacing*opposingFusions[0].GetSize(),
                    stackSpacings,isotropicSpacing)

        # fuse the final stack
        print 'fusing final stack for timepoint %s' %tp
        final = fusion.fuseStacks([alignedViews[0],alignedViews[1]],finalWeights)

        # write final stack to file
        final = sitk.GetArrayFromImage(final)
    else:
        final = sitk.GetArrayFromImage(opposingFusions[0])

    tifffile.imsave(os.path.join(config['fusedDir'],config['filePattern'] %{'t':tp}),final)
    del final

    # fuse additional channels
    if config['processAdditionalTimepoints']:
        for ichannel in range(numberOfAdditionalChannels):
            if numberOfViews > 1:
                tmpFinal = fusion.fuseStacks(addAlignedViewsList[ichannel],finalWeights)
                tmpFinal = sitk.GetArrayFromImage(tmpFinal)
            else:
                tmpFinal = sitk.GetArrayFromImage(addOpposingFusionsList[ichannel][0])
            tifffile.imsave(os.path.join(handler.config['additionalFusedDir%s' %ichannel],config['filePattern'] %{'t':tp}),tmpFinal)
        del tmpFinal

    if not config['keepTemporaryImageFiles'] or config['useFileCaching']:
        for file in tempFiles:
            try:
                os.remove(file)
            except:
                pass
    else:
        for iview,view in enumerate(alignedViews):
            view = sitk.GetArrayFromImage(view)
            if onlyExecutePart not in [2]: tifffile.imsave(os.path.join(config['pipeDir'],'opposingWeightsView%s.tif'%iview),sitk.GetArrayFromImage(sitk.Cast(opposingWeights[iview]*n.power(2,15),3)))
            tifffile.imsave(os.path.join(config['pipeDir'],'finalWeightsView%s.tif'%iview),sitk.GetArrayFromImage(sitk.Cast(finalWeights[iview]*n.power(2,15),3)))

    if config['useFileCaching']:
        if globals().has_key('dataUploads'):
            handler.finalizeProcesses(dataUploads)
            #handler.eliminateTempFiles(origOrFused=True)
        handler.eliminateTempFiles(origOrFused=False)
        dataUploads = handler.uploadDataTimepoint(tp)
        if onlyExecutePart in [0]:
            handler.eliminateDataTimepoint(tp)
        elif onlyExecutePart in [2]:
            handler.eliminateIntermediateDataTimepoint(tp)


    tmpTimeF = time.time()
    timeDict['finalization'].append(tmpTimeF-tmpTimeI)
    tmpTimeI = tmpTimeF

# store obtained parameters
securityString = "# spim preprocessing: obtained parameters\n"
parameterSaveString = securityString + "obtainedParameters = "+str(obtainedParameters)
tmpFile = open(os.path.join(config['pipeDir'],'obtainedParameters.py'),'w')
tmpFile.write(parameterSaveString)

# save times into log file
timeProtocolEntries = ['tp','preparation','elastix','finalization','total']
timeDict['tp'] = config['tps']
timeDict['total'] = list(n.sum([n.array(timeDict[i]) for i in timeProtocolEntries[1:-1]],0))
timeProtocolString = 'Summary over elapsed times:\n\n'
timeProtocolString += '\t\t'.join(timeProtocolEntries)
timeProtocolString += '\n'
for row in zip(*tuple([timeDict[i] for i in timeProtocolEntries])):
    timeProtocolString += '\t\t'.join(['%.2f' %j for j in row])+'\n'
timeProtocolFile = open(os.path.join(config['pipeDir'],'times_%s.log' %handler.timestamp),'w')
timeProtocolFile.write(timeProtocolString)
timeProtocolFile.close()

if config['useFileCaching']:
    if 'dataUploads' in globals():
        handler.finalizeProcesses(dataUploads)
    pipeUpload = handler.copyDir(handler.config['pipeDir'],handler.config['outDir'],isDownload=False)
    handler.finalizeProcesses(pipeUpload)
    handler.removeTempDir()

print '\nFUSIONING end\n'+30*'-'
sys.exit()
