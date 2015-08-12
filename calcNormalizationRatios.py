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
shutil.copy2(configurationFile,os.path.join(config['pipeDir'],configCopyPath))

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
for iviews,views in enumerate(alignParams):
    tmpViews = []
    for imod,mod in enumerate(views):
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
initialParameters = newAlignParams[:,1,1]

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
        if len(config['tps']) > itp+1:
            if onlyExecutePart in [0,1]:
                dataDownloads = handler.downloadDataTimepoint(config['tps'][itp+1])
            else:
                dataDownloads = handler.downloadIntermediateDataTimepoint(config['tps'][itp+1])

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
            tmpims = beads.loadFiles([os.path.join(view,config['filePattern']) %{'t':tp} for view in views])
            if opposingWeights == None: opposingWeights = fusion.calculateOpposingWeights(tmpims[0].GetSize())
            tmpims = fusion.fuseStacks([tmpims[0],tmpims[1]],opposingWeights,[alignParams[iviews][0][0],alignParams[iviews][0][1]])
            tmpims.SetSpacing(stackSpacings[iviews])
            tmpims = filters.constantBackgroundSubtraction(tmpims,config['backgroundLevel'])
            if onlyExecutePart == 1:
                sitk.WriteImage(tmpims,os.path.join(config['alignDirs'][iviews],config['filePattern'] %{'t':tp}))

        # filter image which will serve as input to elastix
        tmpimsFiltered = filters.gaussianBackgroundSubtraction(tmpims,config['gaussianFilterKernelSize'])
        sitk.WriteImage(tmpimsFiltered,tmpFileNameFiltered)
        tempFiles.append(tmpFileNameFiltered)
        tempFiles.append(tmpFileNameFiltered[:-3]+'raw')
        opposingFusions.append(tmpims)

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
        del tmpims,tmpimsFiltered

        # fuse opposing views of additional channels if specified
        if config['processAdditionalTimepoints']:
            for ichannel in range(numberOfAdditionalChannels):
                print 'fusing opposing stacks from additional channel %s' %ichannel
                tmpAddFiles = [view %{'t':tp}.split('.')[0]+'.tif' for view in addDataFilesList[ichannel][iviews]]
                if onlyExecutePart in [2]:
                    tmpims = [sitk.ReadImage(filePath) for filePath in tmpAddFiles]
                else:
                    tmpims = beads.loadFiles([view %{'t':tp} for view in addDataFilesList[ichannel][iviews]])
                tmpims = fusion.fuseStacks([tmpims[0],tmpims[1]],opposingWeights,[alignParams[iviews][0][0],alignParams[iviews][0][1]])
                tmpims = filters.constantBackgroundSubtraction(tmpims,config['backgroundLevel'])
                tmpims.SetSpacing(stackSpacings[iviews])
                if onlyExecutePart in [1]:
                    for iim,im in enumerate(tmpims):
                        sitk.WriteImage(im,tmpAddFiles[iim])
                addOpposingFusionsList[ichannel].append(tmpims)
                del tmpims

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
