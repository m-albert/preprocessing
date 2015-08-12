from imports import *

class AutoIndent(object):
    def __init__(self, stream):
        self.stream = stream
        self.offset = 0
        self.frame_cache = {}

    def indent_level(self):
        i = 0
        base = sys._getframe(2)
        f = base.f_back
        while f:
            if id(f) in self.frame_cache:
                i += 1
            f = f.f_back
        if i == 0:
            # clear out the frame cache
            self.frame_cache = {id(base): True}
        else:
            self.frame_cache[id(base)] = True
        return i

    def write(self, stuff):
        indentation = '  ' * self.indent_level()
        def indent(l):
            if l:
                return indentation + l
            else:
                return l
        stuff = '\n'.join([indent(line) for line in stuff.split('\n')])
        self.stream.write(stuff)

#>>> # Example usage 
#>>>
#>>> def f(x):
#...     print "f(%s)" % x
#...     if x == 0:
#...         return 0
#...     elif x == 1:
#...         return 1
#...     else:
#...         return f(x-1) + f(x-2)
#>>>
#>>> import sys
#>>> sys.stdout = AutoIndent(sys.stdout)

def parseRescaleFactors(filestring):
    lines = filestring.split('\n')
    if not len(lines[-1]): lines.remove(lines[-1])
    resultDict = dict()
    for line in lines:
        line = line.split(' ')
        resultDict[str(int(line[0]))] = float(line[1])
    return resultDict

def createInitialTransformFile(spacing,params,template,outPath):
    spacingString = '\n\n(Spacing %s %s %s)\n' %tuple(spacing)
    paramsString = '\n\n(TransformParameters %s %s %s %s %s %s %s %s %s %s %s %s)\n\n' %tuple(params)
    template = paramsString + spacingString + template
    outFile = open(outPath,'w')
    outFile.write(template)
    outFile.close()
    return

def createParameterFile(spacing,initialTransformFile,template,outPath):
    spacingString = '\n\n(Spacing %s %s %s)\n\n' %tuple(spacing)
    initString = '\n\n(InitialTransformParametersFileName \"%s\")\n\n' %initialTransformFile
    template = initString +spacingString+ template
    outFile = open(outPath,'w')
    outFile.write(template)
    outFile.close()
    return



class handler(object):

    def __init__(self,configDict):
        timestamp = time.localtime()
        self.timestamp = '%02d%02d%02d_%02d%02d%02d' %tuple([timestamp[i] for i in range(6)])
        configDict['elastixDir'] = elastixDir
        self.origConfig = configDict

        # some default configuration entries
        if not self.origConfig.has_key('normalizeRotatedViews'):
            self.origConfig['normalizeRotatedViews'] = False
        if not self.origConfig.has_key('ignoreBeadsForRotatedViews'):
            self.origConfig['ignoreBeadsForRotatedViews'] = False
        if not self.origConfig.has_key('useMask'):
            self.origConfig['useMask'] = True
        if not self.origConfig.has_key('useParametersFromFileWithPath'):
            self.origConfig['useParametersFromFileWithPath'] = ''
        if not self.origConfig.has_key('cropInputImages'):
            self.origConfig['cropInputImages'] = False
        if not self.origConfig.has_key('opposingFusionWithElastix'):
            self.origConfig['opposingFusionWithElastix'] = False
        if not self.origConfig.has_key('useSigmoidalsBetweenViews'):
            self.origConfig['useSigmoidalsBetweenViews'] = True
        if not self.origConfig.has_key('useSigmoidalsBetweenCameras'):
            self.origConfig['useSigmoidalsBetweenCameras'] = True
        if not self.origConfig.has_key('adaptSampleCenterForFusionBetweenViews'):
            self.origConfig['adaptSampleCenterForFusionBetweenViews'] = False
        if not self.origConfig.has_key('h5path'):
            self.origConfig['h5path'] = 't0/channel0'
        if not self.origConfig.has_key('useAverageForFusionBetweenViews'):
            self.origConfig['useAverageForFusionBetweenViews'] = False
        if not self.origConfig.has_key('affineOrRotation'):
            self.origConfig['affineOrRotation'] = True
        if not self.origConfig.has_key('useOnlyBeadParameters'):
            self.origConfig['useOnlyBeadParameters'] = False


        self.method = self.origConfig['fileTransferMethod']
        self.tempFolder = self.origConfig['tempFolder'] %{'timestamp':self.timestamp}
        self.sourceHost = self.origConfig['sourceHost']
        self.downloadProcessID = 0
        self.uploadProcessID = 0
        self.tempFileNameDict = dict()
        self.processCount = 1
        self.processesDict = dict()
        self.processesStatusDict = dict()
        self.lastFolderNumber = 0
        self.prepareFolders()
        self.adapt()
        self.checkForExistingParameters()
        self.temporaryOrigFiles = []
        self.temporaryFusedFiles = []

        return None

    def checkForExistingParameters(self):
        if not len(self.origConfig['useParametersFromFileWithPath']):
            self.useExistingParameters = False
            self.existingParameters = dict()
            return

        # check whether to use existing parameters
        obtainedParameterFile = os.path.join(self.config['pipeDir'],self.config['useParametersFromFileWithPath'])
        if not os.path.exists(obtainedParameterFile):
            raise Exception("the file for existing parameters indicated in config['useParametersFromFileWithPath'] was not found.")

        print '\n'+"#"*10+'\n\n'+'using existing parameters from file:\n\t\t'+obtainedParameterFile+'\n\n'+"#"*10+'\n'
        
        tmpFile = open(obtainedParameterFile)
        tmpString = tmpFile.read()
        tmpFile.close()
        self.useExistingParameters = True
        self.existingParameters = ast.literal_eval(tmpString)
        return

    def checkDir(self,folder,name,ask=False,remote=False):
        #print 'result is %s' %self.fileExists(folder,remote)
        if self.fileExists(folder,remote):
            return folder
        else:
            if ask:
                create = raw_input('The path %(path)s (for configuration entry %(name)s) does not exist.\n Automatically create? y/n (\'n\' means create it yourself in the meantime) ' %{'path':folder,'name':name})
            else:
                create = 'y'
            if create == 'y':
                try:
                    self.makeDirectory(folder,remote)
                    return folder
                except:
                    newfolder = raw_input('Directory %(path)s could not be created (for configuration entry %(name)s). Enter new path: ' %{'path':folder,'name':name})
                    return self.checkDir(newfolder,name,remote)
            else:
                return self.checkDir(folder,name,remote)

    def fileExists(self,path,remote=False):
        returnValue = False
        if self.origConfig['useFileCaching'] and remote:
            if self.origConfig['fileTransferMethod'] == 'scp':
                cmd = ['ssh',self.sourceHost,"""test -e %s && echo 1 || echo 0"""%path]
                test = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=False)
                (test, err) = test.communicate()
                if int(test) == 1: returnValue = True
                else: returnValue = False
            else:
                print 'Error: handler.fileExists not implemented for method'
                sys.exit()
        else:
            returnValue = os.path.exists(path)
        print 'checking if %s exists...%s' %(path,['no','yes'][returnValue])
        return returnValue

    def makeDirectory(self,path,remote=False):
        if self.origConfig['useFileCaching'] and remote:
            if self.origConfig['fileTransferMethod'] == 'scp':
                cmd = ['ssh',self.sourceHost,"""mkdir %s"""%path]
                test = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=False)
                (test, err) = test.communicate()
                if not len(test): return
                else:
                    raise Exception
                #os.system('ssh %s \"mkdir %s\"' %(self.sourceHost,path))
                print 'created path %s on %s machine' %(path,['local (executing)','remote (source)'][int(remote)])
            else:
                print 'Error: handler.makeDir not implemented for method'
                sys.exit()
        else:
            os.mkdir(path)
        return

    def prepareFolders(self):

        if self.origConfig['processAdditionalTimepoints']:
            addDataFiles = list([[self.origConfig['addData_LC'][i],self.origConfig['addData_RC'][i]] for i in range(len(self.origConfig['addData_LC']))])
            dataFiles = list([[self.origConfig['data_LC'][i],self.origConfig['data_RC'][i]] for i in range(len(self.origConfig['data_LC']))])
            self.numberOfAdditionalChannels = int(len(addDataFiles)/len(dataFiles))

        # structure of outDir
        if self.origConfig['useFileCaching']: remote = True
        else: remote = False

        self.origConfig['outDir']                       = self.checkDir(self.origConfig['outDir'],'outDir',remote=remote)

        # pipeline dir for saving results
        self.origConfig['pipeDir']                       = self.checkDir(os.path.join(self.origConfig['outDir'],'pipe'),'pipeDir',remote=remote)

        # directories for intermediate files
        self.origConfig['alignDirs'],self.origConfig['elastixDirs'] = [],[]
        for irot,rot in enumerate(self.origConfig['data_LC']):
            self.origConfig['alignDirs'].append(self.checkDir(os.path.join(self.origConfig['outDir'],'rot'+str(irot)),'alignDirs'+str(irot),remote=remote))
            self.origConfig['elastixDirs'].append(self.checkDir(os.path.join(self.origConfig['alignDirs'][-1],'elastix'),'paramsDir'+str(irot),remote=remote))

        # directory for final fused data
        self.origConfig['fusedDir'] = self.checkDir(os.path.join(self.origConfig['outDir'],'fused'),'fusedDir',remote=remote)

        # folder used for temporary files
        #if self.origConfig['useFileCaching']:
        #    self.checkDir(self.tempFolder,'tempFolder')

        if self.origConfig['processAdditionalTimepoints']:
            self.origConfig['additionalAlignDirs'] = []
            for ichannel in range(self.numberOfAdditionalChannels):
                tmpAddAlignDirs = []
                for irot,rot in enumerate(self.origConfig['data_LC']):
                    tmpAddAlignDirs.append(self.checkDir(os.path.join(self.origConfig['alignDirs'][irot],'additional%s' %ichannel),'addAlignDir'+str(irot),remote=remote))
                self.origConfig['additionalFusedDir%s' %ichannel] = self.checkDir(os.path.join(self.origConfig['fusedDir'],'additional%s'%ichannel),'additional',remote=remote)
                self.origConfig['additionalAlignDirs'].append(tmpAddAlignDirs)
        return True

    def newFolderName(self):
        self.lastFolderNumber += 1
        return os.path.join(self.tempFolder,'folder%04d' % self.lastFolderNumber)

#    def copyPipeDir(self,direction):
#        if direction == 'download':
#            sourceDir = self.origConfig['pipeDir']
#            targetDir = self.config['pipeDir']
#        else:
#            sourceDir = self.config['pipeDir']
#            targetDir = self.origConfig['pipeDir']

#        if sourceDir != targetDir:
#            fileList = os.listdir(sourceDir)
#            sourceFileList = [os.path.join(self.origConfig['pipeDir'],fileName) for fileName in fileList]
#            targetFileList = [os.path.join(self.config['pipeDir'],fileName) for fileName in fileList]
#            self.copyFiles(sourceFileList,targetFileList,isDownload = None)
#        return

    def adapt(self):
        self.config = copy.deepcopy(self.origConfig)
        if self.origConfig['useFileCaching']:
            self.tempFolder = self.checkDir(self.tempFolder,'tempFolder',ask=False)
            confStrings = ['beads_LC','beads_RC','data_LC','data_RC','addData_LC','addData_RC','alignDirs','elastixDirs']
            if self.origConfig['processAdditionalTimepoints']:
                for ichannel in range(self.numberOfAdditionalChannels):
                    #confStrings += ['additionalFusedDir%s' %ichannel]
                    self.config['additionalFusedDir%s' %ichannel] = self.checkDir(self.newFolderName(),'addfused',ask=False)
                    for irot,rot in enumerate(self.origConfig['data_LC']):
                        self.config['additionalAlignDirs'][ichannel][irot] = self.checkDir(self.newFolderName(),'tmp',ask=False)
            for confString in confStrings:
                self.config[confString] = []
                for directory in self.origConfig[confString]:
                    self.config[confString].append(self.checkDir(self.newFolderName(),'tmp',ask=False))
            self.config['pipeDir'] = self.checkDir(self.newFolderName(),'temppipe',ask=False)
            self.config['fusedDir'] = self.checkDir(self.newFolderName(),'tempfused',ask=False)
            self.tempPaths = sum([self.config[confString] for confString in confStrings],[])
            self.resultPath = self.config['fusedDir'] # else set to None
            pipeDownload = self.copyDir(self.origConfig['pipeDir'],self.config['pipeDir'],isDownload = True)
            self.finalizeProcesses(pipeDownload)
            self.config['pipeDir'] = os.path.join(self.config['pipeDir'],os.path.basename(self.origConfig['pipeDir']))
            #self.copyPipeDir('download')
        return None

    def alreadyFused(self,timepoint):

        return self.fileExists(os.path.join(self.origConfig['fusedDir'],self.origConfig['outFilePattern'] %{'t':timepoint}),remote=self.origConfig['useFileCaching'])

    def downloadTimePoint(self):
        return True

    def copyDir(self,source,target,isDownload):
        print 'copying %s to %s' %(source,target)
        if self.method == 'scp':
            if sys.platform[:3] == 'lin':
                print 'scp lin from %s to %s' %(source,target)
                if isDownload:
                    source = self.sourceHost+':'+source
                else:
                    target = self.sourceHost+':'+target
                processOSID = subprocess.Popen(['scp','-r',source,target])
            if sys.platform[:3] == 'win':
                print 'scp win from %s to %s' %(source,target)
        elif self.method == 'mount':
            if sys.platform[:3] == 'lin':
                print 'mount lin from %s to %s' %(source,target)
                processOSID = subprocess.Popen(['cp','-r',source,target])
            if sys.platform[:3] == 'win':
                print 'mount win from %s to %s' %(source,target)
        processID = self.processCount
        self.processCount += 1
        self.processesDict[str(processID)] = processOSID
        return [processID]

    def copyFileInSubprocess(self,source,target,isDownload,markAsTemp=True):
        processID = self.processCount
        self.processCount += 1
        if self.method == 'scp':
            if sys.platform[:3] == 'lin':
                print 'scp lin from %s to %s' %(source,target)
                if isDownload:
                    source = self.sourceHost+':'+source
                else:
                    target = self.sourceHost+':'+target
                processOSID = subprocess.Popen(['scp',source,target])
            if sys.platform[:3] == 'win':
                print 'scp win from %s to %s' %(source,target)
        elif self.method == 'mount':
            if sys.platform[:3] == 'lin':
                print 'mount lin from %s to %s' %(source,target)
                processOSID = subprocess.Popen(['cp',source,target])
            if sys.platform[:3] == 'win':
                print 'mount win from %s to %s' %(source,target)
        print 'copying %s to %s' %(source,target)
        self.processesDict[str(processID)] = processOSID
        if markAsTemp:
            if int(isDownload) == 1: self.temporaryOrigFiles.append(target)
            elif int(isDownload) == 0: self.temporaryFusedFiles.append(source)
        return processID

    def copyFiles(self,sourceFileList,targetFileList,isDownload,markAsTemp=True):
        processIDs = []
        for ii,i in enumerate(sourceFileList):
            processIDs.append(self.copyFileInSubprocess(i,targetFileList[ii],isDownload,markAsTemp))
        return processIDs

    def downloadDataTimepoint(self,tp):
        dataFiles = list([[self.config['data_LC'][i],self.config['data_RC'][i]] for i in range(len(self.config['data_LC']))])
        origDataFiles = list([[self.origConfig['data_LC'][i],self.origConfig['data_RC'][i]] for i in range(len(self.origConfig['data_LC']))])
        targetFiles = [dataDir %{'t':tp} for dataDir in sum(dataFiles,[])]
        sourceFiles = [dataDir %{'t':tp} for dataDir in sum(origDataFiles,[])]
        dataDownloads = self.copyFiles(sourceFiles,targetFiles,isDownload=True)
        if self.config['processAdditionalTimepoints']:
            addDataFiles = list([[self.config['addData_LC'][i],self.config['addData_RC'][i]] for i in range(len(self.config['addData_LC']))])
            origAddDataFiles = list([[self.origConfig['addData_LC'][i],self.origConfig['addData_RC'][i]] for i in range(len(self.config['addData_LC']))])
            for ichannel in range(self.numberOfAdditionalChannels):
                targetFiles = [dataDir %{'t':tp} for dataDir in sum(addDataFiles,[])]
                sourceFiles = [dataDir %{'t':tp} for dataDir in sum(origAddDataFiles,[])]
                dataDownloads += self.copyFiles(sourceFiles,targetFiles,isDownload=True)
        return dataDownloads

    def uploadDataTimepoint(self,tp):
        sources = []
        targets = []
        # fused tp
        sources += [os.path.join(self.config['fusedDir'],self.config['outFilePattern'] %{'t':tp})]
        targets += [os.path.join(self.origConfig['fusedDir'],self.origConfig['outFilePattern'] %{'t':tp})]
        # additional channels
        if self.config['processAdditionalTimepoints']:
            for ichannel in range(self.numberOfAdditionalChannels):
                sources += [os.path.join(self.config['additionalFusedDir%s'%ichannel],self.config['outFilePattern'] %{'t':tp})]
                targets += [os.path.join(self.origConfig['additionalFusedDir%s'%ichannel],self.origConfig['outFilePattern'] %{'t':tp})]
        # perform uploads
        dataUploads = self.copyFiles(sources,targets,isDownload=False)
        return dataUploads

    def uploadIntermediateDataTimepoint(self,tp):
        sources = []
        targets = []
        dataFiles = list([[self.config['data_LC'][i],self.config['data_RC'][i]] for i in range(len(self.config['data_LC']))])
        for iviews,views in enumerate(dataFiles):
            sources += [os.path.join(self.config['alignDirs'][iviews],self.config['outFilePattern'] %{'t':tp})]
            targets += [os.path.join(self.origConfig['alignDirs'][iviews],self.origConfig['outFilePattern'] %{'t':tp})]

        if self.origConfig['processAdditionalTimepoints']:
            for ichannel in range(self.numberOfAdditionalChannels):
                #confStrings += ['additionalFusedDir%s' %ichannel]
                for irot,rot in enumerate(self.origConfig['data_LC']):
                    sources += [os.path.join(self.config['additionalAlignDirs'][ichannel][irot],self.config['outFilePattern'] %{'t':tp})]
                    targets += [os.path.join(self.origConfig['additionalAlignDirs'][ichannel][irot],self.origConfig['outFilePattern'] %{'t':tp})]

        # perform uploads
        dataUploads = self.copyFiles(sources,targets,isDownload=False)
        return dataUploads

    def downloadIntermediateDataTimepoint(self,tp):
        sources = []
        targets = []
        dataFiles = list([[self.config['data_LC'][i],self.config['data_RC'][i]] for i in range(len(self.config['data_LC']))])
        for iviews,views in enumerate(dataFiles):
            targets += [os.path.join(self.config['alignDirs'][iviews],self.config['outFilePattern'] %{'t':tp})]
            sources += [os.path.join(self.origConfig['alignDirs'][iviews],self.origConfig['outFilePattern'] %{'t':tp})]

        if self.origConfig['processAdditionalTimepoints']:
            for ichannel in range(self.numberOfAdditionalChannels):
                for irot,rot in enumerate(self.origConfig['data_LC']):
                    targets += [os.path.join(self.config['additionalAlignDirs'][ichannel][irot],self.config['outFilePattern'] %{'t':tp})]
                    sources += [os.path.join(self.origConfig['additionalAlignDirs'][ichannel][irot],self.config['outFilePattern'] %{'t':tp})]

        # perform uploads
        dataDownloads = self.copyFiles(sources,targets,isDownload=True)
        return dataDownloads

    def uploadDataTimepointOld(self,tp,removeAfter=False):
        dataUploads = []
        dataUploads += self.copyFiles([os.path.join(self.config['fusedDir'],self.config['outFilePattern'] %{'t':tp})],
                                         [os.path.join(self.origConfig['fusedDir'],self.origConfig['outFilePattern'] %{'t':tp})],isDownload=2)
        if self.config['processAdditionalTimepoints']:
            for ichannel in range(self.numberOfAdditionalChannels):
                dataUploads += [self.copyFiles([os.path.join(self.config['additionalFusedDir%s'%ichannel],self.config['outFilePattern'] %{'t':tp})],
                                         [os.path.join(self.origConfig['additionalFusedDir%s'%ichannel],self.origConfig['outFilePattern'] %{'t':tp})],isDownload=2)]
        return dataUploads

    def eliminateDataTimepoint(self,tp):
        fileList = []
        dataFiles = list([[self.config['data_LC'][i],self.config['data_RC'][i]] for i in range(len(self.config['data_LC']))])
        fileList += [dataDir %{'t':tp} for dataDir in sum(dataFiles,[])]
        if self.config['processAdditionalTimepoints']:
            addDataFiles = list([[self.config['addData_LC'][i],self.config['addData_RC'][i]] for i in range(len(self.config['addData_LC']))])
            for ichannel in range(self.numberOfAdditionalChannels):
                fileList += [dataDir %{'t':tp} for dataDir in sum(addDataFiles,[])]
        for filePath in fileList:
            os.remove(filePath)
        return

    def eliminateIntermediateDataTimepoint(self,tp):
        fileList = []
        dataFiles = list([[self.config['data_LC'][i],self.config['data_RC'][i]] for i in range(len(self.config['data_LC']))])
        for iviews,views in enumerate(dataFiles):
            fileList += [os.path.join(self.config['alignDirs'][iviews],self.config['outFilePattern'] %{'t':tp})]

        if self.origConfig['processAdditionalTimepoints']:
            for ichannel in range(self.numberOfAdditionalChannels):
                for irot,rot in enumerate(self.origConfig['data_LC']):
                    fileList += [os.path.join(self.config['additionalAlignDirs'][ichannel][irot],self.config['outFilePattern'] %{'t':tp})]

        for filePath in fileList:
            os.remove(filePath)
        return

    def download(self,sourceFileList,targetFileList):
        targetFileList = [os.path.join(self.tempFolder,os.path.basename(file)) for file in sourceFileList]
        processIDs = []
        for isource,source in enumerate(sourceFileList):
            tmpProcessID = self.copyFileInSubprocess(source,targetFileList[isource])
            processIDs.append(tmpProcessID)
            self.tempFileNameDict[str(tmpProcessID)] = targetFileList[isource]
        return targetFileList

#    def upload(self,sourceFileList,targetFileList):
#        processID = copyFilesInSubprocess(sourceFileList,targetFileList)
#        self.tempFileNameDict[str(processID)] = sourceFileList
#        return True

#    def upload(self,sourceFileList,targetFileList):
#        targetFileList = [os.path.join(self.tempFolder,os.path.basename(file)) for file in sourceFileList]
#        processIDs = []
#        for isource,source in enumerate(sourceFileList):
#            tmpProcessID = copyFilesInSubprocess(source,targetFileList[isource])
#            processIDs.append(tmpProcessID)
#            self.tempFileNameDict[str(tmpProcessID)] = targetFileList[isource]
#        return targetFileList

    def eliminateTempFiles(self,origOrFused=True):
        print 'eliminating temp files:'
        if origOrFused:
            for i in self.temporaryOrigFiles:
                print '   %s' %i
                try:
                    os.remove(i)
                except: pass
            self.temporaryOrigFiles = []
        else:
            for i in self.temporaryFusedFiles:
                print '   %s' %i
                try:
                    os.remove(i)
                except: pass
            self.temporaryFusedFiles = []
        return True

    def waitForProcess(self,processID):
        print 'waiting for process %s' %processID
        processOSID = self.processesDict[str(processID)]
        processOSID.wait()
        return True

    def finalizeProcesses(self,processIDs):
        for processID in processIDs:
            self.finalizeProcess(processID)
        return True

    def finalizeProcess(self,processID):
        self.waitForProcess(processID)
        return True

    def removeTempDir(self):
        print 'removing all temporary files...'
        def recursiveRemove(folder):
            for ii,i in enumerate(os.listdir(folder)):
                path = os.path.join(folder,i)
                if os.path.isdir(path): recursiveRemove(path)
                else: os.remove(path)
            os.rmdir(folder)
            return
        recursiveRemove(self.tempFolder)
        return
