# SPIM preprocessing configuration file
# (don't change the first line)
#
# The configuration file is divided into three sections:
#   1) required user input
#   2) additional options (not necessary for standard run)
#   3) the parameters passed on to elastix (generally remain unchanged)


import os,time
#from extras import checkDir
config = dict()

############################################
############################################
############################################
######### start of required user input

###########
# input
###########

# stack definitions:
#     LC has high contrast at the end of the stack
#     RC has high contrast at the start of the stack
#
#
# file names:
#   %(t)0xd replaces the time with x digits

base = '/media/vip2data/malbert/registration/20120803/Data_scaled'

# h5 hierarchy path
config['h5path'] = 'Data'

config['beads_LC'],config['beads_RC'] = [],[]
config['beads_LC'].append(base+'/LC/Beads/Stack0000/Cam_Left_%(t)05d.h5')
config['beads_LC'].append(base+'/LC/Beads/Stack0001/Cam_Left_%(t)05d.h5')
config['beads_RC'].append(base+'/RC/Beads/Stack0000/Cam_Right_%(t)05d.h5')
config['beads_RC'].append(base+'/RC/Beads/Stack0001/Cam_Right_%(t)05d.h5')

# timepoints
config['data_LC'],config['data_RC'] = [],[]
config['data_LC'].append(base+'/LC/20120803/Stack0000/Cam_Left_%(t)05d.h5')
config['data_LC'].append(base+'/LC/20120803/Stack0001/Cam_Left_%(t)05d.h5')
config['data_RC'].append(base+'/RC/20120803/Stack0000/Cam_Right_%(t)05d.h5')
config['data_RC'].append(base+'/RC/20120803/Stack0001/Cam_Right_%(t)05d.h5')


# additional timepoints to process
# (no separate image based alignment, f.e. other channel)
config['processAdditionalTimepoints'] = False
config['addData_LC'],config['addData_RC'] = [],[]
config['addData_LC'].append(base+'/LC/20120803/Stack0000/Cam_Left_%(t)05d.h5')
config['addData_LC'].append(base+'/LC/20120803/Stack0001/Cam_Left_%(t)05d.h5')
config['addData_RC'].append(base+'/RC/20120803/Stack0000/Cam_Right_%(t)05d.h5')
config['addData_RC'].append(base+'/RC/20120803/Stack0001/Cam_Right_%(t)05d.h5')


#########
# output
#########

config['outDir'] = os.path.join(base,'Result6') # enough to set this variable to an empty folder
config['outFilePattern'] = 'Time%(t)06d_00.tif' # only tif supported currently

#########
# settings
# (generally only timepoints to be processed need to be indicated)
#########

# timepoints to process
# (lists are written as:
#      [1,34,56] or
#      range(4) meaning [0,1,2,3] or
#      range(1,4) meaning [1,2,3] or
#      range(1,6,2) meaning [1,3,5])

config['tps'] = [0]

# debugging

config['keepTemporaryImageFiles']       = True   # referring to: opposing stacks (unaligned and aligned), elastix masks

# image related
config['useSigmoidalsBetweenViews']     = True   # fusion between rotational views: if True, sigmoids are used for only taking into account the optimal views,
                                                 # otherwise a simple sum of the views is performed
config['useSigmoidalsBetweenCameras']   = True   # fusion between opposing camera views: if True, sigmoids are used for only taking into account the optimal views,
                                                 # otherwise a simple sum of the views is performed
config['adaptSampleCenterForFusionBetweenViews'] = False
                                                 # if True, for every plane orthogonal to the rotation axis
                                                 # the sample center is calculated and used as a basis
                                                 # for sigmoidal fusioning
config['useAverageForFusionBetweenViews'] = False # ignores the above options for between views and fusions
                                                 # are the averages of the views
config['outputImageFormat']             = 'tif'  # only format currently supported is tif
config['isotropicOutputSpacing']        = 2.     # unit: distance between pixels in x-y plane (generally 0.26um)
                                                 # (i.e.: if isotropicOutPutSpacing = 2   ==>
                                                 # final spacing is 2*0.26um = 0.52um
                                                 # and original stack of size (2000,1000,200) ==>
                                                 # (1000,500,500) (if spacing_z=5*spacing_y)
config['cropInputImages']               = False  # crop all images when loading in (not for flag 2, there images
                                                 # are assumed to already be cropped)
                                                 # axes are ordered by descending length, so typically (z,y,x)
config['cropOffsetLC']                  = (10,10,10)
config['cropSizeLC']                    = (400,150,30)
config['cropOffsetRC']                  = (10,10,10)
config['cropSizeRC']                    = (400,150,30)

# Optionally use previously obtained parameters.
# At the end of a run a file containing the final fusion parameters is written
# to the pipe folder ('obtainedParameters_<timestamp>.py).
# Specify such a file here to omit image based registration and
# use the previously obtained parameters for fusion.
# For disabling this option leave the string empty.
# (file name is interpreted relative to pipe directory)

config['useParametersFromFileWithPath'] = ''#obtainedParameters_20140224_155453.py'


#########
# file handling over the network
#    useful for:
#       1) running fusion on cluster or
#       2) to precache files
#########

config['useFileCaching'] = False                 # for standard use indicate False
config['fileTransferMethod'] = 'scp'            # possible options:
                                                # 1) 'scp': copy over ssh (key pairs must
                                                #           be installed on source and local host, see README.txt)
                                                # 2) 'mount': standard copy process
config['sourceHost'] = 'hufnagel-vip2.pool'           # host name for ssh connection, file paths above should be
                                                # indicated respective to the host file system

config['tempFolder'] = '/tmp/fusion_%(timestamp)s'    # folder on local host for temporarily storing files


######### end of required user input
###############################################
###############################################
###############################################



###############################################
###############################################
###############################################
# additional options

config['ignoreBeadsForRotatedViews']    = False  # set initial parameters for image based registrat$
                                                 # using beads or unity matrix
config['normalizeRotatedViews']         = False  # normalize images using the sum over
                                                 # all pixels under elastix mask
                                                 # in case of onlyExecutePart!=0 only one mask is used
config['mexicanHatKernelSize']          = 2      # only for segmenting beads (unit: pixels), unit: multiples of x-y spacing
config['backgroundLevel']               = 100    # is subtracted from fusions of opposing stacks
config['gaussianFilterKernelSize']      = 50     # unit: multiples of x-y spacing
config['maskThreshold']                 = 5      # elastix only considers
                                                 # pixels satisfying (im-gaussianfilter(im,gaussianFilterKernelSize))>maskThreshold
config['useMask']                       = False  # whether to use mask for faster image based fusion
config['useNoMaskForFirstIteration']    = True   # get good registration results for first iteration which are then
                                                 # used as starting point for the following iterations (possibly without mask)
                                                 # (option doesn't conflict with useMask = True)
config['opposingFusionWithElastix']     = False  # ignore beads and perform image based registration for opposing views
config['opposingFusionWithElastix_startAndStopZPlane'] = None
                                                 # start and stop indices for planes used for elastix registration
                                                 # example: [80,120] for planes 80 to 120
                                                 #          None     for all planes
config['opposingFusionWithElastix_eliminateZContribution'] = False
                                                 # eliminate z contribution from parameters obtained
                                                 # from elastix for opposing fusion (True or False)
config['affineOrRotation']              = True   # True: searches for affine transformations between views (12 parameters)
                                                 # False: searches for rotations (6 parameters)
                                                 # default is True


config['rescaleBeforeImageBasedFusion'] = False  # rescale fused opposing views before image based fusion

# LIST of files containing factors for rescaling:
# there is one file for each of the n rotational views for each channel
# [ <factors_rot0>,...,<factors_rotn>,<factors_add0_0>, ...,<factores_add0_n>, ..., <factors_addm_0>, ...,<factores_addm_n>]
# filenames are relative to pipe directory (and should be located there)
# the files contain lines organized as:
# absolutetime factor

config['rescaleFactorsFiles']           = ['rot0f.txt','rot1f.txt','rot0f_add0.txt','rot1f_add0.txt']   # list of file names


# end of additional options
###############################################
###############################################
###############################################



###############################################
###############################################
###############################################
# parameters passed on to elastix

# template for initial transformation applied before actual registration
# (the actual affine parameters are added on the run)
config['initialTransformTemplateString'] = """
(Transform "AffineTransform")
(NumberOfParameters 12)

(HowToCombineTransforms "Compose")

(InitialTransformParametersFileName "NoInitialTransform")

// Image specific
(FixedImageDimension 3)
(MovingImageDimension 3)
(FixedInternalImagePixelType "short")
(MovingInternalImagePixelType "short")
//(UseDirectionCosines "false")



(CenterOfRotationPoint 0 0 0)
"""

# template for main parameters regarding the registration itself
# (slightly modified on the run)
config['parameterTemplateStringAffine'] = """
// Description: affine

(GradientMagnitudeTolerance 1e-6)
(NumberOfResolutions 3)

//(ImagePyramidSchedule  8 8 2  4 4 1  2 2 1 )
(ImagePyramidSchedule  30 30 30  10 10 10  4 4 4)

//ImageTypes
(FixedInternalImagePixelType "short")
(FixedImageDimension 3)
(MovingInternalImagePixelType "short")
(MovingImageDimension 3)
(UseDirectionCosines "false")

//Components
(Registration "MultiResolutionRegistration")
(FixedImagePyramid "FixedRecursiveImagePyramid")
//(Registration "MultiMetricMultiResolutionRegistration")
(MovingImagePyramid "MovingRecursiveImagePyramid")
(Interpolator "BSplineInterpolator")
(Metric "AdvancedMattesMutualInformation")
//(Metric "AdvancedKappaStatistic")
//(Metric0Weight 1.0) (Metric1Weight 1.0)
//(Metric "AdvancedMeanSquares")
//(Optimizer "AdaptiveStochasticGradientDescent")


(Optimizer "QuasiNewtonLBFGS")
//(Optimizer ConjugateGradient)

//(StopIfWolfeNotSatisfied "false")

(ResampleInterpolator "FinalBSplineInterpolator")
(Resampler "DefaultResampler")
//(Resampler "CUDAResampler")
(Transform "AffineTransform")
//(CenterOfRotationPoint 0.0 0.0 0.0)

(ErodeMask "false" )

(HowToCombineTransforms "Compose")
(AutomaticTransformInitialization "false")
(Scales 10000 10000000 10000000 10000000 10000 10000000 10000000 10000000 10000 1 1 1)
//(Scales 1000 1000 1000 1 1 1)
(AutomaticScalesEstimation "false")
//(AutomaticTransformInitializationMethod "CenterOfGravity" )

(WriteTransformParametersEachIteration "false")
(WriteResultImage "false")
(CompressResultImage "false")
(WriteResultImageAfterEachResolution "false")
(ShowExactMetricValue "false")

//Maximum number of iterations in each resolution level:
//(MaximumNumberOfIterations 500)

//Number of grey level bins in each resolution level:
(NumberOfHistogramBins 32 )
//(FixedLimitRangeRatio 0.0)
//(MovingLimitRangeRatio 0.0)
//(FixedKernelBSplineOrder 3)
//(MovingKernelBSplineOrder 3)

//Number of spatial samples used to compute the mutual information in each resolution level:
//(ImageSampler "RandomCoordinate")
//(ImageSampler "RandomSparseMask")
(ImageSampler "Full")
//(SampleGridSpacing 1)
//(NumberOfSpatialSamples 5000)
//(NewSamplesEveryIteration "true")
//(CheckNumberOfSamples "true")
//(MaximumNumberOfSamplingAttempts 10)

//Order of B-Spline interpolation used in each resolution level:
(BSplineInterpolationOrder 1)

(FixedImageBSplineInterpolationOrder 1)
(MovingImageBSplineInterpolationOrder 1)

//Order of B-Spline interpolation used for applying the final deformation:
(FinalBSplineInterpolationOrder 1)

//Default pixel value for pixels that come from outside the picture:
(DefaultPixelValue 0)

//(MaximumStepLength 4.0)
(ResultImagePixelType "short")
(ResultImageFormat "mhd")
"""

config['parameterTemplateStringRotation'] = """
// Description: euler

(GradientMagnitudeTolerance 1e-6)
(NumberOfResolutions 3)

//(ImagePyramidSchedule  8 8 2  4 4 1  2 2 1 )
(ImagePyramidSchedule  30 30 30  10 10 10  4 4 4)

//ImageTypes
(FixedInternalImagePixelType "short")
(FixedImageDimension 3)
(MovingInternalImagePixelType "short")
(MovingImageDimension 3)
(UseDirectionCosines "false")

//Components
(Registration "MultiResolutionRegistration")
(FixedImagePyramid "FixedRecursiveImagePyramid")
//(Registration "MultiMetricMultiResolutionRegistration")
(MovingImagePyramid "MovingRecursiveImagePyramid")
(Interpolator "BSplineInterpolator")
(Metric "AdvancedMattesMutualInformation")
//(Metric "AdvancedKappaStatistic")
//(Metric0Weight 1.0) (Metric1Weight 1.0)
//(Metric "AdvancedMeanSquares")
//(Optimizer "AdaptiveStochasticGradientDescent")


(Optimizer "QuasiNewtonLBFGS")
//(Optimizer ConjugateGradient)

//(StopIfWolfeNotSatisfied "false")

(ResampleInterpolator "FinalBSplineInterpolator")
(Resampler "DefaultResampler")
//(Resampler "CUDAResampler")
(Transform "EulerTransform")
//(CenterOfRotationPoint 0.0 0.0 0.0)

(ErodeMask "false" )

(HowToCombineTransforms "Compose")
(AutomaticTransformInitialization "false")
(Scales 1000 1000 1000 1 1 1)
(AutomaticScalesEstimation "false")
//(AutomaticTransformInitializationMethod "CenterOfGravity" )

(WriteTransformParametersEachIteration "false")
(WriteResultImage "false")
(CompressResultImage "false")
(WriteResultImageAfterEachResolution "false")
(ShowExactMetricValue "false")

//Maximum number of iterations in each resolution level:
//(MaximumNumberOfIterations 500)

//Number of grey level bins in each resolution level:
(NumberOfHistogramBins 32 )
//(FixedLimitRangeRatio 0.0)
//(MovingLimitRangeRatio 0.0)
//(FixedKernelBSplineOrder 3)
//(MovingKernelBSplineOrder 3)

//Number of spatial samples used to compute the mutual information in each resolution level:
//(ImageSampler "RandomCoordinate")
//(ImageSampler "RandomSparseMask")
(ImageSampler "Full")
//(SampleGridSpacing 1)
//(NumberOfSpatialSamples 5000)
//(NewSamplesEveryIteration "true")
//(CheckNumberOfSamples "true")
//(MaximumNumberOfSamplingAttempts 10)

//Order of B-Spline interpolation used in each resolution level:
(BSplineInterpolationOrder 1)

(FixedImageBSplineInterpolationOrder 1)
(MovingImageBSplineInterpolationOrder 1)

//Order of B-Spline interpolation used for applying the final deformation:
(FinalBSplineInterpolationOrder 1)

//Default pixel value for pixels that come from outside the picture:
(DefaultPixelValue 0)

//(MaximumStepLength 4.0)
(ResultImagePixelType "short")
(ResultImageFormat "mhd")
"""
