from imports import *

def gaussianBackgroundSubtraction(im,radius=50):
	print 'subtracting background using gaussian filter (radius %s)' %radius
	im = im - sitk.Cast(sitk.SmoothingRecursiveGaussian(im,50),3)
	im = im*sitk.Cast(im<60000,3)
	return im
	
def createMask(im,threshold=10):
    im = im>threshold
    #im = sitk.BinaryErode(im)
    im = beads.transformStack(None,im,
                            outSpacing=im.GetSpacing()*n.array([8.,8.,2]),
                            outShape=n.array(im.GetSize())/n.array([8.,8.,2.]))
    im = sitk.Cast(im,3)
    return im

def constantBackgroundSubtraction(im,background):
    im = sitk.Cast(im,sitk.sitkInt32)
    im = im-sitk.Cast((im>0),sitk.sitkInt32)*background
    im = sitk.Abs(im)
    im = sitk.Cast(im,sitk.sitkUInt16)
    return im