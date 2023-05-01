
"""
Image analysis to find the best threshold/sementation for my 100X
images on O'Toole lab scope. 
"""
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
import random
import plotnine as p9
import joypy
import cv2
import pandas as pd
from numpy.lib.histograms import histogram
import os
import scipy


import skimage.filters
import skimage.io
import skimage.morphology
import skimage.exposure
from skimage import data
from skimage.morphology import disk
from skimage.filters import threshold_otsu, rank
from skimage.filters import threshold_local
from skimage.util import img_as_ubyte
import matplotlib
from scipy import ndimage as ndi
from skimage.color import label2rgb, rgb2gray

from skimage import (
    filters, measure, morphology, segmentation
)

# To allow for inline rendering of plots.
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']
sns.set(style='whitegrid', palette=colors, rc={'axes.labelsize': 16})


import skimage.morphology
import skimage.segmentation
import skimage.io
import skimage.filters

import seaborn as sns
import scipy.stats as st
import scipy.stats
colors = ['#440154FF', '#414487ff', '#2A788EFF', '#22A884FF',
          '#7AD151FF', '#FDE725FF', '#e377c2', '#7f7f7f',
          '#bcbd22', '#17becf']

sns.set(style='whitegrid', palette=colors, rc={'axes.labelsize': 16})

# The following is specific Jupyter notebooks
# %matplotlib inline
# %config InlineBackend.figure_formats = {'png', 'retina'}
# %%
########################################################################
"""
#Functions
"""

########################################################################


def normalize_im(im):

    im_norm = (im - im.min()) / (im.max() - im.min())
    return im_norm

##############################


def local_otsu_im(im, currRadius):
    # convert im to unsigned byte format [0,255]
    img = img_as_ubyte(im)

    # define local neighborhood
    radius = currRadius
    selem = disk(radius)

    # perform local otsu
    local_otsu = rank.otsu(img, selem)

    # select pixels from local otsu
    LO_img = img > local_otsu

    return LO_img, local_otsu


##############################

def median_filter_im(im, currRadius):
    selem = skimage.morphology.square(currRadius)
    im_filt = skimage.filters.median(im, selem)

    return(im_filt)

##############################


def load_images(directory):
    """"
    !!!!!FOR LE ONES I C1 IS PHASE, C2 IS YFP?, and C3 is MKATE so I switched 
        it to c2 for that one if statement 
    
    takes the path as an input and fills arrays w/ full file names 
    for all files of a certain channel and returns them. return(C3_array, C1_array) 
    Assumes C3 is mKate(channel you normalize by) and C1 is YFP (channel being normalized)
    """
    #directory = str(input("Input path to directory with images: "))
    print(directory)

    # fill these arrays and then return them
    C1_array = []
    C3_array = []
    #ShitTimePoints = ["t02", "t06", "t07", "t08"]
    ShitTimePoints = []

    # iterate over files in
    # that directory
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(f):
            if "C3" in f or "c3" in f:
                if any(x in f for x in ShitTimePoints):
                    print("shit be fucky")
                else:
                    C3_array.append(str(f))
            if "C2" in f or "c3" in f:
                if any(x in f for x in ShitTimePoints):
                    print("shit be fucky")
                else:
                    C1_array.append(str(f))
    return(C3_array, C1_array)

###############################


def plt_violin(df):
    (p9.ggplot(data=df,
               mapping=p9.aes(x='timeOrder', y='CellsVolcano', color='TPVolcano'))
     + p9.xlab("Time Point (hr.)")
     + p9.ylab("YFP/mKate")
     + p9.theme_bw()
     + p9.scale_y_continuous(breaks=np.arange(0, .25, .02), limits=[0, .25])
     + p9.theme(axis_text_x=p9.element_text(angle=45, hjust=1))
     + p9.geom_violin()
     + p9.geom_boxplot(width=.1)
     #+ p9.geom_jitter(size=.15)

     + p9.ggtitle("CJG563 0-12hrs, YFP/mKate per Cell; xy=3")
     )

###############################


def medANDnorm(im, currRadius):
    im_filt = median_filter_im(im, currRadius)
    im_normFilt = normalize_im(im_filt)
    return(im_normFilt)


###############################
def cellCheck(im, filt_Thresh):

    
    im_GOtsu_thresh = skimage.filters.threshold_otsu(im)
    if(im_GOtsu_thresh >= filt_Thresh):
        print("there are cells")
        return(True)
    else:
        print("there are no cells in this image")
        return(False)

###############################


def getHists(C3_array, C1_array):
    # Generate histogram of C3 pixel intensities
    C3hist = plt.figure()
    t = 0
    for x in C3_array:
        im_C3 = skimage.io.imread(x)
        hist, bins = skimage.exposure.histogram(im_C3)
        plt.plot(bins, hist, linewidth=1, figure=C3hist)
        t = t + 1

    plt.xlabel('pixel value (a.u.)', figure=C3hist)
    plt.ylabel('count', figure=C3hist)
    plt.title('C3 hists', figure=C3hist)
    # plt.imshow(C3hist)

    # generate histogram of normalized C3 pixel intensities
    normC3hist = plt.figure()
    t = 0
    for x in C3_array:
        im_C3 = skimage.io.imread(x)
        norm_im_C3 = normalize_im(im_C3)
        hist, bins = skimage.exposure.histogram(norm_im_C3)
        plt.plot(bins, hist, linewidth=1, figure=normC3hist)
        t = t + 1

    plt.xlabel('normalized pixel value (a.u.)', figure=normC3hist)
    plt.ylabel('count', figure=normC3hist)
    plt.title('Normalized C3 hists', figure=normC3hist)
    # plt.imshow(normC3hist)

    # Generate histogram of C1 pixel intensities
    C1hist = plt.figure()
    t = 0
    for x in C3_array:
        im_C3 = skimage.io.imread(x)
        hist, bins = skimage.exposure.histogram(im_C3)
        plt.plot(bins, hist, linewidth=1, figure=C1hist)
        t = t + 1

    plt.xlabel('pixel value (a.u.)', figure=C1hist)
    plt.ylabel('count', figure=C1hist)
    plt.title('C1 hists', figure=C1hist)
    # plt.imshow(C1hist)

    return (C3hist, normC3hist, C1hist)

###############################



def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h, h


###############################

def convert(img, target_type_min, target_type_max, target_type):
    imin = img.min()
    imax = img.max()

    a = (target_type_max - target_type_min) / (imax - imin)
    b = target_type_max - a * imax
    new_img = (a * img + b).astype(target_type)
    return new_img

################################


def updateSUMcsv(currDirectory, currMode, aSumDict):
    if(currMode == 'a'):
        print("appending")
        for key in aSumDict:
            tempDF = aSumDict[key]
            tempFilename = "Well_" + str(key) + "_SUM.csv"
            tempFilePath = os.path.join(currDirectory, tempFilename)
            print(tempFilePath)
            tempDF.to_csv(tempFilePath, mode=currMode, header=False)

    if(currMode == 'w'):
        print("writing")
        for key in aSumDict:
            tempDF = aSumDict[key]
            tempFilename = "Well_" + str(key) + "_SUM.csv"
            tempFilePath = os.path.join(currDirectory, tempFilename)
            print(tempFilePath)
            tempDF.to_csv(tempFilePath, mode=currMode, header=True)

##############################


def FIAnalysis(sumInfoDF, cjgANDxys):
    i = 0
    for i in range(len(sumInfoDF)):
        # load info from DF for analysis
        currDirectory = sumInfoDF['paths'][i]
        currBioRep = sumInfoDF['BioRep'][i]
        currKW = sumInfoDF['relInfoKW'][i]

        # Load images from current directory into array
        c3, c1 = load_images(currDirectory)
        C3_array = c3


        # create empty Summary Dictionary
        SumDict = {}
        for k in cjgANDxys.keys():
            SumDict[k] = pd.DataFrame()
            a = SumDict[k]
            a['YFP_ints'] = []
            a['mKate_ints'] = []
            a['Ratios'] = []
            a['cellSizePixel'] = []
            a['timePoint'] = []
            a['BioRep'] = []
            a['XY'] = []

        # Generate histograms of pixel intensity of all images and save to currDir
        C3hist, normC3hist, C1hist = getHists(c3, c1)

        C3hist.savefig(os.path.join(currDirectory, str("C3hist.png")))
        normC3hist.savefig(os.path.join(currDirectory, str("NORM_C3hist.png")))
        C1hist.savefig(os.path.join(currDirectory, str("C1hist.png")))

        ################################################################
        # perform image analysis. Baseline normalized threhold set at .5
        normThresh = .5

        for x in C3_array:

            y = str(x)
            subStr1 = y.split("c3", 1)[0]
            subStr2 = y.split("c3", 1)[1]
            c1_file = subStr1 + "c2" + subStr2

            a = y.split(currKW, 1)[1]
            relInfo = a.split(".tif", 1)[0]
            a = relInfo.split("xy", 1)[1]
            xy = a.split("c", 1)[0]
            b = relInfo.split("t", 1)[1]
            TP = b.split("xy", 1)[0]
            timePoint = (float(TP) - 1)*.25
            searchXY = xy
            foundStrains = [k for k, v in cjgANDxys.items() if searchXY in v]
            corrCJG = "".join(foundStrains)
            currDF = SumDict[corrCJG]

            im_C3 = skimage.io.imread(x)
            normC3filt = medANDnorm(im_C3, 3)

            gauC3 = skimage.filters.gaussian(normC3filt, sigma=3)
            GOthresh = skimage.filters.threshold_otsu(gauC3)
            currThresh = normThresh
            if(GOthresh > currThresh):
                currThresh = GOthresh
            GOgauC3 = gauC3 >= currThresh
            im_cells, num_cells = skimage.measure.label(
                GOgauC3, return_num=True)
            # Make an array where we'll store the cell areas.
            approved_cells = np.zeros_like(im_cells)

            # Loop through each object. Remember that we have to start indexing at 1 in
            # this case!
            for m in range(num_cells):
                cell = (im_cells == m + 1)
                cell_area = np.sum(cell)
                #areas[i] = np.sum(cell)
                if (cell_area > 70.0) & (cell_area < 700):
                    approved_cells += cell

            im_approvedSize, num_approvedSize = skimage.measure.label(
                approved_cells, return_num=True)

            currMask = plt.figure()
            ax = currMask.add_subplot(111)
            ax.imshow(im_approvedSize, cmap=plt.cm.gray)

            TempFileName = "TP" + \
                str(timePoint) + "XY" + str(xy) + "Cells_" + \
                str(num_approvedSize) + "_MASK.png"
            ax.set_title(TempFileName)
            plt.title(TempFileName, figure=currMask)
            currMask.savefig(os.path.join(currDirectory,TempFileName))

            im_C1 = skimage.io.imread(c1_file)
            selem = skimage.morphology.square(3)
            im_c1_filt = skimage.filters.median(im_C1, selem)

            im_C3 = skimage.io.imread(x)
            selem = skimage.morphology.square(3)
            im_c3_filt = skimage.filters.median(im_C3, selem)

            YFP = np.zeros(num_approvedSize)
            C1toC3 = np.zeros(num_approvedSize)
            mKate = np.zeros(num_approvedSize)
            cellSizes = np.zeros(num_approvedSize)
            xyTemp = np.full(shape=num_approvedSize,
                             fill_value=searchXY, dtype=str)
            currTimePoint = np.full(
                shape=num_approvedSize, fill_value=timePoint, dtype=float)
            currBioReps = np.full(shape=num_approvedSize,
                                  fill_value=currBioRep, dtype=str)

            # Loop through each cell.
            for j in range(num_approvedSize):
                # Get the single cell mask.
                cell = (im_approvedSize == j + 1)
                # Store the area.
                cellSizes[j] = np.sum(cell)  # * ip_dist**2
                # Multiply it with the fluorescence image.
                int_im_c3 = cell * im_c3_filt
                curr_c3_int = np.sum(int_im_c3)
                mKate[j] = curr_c3_int

                int_im_c1 = cell * im_c1_filt
                curr_c1_int = np.sum(int_im_c1)
                YFP[j] = curr_c1_int

                ratio = curr_c1_int / curr_c3_int
                C1toC3[j] = ratio

            if(len(C1toC3) > 0):
                mean = statistics.mean(C1toC3)
            else:
                mean = 0

            numCells = len(C1toC3)
            XY = str(xy)

            new_frame = {'YFP_ints': YFP, 'mKate_ints': mKate, 'Ratios': C1toC3, 'cellSizePixel': cellSizes,
                         'timePoint': currTimePoint, 'BioRep': currBioReps, 'XY': XY}
            new_frame_df = pd.DataFrame.from_dict(new_frame)
            SumDict[corrCJG] = SumDict[corrCJG].append(
                new_frame_df, ignore_index=True)

        print("\n\n\n\n\n ####################   UPDATING CSVs    ###################")

        print("current i = " + str(i))

        sumDir = currDirectory.split("594s_LE", 1)[0]

        if(i == 0):
            print("writing")
            updateSUMcsv(sumDir, 'w', SumDict)
        if(i >= 1):
            print("appending")
            updateSUMcsv(sumDir, 'a', SumDict)


####################################################################

def FIAnalysis2(sumInfoDF, cjgANDxys):
    i = 0
    for i in range(len(sumInfoDF)):
        # load info from DF for analysis
        currDirectory = sumInfoDF['paths'][i]
        currBioRep = sumInfoDF['BioRep'][i]
        currKW = sumInfoDF['relInfoKW'][i]

        # Load images from current directory into array
        c3, c2 = load_images(currDirectory)
        C3_array = c3


        # create empty Summary Dictionary
        SumDict = {}
        for k in cjgANDxys.keys():
            SumDict[k] = pd.DataFrame()
            a = SumDict[k]
            a['YFP_ints'] = []
            a['mKate_ints'] = []
            a['Ratios'] = []
            a['cellSizePixel'] = []
            a['timePoint'] = []
            a['BioRep'] = []
            a['XY'] = []

        # Generate histograms of pixel intensity of all images and save to currDir
        C3hist, normC3hist, C1hist = getHists(c3, c2)

        C3hist.savefig(os.path.join(currDirectory, str("C3hist.png")))
        normC3hist.savefig(os.path.join(currDirectory, str("NORM_C3hist.png")))
        C1hist.savefig(os.path.join(currDirectory, str("C2hist.png")))

        ################################################################
        # perform image analysis. Baseline normalized threhold set at .5
        normThresh = .5

        for x in C3_array:

            y = str(x)
            check = y.split(currDirectory + "\\",1)[1]
            
            
            if(check[0] == "T"):
                subStr1 = y.split("C3", 1)[0]
                subStr2 = y.split("C3", 1)[1]
                c1_file = subStr1 + "C2" + subStr2
                a = y.split(currKW, 1)[1]
                relInfo = a.split(".tif", 1)[0]
                a = relInfo.split("XY", 1)[1]
                xy = a.split("_C", 1)[0]
                b = relInfo.split("T", 1)[1]
                TP = b.split("_XY", 1)[0]
                timePoint = (float(TP) - 1)*.16667
                searchXY = xy
                
            if(check[0] == "t"): 
                subStr1 = y.split("c3", 1)[0]
                subStr2 = y.split("c3", 1)[1]
                c1_file = subStr1 + "c2" + subStr2
                # print(c1_file)

                a = y.split(currKW, 1)[1]
                relInfo = a.split(".tif", 1)[0]
                # print(relInfo)
                a = relInfo.split("xy", 1)[1]
                xy = a.split("c", 1)[0]
                b = relInfo.split("t", 1)[1]
                TP = b.split("xy", 1)[0]
                timePoint = round((float(TP) - 1)*.16667, 5)
                searchXY = xy
                
                
                
            foundStrains = [k for k, v in cjgANDxys.items() if searchXY in v]
            corrCJG = "".join(foundStrains)
            # print(corrCJG)
            currDF = SumDict[corrCJG]

            im_C3 = skimage.io.imread(x)
            
            im_mKate = im_C3
            
            # median filter 
            selem = skimage.morphology.square(3)
            medFilt_im = skimage.filters.median(im_mKate, selem)

            # gaussian filter 
            gau_medFilt_im = filters.gaussian(medFilt_im, sigma=1)

            #flat field correction
            ff_im = filters.gaussian(medFilt_im, sigma=70)
            ff_cor = gau_medFilt_im - ff_im


            # Calculate and apply global yen and otsu threshold
            #thresh = threshold_otsu(ff_cor)
            threshYen = skimage.filters.threshold_yen(ff_cor)

            mem_yen = ff_cor > threshYen

            #Local Thresholds
            LT_binary2 = ff_cor > threshold_local(ff_cor, block_size=15, method='gaussian', param=3, offset=0)    

            #Union of global and local thresholds 
            mem_SUM_YENgau = np.logical_and(LT_binary2, mem_yen) 
             
            #Remove Edges
            removeEdges = segmentation.clear_border(mem_SUM_YENgau)

            #binary_opening: erosion followed by dilation to get rid of spots and cracks
            radius = 3
            selem = disk(radius)
            binOpen = skimage.morphology.binary_opening(removeEdges, selem=selem)

            #get rid of objects with a pixel area less than 70
            sizeMask = skimage.morphology.area_opening(binOpen, area_threshold=70)

            
            im_cells_SM, num_cells_SM = skimage.measure.label(sizeMask, return_num=True)
            

            currMask = plt.figure()
            ax = currMask.add_subplot(111)
            ax.imshow(im_cells_SM, cmap=plt.cm.gray)
            

            TempFileName = "TP" + \
                str(timePoint) + "XY" + str(xy) + "Cells_" + \
                str(num_cells_SM) + "_MASK.png"
            ax.set_title(TempFileName)
            plt.title(TempFileName, figure=currMask)
            currMask.savefig(os.path.join(currDirectory,TempFileName))
            # plt.imshow(im_approvedSize)
            # plt.imshow(im_c3_filt)

            im_C1 = skimage.io.imread(c1_file)
            im_c1_filt = im_C1

            im_C3 = skimage.io.imread(x)
            im_c3_filt = im_C3

            YFP = np.zeros(num_cells_SM)
            C1toC3 = np.zeros(num_cells_SM)
            mKate = np.zeros(num_cells_SM)
            cellSizes = np.zeros(num_cells_SM)
            xyTemp = np.full(shape=num_cells_SM,
                             fill_value=searchXY, dtype=str)
            currTimePoint = np.full(
                shape=num_cells_SM, fill_value=timePoint, dtype=float)
            currBioReps = np.full(shape=num_cells_SM,
                                  fill_value=currBioRep, dtype=str)

            # Loop through each cell.
            for j in range(num_cells_SM):
                # Get the single cell mask.
                cell = (im_cells_SM == j + 1)
                # Store the area.
                cellSizes[j] = np.sum(cell)  # * ip_dist**2
                # Multiply it with the fluorescence image.
                int_im_c3 = cell * im_c3_filt
                curr_c3_int = np.sum(int_im_c3)
                mKate[j] = curr_c3_int

                int_im_c1 = cell * im_c1_filt
                curr_c1_int = np.sum(int_im_c1)
                YFP[j] = curr_c1_int

                ratio = curr_c1_int / curr_c3_int
                C1toC3[j] = ratio

            if(len(C1toC3) > 0):
                mean = statistics.mean(C1toC3)
            else:
                mean = 0

            numCells = len(C1toC3)
            XY = str(xy)
#            print("\n \n \n TP: " + str(timePoint) + ", xy: " + XY + ", number of cells: " + str(numCells) +
#                  " , MEAN YFP/mKate: " + str(mean) + ", current Threshold: " + str(currThresh))

            new_frame = {'YFP_ints': YFP, 'mKate_ints': mKate, 'Ratios': C1toC3, 'cellSizePixel': cellSizes,
                         'timePoint': currTimePoint, 'BioRep': currBioReps, 'XY': XY}
            # print(new_frame)
            new_frame_df = pd.DataFrame.from_dict(new_frame)
            print("corrCJG: " + str(corrCJG))
            print("newDF: ")
            print(new_frame_df)
            SumDict[corrCJG] = SumDict[corrCJG].append(
                new_frame_df, ignore_index=True)
            #print("currDF: ")
            #print(SumDict[corrCJG])
            plt.close()

        print("\n\n\n\n\n ####################   UPDATING CSVs    ###################")

        print("current i = " + str(i))

        sumDir = currDirectory.split("594s_LE", 1)[0]

        if(i == 0):
            print("writing")
            updateSUMcsv(sumDir, 'w', SumDict)
        if(i >= 1):
            print("appending")
            updateSUMcsv(sumDir, 'a', SumDict)

####################################################################


def getSUMplots(directories):
    for Dir in directories:
        currDirectory = Dir
        currDF = pd.DataFrame()

        SUM_Dict = {}
        for k in cjgANDxys.keys():
            SUM_Dict[k] = pd.DataFrame()
            a = SUM_Dict[k]
            a['means'] = []
            a['SDs'] = []
            a['SEM'] = []
            a['95CI'] = []
            a['Vars'] = []
            a['TP'] = []
            a['mKate'] = []
            a['YFP'] = []
            a['cellSizePixel'] = []
            a['numCells'] = []
            a['95tInt'] = []

        for key in SUM_Dict:
            tempFilename = "Well_" + str(key) + "_SUM.csv"
            tempFilePath = os.path.join(currDirectory, tempFilename)
            currDF = pd.read_csv(tempFilePath)
            print("\n Number of cells in DF: " + str(len(currDF["timePoint"])))
            print(currDF.head())
            uniqVals = currDF["timePoint"].unique()
          
            for i in uniqVals: 
                print("timePoint pre lookup: " + str(i))
                tempDF = currDF[currDF["timePoint"] == i]
                print("made it to here, timePoint: " + str(i))
                currRatios = tempDF["Ratios"]
                currmKate = np.mean(tempDF["mKate_ints"])
                currYFP = np.mean(tempDF["YFP_ints"])
                currCellSize = np.mean(tempDF["cellSizePixel"])
                currNumCells = len(tempDF["mKate_ints"])
                #currBioRep = tempDF["BioRep"]
                currSEM = scipy.stats.sem(currRatios)
                currMean, CIlow, CIhigh, CI_95 = mean_confidence_interval(
                    currRatios)

                tInt = st.t.interval(alpha=0.95, df=len(currRatios)-1,
                                     loc=np.mean(currRatios), scale=st.sem(currRatios))

                currSD = np.std(currRatios)
                currVar = np.var(currRatios)
                new_frame = {'means': [currMean], 'SDs': [currSD], 'SEM': [currSEM], '95CI': [CI_95],
                             'Vars': [currVar], 'TP': [i], 'mKate': [currmKate], 'YFP': [currYFP],
                             'cellSizePixel': [currCellSize], 'numCells': [currNumCells], '95tInt': [tInt]}
                print("\n new_frame: ", new_frame)
                new_frame_df = pd.DataFrame.from_dict(new_frame)
                print("about to update sumdict")
                SUM_Dict[key] = SUM_Dict[key].append(
                    new_frame_df, ignore_index=True)
                print("updated SUM_Dict")

          

        reduced_Dict = {}
        for key in SUM_Dict:
            reduced_Dict[key] = SUM_Dict[key].dropna()

        return(SUM_Dict, currDF)

##########################################################################

def do_tTests(refDF, testDF, eqVar, alpha): 
    pVals = []
    sigTPs = []
    a = refDF
    b = testDF
    uniqVals = currDf["timePoint"].unique()
    
    for i in uniqVals:
        aTemp = a[a['timePoint'] == i]
        aTemp = aTemp['Ratios']

        bTemp = b[b['timePoint'] == i]
        bTemp = bTemp['Ratios']

        ttests = st.ttest_ind(a=aTemp, b =bTemp, equal_var=eqVar)
        pval = ttests[1]
        pVals.append(ttests[1])
        if(pval < alpha): 
            sigTPs.append(i)
        
    
    return(pVals, sigTPs)

###############################################################################
def sigPlots(WTdf, currDF, sigTPs, label, title): 
    
    point2s = np.full(len(sigTPs), .95)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(sigTPs, point2s, '*', ms=5)
    ax.plot(WTdf["TP"], WTdf["means"], color="r", label = "WT-PA14")
    ax.fill_between(WTdf["TP"], WTdf["means"]-WTdf["95CI"],
                     WTdf["means"]+WTdf["95CI"], color="r", alpha=.2)
    ax.plot(currDF["TP"], currDF["means"], color="b", label = label)
    ax.fill_between(currDF["TP"], currDF["means"]-currDF["95CI"],
                     currDF["means"]+currDF["95CI"], color="b", alpha=.2)
    ax.legend(loc='lower left')
    ax.set_title(title)
    ax.set_ylim(.5,1)
    plt.show()
    return(fig)
    plt.close(fig)
    
    
# %%
# MAIN
###############################################################################
paths = [str(r"Path1"),
         str(r"PATH2")]

BioRep = [1, 2]

relInfoKW = ["T2__", "t2_"]


sumInfoDF = pd.DataFrame({"paths": paths, "BioRep": BioRep, "relInfoKW": relInfoKW},
                         columns=['paths', 'BioRep', 'relInfoKW'])
# create reference dictionary
cjgANDxys = {'strain1': ["01", "02", "03", "04", "05", "06"],
             'strain2': ["07", "08", "09", "10", "11", "12"]}

FIAnalysis2(sumInfoDF, cjgANDxys)

print("\n\n\n\n\n\n\n\n\n #############    DONE    ###########")

# %%
path = [str(r"C:\Users\chris\OneDrive\Desktop\594s_LE")]

reduced_Dict, currDf = getSUMplots(path)


# %%
# 
###############################################################################
a = reduced_Dict['594_WT']['95tInt'][0][0]
b = reduced_Dict['594_WT']['95tInt'][0][1]
c = b - reduced_Dict['594_WT']['means'][0]
print("c: " + str(c))
print("95CI: " + str(reduced_Dict['594_WT']['95CI'][0]))

# %%


cjgANDxys = {'594_WT': ["01", "02", "03", "04", "05", "06"],
             '659_d-pilU': ["07", "08", "09", "10", "11", "12"],
             '651_d-pilT': ["07", "08", "09", "10", "11", "12"],
             # '656_d-cyaAB': ["07", "08", "09", "10", "11", "12"],
             '666_K136A': ["07", "08", "09", "10", "11", "12"],
             '627_E204A': ["07", "08", "09", "10", "11", "12"],
             # '663_d-cpdA': ["07", "08", "09", "10", "11", "12"],
             '630_H222A': ["07", "08", "09", "10", "11", "12"],
             '628_H229A': ["07", "08", "09", "10", "11", "12"],
              # '653_d-pilTU': ["07", "08", "09", "10", "11", "12"],
              '632_H44L': ["07", "08", "09", "10", "11", "12"],
             '625_K58A': ["07", "08", "09", "10", "11", "12"],
             '629_D31K': ["07", "08", "09", "10", "11", "12"]
             }
# %%
########################################################################
#Taking out bad bioreps 

#WT
csvWT_SUM = pd.read_csv(str(r"pathtoSUMMARYcsv.csv"))
WTdf = reduced_Dict['594_WT']

ax1 = csvWT_SUM.plot.scatter(x='timePoint', y='Ratios', c='BioRep', colormap = 'magma')
ax1.set_title('WT_df Raw')

TOTAL_csvWT = csvWT_SUM
csvWT_SUM = TOTAL_csvWT[TOTAL_csvWT["BioRep"] < 5]

#redo graph 
ax1 = csvWT_SUM.plot.scatter(x='timePoint', y='Ratios', c='BioRep', colormap = 'magma')
ax1.set_title('WT_df Corrected')

#save new csv 
csvWT_SUM.to_csv(str(r"pathtoSUMMARY.csv"), mode='w', header=True)

####################


# %%
########################################################################
#re-read in df's and do summary stats

path = [str(r"Path")]

reduced_Dict, currDf = getSUMplots(path)




# %%
########################################################################



csvWT_SUM = pd.read_csv(str(r"PathToCsv.csv"))
csv58_SUM = pd.read_csv(str(r"pathtoSUMMARY.csv"))

WTdf = reduced_Dict['594_WT']
currDF = reduced_Dict['625_K58A']


TOTAL_csvWT = csvWT_SUM
csvWT_SUM = TOTAL_csvWT[TOTAL_csvWT["BioRep"] < 5]



ax2 = csv58_SUM.plot.scatter(x='timePoint', y='Ratios', c='BioRep', colormap = 'magma')


pVals_K58A, sigTPs_K58A = do_tTests(csvWT_SUM, csv58_SUM, True, .05)
pVals_K58A
sigTPs_K58A

currLabel = "pilT-K58A"
currTitle = "WT vs. pilT-K58A"
#when ylim lower is .2
K58A_plot = sigPlots(WTdf, currDF, sigTPs_K58A, currLabel, currTitle)



##########################################################################

# %%

##########################################################################


plt.figure(figsize=(7, 7))
i = 0
for key in reduced_Dict:
    currDF = reduced_Dict[key]
    i = int(i)
    currLabel = str("CJG" + key)
    plt.plot(currDF["TP"], currDF["means"], color=colors[i], label=currLabel)
    plt.fill_between(currDF["TP"], currDF["means"]-currDF["SDs"],
                     currDF["means"]+currDF["SDs"], color=colors[i], alpha=.2)
    i = i+1
plt.legend(loc='lower left')
plt.title("pilT mutants: YFP/mKate 0-7hrs, n=3")
plt.show()





plt.figure(figsize=(15,7))
i = 0
for key in reduced_Dict:
    currDF = reduced_Dict[key]
    currDF = currDF[currDF["TP"] < 6.1]
    
    i = int(i)
    currLabel = str(key)
    plt.plot(currDF["TP"], currDF["means"], color=colors[i], label=currLabel)
    plt.fill_between(currDF["TP"], currDF["means"]-currDF["95CI"],
                     currDF["means"]+currDF["95CI"], color=colors[i], alpha=.2)
    i = i+1
# plt.legend(loc='lower left')
# plt.xlabel("time (hours)")
# plt.ylabel("average YFP/mKate per cell")
plt.show()


###############################################
#for a select few strains

# Number of Cells

plt.figure(figsize=(15, 7))
i = 0
for key in reduced_Dict:
    if key == '594_WT' or key == '630_H222A' or key == '628_H229A' or key == '659_d-pilU' or key == '666_K136A' or key == '627_E204A':
        currDF = reduced_Dict[key]
        i = int(i)
        currLabel = str("CJG" + key)
        plt.plot(currDF["TP"], currDF["numCells"],
                 color=colors[i], label=currLabel)
        # plt.fill_between(currDF["TP"], currDF["means"]-currDF["95CI"],
        # currDF["means"]+currDF["95CI"], color=colors[i], alpha=.2)
    i = i+1
plt.legend(loc='upper left')
plt.title("pilT mutants: Number of Cells 0-8hrs, n=3")
plt.show()

###
# MEANS YFP/mKate

plt.figure(figsize=(15, 7))
i = 0
for key in reduced_Dict:
    if key == '594_WT' or key == '630_H222A' or key == '628_H229A' or key == '659_d-pilU' or key == '666_K136A' or key == '627_E204A':
        currDF = reduced_Dict[key]
        i = int(i)
        currLabel = str("CJG" + key)
        plt.plot(currDF["TP"], currDF["means"],
                 color=colors[i], label=currLabel)
        # plt.fill_between(currDF["TP"], currDF["means"]-currDF["95CI"],
        # currDF["means"]+currDF["95CI"], color=colors[i], alpha=.2)
    i = i+1
plt.legend(loc='lower left')
plt.title("pilT mutants: means YFP/mKate: 0-8hrs, n=3")
plt.show()


###
# YFP / cell

plt.figure(figsize=(15, 7))
i = 0
for key in reduced_Dict:
    if key == '594_WT' or key == '630_H222A' or key == '628_H229A' or key == '659_d-pilU' or key == '666_K136A' or key == '627_E204A':
        currDF = reduced_Dict[key]
        i = int(i)
        currLabel = str("CJG" + key)
        plt.plot(currDF["TP"], currDF["YFP"],
                 color=colors[i], label=currLabel)
        # plt.fill_between(currDF["TP"], currDF["means"]-currDF["95CI"],
        # currDF["means"]+currDF["95CI"], color=colors[i], alpha=.2)
    i = i+1
plt.legend(loc='lower left')
plt.title("pilT mutants: YFP/cell 0-8hrs, n=3")
plt.show()


###
# mKate / cell

plt.figure(figsize=(15, 7))
i = 0
for key in reduced_Dict:
    if key == '594_WT' or key == '630_H222A' or key == '628_H229A' or key == '659_d-pilU' or key == '666_K136A' or key == '627_E204A':
        currDF = reduced_Dict[key]
        i = int(i)
        currLabel = str("CJG" + key)
        plt.plot(currDF["TP"], currDF["mKate"],
                 color=colors[i], label=currLabel)
        # plt.fill_between(currDF["TP"], currDF["means"]-currDF["95CI"],
        # currDF["means"]+currDF["95CI"], color=colors[i], alpha=.2)
    i = i+1
plt.legend(loc='lower left')
plt.title("pilT mutants: mKate/cell 0-8hrs, n=3")
plt.show()

###
# cell size / cell

plt.figure(figsize=(15, 7))
i = 0
for key in reduced_Dict:
    if key == '594_WT' or key == '630_H222A' or key == '628_H229A' or key == '659_d-pilU' or key == '666_K136A' or key == '627_E204A':
        currDF = reduced_Dict[key]
        i = int(i)
        currLabel = str("CJG" + key)
        plt.plot(currDF["TP"], currDF["cellSizePixel"],
                 color=colors[i], label=currLabel)
        # plt.fill_between(currDF["TP"], currDF["means"]-currDF["95CI"],
        # currDF["means"]+currDF["95CI"], color=colors[i], alpha=.2)
    i = i+1
plt.legend(loc='lower left')
plt.title("pilT mutants: cell size (pixels) 0-8hrs, n=3")
plt.show()





plt.figure(figsize=(7, 7))
i = 0
for key in reduced_Dict:
    currDF = reduced_Dict[key]
    i = int(i)
    currLabel = str("CJG" + key)
    plt.plot(currDF["TP"], currDF["cellSizePixel"],
             color=colors[i], label=currLabel)
    # plt.fill_between(currDF["TP"], currDF["means"]-currDF["95CI"],
    # currDF["means"]+currDF["95CI"], color=colors[i], alpha=.2)
    i = i+1
plt.legend(loc='lower left')
plt.title("pilT mutants: average cell size 0-7hrs, n=3")
plt.show()





