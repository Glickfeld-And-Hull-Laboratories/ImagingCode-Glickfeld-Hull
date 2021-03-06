% CALCIUM
%
% ImageJ based TIFF file input / output
%   ijarray2plus       - Converts MATLAB array to ImageJ ImagePlus object
%   ijGetFreeMemory    - prints ImageJ Java Heap usage to standard input
%   ijplus2array       - Converts ImageJ ImagePlus object to MATLAB array
%   ijreadtiff         - Reads TIFF file and returns ImageJ ImagePlus object
%   ijreadtiffinfo     - Reads TIFF file information and returns ImageJ TiffInfo object 
%   ijwritetiff        - Writes ImageJ ImagePlus object to disk
%   readtiff           - Reads TIFF file or sequences of TIFF files into MATLAB array
%   writetiff          - Writes array as TIFF file
%
% Stack analysis (3d images)
%   stackGroupProject     - Downsample stack array by ratio
%   stackSubtract                - Subtract frame from stack array
%   stackCycleAverage     - Trial average
%   stackLocalContrastAdaptation - 
%   stackFixBidirPhase  - (calcium): automatically detect, remove bidir scan offsets
%   stackAnimator                - if 4d, the fourth dimension is the RGB dim
%   stackWriteAvi                - $Id: Contents.m 559 2009-07-22 22:44:47Z kenichi $
%   stackGetTimeCourses          - Calculate pixel intensities inside multiple ROIs
%   stackRegisterJava            - stackRegisterJava (calcium): register a stack (single-threaded)%
%   stack_localcontrastadj       - $Id: Contents.m 559 2009-07-22 22:44:47Z kenichi $
%
% Utility functions
%   imBidirShift        - imBidirShift (calcium): offset every other line; changes bidir scan alignment
%   interpLowpass       - INTERP_LOWPASS (stacks-mh): filter a timecourse using only a subset of frames
%   roiDrawOnImage      - roiDrawOnImage (calcium): compose image with RGB roi mask
%
% Segmentation of cells
%   scriptFindCells       - Runs IMFINDCELLS to segment neurons and glia
%   imFindCells           - Find cells from image
%   imFindCellsParamsSet  - Creates / alters parameters structure for IMFINDCELLS
%   imFindCellsTM         - find cells based on template matching
%
% Time course analysis
%   tcEpochAverage        - calculate average signal value during each epoch
%   tcTrialAverage        - av = tcTrialAverage(timeCourses,trialdur)
%   tcContrast            - 
%   tcRemoveDC            - 
%   tcOffsetPlot          - Plot time courses with offset
%   tcFilter              - (calcium): on tc, do high-pass, dF/F, trial averaging, final smooth
%   tcShuffleCorrect      - Subtracts trial average from time courses
%   tcPlotTcourses        - expt: structure
%   tcLowCut              - Removes low frequencies while keeping DC
%   tcRemoveArtifacts     - Removes artifacts from time courses
%   tcProcess             - Process raw timecourses to get lowcut,normalized, or averaged timecourses and a data table
%   tcStats               - minimal statistics for any stim
%
% Image Processing
%   imShade               - IMSHADECELLS Shade grayscale image with mask
%   imScale               - Rescale image pixel values as desired
%   imFilter2             - 2D filter, modified from filter2, considering boundaries
%   imAdaptHistEqJava            - Wrapper around Filter_RankAK (modification of Filter_Rank IJ plugin)
%   imAdaptHistEqJavaRevBorders  - 
%   adapthisteqRevBorders        - 
%   imCellEditInteractive        - (calcium): add cells to a cell mask interactively
%   imLog                        - Log transform image intensities
%   imColorize            - colorize gray image with a color
%   imHueRotate           - Rotate hue of an image in KO color space (0-1)
%
%
% Interactive 
%   stackZprofiler      - stackZprofiler (calcium): interactively plot third dim (usu. time) of stack
%
%
% Experiment control
%
%   getepochs          - [stims,blanks] = getepochs (off, on, nstim, ntrials);
%   readstimlogfile    - 
%   setupstim          -
%
% Registration
%
%   MOVED TO CORE/REG
%
% Signal Processing
%   filtfilthd         - Zero-phase digital filtering with dfilt objects.
%   dcblock            - determines the filter coefficient a for the dc blocking/high- 
%
% Statistics
%   anova1mult         - Anova with correction for multiple comparisons
%
% Color 
%   hueHSV2KO          - Convert hue values (0-1) in HSV color space to hue vaules in KO color space
%   hueKO2HSV          - Convert hue values (0-1) in KO color space to hue vaules in HSV color space
%   hsv2rgbKO          - In hsv2rgb, red is opposed to cyan. In hsv2rgbKO, red is opposed to green, and blue is opposed to yellow
%
% Misc
%   GetTifFileNames              - get names of tif files in a specified directory
%   GetMatFileNames              - get names of mat files in a specified directory
%   figCellsLabeled              - FIG_CELLS_LABELEDMH (calcium): make an image plot with cell ROIs labeled
%   ezCellMap                    - make a cell map of any "index"
%   readXLSfilelist              - read filelist from an excel file
%   writeXLSfilelist             - write filelist to an excel file
%   getCellCoordinate             - get coordinates of cells from a cell mask
