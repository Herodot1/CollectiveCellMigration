function [UsePIV,met,wins,Overlap,globtrld,loctrld,kernelsize,med,int,...
    UseCellDivDetect,CType,NetName,MakeVideo,MedFiltSize,FiltType,...
    BackThresh,WSize,MaxTimeDiff,UseOrientationAnalysis,MSize,rho] = ParameterFunctionMain()
% Set all relevant parameters for the analysis of velocity and orientation
% fields of dense cell layers
%% Set parameters for PIV:
% Use this module? 1 == yes, 0 == no.
UsePIV = 0;
% Analysis methods. Valid inputs:
% 'single', for single path PIV
% 'multi', for three passes of the PIV using the initial window size and
% ending up with half of it
% 'multin', PIV with n iterations using the n*2 vector "windowsize" for
% each step of iteration. All of the above methods use the "standard cross
% correlation definition"
% 'mqd', single pass PIV using a cross correlation term which does not
% neglect correlation terms between pixels in image 2
% 'norm', refers to a normalized cross correlation
% details to all methods can be found ion the manual of the MatPIV library.
met='multin';
% wins can be a N*2 vector containing the size of the interrogation
% regions (x and y). 
wins=[128 128; 64 64 ;32 32; 32 32]; 
% Overlap of windows:
% Important: Overlap multiplied by the windowsize needs to be an integer!!! 
Overlap = 0.6; 
% threshold for use with globalfiltering of velocity fields
% gives the maximal velocity that is not regarded an outlier in px per
% time between sucessive images. It calculates as follows: mean +-
% globtrld*standard deviation. Outside of these bounds values are
% considered outliers.
globtrld=15;
% threshold for use with local filtering of velocity fields
% a vector is discarded if it is larger/smaller than the mean/median +-
% threshold*standard deviation of its surrounding velocity vectors.
loctrld=2.5; 
% kernelsize of the local filter used:
% defines the number of vectors contributing to the median or mean value of
% each vector as size kernelsize*kernelsize
kernelsize = 3;
% Filtering type used for local filter
% 'median' or 'mean' -> respective type of filtering 
med='median'; 
% Interpolation method used for interpolating outliers
% 'linear' or 'weighted'
int='linear'; 

%% Set parameters for cell division analysis:
% Use this module? 1 == yes, 0 == no.
UseCellDivDetect = 1;
% Cell type to be analysed:
% The approach uses pattern matching to identify candidates of potential
% cell division events. The respective templates are stored along the
% m-files for analysis in the subfolder '\DivisionTemplate\Templates\CType'
CType = 'U138';
% Name of the file containing the neuronal network used for classification 
% of cell division events. Nets are stored along the m-files for analysis 
% in the subfolder '\TrainedANNs'. 
NetName = 'GoogLeNet';
% Make a video of detected cell divisions and associated tracks? 1 = yes, 0
% = no
MakeVideo = 1;

% Minor control parameters:
% Size of median filter used to denoise individual cross correlation maps:
MedFiltSize = [9 9];
% Parameters for peak finding in the final cross correlation map:
% Type of filter used for smooting during peak finding:
FiltType = fspecial('gaussian', 85,5);
% Threshold value in the cross correlation map considered background 
% (values: 0 ... 1):
BackThresh = 0.03;
% Half the diameter of an S-phase cell. Corresponds to the bounding box
% formed around a cell division candidate (center +- window size) to feed
% it to the ANN for classification. Has to correspond to the size used for
% training of the ANN.
WSize = 35;
% Maximal time (in sucessive images) a division can be invisible and still
% be assigned a previous/old detection. Necessary, as sucessive divisions
% can occur at very similar places and might be mismatched.
MaxTimeDiff = 25;
%% Set parameters for cell orientation analysis:
% Use this module? 1 == yes, 0 == no.
UseOrientationAnalysis = 0;
% Meshgrid StepSize. Corresponds to the distance in x and y direction (in
% px) were the next cellular orientation is calculated. Lowering this
% number increases quadratically the amount of RAM needed.
MSize = 10;
% standard deviation (in pixels) of the Gaussian kernel applied
% to the strcture tensor to average the directions of the eigenvectors.
rho=15;


