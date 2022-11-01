function [UsePIV,s,p,r,UseCellDivDetect,CType,NetName,NetNameSeg,MinDivArea, ...
    BlockSize,MakeVideo,MedFiltSize,FiltType,BackThresh,WSize,MaxTimeDiff,...
    UseOrientationAnalysis,MSize,rho,UseBM3D,UseCLAHE,UseDriftCorrect,...
    WinSizeDrift] = ParameterFunctionMainPIVlab()
% Set all relevant parameters for the analysis of velocity and orientation
% fields of dense cell layers
%% Set parameters for PIV:
% Use this module? 1 == yes, 0 == no.
UsePIV = 1;
% Standard PIV Settings
s = cell(15,2); 
%Parameter                          %Setting           %Options
s{1,1}= 'Int. area 1';              s{1,2}=128;        % window size of first pass
s{2,1}= 'Step size 1';              s{2,2}=64;         % step of first pass
s{3,1}= 'Subpix. finder';           s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
s{4,1}= 'Mask';                     s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
s{5,1}= 'ROI';                      s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
s{6,1}= 'Nr. of passes';            s{6,2}=4;          % 1-4 nr. of passes
s{7,1}= 'Int. area 2';              s{7,2}=64;         % second pass window size
s{8,1}= 'Int. area 3';              s{8,2}=32;         % third pass window size
s{9,1}= 'Int. area 4';              s{9,2}=26;         % fourth pass window size
s{10,1}='Window deformation';       s{10,2}='*linear'; % '*spline' is more accurate, but slower
s{11,1}='Repeated Correlation';     s{11,2}=0;         % 0 or 1 : Repeat the correlation four times and multiply the correlation matrices.
s{12,1}='Disable Autocorrelation';  s{12,2}=0;         % 0 or 1 : Disable Autocorrelation in the first pass.
s{13,1}='Correlation style';        s{13,2}=0;         % 0 or 1 : Use circular correlation (0) or linear correlation (1).
s{14,1}='Repeat last pass';         s{14,2}=0;         % 0 or 1 : Repeat the last pass of a multipass analyis
s{15,1}='Last pass quality slope';  s{15,2}=0.025;     % Repetitions of last pass will stop when the average difference to the previous pass is less than this number.

% Standard image preprocessing settings
p = cell(10,1);
%Parameter                       %Setting           %Options
p{1,1}= 'ROI';                   p{1,2}=s{5,2};     % same as in PIV settings
p{2,1}= 'CLAHE';                 p{2,2}=0;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
p{3,1}= 'CLAHE size';            p{3,2}=50;         % CLAHE window size
p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
p{7,1}= 'Wiener';                p{7,2}=0;          % 1 = enable Wiener2 adaptive denoise filter, 0 = disable
p{8,1}= 'Wiener size';           p{8,2}=3;          % Wiener2 window size
p{9,1}= 'Minimum intensity';     p{9,2}=0.0;        % Minimum intensity of input image (0 = no change)
p{10,1}='Maximum intensity';     p{10,2}=1.0;       % Maximum intensity on input image (1 = no change)

% Standard image post processing settings
r = cell(6,1);
%Parameter                                                  %Setting                    % Options
r{1,1}= 'Calibration factor, 1 for uncalibrated data';      r{1,2}=1;                   % Calibration factor for u
r{2,1}= 'Calibration factor, 1 for uncalibrated data';      r{2,2}=1;                   % Calibration factor for v
r{3,1}= 'Valid velocities [u_min; u_max; v_min; v_max]';    r{3,2}=[-50; 50; -50; 50];  % Maximum allowed velocities, for uncalibrated data: maximum displacement in pixels
r{4,1}= 'Stdev check?';                                     r{4,2}=1;                   % 1 = enable global standard deviation test
r{5,1}= 'Stdev threshold';                                  r{5,2}=10;                  % Threshold for the stdev test
r{6,1}= 'Local median check?';                              r{6,2}=1;                   % 1 = enable local median test
r{7,1}= 'Local median threshold';                           r{7,2}=3.82;                % Threshold for the local median test - MAD of 3.82 corresponds to 99% confidence intervall for gaussian distribution


%% Set parameters for cell division analysis:
% Use this module? 1 == yes, 0 == no.
UseCellDivDetect = 1;
% Cell type to be analysed:
% The approach uses pattern matching to identify candidates of potential
% cell division events. The respective templates are stored along the
% m-files for analysis in the subfolder '\DivisionTemplate\Templates\CType'
CType = 'GBM #4';
% CType = 'U138';
% Name of the file containing the neuronal network used for classification 
% of cell division events. Nets are stored along the m-files for analysis 
% in the subfolder '\TrainedANNs'. 
NetName = 'GoogLeNet';
% Name of the network used for semantic segmentation of division events:
NetNameSeg = 'SegmentationUNet';
% Minimal size of a detected cell division figure (used in post processing
% after segmentation). Given in pixel:
MinDivArea = 500;
% Block size used for training the segmentation network:
BlockSize = [256,256];
% Make a video of detected cell divisions and associated tracks? 1 = yes, 0
% = no
MakeVideo = 1;

% Contrast limited adaptive histogram equalization (CLAHE). Recomended for
% division detection if contrast of obtained images is significantly
% different than contrast of templates used for cross-correlation:
% Use this module? 1 == yes, 0 == no.
UseCLAHE = 1; 

% Minor control parameters:
% Size of 2D median filter used to denoise individual cross correlation 
% maps:
MedFiltSize = [9 9];
% Parameters for peak finding in the final cross correlation map:
% Type of filter used for smooting during peak finding:
FiltType = fspecial('gaussian', 85,5);
% Threshold value in the cross correlation map considered background 
% (values: 0 ... 1):
BackThresh = 0.03;
% Approximately half the diameter of an S-phase cell. Corresponds to the 
% bounding box formed around a cell division candidate (center +- window 
% size) to feed it to the ANN for classification. Has to correspond to the 
% size used for training of the ANN.
WSize = 35;
% Maximal time (in sucessive images) a division can be invisible and still
% be assigned a previous/old detection. Necessary, as sucessive divisions
% can occur at very similar places and might be mismatched.
MaxTimeDiff = 25;
%% Set parameters for cell orientation analysis:
% Use this module? 1 == yes, 0 == no.
UseOrientationAnalysis = 1;
% Meshgrid StepSize. Corresponds to the distance in x and y direction (in
% px) were the next cellular orientation is calculated. Lowering this
% number increases quadratically the amount of RAM needed.
MSize = 10;
% standard deviation (in pixels) of the Gaussian kernel applied
% to the strcture tensor to average the directions of the eigenvectors.
rho=15;

%% Other parameters:
% Use BM3D filter to denoise input image? 1 == yes, 0 == no.
UseBM3D = 0;

% Use drift correction for velocity field to adjust stage drift and
% imprecise position recovery when imaging multiple positions:
UseDriftCorrect = 1;
% Parameters for drift correction:
% Windo size used fo pattern matching. Should correspond roughly to image
% size.
WinSizeDrift = 960; 
