% Get order parameter and 4-point susceptibility only. useful for testing
% which cell size (CSize) and overlap threshold (OvThresh) is best.

clear all;
beta = pwd;
beta = addpath(beta);  % Adds the path where the m-files are stored in
StartPath=pwd;
DataPath=uigetdir(StartPath, 'Chose the folder with the images');
alpha = cd(DataPath);   % legt jpg_path2 als aktuellen folder fest
FolderList = dir();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Parameters
% Define points of interest:
s_size = 80;
POIx = [41:s_size:960];
POIy = [41:s_size:1280];
[POIx, POIy] = meshgrid(POIx,POIy);
% POIx = x;
% POIy = y;
% s_size = 40;
% POIx = [8:s_size:960];
% POIy = [8:s_size:1280];
% image size:
im_size = [960,1280];
% Julians measurements:
% im_size = [1040,1392];
% Time between successive images in minutes:
dt = 3; % GBM cells
%dt = 5; % colon carcinoma
% Julians measurements:
% dt = 5;
% Cell Size (in px):
CSize = 80;
% Overlap thresh:
OvThresh = 0.60;
% Window Size:
WSize = 41;
% 6h Intervall for GBM Measurements
WSize = 121;
% % Overlap thresh:
% OvThresh = 0.2;
% % Window Size:
% WSize = 31;

% % Alternativ:
% OvThresh = 0.5;
% WSize = 51;

% Sub-Pixel resolution for the tracking matrix:
SubPxResolution = 0.2; % Colon Carcinoma
SubPxResolution = 1; % GBM Cells
% Take only every x-th image for MSDTemp etc calculation and inerpolate the
% rest:
Sampling = 9; % GBM Cells (Post Calculation)
Sampling = 5; % Colon Carcinoma
Sampling = 61; % Used for time interval Measurements;

% Savename:
SaveName = '6h Interval OvThresh 0.6'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in files and do the velocity field analysis:
% Remove all non directories:
Idx = cat(1,FolderList.isdir);
FolderList(Idx~=1) = [];
% Save start folder:
StartDir = pwd;


for folder_number1 = 3:length(FolderList)
    folder_number1
    % Go to folder
    cd(FolderList(folder_number1).name)    
    
    % Generate object list of current folder:
    FolderListInner = dir();
    % Remove all non directories:
    IdxInner = cat(1,FolderListInner.isdir);
    FolderListInner(IdxInner~=1) = [];
    % Save start folder:
    CurDir = pwd;    
    
    % Allocate variables:
    QAll        = NaN(700,length(FolderListInner)-2);
    QTempAll    = NaN(WSize,500,length(FolderListInner)-2);
    ChiAll      = NaN(700,length(FolderListInner)-2);
    ChiTempAll  = NaN(WSize,500,length(FolderListInner)-2);    
    for folder_number2 = 4:length(FolderListInner)        
        % Go to results folder:
        cd(strcat(FolderListInner(folder_number2).name,'\Results'));
        
        % Load velocity field:
        load('VelocityField.mat')
        m = size(VelField,4);
        % Get tracking matrix from velocity field:
        [TrackMat] = GenTrackMat(VelField,im_size,POIx,POIy,SubPxResolution);
        % Get total MSD, Order Parameter and 4 point susceptibility for one velocity field:
        [MSD, Q, Chi, tau] = MSDChiOrder(TrackMat,dt,CSize,OvThresh);
        % Get time resolved MSD, Order Parameter and 4 point susceptibility for one velocity field:
        [MSDTemp, QTemp, ChiTemp] = MSDChiOrderTimeWindow(VelField,WSize,dt,CSize,OvThresh,im_size,POIx,POIy,SubPxResolution,Sampling);
        
        % Plug parameters into one global construct:
        % Global parameters:
        QAll(1:size(Q,2),folder_number2-2)      	= Q;
        ChiAll(1:size(Chi,2),folder_number2-2)	= Chi;        
        % Time dependent variables:
        ObjNum = size(MSDTemp,2);
        QTempAll(:,1:ObjNum,folder_number2-2)	= QTemp;
        ChiTempAll(:,1:ObjNum,folder_number2-2)  = ChiTemp;
        
        % Go back to top level directory:
        cd(CurDir)
    end
    
    % Save Variables:
    save(sprintf('ChiQ %s.mat', SaveName), 'QAll','QTempAll','ChiAll',...
        'ChiTempAll','tau','WSize', 'dt','im_size',...
        'CSize','OvThresh')
    VelFieldSize = size(VelField);
    VelFieldSize = VelFieldSize(1:2);
    save(sprintf('Settings %s.mat', SaveName),'s_size',...
        'POIx','POIy','im_size','dt','CSize','OvThresh','WSize','SubPxResolution','Sampling','VelFieldSize')   
   % Go back to top level directory:
   cd(StartDir)     
end









