clear all;

mPath = mfilename('fullpath');
Idx = max(strfind(mPath,'\'));
mPath = mPath(1:Idx);
addpath(mPath)
% add path for velocity field analysis:
addpath(strcat(mPath,'\Velocity Fields'))
% add path of the orientation analysis:
addpath(strcat(mPath,'\Orientation Analysis'))
% add path of the circular statistics toolbox:
addpath(strcat(mPath,'\CircStat'))

beta = pwd;
beta = addpath(beta);  % Adds the path where the m-files are stored in
StartPath=mPath;
DataPath=uigetdir(StartPath, 'Chose the folder with the images');
alpha = cd(DataPath);   % legt jpg_path2 als aktuellen folder fest
FolderListOuter = dir();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Parameters
% Define points of interest:
s_size = 80;
POIx = [41:s_size:960];
POIy = [41:s_size:1280];
[POIx, POIy] = meshgrid(POIx,POIy);
% image size:
im_size = [960,1280];
% Julians measurements:
% im_size = [1040,1392];
% Physical dimension of image in µm:
ImPhysSize = [465.80,621.23];
% Julian
% ImPhysSize = [478.18,640.18];
% pixel size of image in µm/px:
PxSize = ImPhysSize/im_size;
% Time between successive images in minutes:
dt = 3; % GBM cells /Ula
%dt = 10; % Scratch Assay Franzi
% dt = 5; % colon carcinoma
% dt = 3;
% Julians measurements:
% dt = 5;
% Cell Size (in px):
CSize = 80;
% Overlap thresh:
OvThresh = 0.1; % Colon Carcinoma
OvThresh = 0.20; % GBM Cells
% OvThresh =0.3;  % A375/ Ulas Messungen/ Alarmine /Sophie
% Window Size:
WSize = 41; % 2h Intervall bei dt = 3min
WSize = 121; % 6h Intervall bei dt = 3min
WSize = 401; % 20h intervall bei dt =3min
WSize = 121; % 6h intervall bei dt =3min
WSize = 73; % 6h Intervall bei dt = 5min
WSize = 241; % 12h intervall bei dt =3min
WSize = 401; % 20h intervall bei dt =3min
WSize = 121; % 20h intervall bei dt =10min
WSize = 401; % 20h intervall bei dt =3min
% Center velocity field if value is 1, do not center otherwise:
% Centering done as follows: Vx = Vx-mean(Vx(:)), etc
CenterSpeed = 0; % Blob Measurements LN229/U138 of Sophie
CenterSpeed = 1; % Astros Scratch Julian, collective migration in dense monolayers;
% Sub-Pixel resolution for the tracking matrix:
SubPxResolution = 0.2; % Colon Carcinoma
SubPxResolution = 1; % GBM Cells/A375/ Ulas Messungen/ Alarmine/Astros Julian
% Take only every x-th image for MSDTemp etc calculation and inerpolate the
% rest:
Sampling = 400; % 20h bei dt=3min
Sampling = 9;  % GBM Cells (Post Calculation)
Sampling = 12; % Colon Carcinoma
Sampling = 40; % GBM/A375 von Ula/ Alarmine/ Ring Test Scratch
Sampling = 240; % 12h bei dt=3min; Astrozytenüberstand, 1% FBS+ Alarmine, 0.1% FBS + Alarmine
Sampling = 40; % 2h bei dt=3min; Astrozytenüberstand 1d, 1% FBS+ Álarmine 1d
Sampling = 120; % 10h bei dt=5min; Astros Scratch von Julian
Sampling = 200; % 10h bei dt=3min; Dense Monolayer 3d measurements,
% Sampling = 60; % 10h bei dt=10min; Scratch Astros Franzi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in files and do the velocity field analysis:

% Remove all non directories:
Idx = cat(1,FolderListOuter.isdir);
FolderListOuter(Idx~=1) = [];
% remove ".." and "." from the folder list:
for i = length(FolderListOuter):-1:1
    if strcmp(FolderListOuter(i).name,'..') || strcmp(FolderListOuter(i).name,'.')
        FolderListOuter(i) = [];
    end
end

% Remove a certain "non-interesting" directory:
for FolderNum = length(FolderListOuter):-1:3
    if strcmp(FolderListOuter(FolderNum).name,'Results Path Analysis') == 1
        FolderListOuter(FolderNum) = [];
    end
end
curr_directory_outer = pwd;

for FolderNum = 1:length(FolderListOuter)
    cd(FolderListOuter(FolderNum).name);
    
    FolderList = dir();
    % Remove all non directories:
    Idx = cat(1,FolderList.isdir);
    FolderList(Idx~=1) = [];
    % Remove a certain "non-interesting" directory:
    for FolderNum2 = length(FolderList):-1:3
        if strcmp(FolderList(FolderNum2).name,'Results Path Analysis') == 1
            FolderList(FolderNum2) = [];
        end
    end
    
    % Save start folder and set save name:
    curr_directory = pwd;
    SaveNameVelField = curr_directory(max(strfind(pwd,'\'))+1:end);
    SaveNameVelField = sprintf('%s.mat',SaveNameVelField);
    SaveNameTrackMat = curr_directory(max(strfind(pwd,'\'))+1:end);
    SaveNameTrackMat = sprintf('%s %s.mat','TrackingMatrixTemp',SaveNameTrackMat);
    SaveNameOrientation = curr_directory(max(strfind(pwd,'\'))+1:end);
    SaveNameOrientation = sprintf('Orientation Analysis %s.mat',SaveNameOrientation);
    SaveNameDivision = curr_directory(max(strfind(pwd,'\'))+1:end);
    SaveNameDivision = sprintf('Division Analysis %s.mat',SaveNameDivision);
    
    % Allocate variables for velocity field analysis:
    pAll        = NaN(1700,length(FolderList)-2);
    VTempAll    = NaN(WSize,1500,length(FolderList)-2);
    AreaAll     = NaN(1700,length(FolderList)-2);
    SpeedAll    = NaN(1700,length(FolderList)-2);
    MeanVelX    = NaN(1700,length(FolderList)-2);
    MeanVelY    = NaN(1700,length(FolderList)-2);
    RMSVelAll   = NaN(1700,length(FolderList)-2);
    AngVelAll   = NaN(1700,length(FolderList)-2);
    MSDAll      = NaN(1700,length(FolderList)-2);
    MSDTempAll  = NaN(WSize,1700,length(FolderList)-2);
    QAll        = NaN(1700,length(FolderList)-2);
    QTempAll    = NaN(WSize,1700,length(FolderList)-2);
    ChiAll      = NaN(1700,length(FolderList)-2);
    ChiTempAll  = NaN(WSize,1700,length(FolderList)-2);
    DelaunayNeighborsTimeReversedAll    = NaN(1700,length(FolderList)-2);
    DelaunayNeighborsAll                = NaN(1700,length(FolderList)-2);
    NumNewNeighborsAll                  = NaN(1700,length(FolderList)-2);
    NewNeighborsPerCellAll              = NaN(1700,length(FolderList)-2);
    NewNeigborsPerDiffusiveCellAll      = NaN(1700,length(FolderList)-2);
    rMeanAll =  NaN(1700,length(FolderList)-2);
    IDStartAll = NaN(1700,length(FolderList)-2);
    NumNewNeighborsTempAll = NaN(1700,length(FolderList)-2);
    TrackMatAll = cell(length(FolderList)-2,1);
    % Allocate variables for orientation analysis:
    AngleVecFields = cell(length(FolderList)-2,1);
    AngleOrientation = cell(length(FolderList)-2,1);
    AngleVelField = cell(length(FolderList)-2,1);
    % Allocate variables for cell division analysis:
    DivTimesAll = NaN(1700,length(FolderList)-2);
    NumDivsAll = NaN(1700,length(FolderList)-2);
    
    for FolderNum2 = 3:length(FolderList)
        FolderNum2
        
        % Velocity field analysis:
        % Check if several folders actually exist:
        folderVelField = fullfile(strcat(FolderList(FolderNum2).name,'\Results'));
        % If it exists do all following analysis:
        if (exist(folderVelField) == 7)
            % Go to results folder:
            cd(strcat(FolderList(FolderNum2).name,'\Results'));
            folder_number = FolderNum2;
            % Load velocity field:
            load('VelocityField.mat')
            load('PositionsX.mat')
            load('PositionsY.mat')
            % Number of images analyzed:
            m = size(VelField,4);
            % Get tracking matrix from velocity field:
            [TrackMat] = GenTrackMat(VelField,im_size,POIx,POIy,SubPxResolution,CenterSpeed);
            % Get total MSD, Order Parameter and 4 point susceptibility for one velocity field:
            [MSD, Q, Chi, tau] = MSDChiOrder(TrackMat,dt,CSize,OvThresh);
            % Get time resolved MSD, Order Parameter and 4 point susceptibility for one velocity field:
            %[MSDTemp, QTemp, ChiTemp] = MSDChiOrderTimeWindow(TrackMat,WSize,dt,CSize,OvThresh);
            [MSDTemp, QTemp, ChiTemp,NumNewNeighborsTemp,TrackMatTemp] = ...
                MSDChiOrderTimeWindow(VelField,WSize,dt,CSize,OvThresh,im_size,POIx,POIy,SubPxResolution,Sampling,CenterSpeed);
            % Get nearest neighbor relations:
            [DelaunayNeighborsTimeReversed,DelaunayNeighbors,NumNewNeighbors,...
                NewNeighborsPerCell,NewNeigborsPerDiffusiveCell] = ...
                NeighborsAndTracks(TrackMat,DataPath,FolderList(FolderNum2).name,PxSize);
            % Get Diffusion Coefficient exponent for one velocity field:
            p = polyfit(log10(tau),log10(nanmean(MSD(1:end-1,:),2)),1);
            yfit = polyval(p,log10(tau));
            % Calculate diffusion exponent and Q for gliding window of size "WSize":
            % Get time dependent diffusion coefficient exponents:
            p = NaN(size(MSDTemp,2),1);
            for i = 1+(WSize-1)/2:size(MSDTemp,2)-1-((WSize-1)/2)
                ptemp = polyfit(log10(tau(1:size(MSDTemp,1)-2)),log10(MSDTemp(1:end-2,i)),1);
                p(i) = ptemp(1);
            end
            % Velocity Order parameter (Szabo et al 2006):
            [VTemp] = VelOrderWindow(VelField,WSize,dt);
            [AutoCorr,Speed,RMSVel,Area,AngVel,MeanVelXTemp,MeanVelYTemp] ...
                = AnalyzeVelField(VelField,im_size,m,POIx,POIy,CenterSpeed);
            
            % Plug parameters into one global construct:
            MSDAll(1:size(MSD,1),folder_number-2)                         = nanmean(MSD,2);
            QAll(1:size(Q,2),folder_number-2)                             = Q;
            ChiAll(1:size(Chi,2),folder_number-2)                         = Chi;
            pAll(1:length(p),folder_number-2)                             = p;
            EndPos                                                        = length(DelaunayNeighbors);
            DelaunayNeighborsTimeReversedAll(1:EndPos,folder_number-2)    = DelaunayNeighborsTimeReversed;
            DelaunayNeighborsAll(1:EndPos,folder_number-2)                = DelaunayNeighbors;
            NumNewNeighborsAll(folder_number-2)                           = NumNewNeighbors;
            NewNeighborsPerCellAll(folder_number-2)                       = NewNeighborsPerCell;
            NewNeigborsPerDiffusiveCellAll(folder_number-2)               = NewNeigborsPerDiffusiveCell;
            NumNewNeighborsTempAll(:,folder_number-2)                     = NumNewNeighborsTemp;
            MeanVelX(1:m,folder_number-2)                                 = MeanVelXTemp;
            MeanVelY(1:m,folder_number-2)                                 = MeanVelYTemp;
            SpeedAll(1:m,folder_number-2)                                 = Speed;
            RMSVelAll(1:m,folder_number-2)                                = RMSVel;
            AreaAll(1:m,folder_number-2)                                  = Area;
            AngVelAll(1:m,folder_number-2)                                = AngVel;
            TrackMatAll{folder_number-2,1}                                = TrackMatTemp;
            
            % Time dependent variables:
            ObjNum = size(MSDTemp,2);
            %ObjNum = size(VTemp,2);
            MSDTempAll(:,1:ObjNum,folder_number-2)              = MSDTemp;
            QTempAll(:,1:ObjNum,folder_number-2)                = QTemp;
            ChiTempAll(:,1:ObjNum,folder_number-2)              = ChiTemp;
            VTempAll(:,1:ObjNum,folder_number-2)                = VTemp;
            AutoCorrAll(1:size(AutoCorr,1),1:m,folder_number-2) = AutoCorr;
            
            % Go one level up in the folder and save the generated tracking matrix
            save('TrackinMatrixFinal.mat','TrackMat','-v7.3') 
            % cd('..');
            cd(curr_directory)
            
            %cd('..')    
            
        end
        
        
        
        % Orientation field analysis:
        % Check if several folders actually exist:
        folderOrientation = fullfile(strcat(FolderList(FolderNum2).name,'\Results Orientation Analysis'));
        % If it exists do all following analysis:
        if (exist(folderVelField,'dir') == 7) && (exist(folderOrientation,'dir') == 7)
            cd(folderOrientation);
            % Load vector field of orientations and parameter for the meshgrid:
            load('CellOrientation.mat')
            load('Settings.mat')
            
            % Create meshgrid of cell orientation field:
            [xCellOrientation,yCellOrientation] = meshgrid(1:MSize:size(CellOrientation{1},2)*MSize,1:MSize:size(CellOrientation{1},1)*MSize);
            % Transform CellOrientation from cell to matrix:
            CellOrientationField = NaN(size(CellOrientation{1},1),size(CellOrientation{1},2),2,m);
            for i = 1:length(CellOrientation)
                if ~isempty(CellOrientation{i})
                    CellOrientationField(:,:,:,i) = CellOrientation{i};
                end
            end
            
            % Get time dependent orientation towards x-axis:
            xAxMat = ones(size(VelField));
            xAxMat(:,:,2,:) = 0;
            [AngleTemp] = RelativeVectorFieldOrientation(CellOrientationField(:,:,:,1:end-1),xAxMat(:,:,:,1:end-1),...
                xCellOrientation,yCellOrientation,x,y,1);
            AngleOrientation{FolderNum2-2} = AngleTemp;
            
            % Get time dependent velocity vector towards x-axis:
            [AngleTemp] = RelativeVectorFieldOrientation(xAxMat(:,:,:,1:end-1),VelField(:,:,:,1:end-1),...
                x,y,x,y,1);
            AngleVelField{FolderNum2-2} = AngleTemp;
            
            % Calculate angles between both fields:
            % Window size for temporal smoothing:
            WinSize = 51;
            % Ignore last entry as this image is often a purely black one
            [AngleTemp] = RelativeVectorFieldOrientation(CellOrientationField(:,:,:,1:end-1),VelField(:,:,:,1:end-1),...
                xCellOrientation,yCellOrientation,x,y,WinSize);
            AngleVecFields{FolderNum2-2} = AngleTemp;
            %cd('..')
            %cd('..')
            cd(curr_directory)
        end
        
        
        % Cell division analysis:        
        % Check if several folders actually exist:
        folderDivisions = fullfile(strcat(FolderList(FolderNum2).name,'\Results Cell Division Analysis'));
        % If it exists do all following analysis:
        if (exist(folderDivisions) == 7) 
            
            cd(folderDivisions);
            % Load vector field of orientations and parameter for the meshgrid:
            load('CellDivisions.mat')
            
            % create another info about the starting point of each division:
            StartPos = NaN(length(CellProps),1);
            for i = 1:length(CellProps)
                StartPos(i) = find(isnan(CellProps(i).centroid(:,1)) == 0,1);
            end
            % Cumulative division events:
            DivTimes= unique(StartPos);
            NumDivs = NaN(length(DivTimes),1);
            for i = 1:length(DivTimes)
                NumDivs(i) =  numel(find(StartPos == DivTimes(i)));
            end
            NumDivs = cumsum(NumDivs);
            DivTimesAll(1:size(DivTimes,1),folder_number-2) = dt*(DivTimes-1);
            NumDivsAll(1:size(NumDivs,1),folder_number-2) = NumDivs;
            %cd('..')
        end
         
        
        
         % Go back to top level directory:
         cd(curr_directory)
        
    end
    
    % Save Variables:
    % Velocity field associated data:
    if (exist(folderVelField) == 7)
        save(SaveNameVelField, 'MSDAll','MSDTempAll','QAll','QTempAll','ChiAll',...
            'ChiTempAll','pAll','tau', 'SpeedAll', 'VTempAll', 'AngVelAll', ...
            'RMSVelAll','AutoCorrAll', 'AreaAll', 'WSize', 'dt','im_size',...
            'CSize','OvThresh','DelaunayNeighborsTimeReversedAll',...
            'DelaunayNeighborsAll','NumNewNeighborsAll',...
            'NewNeighborsPerCellAll','NewNeigborsPerDiffusiveCellAll',...
            'NumNewNeighborsTempAll','MeanVelX','MeanVelY')
        VelFieldSize = size(VelField);
        VelFieldSize = VelFieldSize(1:2);
        save('Settings FlowFieldAnalysis.mat','s_size','POIx','POIy','im_size','dt','CSize','OvThresh','WSize','SubPxResolution','Sampling','VelFieldSize','PxSize','CenterSpeed')
        save(SaveNameTrackMat,'WSize', 'dt','TrackMatAll','-v7.3')     
    end
    % Orientation field associated data:
    if (exist(folderVelField) == 7) && (exist(folderOrientation) == 7)
        % Save Variables:
        save(SaveNameOrientation,'AngleVecFields','AngleOrientation','AngleVelField','-v7.3')
        VelFieldSize = size(VelField);
        VelFieldSize = VelFieldSize(1:2);
        save('Settings VectorFieldAnalysis.mat','im_size','dt','PxSize')
    end
    
    if (exist(folderDivisions) == 7) 
       % Save Variables:
       save(SaveNameDivision,'DivTimesAll','NumDivsAll','-v7.3')
    end
    
    cd(curr_directory_outer)
    
end









