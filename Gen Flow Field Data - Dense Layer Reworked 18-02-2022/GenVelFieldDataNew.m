%% Set folder:
clear all
beta = pwd;
beta = addpath(beta);  % Adds the path where the m-files are stored in
StartPath=pwd;
DataPath=uigetdir(StartPath, 'Chose the folder with the images');
alpha = cd(DataPath);   % legt jpg_path2 als aktuellen folder fest
FolderList = dir();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up storage structure and some parameters:
VelFieldData = struct('Name',[],'pAll',[],'VTempAll',[],'AreaAll',[],...
    'RMSVelAll',[],'AngVelAll',[],'MSDAll',[],'MSDTempAll',[],...
    'QAll',[],'QTempAll',[],'ChiAll',[],'ChiTempAll',[],'tau',[],...
    'Speed',[],'AutoCorrAll',[],'DelaunayNeighborsTimeReversedAll',[],...
    'DelaunayNeighborsAll',[],'NumNewNeighborsAll',[],...
    'NewNeighborsPerCellAll',[],'NewNeigborsPerDiffusiveCellAll',[],...
    'NumNewNeighborsTempAll',[],'MeanVelX',[],'MeanVelY',[],...
    'rMeanAll',[],'IDStartAll',[],'WSize',[]);

VecFieldData = struct('Name',[],'AngleVecFields',[],'AngleOrientation',[],'AngleVelField',[]);

CellDivData = struct('Name',[],'DivTimesAll',[],'NumDivsAll',[]);

TrackMat = struct('Name',[],'TrackMat',[],'WSize',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in files and transfer data into struct:
% Remove all non directories:
Idx = cat(1,FolderList.isdir);
FolderList(Idx~=1) = [];
% remove ".." and "." from the folder list:
for i = length(FolderList):-1:1
    if strcmp(FolderList(i).name,'..') || strcmp(FolderList(i).name,'.')
        FolderList(i) = [];
    end
end
% Save start folder:
curr_directory = pwd;

for folder_number2 = 1:length(FolderList)
    folder_number = folder_number2
    % Go to results folder:
    cd(strcat(FolderList(folder_number).name));
    % List files:
    FileList = dir(pwd);
    % Remove all non file entries:
    Idx = cat(1,FileList.isdir);
    FileList(Idx==1) = [];
    % Get .mat file for all vector fields:
    SaveNameVelField = FolderList(folder_number).name;
    SaveNameVelField = sprintf('%s.mat',SaveNameVelField);
    SaveNameTrackMat = FolderList(folder_number).name;
    SaveNameTrackMat = sprintf('%s %s.mat','TrackingMatrixTemp',SaveNameTrackMat);
    SaveNameOrientation = FolderList(folder_number).name;
    SaveNameOrientation = sprintf('Orientation Analysis %s.mat',SaveNameOrientation);
    SaveNameDivision = FolderList(folder_number).name;
    SaveNameDivision = sprintf('Division Analysis %s.mat',SaveNameDivision);
    
    %% Velocity field data:
    if (exist(SaveNameVelField) == 2)
        % Load data velocity field data:
        load(SaveNameVelField)
        % Get data into struct:
        VelFieldData(folder_number).Name                              = SaveNameVelField(1:end-4);
        MeanVelX(MeanVelX == 0) = NaN;
        VelFieldData(folder_number).MeanVelX                          = MeanVelX;
        MeanVelY(MeanVelY == 0) = NaN;
        VelFieldData(folder_number).MeanVelY                          = MeanVelY;
        pAll(pAll == 0) = NaN;
        VelFieldData(folder_number).pAll                              = pAll;
        VTempAll(VTempAll == 0) = NaN;
        VelFieldData(folder_number).VTempAll                          = VTempAll;
        AreaAll(AreaAll == 0) = NaN;
        VelFieldData(folder_number).AreaAll                           = AreaAll;
        RMSVelAll(RMSVelAll == 0) = NaN;
        VelFieldData(folder_number).RMSVelAll                         = RMSVelAll;
        AngVelAll(AngVelAll == 0) = NaN;
        VelFieldData(folder_number).AngVelAll                         = AngVelAll;
        MSDAll(MSDAll == 0) = NaN;
        VelFieldData(folder_number).MSDAll                            = MSDAll;
        MSDTempAll(MSDTempAll == 0) = NaN;
        VelFieldData(folder_number).MSDTempAll                        = MSDTempAll;
        QAll(QAll == 0) = NaN;
        VelFieldData(folder_number).QAll                              = QAll;
        QTempAll(QTempAll == 0) = NaN;
        VelFieldData(folder_number).QTempAll                          = QTempAll;
        ChiAll(ChiAll == 0) = NaN;
        VelFieldData(folder_number).ChiAll                            = ChiAll;
        ChiTempAll(ChiTempAll == 0) = NaN;
        VelFieldData(folder_number).ChiTempAll                        = ChiTempAll;
        tau(tau == 0) = NaN;
        VelFieldData(folder_number).tau                               = tau;
        SpeedAll(SpeedAll == 0) = NaN;
        VelFieldData(folder_number).Speed                             = SpeedAll;
        AutoCorrAll(AutoCorrAll == 0) = NaN;
        VelFieldData(folder_number).AutoCorrAll                       = AutoCorrAll;
        DelaunayNeighborsTimeReversedAll(DelaunayNeighborsTimeReversedAll == 0) = NaN;
        VelFieldData(folder_number).DelaunayNeighborsTimeReversedAll  = DelaunayNeighborsTimeReversedAll;
        DelaunayNeighborsAll(DelaunayNeighborsAll == 0) = NaN;
        VelFieldData(folder_number).DelaunayNeighborsAll              = DelaunayNeighborsAll;
        NumNewNeighborsAll(NumNewNeighborsAll == 0) = NaN;
        VelFieldData(folder_number).NumNewNeighborsAll                = NumNewNeighborsAll;
        NewNeighborsPerCellAll(NewNeighborsPerCellAll == 0) = NaN;
        VelFieldData(folder_number).NewNeighborsPerCellAll            = NewNeighborsPerCellAll;
        NewNeigborsPerDiffusiveCellAll(NewNeigborsPerDiffusiveCellAll == 0) = NaN;
        VelFieldData(folder_number).NewNeigborsPerDiffusiveCellAll    = NewNeigborsPerDiffusiveCellAll;
        NumNewNeighborsTempAll(NumNewNeighborsTempAll == 0)             = NaN;
        VelFieldData(folder_number).NumNewNeighborsTempAll            = NumNewNeighborsTempAll;
        VelFieldData(folder_number).WSize                             = WSize;
    end
    %% Cell orientation data:
    if (exist(SaveNameOrientation) == 2)
        % Load data orientation field data:
        load(SaveNameOrientation)
        % Get data into struct:
        VecFieldData(folder_number).Name              = SaveNameVelField(1:end-4);
        VecFieldData(folder_number).AngleVecFields    = AngleVecFields;
        VecFieldData(folder_number).AngleOrientation  = AngleOrientation;
        VecFieldData(folder_number).AngleVelField     = AngleVelField;
    end
    %% Cell division data:
    if (exist(SaveNameDivision) == 2)
        % Load cell division data:
        load(SaveNameDivision)
        CellDivData(folder_number).Name        = SaveNameVelField(1:end-4);
        CellDivData(folder_number).DivTimesAll = DivTimesAll;
        CellDivData(folder_number).NumDivsAll  = NumDivsAll;
    end
    
    %% Tracking matrix data:
    if (exist(SaveNameTrackMat) == 2)
        % Load data orientation field data:
        load(SaveNameTrackMat)
        % Get data into struct:
        TrackMat(folder_number).TrackMat = TrackMatAll;
        TrackMat(folder_number).Name     = SaveNameVelField(1:end-4);
        TrackMat(folder_number).WSize    = WSize;
    end
    
    % Go back to top level directory:
    cd(curr_directory)
end
% Save all data:
save('VelFieldData.mat','VelFieldData','-v7.3')
save('VecFieldData.mat','VecFieldData','-v7.3')
save('CellDivData.mat','CellDivData','-v7.3')
save('TrackingMatrix All.mat','TrackMat','-v7.3')
