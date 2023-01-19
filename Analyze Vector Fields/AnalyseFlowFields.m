% Analyze the obtained orientation and velocity fields and store the
% resulting data for further processing.


clear all;

mPath = mfilename('fullpath');
Idx = max(strfind(mPath,filesep));
mPath = mPath(1:Idx);
addpath(mPath)
% add path for velocity field analysis:
addpath(strcat(mPath,filesep,'Velocity Fields'))
% add path of the orientation analysis:
addpath(strcat(mPath,filesep,'Orientation Analysis'))
% add path of the circular statistics toolbox:
% addpath(strcat(mPath,'\CircStat'))

beta = pwd;
beta = addpath(beta);  % Adds the path where the m-files are stored in
StartPath=mPath;
DataPath=uigetdir(StartPath, 'Chose the folder with the images');
alpha = cd(DataPath);   % legt jpg_path2 als aktuellen folder fest
FolderListOuter = dir();

% load parameters:
[s_size,im_size,ImPhysSize,dt,CSize,OvThresh,WSize,...
    CenterSpeed,SubPxResolution,Sampling,MaxVisibleTime] = ParameterFunctionMain();
% Check if inputs are valid:
CheckForValidInputs(s_size,im_size,ImPhysSize,dt,CSize,OvThresh,WSize,...
    CenterSpeed,SubPxResolution,Sampling,MaxVisibleTime)

% pixel size of image in µm/px:
PxSize = ImPhysSize/im_size;
% calculate meshgrid based on the points of interest difference:
s_size = ceil(s_size(1));
POIx = [s_size/2+1:s_size:im_size(1)];
POIy = [s_size/2+1:s_size:im_size(2)];
[POIx, POIy] = meshgrid(POIx,POIy);
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
    
    % Save start folder and set save names:
    curr_directory = pwd;
    SaveNameVelField = curr_directory(max(strfind(pwd,filesep))+1:end);
    SaveNameVelField = sprintf('%s.mat',SaveNameVelField);
    SaveNameTrackMat = curr_directory(max(strfind(pwd,filesep))+1:end);
    SaveNameTrackMat = sprintf('%s %s.mat','TrackingMatrixTemp',SaveNameTrackMat);
    SaveNameOrientation = curr_directory(max(strfind(pwd,filesep))+1:end);
    SaveNameOrientation = sprintf('Orientation Analysis %s.mat',SaveNameOrientation);
    SaveNameDivision = curr_directory(max(strfind(pwd,filesep))+1:end);
    SaveNameDivision = sprintf('Division Analysis %s.mat',SaveNameDivision);
    
    % Allocate variables for velocity field analysis:    
    SpeedAll    = NaN(1700,length(FolderList)-2);   
    RMSVelAll   = NaN(1700,length(FolderList)-2);   
    MSDAll      = NaN(1700,length(FolderList)-2);
    MSDTempAll  = NaN(WSize,1700,length(FolderList)-2);
    MSDCagedVarAll  = NaN(1700,10,length(FolderList)-2);
    MSDCagedTempVar = NaN(WSize,1700,10,length(FolderList)-2);
    QAll        = NaN(1700,length(FolderList)-2);
    QTempAll    = NaN(WSize,1700,length(FolderList)-2);
    ChiAll      = NaN(1700,length(FolderList)-2);
    ChiTempAll  = NaN(WSize,1700,length(FolderList)-2);       
    NumNewNeighborsAll     = NaN(1700,length(FolderList)-2);
    NumNewNeighborsTempAll = NaN(1700,length(FolderList)-2); 
    TrackMatAll = cell(length(FolderList)-2,1);
    % Allocate variables for orientation analysis:
    AngleOrientation = cell(length(FolderList)-2,1);
    % Allocate variables for cell division analysis:
    DivTimesAll = NaN(1700,length(FolderList)-2);
    MeanSpeedDivAll = NaN(1700,length(FolderList)-2);
    MeanSpeedNonDivAll = NaN(1700,length(FolderList)-2);
    NumDivsAll = NaN(1700,length(FolderList)-2);
    StartDivsion = NaN(1700,length(FolderList)-2);
    EndDivsion = NaN(1700,length(FolderList)-2);
    Visibility = NaN(1700,length(FolderList)-2);
    DivisionCenters = NaN(ceil(MaxVisibleTime),2,1000,length(FolderList)-2);
    
    for FolderNum2 = 3:length(FolderList)
        folder_number = FolderNum2;
        % Velocity field analysis:
        % Check if several folders actually exist:
        folderVelField = fullfile(strcat(curr_directory ,filesep,FolderList(FolderNum2).name,filesep,'Results'));
        % If it exists do all following analysis:
        if (exist(folderVelField) == 7)
            % Go to results folder:
            cd(strcat(FolderList(FolderNum2).name,filesep,'Results'));            
            % Load velocity field:
            load('VelocityField.mat')
            load('PositionsX.mat')
            load('PositionsY.mat')            
            % remove drift in velocity field            
            for i = 1:size(VelField,4)
                VelField(:,:,1,i) = VelField(:,:,1,i) - VelFieldDrift(1,i);
                VelField(:,:,2,i) = VelField(:,:,2,i) - VelFieldDrift(2,i);
            end
            
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
            % Get caged MSD:
            [CagedMSD] = MSDCaged(TrackMat,dt);
            [CagedMSDTemp] = MSDCagedTemp(VelField,WSize,dt,im_size,POIx,POIy,SubPxResolution,Sampling,CenterSpeed);
            % Get nearest neighbor relations:
            [NumNewNeighbors] = NearestNeighbors(TrackMat);
            % Calculate speeds:
            [Speed,RMSVel] = AnalyzeVelField(VelField,m,CenterSpeed);
            
            % Plug parameters into one global construct:
            MSDAll(1:size(MSD,1),folder_number-2)                         = nanmean(MSD,2);
            MSDCagedVarAll(1:size(MSD,1),:,folder_number-2)               = CagedMSD;
            MSDCagedTempVar(:,1:size(CagedMSDTemp,2),:,folder_number-2)   = CagedMSDTemp;
            QAll(1:size(Q,2),folder_number-2)                             = Q;
            ChiAll(1:size(Chi,2),folder_number-2)                         = Chi;
            NumNewNeighborsAll(folder_number-2)                           = NumNewNeighbors;            
            NumNewNeighborsTempAll(1:length(NumNewNeighborsTemp),folder_number-2) = NumNewNeighborsTemp;            
            SpeedAll(1:m,folder_number-2)                                 = Speed;
            RMSVelAll(1:m,folder_number-2)                                = RMSVel;            
            TrackMatAll{folder_number-2,1}                                = TrackMatTemp;            
            % Time dependent variables:
            ObjNum = size(MSDTemp,2);
            MSDTempAll(:,1:ObjNum,folder_number-2)              = MSDTemp;
            QTempAll(:,1:ObjNum,folder_number-2)                = QTemp;
            ChiTempAll(:,1:ObjNum,folder_number-2)              = ChiTemp;     
            % Go one level up in the folder and save the generated tracking matrix
            save('TrackinMatrixFinal.mat','TrackMat','-v7.3') 
            cd(curr_directory)           
        end
        
        % Orientation field analysis:
        % Check if several folders actually exist:
        folderOrientation = fullfile(strcat(curr_directory ,filesep,FolderList(FolderNum2).name,filesep,'Results Orientation Analysis'));
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
            cd(curr_directory)
        end
        
        % Cell division analysis:        
        % Check if several folders actually exist:
        folderDivisions = fullfile(strcat(curr_directory ,filesep,FolderList(FolderNum2).name,filesep,'Results Cell Division Analysis'));
        % If it exists do all following analysis:
        if (exist(folderDivisions,'file') == 7) 
            
            cd(folderDivisions);
            % Load vector field of orientations and parameter for the meshgrid:
            load('CellDivisions.mat')
            
            % create another info about the starting point of each division:
            StartPos = NaN(length(CellProps),1);
            EndPos = NaN(length(CellProps),1);
            VisibleTime = StartPos;
            for i = 1:length(CellProps)
                StartPos(i) = find(isnan(CellProps(i).centroid(:,1)) == 0,1);
                EndPos(i) = find(isnan(CellProps(i).centroid(:,1)) == 0,1,'last');
                % How long was cell division figure visible?
                VisibleTime(i) = CellProps(i).totalVisibleCount;                
            end
            % Exclude all divisionevents persistent longer than 50 images,
            % as these are likely cell extrusions or similar things. 50
            % images corresponds to 150min.                  
            IdxDivTooLong = find(VisibleTime>MaxVisibleTime);
            StartPos(IdxDivTooLong) = [];
            EndPos(IdxDivTooLong) = [];
            VisibleTime(IdxDivTooLong) = [];
            Idx = 1:size(CellProps,2);
            Idx(IdxDivTooLong) = [];
            % Cumulative division events:
            DivTimes= unique(StartPos);
            NumDivs = NaN(length(DivTimes),1);
            for i = 1:length(DivTimes)
                NumDivs(i) =  numel(find(StartPos == DivTimes(i)));
            end
            NumDivs = cumsum(NumDivs);
            DivTimesAll(1:size(DivTimes,1),folder_number-2) = dt*(DivTimes-1);
            NumDivsAll(1:size(NumDivs,1),folder_number-2) = NumDivs;
            % Center start and length of divisions:
            StartDivsion(1:size(StartPos,1),folder_number-2) = StartPos;
            EndDivsion(1:size(StartPos,1),folder_number-2) = EndPos;
            Visibility(1:size(StartPos,1),folder_number-2) = VisibleTime;
            count = 0;
            for j = Idx
                count = count+1;
                DivisionCenters(1:EndPos(count)+1-StartPos(count),1:2,count,folder_number-2) = CellProps(j).centroid(StartPos(count):EndPos(count),:);
            end
            % Go back to top level directory:
            cd(curr_directory)
        end   
        
        % Calculate local velocities around divisions and in the rest of
        % the image for each single image:
        if (exist(folderVelField) == 7) && (exist(folderDivisions,'file') == 7)
            cd(strcat(FolderList(FolderNum2).name,filesep,'Results'));            
            [MeanSpeedDiv,MeanSpeedNonDiv] = VelFieldAroundDivisions(CellProps,VelField,VelFieldDrift,CenterSpeed,x,y,StartPos,EndPos,Idx);           
            MeanSpeedNonDivAll(1:size(MeanSpeedNonDiv,1),folder_number-2) = MeanSpeedNonDiv;
            MeanSpeedDivAll(1:size(MeanSpeedNonDiv,1),folder_number-2) = MeanSpeedDiv;
            %figure; plot(movmean(MeanSpeedDiv./MeanSpeedNonDiv,11,'omitnan'))
            %nanmean(MeanSpeedDiv./MeanSpeedNonDiv)
            %min(MeanSpeedDiv./MeanSpeedNonDiv)
            %nanstd(MeanSpeedDiv./MeanSpeedNonDiv)
        end
        % Go back to top level directory:
        cd(curr_directory) 
        
    end
    
    if (exist(folderVelField) == 7)
        save(SaveNameVelField, 'MSDAll','MSDTempAll','QAll','QTempAll','ChiAll',...
            'ChiTempAll','tau', 'SpeedAll', 'RMSVelAll','NumNewNeighborsAll','NumNewNeighborsTempAll',...
            'WSize','dt','im_size','CSize','OvThresh', 'MSDCagedVarAll', 'MSDCagedTempVar')
        VelFieldSize = size(VelField);
        VelFieldSize = VelFieldSize(1:2);
        save('Settings FlowFieldAnalysis.mat','s_size','POIx','POIy','im_size','dt','CSize','OvThresh','WSize','SubPxResolution','Sampling','VelFieldSize','PxSize','CenterSpeed')
        save(SaveNameTrackMat,'WSize', 'dt','TrackMatAll','-v7.3')     
    end
    % Orientation field associated data:
    if (exist(folderVelField) == 7) && (exist(folderOrientation) == 7)
        % Save Variables:
        save(SaveNameOrientation,'AngleOrientation','-v7.3')
        VelFieldSize = size(VelField);
        VelFieldSize = VelFieldSize(1:2);
        save('Settings VectorFieldAnalysis.mat','im_size','dt','PxSize')
    end
    % Cell division data:
    if (exist(folderVelField) == 7) && (exist(folderDivisions,'file') == 7)
        save(SaveNameDivision,'DivTimesAll','NumDivsAll','MaxVisibleTime','MeanSpeedDivAll','MeanSpeedNonDivAll','-v7.3')
    elseif (exist(folderDivisions) == 7)
        % Save Variables:
        save(SaveNameDivision,'DivTimesAll','NumDivsAll','MaxVisibleTime','-v7.3')
    end    
    cd(curr_directory_outer)    
end









