% to do:
% 1) add all input parameters for the cell division detection
% 2) save settings for both cell division detection and cell orientation
% analysis
% 3) make an automatic attempt at predicting the size of the velocity field
% containing matrix



%Use imhistmatch for subsequent images for better results?

clear all;
mPath = mfilename('fullpath');
Idx = max(strfind(mPath,'\'));
mPath = mPath(1:Idx);
addpath(mPath)


% add path for BM3D filter:
addpath(strcat(mPath,'\BM3D Filter'))
% add path of the circular statistics toolbox:
addpath(strcat(mPath,'\CircStat'))
% add path for MatPIV:
addpath(strcat(mPath,'\MatPIV 1.7\src'))
addpath(strcat(mPath,'\MatPIV 1.7\masking'))
addpath(strcat(mPath,'\MatPIV 1.7\filters'))
addpath(strcat(mPath,'\MatPIV 1.7\postprocessing'))
addpath(strcat(mPath,'\MatPIV 1.7\PTV'))
addpath(strcat(mPath,'\MatPIV 1.7\Demo3'))
% add path for cell division detection:
addpath(strcat(mPath,'\Cell Division Analysis'))
% add path for orientation analysis:
addpath(strcat(mPath,'\Orientation Analysis'))

% load parameters:
[UsePIV,met,wins,Overlap,globtrld,loctrld,kernelsize,med,int,...
    UseCellDivDetect,CType,NetName,MakeVideo,MedFiltSize,FiltType,...
    BackThresh,WSize,MaxTimeDiff,UseOrientationAnalysis,MSize,rho] = ParameterFunctionMain;


% Load templates used for cell division detection and the ANN:
if UseCellDivDetect == 1
    % load templates used in cross correlation analysis:
    cd(strcat(mPath,'\Cell Division Analysis\DivisionTemplate\Templates\',CType))
    FileList = dir();
    FileList([FileList.isdir]==1) = [];
    Template = cell(length(FileList),1);
    for i = 1:length(FileList)
        tmp = imread(FileList(i).name);
        if size(tmp,3) == 3
            tmp = rgb2gray(tmp);
        elseif size(tmp,3) == 4
            tmp = rgb2gray(tmp(:,:,1:3));
        end
        Template{i} = double(tmp);
    end
    
    % load the trained ANN:
    cd(strcat(mPath,'\Cell Division Analysis\TrainedANNs'))
    load(strcat(NetName,'.mat'))
end


% Get image path:
Startpath=mPath;
jpg_path2=uigetdir(Startpath, 'Chose the folder with the images');
alpha = cd(jpg_path2);   % legt jpg_path2 als aktuellen folder fest
folderList = dir();

% Go through all folders
for folder_number = 3:length(folderList)
    
    curr_directory = pwd;
    if strcmp(folderList(folder_number).name(1:end),'Thumbs.db') == 1
        folderList(folder_number) = [];
    elseif strcmp(folderList(folder_number).name(1:end),'SavedData.mat') == 1
        folderList(folder_number) = [];
    end
    cd(folderList(folder_number).name);
    
    if UsePIV
        % Check if several folders actually exist:
        folderVelField = fullfile(strcat(jpg_path2,'\',folderList(folder_number).name),'Results');
        % Checks the existence of the folder with the name "Results". If it
        % does not exist it is created.
        if (exist(folderVelField) == 0)
            mkdir('Results');
        end
    end
    
    if UseCellDivDetect
        folderDivision = fullfile(strcat(jpg_path2,'\',folderList(folder_number).name),'Results Cell Division Analysis');
        % Checks the existence of the folder with the name "Results". If it
        % does not exist it is created.
        if (exist(folderDivision) == 0)
            mkdir('Results Cell Division Analysis');
        end
    end
    
    if UseOrientationAnalysis
        % Check if several folders actually exist:
        folderOrientation = fullfile(strcat(jpg_path2,'\',folderList(folder_number).name),'Results Orientation Analysis');
        % Checks the existence of the folder with the name "Results". If it
        % does not exist it is created.
        if (exist(folderOrientation) == 0)
            mkdir('Results Orientation Analysis');
        end
    end
    
    % set first image:
    [im1,~,m] = Reader(strcat(jpg_path2,'\',folderList(folder_number).name),1);
    if ~isa(im1,'uint8')
        im1 = im2uint8(im1);
    end
    
    % Denoising using the BM3D filter:
    [~, im1] = BM3D(1, double(im1), 0.5*std(double(im1(:))));
    
    %     % Velocity field:
    %     VelField = NaN(119,159,2,m);
    %     VelField(:,:,1,1) = zeros(119, 159);
    %     VelField(:,:,2,1) = zeros(119, 159);
    %     % Size for Final images with 16Px cross Correlation with 0.25 overlap
    %     VelField = NaN(79,106,2,m);
    %     VelField(:,:,1,1) = zeros(79, 106);
    %     VelField(:,:,2,1) = zeros(79, 106);
    %     % Size for Final images with 32Px cross Correlation with 0.6 overlap
    %     VelField = NaN(72,97,2,m);
    %     VelField(:,:,1,1) = zeros(72, 97);
    %     VelField(:,:,2,1) = zeros(72, 97);
    
    % Velocity field - pre allocation:
    [x,y,u,v,snr] = matpiv(im1,im1,wins,1,Overlap,met);
    VelField = NaN(size(x,1),size(x,2),2,m);
    VelField(:,:,1,1) = 0;
    VelField(:,:,2,1) = 0;
    % Centers of cell divisions:
    CenterDivisions = NaN(250,2,m);
    % Direction of cells
    CellOrientation = cell(m,1);
    for i=1:m-1
        
        Percentage =((folder_number-3)*(m-1) + i)/((length(folderList)-2)*m-1)
        [images]=Reader(strcat(jpg_path2,'\',folderList(folder_number).name),[i+1]);
        %images = imhistmatch(double(images),im1,256);
        if ~isa(images,'uint8')
            images = im2uint8(images);
        end
        im2 = images;
        [~, im2] = BM3D(1, double(im2), 0.5*std(double(im2(:))));
        
        % PIV calculation
        if UsePIV
            % The result consists of four matrices x, y, u and v which are measured
            % (when no world coordinates are given) in pixels and pixels/second.
            % With T = 1 (as here) it gives pixels change between two sucessive
            % images (whatever there time difference may be)
            
            [x,y,u,v,snr] = matpiv(im1,im2,wins,1,Overlap,met);
            % For final windowsize of 64px to receive a similar spatial
            % resolution as for the 16px window with 25% overlap
            %[x,y,u,v,snr] = matpiv(im1,im2,wins,1,0.83,met);
            
            %     [HsRt] = DriftCorrection(images(:,:,1),images(:,:,2));
            %     du = HsRt(3,1)/T;
            %     dv = HsRt(3,2)/T;
            %     u = u + du;
            %     v = v + dv;
            
            [gu,gv]=globfilt(x,y,u,v,globtrld);
            % local median filter:
            [lu,lv]=localfilt(x,y,gu,gv,loctrld,med,kernelsize);
            % interpolate outliers
            [fu,fv]=naninterp(lu,lv,int);
            
            % Set velocity field:
            VelField(:,:,1,i+1) = fu;
            VelField(:,:,2,i+1) = fv;
            nanmean(sqrt((fu(:)).^2+(fv(:)).^2))
        end
        
        % Get cell division events:
        if UseCellDivDetect
            CenterDivisionsTemp = GetCellDivisions(im1,Template,net,MedFiltSize,BackThresh,FiltType,WSize);
            CenterDivisions(1:size(CenterDivisionsTemp,1),:,i) = CenterDivisionsTemp;
        end
        
        % Get cellular orientation:
        if UseOrientationAnalysis
            % Compute eigenvectors:
            EigInfo = coherence_orientation(double(im1),rho);
            % Largest eigen-vector:
            EigVec = EigInfo.w1;
            % Take only every 'MSize' value in all dimensions
            EigVec = EigVec(1:MSize:end,1:MSize:end,:);
            CellOrientation{i} = EigVec;
            
%                     % Make Plot:
%                     [images]=Reader(strcat(jpg_path2,'\',folderList(folder_number).name),[i+1]);
%                     %images = imhistmatch(double(images),im1,256);
%                     if ~isa(images,'uint8')
%                         images = im2uint8(images);
%                     end
%                     im2 = images;                    
%                     [~, im2] = BM3D(1, double(im2), 0.5*std(double(im2(:))));
%                     
%                     para.Step  = 35; %%% intensity of orientation
%                     para.scl   = 10; %%% length of orientation
%                     rho=15;
%                     EigInfo = coherence_orientation(double(im2),rho);
%                     ConvInfo.imconv = ones(size(im2));
%                     DisplayImage(im2,EigInfo,ConvInfo,para)
%                     %cd(folder)
%                     saveas(gcf, sprintf('Orientation Image #10 V1 01 im%d.png',i));
%                     close all
%                     cd('..');
%             
%                     [x,y] = meshgrid(1:size(im1,2),1:size(im1,1));
%                     figure;
%                     imshow(im1,[])
%                     hold on
%                     quiver(x(1:MSize:end,1:MSize:end),y(1:MSize:end,1:MSize:end),EigVec(1:1:end,1:1:end,2),EigVec(1:1:end,1:1:end,1))
%             
            
        end
        
        
        % Set new im1 as old im2:
        im1 = im2;
    end
    
    
    
    % save velocity field in results folder:
    if UsePIV
        cd(folderVelField)
        save('VelocityField.mat','VelField','-v7.3' )
        save('PositionsX.mat','x')
        save('PositionsY.mat','y')
        save('Settings.mat','met', 'wins', 'globtrld', 'loctrld','kernelsize','med','int','Overlap')
        cd(curr_directory);
    end
    
    % save cell division events and create video if wanted:
    if UseCellDivDetect
        % Multiobject tracking assigning positions of cell divisions to tracks:
        CellProps = MultiObjectTracking(CenterDivisions,MaxTimeDiff);
        % add another info about the starting point of each division:
        StartPos = NaN(length(CellProps),1);
        for i = 1:length(CellProps)
            StartPos(i) = find(isnan(CellProps(i).centroid(:,1)) == 0,1);
        end
        
        % % Cumulative division events:
        % DivTimes= unique(StartPos);
        % NumDivs = NaN(length(DivTimes),1);
        % for i = 1:length(DivTimes)
        %     NumDivs(i) =  numel(find(StartPos == DivTimes(i)));
        % end
        % NumDivs = cumsum(NumDivs);
        %figure; plot(DivTimes,NumDivs)
        
        % go to save folder:
        cd(folderDivision)
        save('CellDivisions.mat','CellProps','-v7.3')
        save('Settings.mat','CType','NetName', 'MedFiltSize','FiltType',...
            'BackThresh','WSize','MaxTimeDiff')
        if MakeVideo
            MakeDivisionVideo(CellProps,StartPos,jpg_path2,folderList,folder_number,m)
        end
        % go back to top level folder:
        cd(curr_directory);
    end
    
    if UseOrientationAnalysis
        % save orientation analysis:
        cd(folderOrientation)
        save('CellOrientation.mat','CellOrientation','-v7.3' )
        save('Settings.mat','rho','MSize')
        cd(curr_directory);
    end
    
    % Removes all unnescessary variables:
    clearvars -EXCEPT folder_number jpg_path2 folderList met wins ...
        globtrld loctrld kernelsize med int x y Overlap net Template ...
        MSize rho MakeVideo MedFiltSize FiltType BackThresh WSize ...
        MaxTimeDiff UseOrientationAnalysis UsePIV UseCellDivDetect ...
        curr_directory CType NetName
end