clear all;
mPath = mfilename('fullpath');
Idx = max(strfind(mPath,filesep));
mPath = mPath(1:Idx);
addpath(mPath)

% add path for BM3D filter:
addpath(strcat(mPath,filesep,'BM3D Filter'))
% add path of the circular statistics toolbox:
addpath(strcat(mPath,filesep,'CircStat'))
% add path for PIVlab:
addpath(strcat(mPath,filesep,'PIVlab Analysis'))
% add path for cell division detection:
addpath(strcat(mPath,filesep,'Cell Division Analysis'))
% add path for orientation analysis:
addpath(strcat(mPath,filesep,'Orientation Analysis'))

% Load parameters and check validity:
[UsePIV,s,p,r,UseCellDivDetect,CType,NetName,NetNameSeg,MinDivArea, ...
    BlockSize,MakeVideo,MedFiltSize,FiltType,BackThresh,WSize,MaxTimeDiff,...
    UseOrientationAnalysis,MSize,rho,UseBM3D,UseCLAHE,UseDriftCorrect,...
    WinSizeDrift] = ParameterFunctionMainPIVlab;
CheckForValidInputsPIVlab(UsePIV,s,p,r,UseCellDivDetect,...
    CType,NetName,NetNameSeg,MinDivArea,BlockSize,MakeVideo,MedFiltSize,...
    FiltType,BackThresh,WSize,MaxTimeDiff,UseOrientationAnalysis,MSize,...
    rho,UseBM3D,UseDriftCorrect,WinSizeDrift)

% Load templates used for cell division detection and the ANN:
if UseCellDivDetect == 1
    % load templates used in cross correlation analysis:
    cd(strcat(mPath,filesep,'Cell Division Analysis',filesep,'DivisionTemplate',filesep,'Templates',filesep,CType))
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
    cd(strcat(mPath,filesep,'Cell Division Analysis',filesep,'TrainedANNs'))
    % As the name of the variable is not clear, use field names, as the
    % model is the only constituent (I hope ...):
    net = load(strcat(NetName,'.mat'));
    VariableNames = fieldnames(net);
    net = net.(VariableNames{1});    
    
    % load the trained ANN for segmentation:
    cd(strcat(mPath,filesep,'Cell Division Analysis',filesep,'Trained Segmentation Networks'))
    % As the name of the variable is not clear, use field names, as the
    % model is the only constituent (I hope ...):
    SegNet = load(strcat(NetNameSeg,'.mat'));
    VariableNames = fieldnames(SegNet);
    SegNet = SegNet.(VariableNames{1});
    
    % Check if the input layer has the right size, as determined by the "WSize"
    % parameter:
    ImSize = net.Layers(1,1).InputSize;
    if (WSize * 2 +1) ~= ImSize(1)
        error(['Error in ParameterFunctionMain.m. \nInput of ''WSize'' not' ...
            ' correct. 2*WSize + 1 must equal the image size in the neuronal' ...
            ' network used for classification of cell division events.\n' ...
            'Currently: WSize == %g and ImSize = %g'],WSize,ImSize(1))
    end
end

% Get image path:
Startpath=mPath;
FolderPath=uigetdir(Startpath, 'Chose the folder with the images');
cd(FolderPath);   
FolderList = dir();
FolderList([FolderList.isdir]~=1) = [];

% get all relevant folders for saving outside of the parallelized loop:
folderVelField=cell(0);
folderDivision=cell(0);
folderOrientation=cell(0);
count = 1;
for FolderNumber=3:length(FolderList)
	folderVelField{count} = fullfile(strcat(FolderPath,filesep,FolderList(FolderNumber).name),'Results');
	folderDivision{count} = fullfile(strcat(FolderPath,filesep,FolderList(FolderNumber).name),'Results Cell Division Analysis');
    folderOrientation{count} = fullfile(strcat(FolderPath,filesep,FolderList(FolderNumber).name),'Results Orientation Analysis');
	count = count+1;
end

CurrDirectory = pwd;
% Go through all folders
parfor FolderNumber = 3:length(FolderList)
    cd(FolderList(FolderNumber).name);
    
    if UsePIV
        % Check if several folders actually exist:
        % folderVelField = fullfile(strcat(FolderPath,filesep,FolderList(FolderNumber).name),'Results');
        % Checks the existence of the folder with the name "Results". If it
        % does not exist it is created.
        if (exist(folderVelField{FolderNumber-2},'dir') == 0)
            mkdir('Results');
        end
    end
    
    if UseCellDivDetect
        % folderDivision = fullfile(strcat(FolderPath,filesep,FolderList(FolderNumber).name),'Results Cell Division Analysis');
        % Checks the existence of the folder with the name "Results Cell Division Analysis". If it
        % does not exist it is created.
        if (exist(folderDivision{FolderNumber-2},'dir') == 0)
            mkdir('Results Cell Division Analysis');
        end
    end
    
    if UseOrientationAnalysis
        % Check if several folders actually exist:
        % folderOrientation = fullfile(strcat(FolderPath,filesep,FolderList(FolderNumber).name),'Results Orientation Analysis');
        % Checks the existence of the folder with the name "Results Orientation Analysis". If it
        % does not exist it is created.
        if (exist(folderOrientation{FolderNumber-2},'dir') == 0)
            mkdir('Results Orientation Analysis');
        end
    end
    
    % set first image:
    [im1,~,m] = Reader(strcat(FolderPath,'\',FolderList(FolderNumber).name),1);
    if ~isa(im1,'uint8')
        im1 = im2uint8(im1);
    end
    
    if UseBM3D % denoise image
    % Denoising using the BM3D filter:
        [~, im1] = BM3D(1, double(im1), 0.5*std(double(im1(:))));
    else % If not denoising, just rescale and transform to double so the 
         % rest works as well
        im1 = double(im1)./max(double(im1(:)));
    end
    
    % Allocate variables
    % Velocity field - pre allocation of size:   
    if UsePIV
        [x, ~, ~, ~, ~,~] = ...
            piv_analysis(im1,im1,p,s,1,false);
        VelField = NaN(size(x,1),size(x,2),2,m);
        % Set speed at first image to zero:
        VelField(:,:,1,1) = 0;
        VelField(:,:,2,1) = 0;
	% Set drift speed:
        sDrift = s;
        sDrift{1,2} = WinSizeDrift;
        VelFieldDrift =NaN(2,m);
        VelFieldDrift(:,1) = 0;
    end
    % Centers of cell divisions:
    CenterDivisions = NaN(250,2,m);
    % Direction of cells
    CellOrientation = cell(m,1);
    for i=1:m-1
        
        % Percentage =((FolderNumber-3)*(m-1) + i)/((length(folderList)-2)*m-1)        
        [images]=Reader(strcat(FolderPath,'\',FolderList(FolderNumber).name),i+1);
        %images = imhistmatch(double(images),im1,256);
        if ~isa(images,'uint8')
            images = im2uint8(images);
        end
        im2 = images;
        
        if UseBM3D % denoise image
            [~, im2] = BM3D(1, double(im2), 0.5*std(double(im2(:))));
        else % Just rescale and transform to double so the rest works as well
            im2 = double(im2)./max(double(im2(:)));
        end
        % display current folder and image:
        sprintf('Folder: %s  \n Image: %05d', ...
                FolderList(FolderNumber).name, i+1)
            
        % PIV calculation
        if UsePIV            
            % PIV
            [x, y, u, v, typevector,correlation_map] = ...
            piv_analysis(im1,im2,p,s,1,false);
            % Post processing:
            [u_filt, v_filt,typevector_filt]= ...
            post_proc_wrapper(u,v,typevector,r,true);  
            % Set velocity field:
            VelField(:,:,1,i+1) = u_filt;
            VelField(:,:,2,i+1) = v_filt;      
	    
	    % Calculate drift:
            if UseDriftCorrect
                % Drift estimate
                [x, y, u, v, typevector,correlation_map] = ...
                    piv_analysis(im1,im2,p,sDrift,1,false);
                % Post processing:
                [u_filt_drift, v_filt_drift,typevector_filt]= ...
                    post_proc_wrapper(u,v,typevector,r,true);
                % Set drift speed:
                u_filt_drift = nanmean(u_filt_drift(:));
                v_filt_drift = nanmean(v_filt_drift(:));
                VelFieldDrift(1,i+1) = u_filt_drift;
                VelFieldDrift(2,i+1) = v_filt_drift;
            else
                VelFieldDrift(:,i+1) = 0;
            end

            % display current average speed in px/ whatever time between
            % sucessive images is ...
            sprintf('MeanSpeed: %2.3f', mean(sqrt((u_filt(:)).^2+(v_filt(:)).^2),'omitnan'))
        end
        
        % Get cell division events:
        if UseCellDivDetect
            tmp = Reader(strcat(FolderPath,'\',FolderList(FolderNumber).name),i);
            if ~isa(tmp,'uint8')
                tmp = im2uint8(tmp);
            end
            CenterDivisionsTemp = GetCellDivisionsSegmentation(tmp,Template,SegNet,net,MinDivArea,BlockSize,UseCLAHE,MedFiltSize,BackThresh,FiltType,WSize);       
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
        end        
        % Set new im1 as old im2:
        im1 = im2;
    end
    
    % save velocity field in results folder:
    if UsePIV
        cd(folderVelField{FolderNumber-2})
        SavePIVlab(VelField,VelFieldDrift,x,y,p,s,r,WinSizeDrift)
        cd(CurrDirectory);
    end
    
    % save cell division events and create video if wanted:
    if UseCellDivDetect
        % Multiobject tracking assigning positions of cell divisions to tracks:
        CellProps = MultiObjectTracking(CenterDivisions,MaxTimeDiff);
        % starting point of each division event:
        StartPos = NaN(length(CellProps),1);
        for i = 1:length(CellProps)
            StartPos(i) = find(isnan(CellProps(i).centroid(:,1)) == 0,1);
        end

        % go to save folder and save data:
        cd(folderDivision{FolderNumber-2})
        SaveDivision(CellProps,CType,NetName,NetNameSeg,MedFiltSize,FiltType,...
        BackThresh,WSize,MaxTimeDiff,BlockSize,MinDivArea)
        if MakeVideo           
            MakeDivisionVideo(CellProps,StartPos,FolderPath,FolderList,FolderNumber,m)
        end
        % go back to top level folder:
        cd(CurrDirectory);
    end
    
    if UseOrientationAnalysis
        % save orientation analysis:
        cd(folderOrientation{FolderNumber-2})
        SaveOrientation(CellOrientation,rho,MSize)
        cd(CurrDirectory);
    end
end
% Remove PIVlab path as it contains several changes to default-functions 
% such as nanstd, nanmean, etc.
rmpath(strcat(mPath,filesep,'PIVlab Analysis'));
