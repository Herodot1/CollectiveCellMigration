function [TrackMat] = GenTrackMat(VelField,im_size,POIx,POIy,SubPxResolution,CenterSpeed)
% Creates tracking matrix for defined spots.
% Input: 
% VelField Velocity field of size MxNx2xT. M,N spatial coordinates,
% dimension3 = x,y component of velocity. T time component
% POIx,POIy = meshgrid of points of interst, created via the spacing
% defined by s_size in ParameterFunctionMain.m
% remaining parameters taken from ParameterFunctionMain.m
% Output:
% TrackMat = Tracking matrix of individual points as a function of time,
% ordered the same way as the velocity field but instead of local
% velocities it gives the position of the points defined by TrackMat(:,:,:,1)


% SubPxResolution = value between 0 and 1 giving the sub-pixel resolution
SubPxResolution = 1/SubPxResolution;

% Tracking matrix
m = size(VelField);
m = m(end);
TrackMat = NaN(size(POIx,1),size(POIx,2),2,m);
TrackMat(:,:,1,1) = POIx;
TrackMat(:,:,2,1) = POIy;
% Velocity field:
% VelField = NaN(im_size(1),im_size(2),2,size(fu,3));
% VelField(:,:,1,1) = zeros(im_size(1), im_size(2));
% VelField(:,:,2,1) = zeros(im_size(1), im_size(2));

% Create Meshgrid:
[X,Y] = meshgrid(1:im_size(2)/size(VelField,2):im_size(2),1:im_size(1)/size(VelField,1):im_size(1));
[Xq,Yq] = meshgrid(1:1/SubPxResolution:im_size(2),1:1/SubPxResolution:im_size(1));
for i=1:m-1
    
    %     Percentage  = i/size(fu,3)
    %     if rem(round(Percentage*100),10) == 0
    %        %disp(sprintf('%d%% done',Percentage))
    %     end
    
    % Correct velocity field for drifts (because motion is supposed to have zero mean):
    VelX = VelField(:,:,1,i+1);
    VelX = imgaussfilt(VelX,0.71);
    if CenterSpeed == 1
        VelX = VelX - nanmean(VelX(:));
    end
    VelY = VelField(:,:,2,i+1);
    VelY = imgaussfilt(VelY,0.71);
    if CenterSpeed == 1
        VelY = VelY - nanmean(VelY(:));
    end
    VelInterpX = interp2(gpuArray(X),gpuArray(Y),gpuArray(squeeze(VelX)),gpuArray(Xq),gpuArray(Yq),'linear');
    VelInterpY = interp2(gpuArray(X),gpuArray(Y),gpuArray(squeeze(VelY)),gpuArray(Xq),gpuArray(Yq),'linear');
  
    % VelInterpX = interp2(X,Y,squeeze(VelField(:,:,1,i+1)),Xq,Yq,'spline');
    % VelInterpY = interp2(X,Y,squeeze(VelField(:,:,2,i+1)),Xq,Yq,'spline');
    % Ab version 2018a
    % VelInterpX = interp2(X,Y,fu(:,:,i),Xq,Yq,'makima');
    % VelInterpY = interp2(X,Y,fv(:,:,i),Xq,Yq,'makima');
    

    % Get linear indices: (round to 0.x)
    PosY = round(SubPxResolution*TrackMat(:,:,1,i))/SubPxResolution;
    PosY = PosY(:);    
    PosX = round(SubPxResolution*TrackMat(:,:,2,i))/SubPxResolution;
    PosX = PosX(:);
    
    % Set spots that moved out of the image to x=y=1:
    PosX(PosY<1) = 1;
    PosY(PosY<1) = 1;
    PosY(PosX<1) = 1;
    PosX(PosX<1) = 1;
    
    PosX(PosY>im_size(1)) = 1;
    PosY(PosY>im_size(1)) = 1;
    PosY(PosX>im_size(2)) = 1;
    PosX(PosX>im_size(2)) = 1;
     
    % Adjust rounding errors from the interp step:
    PosX(SubPxResolution*PosY>size(Xq,1)) = 1;
    PosY(SubPxResolution*PosY>size(Xq,1)) = 1;
    PosY(SubPxResolution*PosX>size(Xq,2)) = 1;
    PosX(SubPxResolution*PosX>size(Xq,2)) = 1;
    
    idx = find((PosY) == 1);
    
    LinIdx = sub2ind(size(VelInterpX), SubPxResolution*PosY,SubPxResolution*PosX);
    % LinIdx(isnan(LinIdx) == 1) = [];
    % Actualize tracking matrix:
    if sum(isnan(VelInterpY(LinIdx))) ~=0 || sum(isnan(VelInterpX(LinIdx))) ~=0
        VelInterpY(isnan(VelInterpY) == 1) = 0;
        VelInterpX(isnan(VelInterpX) == 1) = 0;
    end

    % Actualize the tracking matrix and transform gpuarray into regular
    % one:
    VelX = gather(reshape(VelInterpX(LinIdx),size(TrackMat(:,:,1,i))));
    VelY = gather(reshape(VelInterpY(LinIdx),size(TrackMat(:,:,1,i))));
    
    % if idx is not empty (NaNs were removed, objects out of bound) set
    % velocity to zero:
    if ~isempty(idx)
        % find the assigened values and replace them:
        VelX(VelX == VelInterpX(1)) = 0;
        VelY(VelY == VelInterpY(1)) = 0;
    end

    TrackMat(:,:,1,i+1) = round(SubPxResolution*TrackMat(:,:,1,i))/SubPxResolution + VelY;
    TrackMat(:,:,2,i+1) = round(SubPxResolution*TrackMat(:,:,2,i))/SubPxResolution + VelX;    
end