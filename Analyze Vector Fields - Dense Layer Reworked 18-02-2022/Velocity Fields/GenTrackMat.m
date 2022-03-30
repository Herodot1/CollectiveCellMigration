function [TrackMat] = GenTrackMat(VelField,im_size,POIx,POIy,SubPxResolution,CenterSpeed) %[TrackMat,VelField] = GenTrackMat(fu,fv,im_size,POIx,POIy)
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
    
%     % Set velocity field:
%     VelField(:,:,1,i+1) = VelInterpX;
%     VelField(:,:,2,i+1) = VelInterpY;
    
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
    
%     % Create figure showing tracks:
%     fig1 = figure('name','fig1');
%     imshow(im2)
%     hold all
%     % Plots the detected boundaries.
%     for o = 1:i
%         
%         PosX2 = (TrackMat(:,:,2,o+1));
%         PosX2 = PosX2(:);
%         PosY2 = (TrackMat(:,:,1,o+1));
%         PosY2 = PosY2(:);
%         
%         PosX1 = round(TrackMat(:,:,2,o));
%         PosX1 = PosX1(:);
%         PosY1 = round(TrackMat(:,:,1,o));
%         PosY1 = PosY1(:);
%         
%         Pos1  = [PosX2,PosX1]';
%         Pos2  = [PosY2,PosY1]';
%         
%         % Plots the travelled path.
%         line(Pos1, Pos2,'LineWidth',1.5,'Color','r')
%         if o == i
%             plot(PosX2, PosY2,'r.', 'LineWidth', 2, 'MarkerSize', 15)
%         end
%         
%         %line([stats{start2(k,i) +n-1}(k,1).Centroid(1),stats{start2(k,i)+n}(k,1).Centroid(1)], ...
%         %    [stats{start2(k,i) +n -1}(k,1).Centroid(2),stats{start2(k,i)+n}(k,1).Centroid(2)],'LineWidth',1.5,'Color',colormap(k,:));
%         
%     end
%     
%     hold off
%     
%     % Gives the saved images names according to the image number
%     saveas(gcf,sprintf('img%04d.tif',i));
%     %close curr_frame;
%     close fig1;    
%     
%     % Make an image with original and heatmap image side by side
%     % Interpolate data for the heatmap to the rest of the image
%     var = (fu.^2 + fv.^2).^(0.5);
%     test = [];
%     [y_new,x_new] = meshgrid(1:size(images,1),1:size(images,2));
%     test = interp2(x,y,var,x_new,y_new,'linear');
%     test = test';
%     I = double(test);
%     %I = exp(I);
%     %I = I./max(I(:));
%     I(I<0) = 0;
%     %figure;
%     %imshow(I)
%     test = uint8(255*(test - min(test(:)))./(max(test(:))-min(test(:))));
%     colorMap = jet(256);
%     test = ind2rgb(test,colorMap);
    
    
end