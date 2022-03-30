function [AutoCorrAll,SpeedAll,RMSVelAll,AreaAll,AngVelAll,MeanVelX,MeanVelY] = AnalyzeVelField(VelField,im_size,m,POIx,POIy,CenterSpeed)
% Analyzes the velocity field generated via Mat PIV. 

% Set meshgrid:
[X,Y] = (meshgrid(1:im_size(2)/size(VelField,2):im_size(2),1:im_size(1)/size(VelField,1):im_size(1)));
[Xq,Yq] = meshgrid(1:im_size(2),1:im_size(1));
for i = 1:m
    % Get absolute speed values:
    % Correct velocity field for drifts (because motion is supposed to have zero mean):
    VelX = VelField(:,:,1,i);
    VelX = imgaussfilt(VelX,0.71);
    VelY = VelField(:,:,2,i);
    VelY = imgaussfilt(VelY,0.71);
    MeanVelX(i) = nanmean(VelX(:));
    MeanVelY(i) = nanmean(VelY(:));
    if CenterSpeed == 1
        VelX = VelX - nanmean(VelX(:));
    end
    
    if CenterSpeed == 1
        VelY = VelY - nanmean(VelY(:));
    end
   
    Speed = sqrt(VelX.^2+VelY.^2);
    SpeedAll(i) = nanmean(Speed(:));
    % Root mean squared speed:
    RMSVelAll(i) = rms(Speed(:));
    % Calculate heterogeneity parameter (Angelini et al, 2011):
    % Get all elements of the positions of cells belonging to the
    % fastest 20% of all cells:
    % Use cell size defined by POIx and POIy:
    VelInterpX = interp2(X,Y,squeeze(VelX),POIx,POIy,'spline');
    VelInterpY = interp2(X,Y,squeeze(VelY),POIx,POIy,'spline');
    q20 = quantile(Speed(:),0.8);
    pos = find(Speed(:)>q20);
    FastCells = false(size(Speed,1),size(Speed,2));
    FastCells(pos) = true;
    % Connected components:
    CC = bwconncomp(FastCells);
    Area = regionprops(FastCells, 'Area');
    Area = cat(1,Area.Area);
    % Remove Single Pixels:
    % Area(Area<2) = [];
    % Transfer to main variable:
    AreaAll(i) = nanmean(Area);
    % Get angular velocity:
    [curlz,cav]= curl(X,Y,VelX,VelY);
    AngVelAll(i) = nanmean(cav(:));
    
    % Get velocity autocorrelation:
    Speed = sqrt(VelX.^2+VelY.^2);
    %D = autocorr2d(Speed-nanmean(Speed(:)));
    % Alternative Autocorrelationsberechnung:
    D = xcorr2(Speed-nanmean(Speed(:)));
    if max(D(:)) ~= 0
        D = D./max(D(:));
    end
    [r,c] =find(D == max(D(:)),1);
    y = 1:size(D,1);
    x = 1:size(D,2);
    Dist = round(1*sqrt(bsxfun(@plus, (y.' - r) .^ 2, (x - c) .^ 2)))/1;
    count = 0;
    for k = 0:max(Dist(:))
        pos = find(Dist == k);
        %test(3,k+1) = nanmean(D(pos));
        AutoCorrAll(count+1,i) = nanmean(D(pos));
        count = count+1;
    end
end

%     % Is slower :-(
%     Speed = (squeeze((VelField(:,:,1,:).^2+VelField(:,:,2,:).^2))).^0.5;
%     Speed = reshape(Speed,[],size(Speed,3));
%     Speed = nanmean(Speed,1)';
%     SpeedAll(1:length(Speed),folder_number-2) = Speed;
