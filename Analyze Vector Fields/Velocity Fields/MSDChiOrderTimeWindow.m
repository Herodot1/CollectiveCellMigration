function [MSDTemp, QTemp, ChiTemp,NumNewNeighbors,TrackMatTemp] = MSDChiOrderTimeWindow(VelField,WSize,dt,CSize,OvThresh,im_size,POIx,POIy,SubPxResolution,Sampling,CenterSpeed)
% Get time dependent MSD, order parameters and 4 point susceptibility and
% neighbor statistics:
Trajectory = [];
NumObj = [];
QTemp = NaN(WSize,size(VelField,4));
MSDTemp = NaN(WSize,size(VelField,4));
ChiTemp = NaN(WSize,size(VelField,4));
NumNewNeighbors = NaN(length(1+(WSize-1)/2:Sampling:size(VelField,4)-1-((WSize-1)/2)),1);
TrackMatTemp = cell(1,numel(1+(WSize-1)/2:Sampling:size(VelField,4)-1-((WSize-1)/2)));
count = 0;
for k = 1+(WSize-1)/2:Sampling:size(VelField,4)-1-((WSize-1)/2)
    count = count + 1;
    idx = k - (WSize-1)/2;    
    %Percentage = (k-(WSize-1)/2)/(size(VelField,4)-1-(WSize-1)/2)    
    interval = k-(WSize-1)/2:k+(WSize-1)/2;   
    
    [TrackMat] = GenTrackMat(VelField(:,:,:,interval),im_size,POIx,POIy,SubPxResolution,CenterSpeed);
    TrackMatTemp{count} = TrackMat;
    
    % Do the nearest neighbor calculations:
    [Neighbors] = NearestNeighbors(TrackMat);
    NumNewNeighbors(idx) = Neighbors;
    
    % Allocate variables:
    RelDistTemp = NaN(WSize,size(TrackMat,1)*size(TrackMat,2));
    MSD = NaN(WSize,size(TrackMat,1)*size(TrackMat,2));
    QTempVar = [];
    for i = 1:size(TrackMat,1)
        for j = 1:size(TrackMat,2)
            % Generate trajecteroy(x,y,t)            
            Trajectory(:,1) = TrackMat(i,j,1,:);
            Trajectory(:,2) = TrackMat(i,j,2,:);
            Trajectory(:,3) = [0:WSize-1]*dt;
            % Get points were objects move out of image (<=0) or detection has
            % been abandoned (=NaN)
            ende = find(isnan(Trajectory) == 1,1);
            ende2 = find(Trajectory(:,1) <= 0,1);
            ende3 = find(Trajectory(:,2) <= 0,1);
            ende = min([ende;ende2;ende3;WSize]);
            % Calculate relative distance for Order Parameter and MSD:
            pos = size(TrackMat,2)*(i-1) + j;                     
            if ende > 2
                [temp,tau] = Kehl(Trajectory(1:ende-1,:));  
                if ~isempty(temp)
                    MSD(1:ende-2,pos) = temp;
                end
                [RelDistTemp(1:ende-1,pos),~,DD,DT] = RelativeDistance(Trajectory(1:ende,:));                
            end
            
            % Get order parameter for each object and timepoint and time difference:
            pos1 = find(DD>OvThresh*CSize);
            pos2 = find(DD<=OvThresh*CSize);
            DD(pos1) = 0;
            DD(pos2) = 1;
            DD(DT == 0) = NaN;
            QTempVar(1:size(DD,1),1:size(DD,2),i,j) = DD;
        end
    end
    
    % Get MSD:
    MSDTemp(:,idx) = nanmean(MSD,2);
    
    % Number of Objects:
    NumObj = [];
    for i = 1:size(RelDistTemp,1)
        NumObj(i) = size(RelDistTemp,2)-sum(isnan(RelDistTemp(i,:)));
    end    
    % Get order parameter Q:
    Q = squeeze(nansum(nansum(QTempVar,4),3));
    % Normalize QSquared:
    Q = Q./NumObj(1:end)';
    % Make time average:
    for l = 1:size(Q,1)
        TempVar = circshift(Q(l,:),-l,2);
        Q(l,:) = TempVar;
    end
    Q(Q == 0) = NaN;
    Q = nanmean(Q,1);
    QTemp(1:length(Q),idx) = Q;
    
    
    % Get QSquared:
    QSquared = squeeze(nansum(nansum(QTempVar,4),3));
    % Normalize QSquared and actually square it:
    QSquared = QSquared./NumObj(1:end)';
    QSquared = QSquared.^2;
    % Make time average:
    for l = 1:size(QSquared,1)
        TempVar = circshift(QSquared(l,:),-l,2);
        QSquared(l,:) = TempVar;
    end
    QSquared(QSquared == 0) = NaN;
    QSquared = nanmean(QSquared,1);
    % Set up 4 point susceptibility
    ChiTemp(1:length(Q),idx) = (size(TrackMat,1)*size(TrackMat,2))*(QSquared-Q.^2);
    
end

% % If sampling is used interpolate the data to fit the actual values,
% % using a linear 1d interpolation:
% if Sampling ~= 1
%     t = [0:size(MSDTemp,2)-1]*dt;
%     for i = 1:WSize-2
%         MSDTemp(i,:) = interp1( t(~isnan(MSDTemp(i,:))), MSDTemp(i,~isnan(MSDTemp(i,:))) , t, 'linear');
%         ChiTemp(i,:) = interp1( t(~isnan(ChiTemp(i,:))), ChiTemp(i,~isnan(ChiTemp(i,:))) , t, 'linear');
%         QTemp(i,:) = interp1( t(~isnan(QTemp(i,:))), QTemp(i,~isnan(QTemp(i,:))) , t, 'linear');
%     end
%     QTemp(WSize-1,:) = interp1( t(~isnan(QTemp(WSize-1,:))), QTemp(WSize-1,~isnan(QTemp(WSize-1,:))) , t, 'linear');
% end