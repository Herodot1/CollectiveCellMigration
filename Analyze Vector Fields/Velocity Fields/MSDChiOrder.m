function [MSD, Q, Chi,tau] = MSDChiOrder(TrackMat,dt,CSize,OvThresh)

% Calculate MSD:
MSD = NaN(size(TrackMat,4),size(TrackMat,1)*size(TrackMat,2));
RelDist = NaN(size(TrackMat,4),size(TrackMat,1)*size(TrackMat,2));
Trajectory = [];
QTemp = [];
for i = 1:size(TrackMat,1)
    for j = 1:size(TrackMat,2)
        % Generate trajecteroy(x,y,t)
        Trajectory(:,1) = TrackMat(i,j,1,:);
        Trajectory(:,2) = TrackMat(i,j,2,:);
        Trajectory(:,3) = [0:size(TrackMat,4)-1]*dt;
        % Get points were objects move out of image (<=0) or detection has
        % been abandoned (=NaN)
        ende = find(isnan(Trajectory) == 1,1);
        ende2 = find(Trajectory(:,1) <= 0,1);
        ende3 = find(Trajectory(:,2) <= 0,1);
        ende = min([ende;ende2;ende3;size(TrackMat,4)+1]);
        % Get MSD:
        [temp,tau] = Kehl(Trajectory(1:ende-1,:));
        pos = size(TrackMat,2)*(i-1) + j;
        MSD(1:ende-2,pos) = temp;
        % Calculate relative distance for Order Parameter:
        [RelDist(1:ende-2,pos),~,DD,DT] = RelativeDistance(Trajectory(1:ende-1,:));
        RelDist(RelDist == 0) = NaN;
        
        % Get squared order parameter for each object:
        pos1 = find(DD>OvThresh*CSize);
        pos2 = find(DD<=OvThresh*CSize);
        DD(pos1) = 0;
        DD(pos2) = 1;
        DD(DT == 0) = NaN;        
        QTemp(1:size(DD,1),1:size(DD,2),i,j) = DD;
    end
end
MSD(MSD == 0) = NaN;
RelDist(RelDist == 0) = NaN;

% Number of Objects:
NumObj = [];
for i = 1:size(RelDist,1)
    NumObj(i) = size(RelDist,2)-sum(isnan(RelDist(i,:)));
end

% Get order parameter Q:
Q = squeeze(nansum(nansum(QTemp,4),3));
% Normalize QSquared:
Q = Q./NumObj(1:end)';
% Make time average:
for k = 1:size(Q,1)
    TempVar = circshift(Q(k,:),-k,2);
    Q(k,:) = TempVar;
end
Q(Q == 0) = NaN;
Q = nanmean(Q,1);

% Get QSquared:
QSquared = squeeze(nansum(nansum(QTemp,4),3));
% Normalize QSquared and actually square it:
QSquared = QSquared./NumObj(1:end)';
QSquared = QSquared.^2;
% Make time average:
for k = 1:size(QSquared,1)
    TempVar = circshift(QSquared(k,:),-k,2);
    QSquared(k,:) = TempVar;
end
QSquared(QSquared == 0) = NaN;
QSquared = nanmean(QSquared,1);

% Get four point susceptibility Chi:
Chi = (size(TrackMat,1)*size(TrackMat,2))*(QSquared-Q.^2);