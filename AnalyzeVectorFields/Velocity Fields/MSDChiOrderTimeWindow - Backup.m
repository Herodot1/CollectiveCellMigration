function [MSDTemp, QTemp, ChiTemp] = MSDChiOrderTimeWindow(TrackMat,WSize,dt,CSize,OvThresh)

% Get time dependent order parameters and 4 point susceptibility:
Trajectory = [];
NumObj = [];
QTemp = NaN(WSize,size(TrackMat,4));
MSDTemp = NaN(WSize,size(TrackMat,4));
ChiTemp = NaN(WSize,size(TrackMat,4));
for k = 1+(WSize-1)/2:size(TrackMat,4)-1-((WSize-1)/2)
    interval = k-(WSize-1)/2:k+(WSize-1)/2;
    RelDistTemp = NaN(WSize,size(TrackMat,1)*size(TrackMat,2));
    MSD = NaN(WSize,size(TrackMat,1)*size(TrackMat,2));
    for i = 1:size(TrackMat,1)
        for j = 1:size(TrackMat,2)
            % Generate trajecteroy(x,y,t)
            Trajectory(:,1) = TrackMat(i,j,1,interval);
            Trajectory(:,2) = TrackMat(i,j,2,interval);
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
                RelDistTemp(1:ende-1,pos) = RelativeDistance(Trajectory(1:ende,:));
            end
        end
    end
    
    % Get MSD:
    MSDTemp(:,k) = nanmean(MSD,2);
    % Set up order parameter:
    pos1 = find(RelDistTemp>OvThresh*CSize);
    pos2 = find(RelDistTemp<=OvThresh*CSize);
    TempVar = RelDistTemp;
    TempVar(pos1) = 0;
    TempVar(pos2) = 1;
    for i = 1:size(TempVar,1)
        NumObj(i) = size(TempVar,2)-sum(isnan(TempVar(i,:)));
    end
    TempVar = nansum(TempVar,2)./NumObj';
    QTemp(:,k) = TempVar;
    % Set up 4 point susceptibility
    Trajectory = [];
    Trajectory(:,1) = QTemp(:,k);
    Trajectory(:,2) = [0:WSize-1]*dt;
    var1 = RelativeDistance(Trajectory(1:end-1,:));
    var2 = Kehl(Trajectory(1:end-1,:));
    ChiTemp(1:length(var1),k) = (size(TrackMat,1)*size(TrackMat,2))*(var2-var1.^2);
end
