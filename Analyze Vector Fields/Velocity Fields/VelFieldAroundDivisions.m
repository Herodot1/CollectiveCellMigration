function [MeanSpeedDiv,MeanSpeedNonDiv] = VelFieldAroundDivisions(CellProps,VelField,x,y,StartPos,EndPos,Idx)
% Calculate local speed around a determine division center:

% Distance in pixels of the velocity field from the center of the division 
% event that is taken to calculate the speed around the division.
Neighborhood = 1;
% Offset:
% Modify StartPos and EndPos to account for the fact that only
% mitotically rounded cells are detected and the main effect on
% local speed is expected to occur before or after the cell
% attained said shape:
Offset = 7;

% get possible x and y coordinates:
xPos = unique(x);
yPos = unique(y);

% Modify StartPos and EndPos to account for the fact that only
% mitotically rounded cells are detected and the main effect on
% local speed is expected to occur before or after the cell
% attained said shape:
StartPos = StartPos-Offset;
EndPos = EndPos+Offset;
StartPos(StartPos<1) = 1;
EndPos(EndPos>size(VelField,4)) = size(VelField,4);
% Create list of indices containing division events present in
% each single image:
DivList = cell(size(VelField,4),1);
for i = 1:size(VelField,4)
    for j = 1:length(Idx)
        if ismember(i,[StartPos(j):EndPos(j)])
            DivList{i} = [DivList{i},j];
        end
    end
end
MeanSpeedNonDiv = NaN(size(VelField,4),1);
MeanSpeedDiv = NaN(size(VelField,4),1);
for i = 2:size(VelField,4) % start at i = 2 becaus velocities are set to 0 at first image (no difference can be taken).
    % Get absolute speed values:
    % Correct velocity field for drifts (because motion is supposed to have zero mean):
    VelX = VelField(:,:,1,i);
    VelX = imgaussfilt(VelX,2.71);
    VelY = VelField(:,:,2,i);
    VelY = imgaussfilt(VelY,2.71);
    %if CenterSpeed == 1
    VelX = VelX - nanmean(VelX(:));
    %end
    %if CenterSpeed == 1
    VelY = VelY - nanmean(VelY(:));
    %end
    Speed = sqrt(VelX.^2+VelY.^2);
    
    % Get speed values for those pixels that are closest to the
    % division event +- Neighborhood in each direction:
    DivEvents = DivList{i};
    DivCoordinates = NaN(length(DivEvents),2);
    Coords = NaN(length(DivEvents),2);
    CoordsPlus = Coords;
    CoordsMinus = Coords;
    MeanSpeedDivTmp = NaN(length(DivEvents),1);
    if ~isempty(DivEvents)
        for j = 1:length(DivEvents)
            if ~isnan(CellProps(Idx(DivEvents(j))).centroid(i,1))
                DivCoordinates(j,1:2) = CellProps(Idx(DivEvents(j))).centroid(i,:);
            else % if not detected take last/next identified detection as approximation:
                tmp = find(~isnan(CellProps(Idx(DivEvents(j))).centroid(:,1)) == 1);
                if length(tmp)>=3
                    % 1D interpolation
                    tmpX = interp1(tmp,CellProps(Idx(DivEvents(j))).centroid(tmp,1),i,'linear','extrap');
                    tmpY = interp1(tmp,CellProps(Idx(DivEvents(j))).centroid(tmp,2),i,'linear','extrap');
                    DivCoordinates(j,1:2) = [tmpX,tmpY];
                else
                    % nearest neighbor (temporal sense):
                    [~,tmpIdx] = min(abs(i-tmp));
                    DivCoordinates(j,1:2) = CellProps(Idx(DivEvents(j))).centroid(tmp(tmpIdx),:);
                end
            end
            % Closest velocity field data point:
            % [~,MinIdxX] = min(abs(xPos-CellProps(Idx(DivEvents(j))).centroid(i,1)));
            % [~,MinIdxY] = min(abs(yPos-CellProps(Idx(DivEvents(j))).centroid(i,2)));
            [~,MinIdxX] = min(abs(xPos-DivCoordinates(j,1)));
            [~,MinIdxY] = min(abs(yPos-DivCoordinates(j,2)));
            % Merged Coordinates:
            Coords(j,1) = MinIdxX;
            Coords(j,2) = MinIdxY;
            % Get expanded coordinates around the specific
            % event:            
            CoordsPlus(j,:) = Coords(j,:) + Neighborhood;
            CoordsMinus(j,:) = Coords(j,:) - Neighborhood;
            % Check if anything is out of bounds, if so set it
            % to bound value:
            if CoordsPlus(j,1) > length(xPos)
                CoordsPlus(j,1) = length(xPos);
            end
            if CoordsPlus(j,2) > length(yPos)
                CoordsPlus(j,2) = length(yPos);
            end
            if CoordsMinus(j,1) < 1
                CoordsMinus(j,1) = 1;
            end
            if CoordsMinus(j,2) <1
                CoordsMinus(j,2) = 1;
            end
            % Get mean velocity around divisions:
            tmp = Speed(CoordsMinus(j,2):CoordsPlus(j,2),CoordsMinus(j,1):CoordsPlus(j,1));
            MeanSpeedDivTmp(j) = nanmean(tmp(:));
        end        
        MeanSpeedDiv(i) = nanmean(MeanSpeedDivTmp);
        
        % Get mean speed of the remaining velocity field
        % (without divisions):
        tmpSpeed = Speed;
        for j = 1:length(DivEvents)
            tmpSpeed(CoordsMinus(j,2):CoordsPlus(j,2),CoordsMinus(j,1):CoordsPlus(j,1)) = NaN;
        end
        MeanSpeedNonDiv(i) = nanmean(tmpSpeed(:));
        %nanmean(tmpSpeed(:));
        %nanmean(MeanSpeedDivTmp(:));
        %nanmean(tmpSpeed(:))./nanmean(MeanSpeedDivTmp(:));
    end
    
end
