function [MSDCaged, CagedMSDTemp] = MSDCaged(TrackMat,dt)

TrackMat = TrackMat(2:end-1,2:end-1,:,:);
MSDCaged = NaN(size(TrackMat,4),(size(TrackMat,1)-1)*(size(TrackMat,2)-1),10);
Trajectory = [];
% Nearest Neighbors of each non-boundary point:
TempVarXAll = squeeze(TrackMat(1:end,1:end,2,1));
TempVarYAll = squeeze(TrackMat(1:end,1:end,1,1));
XStartAll = [TempVarXAll(:),TempVarYAll(:)];

   
countNeighbors = 0;
for NumNeighborsKNN = [1:10].^2%[1,9,25,49,81,121]
    countNeighbors = countNeighbors+1;
    % size constraints for usage of points in tracking matrix:
    Constraint = 1+ceil((sqrt(NumNeighborsKNN)-1)/2);
    
    % get constraint points of interest:
    TempVarXConstraint = squeeze(TrackMat(Constraint:end-Constraint,Constraint:end-Constraint,2,1));
    TempVarYConstraint = squeeze(TrackMat(Constraint:end-Constraint,Constraint:end-Constraint,1,1));
    XStartConstraint = [TempVarXConstraint(:),TempVarYConstraint(:)];
    % Nearest Neighbors of each point:
    [IdxStart,DStart] = knnsearch(XStartAll,XStartAll,'K',NumNeighborsKNN);
    % Identify points of interest:
    IdxPOI = find(ismember(XStartAll,XStartConstraint,'rows') == 1);
    % Get the associated coordinates in the tracking matrix:
    [rowPOI,colPOI] = ind2sub(size(squeeze(TrackMat(:,:,1,1))), IdxPOI);
    for j = 1:length(IdxPOI)
        % get rows and columns of relevant data. denote that the point of
        % interest is also included in the list (as it is one of its
        % nearest neighbors ...):
        [rowStart,colStart] = ind2sub(size(TempVarXAll), IdxStart(IdxPOI(j),1:end));
        % get all coordinates of neighbors:
        Coord = NaN(2,length(rowStart),size(TrackMat,4));
        for k = 1:length(rowStart)
            Coord(1,k,:) = squeeze(TrackMat(rowStart(k),colStart(k),1,:));
            Coord(2,k,:) = squeeze(TrackMat(rowStart(k),colStart(k),2,:));
        end
        % Generate trajecteroy(x,y,t) of cage by substracting the mean
        % position of the neighbors from the position of interest.
        % denote the substraction of the point of interest before
        % calculating the mean, as it is also included in the list of
        % neighbors.
        % the initial
        
        Trajectory(:,1) = squeeze(TrackMat(rowPOI(j),colPOI(j),1,:)) - (squeeze(sum(Coord(1,:,:),2)) - squeeze(TrackMat(rowPOI(j),colPOI(j),1,:)))./max([length(rowStart)-1,1]);                                                                               
        Trajectory(:,2) = squeeze(TrackMat(rowPOI(j),colPOI(j),2,:)) - (squeeze(sum(Coord(2,:,:),2)) - squeeze(TrackMat(rowPOI(j),colPOI(j),2,:)))./max([length(rowStart)-1,1]);
        Trajectory(:,3) = [0:size(TrackMat,4)-1]*dt;
        % Get points were objects move out of image (change ==0) or detection has
        % been abandoned (=NaN)
        ende = find(isnan(Trajectory) == 1,1);
        ende2 = find(diff(TrackMat(rowPOI(j),colPOI(j),1,:)) == 0,1);
        ende3 = find(diff(TrackMat(rowPOI(j),colPOI(j),2,:)) == 0,1);
        ende = min([ende;ende2;ende3;size(TrackMat,4)+1]);
        % Get MSD:
        [temp,tau] = Kehl(Trajectory(1:ende-1,:));
        %pos = size(TrackMat,2)*(i-1) + j;
        %MSDCaged(1:ende-2,pos,countNeighbors) = temp;
        MSDCaged(1:ende-2,j,countNeighbors) = temp;
    end
    MSDCaged(MSDCaged == 0) = NaN;
end
% Average over all virtual particles of one cage size:
 MSDCaged = nanmean(MSDCaged,2);
