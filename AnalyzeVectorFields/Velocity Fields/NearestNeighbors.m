function [NumNewNeighbors] = NearestNeighbors(TrackMat)
% Function that gives the relative amount of cells that made new neighbors
% during the measurements.
% Input: TrackMat = 4d matrix of x,y coordinates for multiple objects
% tracked over time.  dimension 1 and 2: initial locations of objects as
% defined via a meshgrid. dimension 3: y/x coordinate. dimension 3: time.  
% Chose 9 nearest neighbors to include the query point itself and its 8
% nearest neighbors. 
NumNeighborsKNN = 9;

%% Alterntive nearest neighbors estimate using a simple euclidean metric for estimation of new neighbors
% Use a simple euclidean metric to determine nearest neighbors
% find nearest neighbors in the beginning and end of TrackMat:
% Start with the beginning:
TempVarX = squeeze(TrackMat(2:end-2,2:end-2,2,1));
TempVarY = squeeze(TrackMat(2:end-2,2:end-2,1,1));
XStart = [TempVarX(:),TempVarY(:)];
[IdxStart,DStart] = knnsearch(XStart,XStart,'K',NumNeighborsKNN);

% Continue with track ends:
TempVarX = squeeze(TrackMat(2:end-2,2:end-2,2,end));
TempVarY = squeeze(TrackMat(2:end-2,2:end-2,1,end));
XEnd = [TempVarX(:),TempVarY(:)];
[IdxEnd,DEnd] = knnsearch(XEnd,XEnd,'K',NumNeighborsKNN);

% Get new neighbors:
for i = 1:size(IdxStart,1)
    TempVar = [IdxStart(i,:)',IdxEnd(i,:)'];
    TempVar = TempVar(:);
    [tmp,Idx] = sort(TempVar);
    idp = diff(tmp)>0;
    TempVar = TempVar(Idx([true;idp]&[idp;true]));
    if isempty(TempVar)
        NewNeighbors(i) = 0;
    else
        % Division by 2 is necessary as this approach assignes each changed
        % neigbor 2 new ones. The "old" object is found and the "new" one
        % substituting the "old" one, effectively doubling the number of
        % new neighbors. Thus, divide by two.
        NewNeighbors(i) = numel(TempVar)/2;
    end
end
NewNeighbors = NewNeighbors./(NumNeighborsKNN-1);

% Relative number of cells with new neighbors;
NumNewNeighbors = numel(find(NewNeighbors>0))./length(NewNeighbors);

% % % average number of new neighbors per cell:
% NewNeighborsPerCell = mean(NewNeighbors).*(NumNeighborsKNN-1);
% % % average number of new neighbors per cell only of those cells that made
% % % new neighbors:
% NewNeigborsPerDiffusiveCell = mean(NewNeighbors(NewNeighbors>0)).*(NumNeighborsKNN-1);