function [tracks, rMean, IDStart]= AnalyzeCellDivisions(centers,radii)
% Uses multiple object tracking to match potential cell divisions to true
% cell divisions, based on multiple object tracking:

m = size(centers,3);

% Object tracking:
tracks = MultiObjectTracking(centers(:,:,:),radii(:,:));
% Allocate Variables:
IDList = [];
IDStart = [];
for j = 1:length(tracks)
    ID(j) = tracks(j).id;
end

% Get SAtarting positions of each cell division:
count = 0;
for i = 1:m-1
    centersTemp = [];
    % Get cell centers at current time point:
    for j = 1:length(tracks)
        centersTemp(j,1:2) = tracks(j).centroid(i,:);
    end    
    % Get for each cell when the division started:
    for j = 1:length(tracks)
        if ~isnan(centersTemp(j,1)) && sum(ismember(IDList,tracks(j).id)) == 0
            count = count+1;
            IDList(count) = tracks(j).id;
            IDStart(count) = i;
        end
    end    
end

% Average radius of the s-phase cell:
rMean = NaN(length(tracks),1);
for j = 1:length(tracks)
    rMean(j) = nanmean(tracks(j).radii);
end