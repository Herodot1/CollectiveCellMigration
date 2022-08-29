function tracks = MultiObjectTracking(Center,MaxTimeDiff)
% Matches multiple detections to tracks:
% Input:
% Center = center coordinates (output of the "GetCellDivisions" function)
% MaxTimeDiff = maximal time an object remains invisible before the track
% is considered to be terminated.
% Output:
% tracks = ordered tracks matched to the detections

tracks = initializeTracks(); % Create an empty array of tracks.
nextId = 1; % ID of the next track
predictedCentroid  = []; % Assign list of predicted Centroid
for j = 1:size(Center,3)  
    centroids       = Center(:,:,j);  
    pos_in = find(isnan(centroids(:,1)) == 1,1);
    if ~isempty(pos_in)        
        centroids = centroids(1:pos_in-1,:);        
    end

    predictNewLocationsOfTracks();
    [assignments, unassignedTracks, unassignedDetections] = ...
        detectionToTrackAssignment();
    updateAssignedTracks();
    updateUnassignedTracks();
    %deleteLostTracks();
    createNewTracks();
    %displayTrackingResults();
end
deleteLostTracks();

function tracks = initializeTracks()
% create an empty array of tracks
tracks = struct(...
    'id', {}, ...
    'kalmanFilter', {}, ...
    'age', {}, ...
    'totalVisibleCount', {}, ...
    'consecutiveInvisibleCount', {}, ...
    'centroid',{}, ...
    'centroidPredicted',{}, ...
    'TrackTerminated', {});
end

% Use the Kalman filter to predict the centroid of each track in the 
% current frame, and update its bounding box accordingly.
function predictNewLocationsOfTracks()
for i = 1:length(tracks)
        
    % Predict the current location of the track.
    predictedCentroid = predict(tracks(i).kalmanFilter);       
end
end

% Assign detections to tracks:
function [assignments, unassignedTracks, unassignedDetections] = ...
    detectionToTrackAssignment()

nTracks = length(tracks);
nDetections = size(centroids, 1);

% Compute the cost of assigning each detection to each track.
cost = zeros(nTracks, nDetections);
% Scale cost function with amount of invisible time:
ages = [tracks(:).age];
totalVisibleCounts = [tracks(:).totalVisibleCount];
invisibility = 1-totalVisibleCounts ./ ages;
% if total invisibility count is larger than 25 images, set invisibility to
% a very high value so a new track is created:
invisibility([tracks(:).consecutiveInvisibleCount]>MaxTimeDiff) = 10^10;

for i = 1:nTracks
    %cost(i, :) = (1+5*invisibility(i))*distance(tracks(i).kalmanFilter, centroids);
    % Get centroid of last detection:
    LastIdx = find(~isnan(tracks(i).centroid(:,1)),1,'last');
    % calculate cost function as euclidean distance between cell centers.
    % Take into account the time the track was invisible.
    cost(i, :) = (1+5*invisibility(i))*hypot(tracks(i).centroid(LastIdx,1)-centroids(:,1), tracks(i).centroid(LastIdx,2)-centroids(:,2));
    
end

% Solve the assignment problem.
costOfNonAssignment = 150;
costOfNonAssignment = 50; % eucledean metric
[assignments, unassignedTracks, unassignedDetections] = ...
    assignmunkres(cost, costOfNonAssignment);
end

% The updateAssignedTracks function updates each assigned track with the 
% corresponding detection. It calls the correct method of 
% vision.KalmanFilter to correct the location estimate. Next, it stores 
% the new bounding box, and increases the age of the track and the total 
% visible count by 1. Finally, the function sets the invisible count to 0.
function updateAssignedTracks()
numAssignedTracks = size(assignments, 1);
for i = 1:numAssignedTracks
    trackIdx = assignments(i, 1);
    detectionIdx = assignments(i, 2);
    centroid = centroids(detectionIdx, :);   
    
    center = tracks(trackIdx).centroid;
          
    % Correct the estimate of the object's location
    % using the new detection.
    tracked_location = correct(tracks(trackIdx).kalmanFilter, centroid);
    center(end+1,1:2) = tracked_location(:,1:2);   
    % Update center positions over time:
    tracks(trackIdx).centroid = center;
    % Update other parameters:
    tracks(trackIdx).centroidPredicted = [tracks(trackIdx).centroidPredicted;predictedCentroid];
    
    % Update track's age.
    tracks(trackIdx).age = tracks(trackIdx).age + 1;
    
    % Update visibility.
    tracks(trackIdx).totalVisibleCount = ...
        tracks(trackIdx).totalVisibleCount + 1;
    tracks(trackIdx).consecutiveInvisibleCount = 0;
end
end

% Mark each unassigned track as invisible, and increase its age by 1.
function updateUnassignedTracks()
for i = 1:length(unassignedTracks)
    ind = unassignedTracks(i);
    tracks(ind).age = tracks(ind).age + 1;
    tracks(ind).consecutiveInvisibleCount = ...
    tracks(ind).consecutiveInvisibleCount + 1;
    tracks(ind).centroid(end+1,1:2) = NaN;
    tracks(ind).centroidPredicted(end+1,1:2) = NaN;
end
end

% The deleteLostTracks function deletes tracks that have been invisible for 
% too many consecutive frames. It also deletes recently created tracks that 
% have been invisible for too many frames overall.
function deleteLostTracks()
if isempty(tracks)
    return;
end

% Set parameters in such a way, that each detection - even for one image -
% is valid and all detections are kept:
invisibleForTooLong = 10^7;
ageThreshold = 1;

% Compute the fraction of the track's age for which it was visible.
ages = [tracks(:).age];
totalVisibleCounts = [tracks(:).totalVisibleCount];
visibility = totalVisibleCounts ./ ages;

% Find the indices of 'lost' tracks.
% lostInds = (ages < ageThreshold & visibility < 0.8) | ...
%     [tracks(:).consecutiveInvisibleCount] >= invisibleForTooLong;
lostInds = (ages < ageThreshold) | totalVisibleCounts < ageThreshold |...
    [tracks(:).consecutiveInvisibleCount] >= invisibleForTooLong;
% Delete lost tracks.
tracks = tracks(~lostInds);
end

% Create new tracks from unassigned detections. Assume that any unassigned 
% detection is a start of a new track. In practice, you can use other cues 
% to eliminate noisy detections, such as size, location, or appearance.
function createNewTracks()
centroids = centroids(unassignedDetections, :);
% From Huth et al:
param.motionModel           = 'ConstantVelocity';
param.initialEstimateError  = [5,sqrt(5)]; % Estimates from Huth2010 Supplement
param.motionNoise           = [5,sqrt(5)];
param.measurementNoise      = 5;

for i = 1:size(centroids, 1)
    centroid                = centroids(i,:);
    param.initialLocation	= centroid;
    
    % Add NaNs to the beginning to account for images that object was not
    % yet visible:
    centroid = [NaN(j-1,2);centroid];  
    % kalmanFilter = 'NotUsed';
    kalmanFilter = configureKalmanFilter(param.motionModel, ...
    	param.initialLocation, param.initialEstimateError, param.motionNoise , param.measurementNoise);
    % Create a new track.
    newTrack = struct(...
        'id', nextId, ...
        'kalmanFilter', kalmanFilter, ...
        'age', 1, ...
        'totalVisibleCount', 1, ...
        'consecutiveInvisibleCount', 0, ...
        'centroid', centroid, ...
        'centroidPredicted',centroid, ...
        'TrackTerminated', 0);
    
    % Add it to the array of tracks.
    tracks(end + 1) = newTrack;
    
    % Increment the next id.
    nextId = nextId + 1;
end
end

end







