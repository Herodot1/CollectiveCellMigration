function tracks = MultiObjectTracking(centroids_start,radii_start)
tracks = initializeTracks(); % Create an empty array of tracks.
nextId = 1; % ID of the next track
for j = 1:size(centroids_start,3)
    
    centroids  = centroids_start(:,:,j);
    pos_in = find(isnan(centroids(:,1)) == 1,1);
    if ~isempty(pos_in)
        centroids = centroids(1:pos_in-1,:);
    end
    
    radiiTemp  = radii_start(:,j);
    if ~isempty(pos_in)
        radiiTemp = radiiTemp(1:pos_in-1,:);
    end
    
    % bboxes = M-by-4 matrix of [x y width height] bounding box
    % coordinates, where M represents the number of blobs and [x y]
    % represents the upper left corner of the bounding box.
    bboxes = zeros(size(centroids,1),4);
    bboxes(:,1) = centroids(:,1);
    bboxes(:,2) = centroids(:,2);
    bboxes(:,3) = radiiTemp;
    bboxes(:,4) = radiiTemp;
    
    predictNewLocationsOfTracks();
    [assignments, unassignedTracks, unassignedDetections] = ...
        detectionToTrackAssignment();
    updateAssignedTracks();
    updateUnassignedTracks();
    deleteLostTracks();
    createNewTracks();
    %displayTrackingResults();
end

function tracks = initializeTracks()
% create an empty array of tracks
tracks = struct(...
    'id', {}, ...
    'bbox', {}, ...
    'kalmanFilter', {}, ...
    'age', {}, ...
    'totalVisibleCount', {}, ...
    'consecutiveInvisibleCount', {}, ...
    'centroid',{}, ...
    'radii',{}, ...
    'Boarder', {}, ...
    'Mitosis', {}, ...
    'TrackTerminated', {});
end

% Use the Kalman filter to predict the centroid of each track in the 
% current frame, and update its bounding box accordingly.
function predictNewLocationsOfTracks()
for i = 1:length(tracks)
    bbox = tracks(i).bbox;
    
    % Predict the current location of the track.
    predictedCentroid = predict(tracks(i).kalmanFilter);
    
    % Shift the bounding box so that its center is at
    % the predicted location.
    % predictedCentroid = int32(predictedCentroid) - bbox(3:4) / 2;
    % tracks(i).bbox = [predictedCentroid, bbox(3:4)];
    predictedCentroid = int32(predictedCentroid);
    tracks(i).bbox = [predictedCentroid, predictedCentroid];
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
for i = 1:nTracks
    cost(i, :) = (1+1*invisibility(i))*distance(tracks(i).kalmanFilter, centroids);
end
% Scale cost function with amount of invisible time:



% Solve the assignment problem.
costOfNonAssignment = 10;
[assignments, unassignedTracks, unassignedDetections] = ...
    assignDetectionsToTracks(cost, costOfNonAssignment);
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
    bbox = bboxes(detectionIdx, :);
    center = tracks(trackIdx).centroid;
    radius = tracks(trackIdx).radii; 
    
    % Correct the estimate of the object's location
    % using the new detection.
    tracked_location = correct(tracks(trackIdx).kalmanFilter, centroid);
    
    % redefine the bounding box: 
    bbox(:,1) = tracked_location(:,1);
    bbox(:,2) = tracked_location(:,2);
    bbox(:,3) = tracked_location(:,1);
    bbox(:,4) = tracked_location(:,2);
    center(end+1,1:2) = tracked_location(:,1:2);
    radius(end+1,1) = radiiTemp(detectionIdx); 
    
    % Replace predicted bounding box with detected
    % bounding box.
    tracks(trackIdx).bbox = bbox;
    
    % Update center positions over time:
    tracks(trackIdx).centroid = center;
    
    % Update radii over time:
    tracks(trackIdx).radii = radius;
    
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
    tracks(ind).radii(end+1) = NaN;
end
end

% The deleteLostTracks function deletes tracks that have been invisible for 
% too many consecutive frames. It also deletes recently created tracks that 
% have been invisible for too many frames overall.
function deleteLostTracks()
if isempty(tracks)
    return;
end

invisibleForTooLong = 10000;
ageThreshold = 3;

% Compute the fraction of the track's age for which it was visible.
ages = [tracks(:).age];
totalVisibleCounts = [tracks(:).totalVisibleCount];
visibility = totalVisibleCounts ./ ages;

% Find the indices of 'lost' tracks.
lostInds = (ages < ageThreshold & visibility < 0.6) | ...
    [tracks(:).consecutiveInvisibleCount] >= invisibleForTooLong;

% Delete lost tracks.
tracks = tracks(~lostInds);
end

% Create new tracks from unassigned detections. Assume that any unassigned 
% detection is a start of a new track. In practice, you can use other cues 
% to eliminate noisy detections, such as size, location, or appearance.
function createNewTracks()
centroids = centroids(unassignedDetections, :);
radius = radiiTemp(unassignedDetections);
%bboxes: M-by-4 matrix of [x y width height] bounding box coordinates, 
% where M represents the number of blobs and [x y] represents the upper 
% left corner of the bounding box.
bboxes = bboxes(unassignedDetections, :);

% From Huth et al:
param.motionModel           = 'ConstantVelocity';
param.initialEstimateError  = [5,sqrt(5)]; % Estimates from Huth2010 Supplement
param.motionNoise           = [5,sqrt(5)];
param.measurementNoise      = 5;

for i = 1:size(centroids, 1)
    
    centroid = centroids(i,:);
    bbox = bboxes(i, :);
    radius = radiiTemp(i);
    param.initialLocation       = centroid;
    
    % Add NaNs to the beginning to account for images that object was not
    % yet visible:
    centroid = [NaN(j-1,2);centroid];  
    radius = [NaN(j-1,1);radius]; 
    % Create a Kalman filter object.
    %kalmanFilter = configureKalmanFilter('ConstantVelocity', ...
    %    centroid, [200, 50], [100, 25], 100);
    
    kalmanFilter = configureKalmanFilter(param.motionModel, param.initialLocation,...
        param.initialEstimateError, param.motionNoise, param.measurementNoise);
    
    % Create a new track.
    newTrack = struct(...
        'id', nextId, ...
        'bbox', bbox, ...
        'kalmanFilter', kalmanFilter, ...
        'age', 1, ...
        'totalVisibleCount', 1, ...
        'consecutiveInvisibleCount', 0, ...
        'centroid', centroid, ...
        'radii', radius,...
        'Boarder', 0, ...
        'Mitosis', 0, ...
        'TrackTerminated', 0);
    %'bbox', bbox, ...
    
    %newTrack(1,1).id
    
    % Add it to the array of tracks.
    tracks(end + 1) = newTrack;
    
    % Increment the next id.
    nextId = nextId + 1;
end
end

end







