function [VTemp] = VelOrderWindow(VelField,WSize,dt)

% Get time dependent velocity order parameter:
VTemp = NaN(WSize,size(VelField,4));
N = size(VelField,1)*size(VelField,2);
for k = 1+(WSize-1)/2:size(VelField,4)-1-((WSize-1)/2)
    % Get time window:
    interval = k-(WSize-1)/2:k+(WSize-1)/2;       
    % Set velocity Field in x and y direction. Correct for drift as motion
    % is supposed to have zero mean:
    VelTempX = squeeze(VelField(:,:,1,interval));
    VelTempX = VelTempX - nanmean(VelTempX(:));
    VelTempY = squeeze(VelField(:,:,2,interval));
    VelTempY = VelTempY - nanmean(VelTempY(:));
    % Calculate absolute speed values:
    AbsVel = sqrt(VelTempX.^2+VelTempY.^2);
    % Normalize velocities in x and y direction:
    VelTempX = VelTempX./AbsVel;
    VelTempY = VelTempY./AbsVel;
    % Calculate order parameter:
    VelTemp = abs(squeeze((sum(sum(VelTempX,1),2) + sum(sum(VelTempY,1),2))))/N;
    % Remove NaNs:  
    offset = sum(isnan(VelTemp));
    VelTemp(isnan(VelTemp)) = [];
    % Set Trajectory and get the ensemble mean:
    Trajectory = [];
    Trajectory(:,1) = VelTemp;
    Trajectory(:,2) = [0:WSize-1-offset]*dt;
    VTemp(1:end-1-offset,k-(WSize-1)/2) = RelativeDistance(Trajectory);

end
