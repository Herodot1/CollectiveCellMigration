function [SpeedAll,RMSVelAll] = AnalyzeVelField(VelField,m,CenterSpeed)

for i = 1:m
    % Get absolute speed values:
    % Correct velocity field for drifts (because motion is - often - supposed to have zero mean):
    VelX = VelField(:,:,1,i);
    VelX = imgaussfilt(VelX,0.71);
    VelY = VelField(:,:,2,i);
    VelY = imgaussfilt(VelY,0.71);
    MeanVelX(i) = nanmean(VelX(:));
    MeanVelY(i) = nanmean(VelY(:));
    if CenterSpeed == 1
        VelX = VelX - nanmean(VelX(:));
    end    
    if CenterSpeed == 1
        VelY = VelY - nanmean(VelY(:));
    end   
    % Calculate speed:
    Speed = sqrt(VelX.^2+VelY.^2);
    SpeedAll(i) = nanmean(Speed(:));
    % Root mean squared speed:
    RMSVelAll(i) = rms(Speed(:));
end

