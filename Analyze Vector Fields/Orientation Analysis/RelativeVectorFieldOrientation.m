% Analyze the relative directions between cell orientation and the velocity
% fields. As both vectorfields are not necessarily of equal size the cell
% orientation field will be interpolated to match the velocity field
% coordinates, as the cell orientation vector field is likely set to a
% better resolution

% Input:
% CellOrientation: m x n x 2 matrix of cell orientation vectors (2D)
% VelField: m x n x 2 matrix of cell velocity vectors (2D)
% xVelField: x coordinates of the velocity field <- output of meshgrid
% yVelField: y coordinates of the velocity field <- output of meshgrid
% xCellOrientation: x coordinates of the orientation vector field <- output of meshgrid
% yCellOrientation: y coordinates of the orientation vector field <- output of meshgrid
% WinSize: moving average used for calculation of the mean of both vector
% fields in temporal dimension
% Binary: Binary images seperating fore and background.

% Output:
% Angle: Relative angle (0-90 degree) of the cellular orientation relative
% the velocity field


function [Angle] = RelativeVectorFieldOrientation(CellOrientationField,VelField,...
    xCellOrientation,yCellOrientation,xVelField,yVelField,WinSize,Binary)

% CorseGraining:
CG = 1;

% Allocate variables:
Angle = NaN(size(VelField,1)*size(VelField,2),size(VelField,4));
for i = 1:size(VelField,4)
    % Interpolate CellOrientation to the respective values of the
    % velocityfield:
    CellOrientationInterp = [];
    % x values:
    CellOrientationInterp(:,:,1) = interp2(xCellOrientation,yCellOrientation,CellOrientationField(:,:,2,i),xVelField(1:CG:end,1:CG:end),yVelField(1:CG:end,1:CG:end));
    % y values -> denote the change in x and y coordinates ... dont ask me
    % why the original algorithm switched this up ...
    CellOrientationInterp(:,:,2) = interp2(xCellOrientation,yCellOrientation,CellOrientationField(:,:,1,i),xVelField(1:CG:end,1:CG:end),yVelField(1:CG:end,1:CG:end)); 
    
    % restrict calculation to region of interest:
    if size(Binary,3) ~= 1
        tmp = Binary(:,:,i);
    else
        tmp = Binary;
    end
    tmp = tmp(yVelField(1:CG:end,1),xVelField(1,:)');
    CellOrientationInterp(:,:,1) = CellOrientationInterp(:,:,1).*tmp;
    CellOrientationInterp(:,:,2) = CellOrientationInterp(:,:,2).*tmp;
    % Calculate angle of intersecting vectors:
    AngleTemp = atan2d(CellOrientationInterp(:,:,1).*VelField(1:CG:end,1:CG:end,2,i) - CellOrientationInterp(:,:,2).*VelField(1:CG:end,1:CG:end,1,i), ...
        CellOrientationInterp(:,:,1).*VelField(1:CG:end,1:CG:end,1,i) + CellOrientationInterp(:,:,2).*VelField(1:CG:end,1:CG:end,2,i));
    % Remove points outside region of interest
    AngleTemp(CellOrientationInterp(:,:,1) == 0) = NaN;
    Angle(:,i) = AngleTemp(:);
end
% Transform to 0-180 degree scale:
Angle = abs(Angle);
% Is this step necessary?
Angle(Angle>90) = 180 - Angle(Angle>90);








