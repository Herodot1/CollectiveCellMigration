function [CenterDivisions] = GetCellDivisions(im,Template,net,MedFiltSize,BackThresh,FiltType,WSize)
% Input:
% im = input image
% Template = templates (cell) of cell divisions used for pattern matching
% net = trained neuronal network used for classification of cell division
% events.
% Output:
% CenterDivisions = x/y coordinates of detected cell division events
%%
% Calculate  2d cross-correlation associated with division template:
% Maybe the sum of correlation maps after usage of multiple
% templates:
CCorr = ones(size(im));
for j = 1:size(Template,1)
    %tmp = xcorr2e(Template{j},im2./max(im2(:)),'same');
    tmp = xcorr2e(Template{j},im,'same');
    tmp = medfilt2(tmp,MedFiltSize);
    tmp = flip(tmp,1);
    tmp = flip(tmp,2);
    tmp = (tmp-min(tmp(:))) ./ (max(tmp(:))-min(tmp(:)));
    CCorr = CCorr .* tmp;
end
CCorr = CCorr./max(CCorr(:));
%figure; imshow(CCorr,[])
% Get peak centers:
[Cent]=FastPeakFind(CCorr, BackThresh*2^16, FiltType ,5, 1);
% x and y coordinates:
x = Cent(1:2:end);
y = Cent(2:2:end);

% create patches for ANN classification, assuming a window of
% +-35 px around detected cell division figure is enough to capture
% the whole dividing cell. Afterwards classify it.
% WSize = 35;
ImClass = NaN(size(x,1),1);
if ~isempty(x)
    for j = 1:size(x,1)
        if x(j)-WSize>0 && x(j)+WSize<size(im,2) && y(j)-WSize>0 && y(j)+WSize<size(im,1)
            tmp = im(y(j)-WSize:y(j)+WSize,x(j)-WSize:x(j)+WSize);
            %figure; imshow(tmp)
            % Classify and assign:
            [ImClass(j),prob(j,1:2)] = classify(net,tmp);
        end
    end
end
% Remove detections not identified as divisions:
x(ImClass == 0) = [];
y(ImClass == 0) = [];
CenterDivisions(1:length(x),1) = x;
CenterDivisions(1:length(x),2) = y;