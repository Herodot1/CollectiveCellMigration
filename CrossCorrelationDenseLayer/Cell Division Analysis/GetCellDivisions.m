function [CenterDivisions] = GetCellDivisions(im,Template,net,UseCLAHE,MedFiltSize,BackThresh,FiltType,WSize)
% Get the center of cell division events. The algorithm uses multiple 
% previously % segmented cell divisions events as templates for cross
% correlation. Individual cross correlation maps are multiplied and regions
% around peaks in the cross correlation matches are than fed into a ANN for
% classification (true cell division or not).

% Input:
% im = input image
% Template = templates (cell) of cell divisions used for pattern matching
% net = trained neuronal network used for classification of cell division
% events.
% MedFiltSize = size of median filter used for denoising 2d cross
% correlation mao
% BackThresh = background threshold -> minimal value (range: 0-1) of a peak
% in the cross correlation map
% FiltType =used filter type for denoising the joint cross correlation map
% WSize = window size -> create a square of diameter 2*WSize + 1 pixel
% around each suggested cell division for classification using the trained
% neuronal network (net).
% Output:
% CenterDivisions = x/y coordinates of detected cell division events
%%

% Use CLAHE:
% Use filter size of the same dimension as the images fed to the ANN:

if UseCLAHE
    FiltSize = net.Layers(1,1).InputSize(1);
    if max(im(:))>1
        im = im./max(im(:));
    end
    im = adapthisteq(im,'NumTiles',[round(size(im,1)/FiltSize),round(size(im,2)/FiltSize)],'ClipLimit', 0.003, 'NBins', 256, 'Range', 'full', 'Distribution', 'exponential'); 
end
% Calculate  joint 2d cross-correlation associated with division templates,
% as the product of multiple templates. 
CCorr = ones(size(im));
for j = 1:size(Template,1)
    %tmp = xcorr2e(Template{j},im2./max(im2(:)),'same');
    tmp = xcorr2e(Template{j},im,'same');
    tmp = medfilt2(tmp,MedFiltSize);
    tmp = flip(tmp,1);
    tmp = flip(tmp,2);
    tmp = (tmp-min(tmp(:))) ./ (max(tmp(:))-min(tmp(:)));
    CCorr = CCorr .* tmp;
    % write as image:
    % imwrite(CCorr./max(CCorr(:)),sprintf('CrossCorrImage %02d.png',j),'png')
end
% transform CCorr back to normal array (from gpuArray):
CCorr = gather(CCorr);
CCorr = CCorr./max(CCorr(:));
%figure; imshow(CCorr,[])
% Get peak centers:
[Cent]=FastPeakFind(CCorr, BackThresh*2^16, FiltType ,5, 1);
% x and y coordinates:
x = Cent(1:2:end);
y = Cent(2:2:end);

% create patches for ANN classification, assuming a window of
% +-WSize px around detected cell division figure is enough to capture
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