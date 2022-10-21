function [CenterDivisions] = GetCellDivisionsSegmentation(Im,Template,SegNet,net,MinArea,BlockSize,UseCLAHE,MedFiltSize,BackThresh,FiltType,WSize)
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

% Block size used for training:
%BlockSize = [256,256];
% Mininmal object size:
%MinArea = 500;

% Use CLAHE:
% Use filter size of the same dimension as the images fed to the ANN:
% clahe currently destroys classification of the classification ANN
if UseCLAHE
    FiltSize = net.Layers(1,1).InputSize(1);
    if max(Im(:))>1
        Im = Im./max(Im(:));
    end
    Im = adapthisteq(Im,'NumTiles',[round(size(Im,1)/FiltSize),round(size(Im,2)/FiltSize)],'ClipLimit', 0.003, 'NBins', 256, 'Range', 'full', 'Distribution', 'exponential'); 
end


out = segmentImage(Im, SegNet, BlockSize);
% make it binary:
out = out >1.5;
% Reverse image as positive signals are currently set to false ... stupid
% me!
out = ~out;
% morphological opening:
out = imopen(out,strel('disk',5));
out = imclose(out,strel('disk',5));
% Remove small objects:
out = bwareaopen(out,MinArea);
% get centers for the respective objects:
stats = regionprops(out,'Centroid','BoundingBox');

CenterDivisions = [];
if ~isempty(stats)
    for j = 1:size(stats,1)
        xMin = ceil(stats(j).BoundingBox(1))-10;
        xMax = xMin + stats(j).BoundingBox(3) - 1+2*10;
        yMin = ceil(stats(j).BoundingBox(2))-10;
        yMax = yMin + stats(j).BoundingBox(4) - 1+2*10;
        xMin(xMin<=0) = 1;
        yMin(yMin<=0) = 1;
        xMax(xMax>=size(Im,2)) = size(Im,2);
        yMax(yMax>=size(Im,1)) = size(Im,1);
        %j
        tmp = Im(yMin:yMax,xMin:xMax);
        CCorr = ones(size(tmp));
        %figure; imshow(tmp)
        for k = 1:min([size(Template,1),3])
            tmpCC = xcorr2e(Template{k},double(tmp),'same');
            tmpCC = medfilt2(tmpCC,[5,5]);
            tmpCC = flip(tmpCC,1);
            tmpCC = flip(tmpCC,2);
            %tmpCC = (tmpCC-NormFactor(k,2)-min(tmpCC)) ./ (NormFactor(k,1)-NormFactor(k,2));
            % normalization factors:
            tmpCC = (tmpCC-min(tmpCC(:))) ./ (max(tmpCC(:))-min(tmpCC(:)));
            CCorr = CCorr .* tmpCC;
        end
        %figure; imshow(CCorr,[])
        % figure; imshow(tmp,[])
        % Classify and assign:
        %[ImClass(j),prob(j,1:2)] = classify(net,tmp);
        
        % get cross correlation map and find all potential peaks to feed to
        % the neuronal network. Eliminate candidates with too low
        CCorr = gather(CCorr);
        CCorr = CCorr./max(CCorr(:));
        %figure; imshow(CCorr,[])
        % Get peak centers:
        [Cent]=FastPeakFind(CCorr, BackThresh*2^16, FiltType ,5, 1);
        % x and y coordinates:
        x = Cent(1:2:end);
        y = Cent(2:2:end);
        % Adjust coordinates via xMin/yMin:
        x = x +xMin;
        y = y +yMin;
        
        
        % Create window of given size to classify the given object
        % create patches for ANN classification, assuming a window of
        % +-WSize px around detected cell division figure is enough to capture
        % the whole dividing cell. Afterwards classify it.
        % WSize = 35;
        ImClass = NaN(size(x,1),1);
        if ~isempty(x)
            for k = 1:size(x,1)
                if x(k)-WSize>0 && x(k)+WSize<size(Im,2) && y(k)-WSize>0 && y(k)+WSize<size(Im,1)
                    tmp = Im(y(k)-WSize:y(k)+WSize,x(k)-WSize:x(k)+WSize);
                    %figure; imshow(tmp)                    
                else
                    %occ = [j,k]
                    % Adjust window to fit into image and match size of
                    % used template:
                    yS = y(k)-WSize;
                    yE = y(k)+WSize;
                    xS = x(k)-WSize;
                    xE = x(k)+WSize;
                    if xS <= 0
                       xS = 1;
                       xE = 1+2*WSize;
                    end
                    if yS <= 0
                       yS = 1;
                       yE = 1+2*WSize;
                    end
                    if xE > size(Im,2)
                       xS = size(Im,2)-2*WSize;
                       xE = size(Im,2);
                    end
                    if yE > size(Im,1)
                       yS = size(Im,1)-2*WSize;
                       yE = size(Im,1);
                    end
                    tmp = Im(yS:yE,xS:xE);
                    
                end
                % While formally correct and usefull, when working with the
                % ANN, the current goodness of the classificator is too poor to
                % be usefull. Thus the approach as used above is used even
                % though no class of value 0 (see below) will ever be used due
                % to the casting from categorical to double.
                % Classify and assign:
                [ImClassTmp,prob(k,1:2)] = classify(net,tmp);
                % For gods sake why is the categorical [0,1] transformed
                % into [1,2] when casting to double? That is why I have
                % to do the ruckus with string command inbetween.
                ImClass(k) = double(string(ImClassTmp));
                
                % Classify and assign:
                %[ImClass(k),prob(k,1:2)] = classify(net,tmp);
            end
        end
        % Remove detections not identified as divisions:
        x(ImClass == 0) = [];
        y(ImClass == 0) = [];
        Start = length(CenterDivisions)+1;
        Ende = length(CenterDivisions)+length(x);
        CenterDivisions(Start:Ende,1) = x;
        CenterDivisions(Start:Ende,2) = y;
        
    end
end






