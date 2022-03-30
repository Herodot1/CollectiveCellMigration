function [CorrLength,AutoCorrMean,AutoCorrStd] = FitAutoCorrelation(VelFieldData,VelFieldPxSize)

% Fit autocorrelation data to exponential decay function in µm:
x = (0:size(VelFieldData(2).AutoCorrAll,1)-1)*VelFieldPxSize;
x = x';
ft=fittype('exp1');
CorrLength = NaN(size(VelFieldData(1).AutoCorrAll,2), length(VelFieldData));
AutoCorrMean = NaN(size(VelFieldData(1).AutoCorrAll,2),size(VelFieldData(1).AutoCorrAll,1),length(VelFieldData));
AutoCorrStd = NaN(size(VelFieldData(1).AutoCorrAll,2),size(VelFieldData(1).AutoCorrAll,1),length(VelFieldData));
for i = 1:length(VelFieldData)
    for j = 2:size(VelFieldData(i).AutoCorrAll,2)-1        
        y = nanmean(VelFieldData(1,i).AutoCorrAll(:,j,:),3);
        AutoCorrMean(j,1:length(y),i) = y;         
        AutoCorrStd(j,1:length(y),i) = nanstd(squeeze(VelFieldData(1,i).AutoCorrAll(:,j,:)),[],2); 
        if nansum(y(:)) ~=0
            % Old, changed with 06-04-2021
            %pos = find(y <= 0,1)-1;
            pos = find(isnan(y) == 1,1)-1;
            if isempty(pos)
                pos = length(y);
            end
            [cf,G]=fit(x(1:pos),y(1:pos),ft,'Start', [1.0, -0.08]);
            CorrLength(j,i) = cf.b;
        end
        
%         figure;
%         hold on
%         plot(x,y)
%         plot(x, cf.a*exp(cf.b*x))
    end
end

% % Generiert keinen Unterschied
% CorrLength = NaN(size(VelFieldData(1).AutoCorrAll,2), length(VelFieldData));
% for i = 1:length(VelFieldData)
%     i
%     for j = 2:size(VelFieldData(1).AutoCorrAll,2) 
%         TempCorrLength = [];
%         for k = 1:size(VelFieldData(1,i).AutoCorrAll,3)
%             y = VelFieldData(1,i).AutoCorrAll(:,j,k);
%             if sum(y) ~= 0
%                 pos = find(y == 0,1)-1;
%                 cf=fit(x(1:pos),y(1:pos),ft,'Start', [1.1, -0.8]);
%                 TempCorrLength(k) = cf.b;
%             end
%             
%         end        
%         CorrLength(j,i) = nanmean(TempCorrLength);        
%     end
% end

CorrLength = abs(1./CorrLength);
