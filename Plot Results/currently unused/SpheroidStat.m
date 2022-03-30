function [stats,DataIn] = SpheroidStat(Data,TimePoint,FiltSize)

for i = 1:size(Data,2)
    NumEl(i) = size(Data{i},2);
end
DataIn = NaN(max(NumEl),size(Data,2));

for i = 1:size(Data,2)
    %DataIn(1:NumEl(i),i) = nanmean(Data{i}(TimePoint-((FiltSize-1)/2):TimePoint+((FiltSize-1)/2),:),1);
    DataIn(1:NumEl(i),i) = Data{i}(TimePoint,:);
end

[p,tbl,stats] = anova1(DataIn);