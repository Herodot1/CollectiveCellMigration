function [MeanDivs,StdDivs,DivTimes] = CumCellDivs(CellDivData,t_max,dt,groups,ImPhysArea)
% Calculates cumulative cell division events

% Output: 
% MeanDivs/StdDivsMean: number/standard devision of divisions at a specific time
% DivTimes: Time in minutes 

% Allocate variables:
DivTimes = [0:dt:t_max-dt]'.*ones(t_max./dt,1);
MeanDivs = NaN(t_max./dt,length(groups));
StdDivs = NaN(t_max./dt,length(groups));
for i = groups
    % Divisions at specific time:
    % Just fill in times with no divisions as identical to the values
    % before:
    tmpDivs = CellDivData(i).NumDivsAll;
    tmpTime = CellDivData(i).DivTimesAll;
    % Fill missing values with the ones from the previos time step as the
    % number of cell divisions did not change:
    DivTmp = NaN(length(DivTimes),size(tmpDivs,2));
    for j = 1:size(tmpTime,2)
        tmp = tmpTime(:,j);
        IdxIn = ismember(DivTimes,tmp);
        IdxIn = find(IdxIn);
        IdxOut = ~ismember(DivTimes,tmp);
        IdxOut = find(IdxOut);
        DivTmp(IdxIn,j) = tmpDivs(~isnan(tmpDivs(1:length(IdxIn),j)),j);
        for k = 1:length(IdxOut)
            % Set all NaNs in DivTmp to the first previously measured value:
            Difference = IdxIn-IdxOut(k);
            Difference(Difference>0) = [];            
            if ~isempty(Difference)
                % The last element is the point of interest (due to the ordering
                % of IdxIn), if it is not empty, otherwise 0 divisions occurred
                % so far:
                DivTmp(IdxOut(k),j) =  DivTmp(IdxIn(length(Difference)),j);
            else
                DivTmp(IdxOut(k),j) = 0;
            end
        end
    end  
    MeanDivs(:,i) = nanmean(DivTmp,2);
    StdDivs(:,i) = nanstd(DivTmp,[],2);
end   
% Rescale number of divisions to divisions per mm²;
MeanDivs = MeanDivs.* 1/ImPhysArea;
StdDivs = StdDivs.* 1/ImPhysArea;