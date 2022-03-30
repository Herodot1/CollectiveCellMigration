function [MeanPeakTime,StdPeakTime,MeanCellsMoving,StdCellsMoving,...
    MeanPeakTimeTemp,StdPeakTimeTemp,MeanCellsMovingTemp,StdCellsMovingTemp] ...
    = CellClusterSizes(VelFieldData,dt,Sampling)

% Analyze the order parameter and 4-poin susceptibility in order to obtain
% the number of cells moving as a pack and the pack life time:

% ChiPeakHeight       = NaN(900,length(VelFieldData));
% ChiPeakStd          = NaN(900,length(VelFieldData));
% ChiPeakPos          = NaN(900,length(VelFieldData));
% QDropPos            = NaN(900,length(VelFieldData));
% for i = 1:length(VelFieldData)
%     % Peak position and value of 4 point susceptibility over time:
%     [PeakHeight,PeakPos]                    = max(nanmean(VelFieldData(i).ChiTempAll,3),[],1);
%     PeakStd                                 = nanstd(max(VelFieldData(i).ChiTempAll,[],1),[],3);
%     ChiPeakStd(1:length(PeakStd),i)         = PeakStd;
%     ChiPeakHeight(1:length(PeakHeight),i)   = PeakHeight;
%     ChiPeakPos(1:length(PeakPos),i)         = PeakPos;
%     % Find Point where Q reaches 0.9 first:
%     [~,PeakPos]                             = min(abs(nanmean(VelFieldData(i).QTempAll,3)-0.9),[],1);
%     QDropPos(1:length(PeakPos),i)           = PeakPos;
%      % Find Point where autocorrelation reaches 1/e  ~ 0.3679 first:
%     [~,PeakPos]             = min(abs(nanmean(VelFieldData(i).AutoCorrAll,3)-1/exp(1)),[],1);
%     AutoCorrDropPos(1:length(PeakPos),i)    = PeakPos;
%     % Velocity order parameter:
%     VTemp(1:size(VelFieldData(i).VTempAll,2),i) = nanmean(nanmean(VelFieldData(i).VTempAll,3),1);
%     
% end

% Pack size of collectively moving cells:
for i = 1:length(VelFieldData)
    % Peak position and value of 4 point susceptibility over time:
    [PeakHeight,PeakPos]                    = max(nanmean(VelFieldData(i).ChiAll,2),[],1);
    PeakStd                                 = nanstd(max(VelFieldData(i).ChiAll,[],1),[],3);
    ChiPeakStd(1:length(PeakStd),i)         = PeakStd;
    ChiPeakHeight(1:length(PeakHeight),i)   = PeakHeight;
    ChiPeakPos(1:length(PeakPos),i)         = PeakPos;
    Temp = nanmean(VelFieldData(i).QAll,2);
    CellsMoving(i) = PeakHeight./(1-Temp(PeakPos));
    PeakTime(i) = PeakPos*dt;
end




% Pack size of collectively moving cells extracted from global parameters:
for i = 1:length(VelFieldData)
    % Peak position and value of 4 point susceptibility over time:
    for j = 1:size(VelFieldData(i).ChiAll,2)
        [PeakHeight,PeakPos]                    = max(VelFieldData(i).ChiAll(:,j),[],1);
        ChiPeakHeight(j,i)   = PeakHeight;
        ChiPeakPos(j,i)      = PeakPos;
        Temp = VelFieldData(i).QAll(:,j);
        CellsMoving(j,i) = PeakHeight./(1-Temp(PeakPos));
        PeakTime(j,i) = PeakPos*dt;
    end
    %[PeakHeight,PeakPos]                    = max(nanmean(ChiQData(i).ChiAll,2),[],1);
    %PeakStd                                 = nanstd(max(ChiQData(i).ChiAll,[],1),[],3);
    %ChiPeakStd(1:length(PeakStd),i)         = PeakStd;
    %ChiPeakHeight(1:length(PeakHeight),i)   = PeakHeight;
    %ChiPeakPos(1:length(PeakPos),i)         = PeakPos;
    %Temp = nanmean(ChiQData(i).QAll,2);
    %CellsMoving(i) = PeakHeight./(1-Temp(PeakPos));
    %PeakTime(i) = PeakPos*dt;
end
% Remove artificial values:
Idx = [find(ChiPeakPos == 0); find(ChiPeakPos == 1)];
ChiPeakHeight(Idx) = NaN;
%nanmean(ChiPeakHeight,1)
%nanstd(ChiPeakHeight,[],1)
ChiPeakPos(Idx) = NaN;
%nanmean(ChiPeakPos,1)
%nanstd(ChiPeakPos,[],1)
PeakTime(Idx) = NaN;
MeanPeakTime = nanmean(PeakTime,1);
StdPeakTime = nanstd(PeakTime,[],1);
CellsMoving(Idx) = NaN;
MeanCellsMoving = nanmean(CellsMoving,1);
StdCellsMoving = nanstd(CellsMoving,[],1);

%%%%%
% Pack size of collectively moving cells extracted from local parameters:
for i = 1:length(VelFieldData)
    
    % Peak position and value of 4 point susceptibility over time:
    for j = 1:size(VelFieldData(i).ChiTempAll,3)
        %Idx = find(isnan(VelFieldData(i).ChiTempAll(3,:,j)) == 0);
        Idx = 1:Sampling:size(VelFieldData(i).ChiTempAll,2);
        count = 0;
        for k = Idx
            count = count + 1;
            [PeakHeight,PeakPos]                    = max(VelFieldData(i).ChiTempAll(:,k,j),[],1);
            ChiPeakHeightTemp(count,j,i)   = PeakHeight;
            ChiPeakPosTemp(count,j,i)      = PeakPos;
            Temp = VelFieldData(i).QTempAll(:,k,j);
            CellsMovingTemp(count,j,i) = PeakHeight./(1-Temp(PeakPos));
            PeakTimeTemp(count,j,i) = PeakPos*dt;
        end
    end
    %[PeakHeight,PeakPos]                    = max(nanmean(ChiQData(i).ChiAll,2),[],1);
    %PeakStd                                 = nanstd(max(ChiQData(i).ChiAll,[],1),[],3);
    %ChiPeakStd(1:length(PeakStd),i)         = PeakStd;
    %ChiPeakHeight(1:length(PeakHeight),i)   = PeakHeight;
    %ChiPeakPos(1:length(PeakPos),i)         = PeakPos;
    %Temp = nanmean(ChiQData(i).QAll,2);
    %CellsMoving(i) = PeakHeight./(1-Temp(PeakPos));
    %PeakTime(i) = PeakPos*dt;
end
% Remove artificial values:
Idx = [find(ChiPeakPosTemp == 0); find(ChiPeakPosTemp == 1)];
ChiPeakHeightTemp(Idx) = NaN;
%nanmean(ChiPeakHeightTemp,2)
%nanstd(ChiPeakHeightTemp,[],2)
ChiPeakPosTemp(Idx) = NaN;
%nanmean(ChiPeakPosTemp,2)
%nanstd(ChiPeakPosTemp,[],2)
PeakTimeTemp(Idx) = NaN;
MeanPeakTimeTemp = squeeze(nanmean(PeakTimeTemp,2));
StdPeakTimeTemp = squeeze(nanstd(PeakTimeTemp,[],2));
CellsMovingTemp(Idx) = NaN;
MeanCellsMovingTemp = squeeze(nanmean(CellsMovingTemp,2));
StdCellsMovingTemp = squeeze(nanstd(CellsMovingTemp,[],2));

save('CellsMovingPeakTime.mat','MeanPeakTime','StdPeakTime','MeanCellsMoving'...
    ,'StdCellsMoving','MeanPeakTimeTemp','StdPeakTimeTemp','MeanCellsMovingTemp','StdCellsMovingTemp')
