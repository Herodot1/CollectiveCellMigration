% Remarks: When applying all this to scratch normalize divisions to rhe
% region of interest, such as the cell covered area

%e.g.:
% % Calculate cell divisions over time:
% [SumDivAll,SumDivStdAll] = CumSumDivisions(VelFieldData);
% % F�r die Statistik:
% StatData = squeeze(max(SumDivAll,[],1));
% % Normalize data to mean of all experiments of the temporal median values 
% % of the covered area (complex wording ...):
% RelArea = [];
% for i = 1:length(VelFieldData)
% 	RelArea(i) = nanmean(nanmedian(VelFieldData(i).CoveredAreaAll));
% end
% % Give relative area in units of the control conditions:
% RelArea = RelArea./RelArea(1);
% % Normalize data:
% StatData = StatData./RelArea';
% % Save Data as new xls sheet:
% writetable(array2table(StatData),'CellDivisionsCenter.xlsx')
% 
% [p,tbl,stats] = anova1(StatData');
% [c,m,h,gnames] = multcompare(stats,'CType','hsd');
% 
% SumDivAll = nanmean(SumDivAll,3);










% Plot comparison plots:
a
addpath(pwd)
addpath(strcat(pwd,'\Vectorfield Visualization'))
clearvars -except VelFieldData VecFieldData TrackMat CellDivData;
%% Set constants:
% Groups to compare:
groups = [2,5,10];
% time difference between cell density estimates in minutes:
DensDt = 240;
% time difference between sucessive images in minutes:
dt = 15;
%dt = 5;
dt = 3;
% max time to plot in min:
t_max = 3*1200; %4000;
t_max = 3*1400; %4000;
%t_max = 1100; % Colon Carcinoma;
% Sampling of MSD, Chi and Q;
Sampling = 40; %240;
Sampling = 200; %240;
% Julian:
% dt = 5;
% Set offset corresponding to the number of images before and after the set
% timepoint that are taken as reference:
offset = 10;
% Area proportion where cells have been counted:
dA = 0.3^2;
% Image size in pixels:
ImSize = [1280,960];
% julian
% ImSize = [1040,1392];
% Physical dimension of image in �m:
ImPhysSize = [621.23,465.80];
% Julian
% ImPhysSize = [478.18,640.18];
% pixel size of image in �m/px:
PxSize = ImPhysSize/ImSize;
% Area of image in mm^2:
ImPhysArea = ImPhysSize(1)*ImPhysSize(2)./10^6;
% Factor for absolute density (cells/mm^2);
CellDensFac = 1/(ImPhysArea*dA);
% Size Velocity Field in px:
% GBM Cells
VelFieldSize = [106,79];
% Colon Carcinoma/ Ula Messungen:
VelFieldSize = [97,72];
% Julian:
% VelFieldSize = [86,115];
% Pixelsize of velocity field in �m/px:
VelFieldPxSize = mean((ImSize./VelFieldSize))*PxSize;
% Text font size:
text_font = 24;
% Windowsize:
WSize = VelFieldData(1).WSize;

% Times of sampling:
Idx = [1:DensDt/dt:length(VelFieldData(1).tau)+1];
Idx(1) = 2;

% Supress k-means warnings for NaN in data-set
warning('off','stats:kmeans:MissingDataRemoved')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time evolution of some critical parameters:
[MeanPeakTime,StdPeakTime,MeanCellsMoving,StdCellsMoving,...
    MeanPeakTimeTemp,StdPeakTimeTemp,MeanCellsMovingTemp,StdCellsMovingTemp] ...
    = CellClusterSizes(VelFieldData,dt,Sampling);

save('CellsMovingPeakTime.mat','MeanPeakTime','StdPeakTime','MeanCellsMoving'...
    ,'StdCellsMoving','MeanPeakTimeTemp','StdPeakTimeTemp','MeanCellsMovingTemp','StdCellsMovingTemp')

% Create scaling coefficient as function of used lag-time for the whole
% time window:
[pMeanAll,pStdAll] = MSDScalingCoeff(VelFieldData);

% Create scaling coefficient as function of used lag-time for different
% time windows only:
[pMeanAllTemp,pStdAllTemp] = MSDScalingCoeffTemp(VelFieldData,Sampling);

% Calculate Autocorrelation length:
tic
[CorrLength,AutoCorrMean,AutoCorrStd] = FitAutoCorrelation(VelFieldData,VelFieldPxSize);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot means of VelFieldData for different cell densities. Global parameters.
% Groups to compare:
groups = [1:5];
groups = [1:8];
% Savename:

savename = 'MultiCellStreaming Colagenase';
savename = 'Astros + LN + U138 3d';
savename = 'LN229 Center Part 2';
savename = 'LN229 + NM + HY';
savename = 'LN229 + NM + Low FBS';
savename = 'LN229 + ACM';
savename = 'LN229 + Alarmines 1d';
savename = 'LN229 + MGs';
savename = 'LN229 + CBs';
savename = 'A375 + ACM/BV2CM';
savename = 'SW480';
savename = 'SW620';
savename = 'Astros ROCK + Myosin II Inhibitor';
savename = 'U138 Monolayer';
savename = 'Astros+LN+U138 1d';
savename = 'Astros + MG Scratch';
savename = 'Astros Monolayer';
savename = 'LN229 ROCK Inhibitor';
savename = 'LN229 Myosin II Inhibitor';
savename = 'LN229 ROCK + Myosin II Inhibitor';
savename = 'U138 ROCK + Myosin II Inhibitor';
LegendNames = [];
count = 0;
for i = groups
    count = count + 1;
    LegendNames{count} = VelFieldData(i).Name;
end
% 5 groups:
color_map = [1,0,0; 0,0,1; 0,0,0; 0.5,0,0; 0.5,0.5,0.5;];
% 9 groups
color_map = [0,0,0; 0,0,1; 0,0,0.5; 1,0,0; 0.5,0,0; 0,1,0; 0,0.5,0; 0,1,1; 0,0.5,0.5];
% 6 groups
color_map = [0,0,0; 0,0,1; 0,0,0.5; 1,0,0; 0.66,0,0; 0.33,0,0;];
% 3 groups
color_map = [1,0,0; 0,0,1; 0,0,0;];
% 4 groups
color_map = [1,0,0; 0,0,1; 0,1,0; 0,0,0;];
% 5 groups
color_map = [1,0,0; 0.75,0,0; 0.5,0,0; 0.25,0,0; 0,0,0;];
% 6 groups
color_map = [0,0,0; 1,0,0; 0.75,0,0; 0.5,0,0; 0,0,1; 0,0,0.75;];
% 6 groups
color_map = [0,0,0; 0,0,1; 0,0,0.75; 1,0,0; 0.75,0,0; 0.5,0,0;];
% 8 groups
color_map = [0,0,0; 1,0,0; 0.75,0,0; 0.5,0,0; 0.25,0,0; 0,0,1; 0,0,0.75; 0,0,0.5];
% 9 groups
color_map = [1,0,0; 0.66,0,0; 0.33,0,0; 0,0,1; 0,0,0.66; 0,0,0.33; 0,0,0; 0.66,0.66,0.66; 0.33,0.33,0.33];
% 9 groups
color_map = [1,0,0;        0,0,1;          0,0,0; ...
    0.66,0,0;     0,0,0.66;       0.66,0.66,0.66; ...
    0.33,0,0;     0,0,0.33;       0.33,0.33,0.33];
% 12 groups
color_map = [0,0,0; 1/4,1/4,1/4; 2/4,2/4,2/4; 3/4,3/4,3/4; ...
    1,0,0; 3/4,0,0; 2/4,0,0; 1/4,0,0;  ...
    0,0,1; 0,0,3/4; 0,0,2/4; 0,0,1/4];
% 14 groups
color_map = [0,0,0; 1,0,0;   2/4,2/4,2/4; 3/4,3/4,3/4; ...
            0,0,1; 0,0,2/4; 0,1,1; 0,2/4,2/4
            1,0,1; 2/4,0,2/4; 1,1,0; 2/4,2/4,0
            2/4,0,0; 1/4,0,0;];
% 8 groups:
color_map = [0,0,0; 1,0,0;   2/4,2/4,2/4; 3/4,3/4,3/4; ...
            0,0,1; 0,0,2/4; 0,1,1; 0,2/4,2/4];        
% 8 groups:        
color_map = [1,0,0; 3/4,0,0;   2/4,0,0; 1/4,0,0; ...
            0,0,1; 0,0,3/4; 0,0,2/4; 0,0,1/4];
% 4 Groups:
color_map = [0,0,0; 1,0,0; 0,0,1; 1,0,1];

% 5 Groups:
color_map = [0.5,0.5,0.5; 0,0,0; 1,0,0; 0,0,1; 1,0,1];

%%
% Actual plots:

% Plot nearest neighbor stuff:
% Agglomerate Variable:
PlotVar = NaN(length(groups),1700);
count = 0;
for i = groups
    count = count + 1;    
    PlotVar(count,1:size(VelFieldData(i).NumNewNeighborsAll,1)) = 1-VelFieldData(i).NumNewNeighborsAll(:,1);    
end
% shape factor over density:
h=figure('units','normalized','outerposition',[0 0 1 1]);
set(h,'DefaultTextFontSize',24)
hold all
%h = boxplot(flip(Young,2),{'10�M Blebbistatin', '40�M Y-27632','U138','5�M Blebbistatin','20�M Y-27632','LN229'},'BoxStyle','outline','Width',0.5,'LabelOrientation','inline','LabelVerbosity','majorminor','Orientation','horizontal');
h = boxplot(PlotVar',LegendNames,'BoxStyle','outline','Width',0.5,'LabelOrientation','inline','LabelVerbosity','majorminor','Orientation','horizontal');
set(h,{'linewidth'},{2})
xlabel('Relative Number of Cells with no Changed Neighbors','FontSize',24)
set(gca,'fontsize', 24);
title('Preserved Neighborhood','FontSize',30);
saveas(gcf, sprintf('%s Relative Number of Cells with no Changed Neighbors.png',savename));
savefig(sprintf('%s Relative Number of Cells with no Changed Neighbors.fig',savename))
close all


% Plot nearest neighbor stuff as function of time:
CellType = 4;
Names = {'0-20 h','10-30 h','20-40 h','30-50 h','40-60 h'};
% Agglomerate Variable:
PlotVar = NaN(100,50);
count = 0;
for i = 1:Sampling:size(VelFieldData(CellType).NumNewNeighborsTempAll,1)
    count = count + 1;    
    PlotVar(count,1:size(VelFieldData(CellType).NumNewNeighborsTempAll,2)) = 1-VelFieldData(CellType).NumNewNeighborsTempAll(i,:);    
end
% shape factor over density:
h=figure('units','normalized','outerposition',[0 0 1 1]);
set(h,'DefaultTextFontSize',24)
hold all
%h = boxplot(flip(Young,2),{'10�M Blebbistatin', '40�M Y-27632','U138','5�M Blebbistatin','20�M Y-27632','LN229'},'BoxStyle','outline','Width',0.5,'LabelOrientation','inline','LabelVerbosity','majorminor','Orientation','horizontal');
h = boxplot(PlotVar(1:length(Names),:)',Names,'BoxStyle','outline','Width',0.5,'LabelOrientation','inline','LabelVerbosity','majorminor','Orientation','horizontal');
set(h,{'linewidth'},{2})
xlabel('Relative Number of Cells with no Changed Neighbors','FontSize',24)
set(gca,'fontsize', 24);
title(sprintf('Preserved Neighborhood - %s',char(LegendNames(CellType))),'FontSize',30);
saveas(gcf, sprintf('%s Relative Number of Cells with no Changed Neighbors over Time %s.png',savename,char(LegendNames(CellType))));
savefig(sprintf('%s Relative Number of Cells with no Changed Neighbors over Time %s.fig',savename,char(LegendNames(CellType))))
close all



% Cumulative cell divisions
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau,SumDivAll(1:length(VelFieldData(groups(1)).tau),1),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(VelFieldData(i).tau,SumDivAll(1:length(VelFieldData(i).tau),i),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count + 1;
    SEM = SumDivStdAll(1:length(VelFieldData(i).tau),i)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    %errorbar(VelFieldData(i).tau((WSize-1)/2:end),nanmean(VelFieldData(i).pAll((WSize-1)/2:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
    shadedErrorBar(VelFieldData(i).tau,SumDivAll(1:length(VelFieldData(i).tau),i),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('Time in min','FontSize',text_font)
ylabel('Number of Divisions','FontSize',text_font)
%ylim([0.90 1.47]);
xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('Number of Cell Divisions','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Divisions %s.png',savename));
savefig(sprintf('Cell Divisions %s.fig',savename))
close all

% Number of cell divisions normalized to the area covered by cells (in
% units of the control):
data_in = max(SumDivAll,[],1);
data_in = data_in(groups);
% Normalize data to mean of all experiments of the temporal median values 
% of the covered area (complex wording ...):
count = 0;
RelArea = [];
for i = groups
    count = count+1;
    RelArea(count) = nanmean(nanmedian(VelFieldData(i).CoveredAreaAll));
end
% Give relative area in units of the control conditions:
RelArea = RelArea./RelArea(1);
% Normalize data:
data_in = data_in./RelArea;
% SEM calculation:
count = 0;
SEM = [];
for i = groups
    count = count+1;
    SEM(count) = sqrt(size(VelFieldData(i).MSDAll(1,:),2));
end
SEM = max(SumDivStdAll(:,groups),[],1)./SEM;
standard_in = SEM;
Name = LegendNames;
%Name = {'A375', '5% MG','10% MG','15% MG','30% MG'};
Title = 'Number of Cell Divisions';
%Title = 'Spheroid Size A375 70h';
significance = {'','',''};
group_name = {'','',''};
significance = {'','','','','','','',''};
group_name =   {'','','','','','','',''};
significance = {'','','','','','','','',''};
group_name =   {'','','','','','','','',''};
significance = {'','','','','',''};
group_name = {'','','','','',''};
significance = {'','','','','','','','','','','','','',''};
group_name =   {'','','','','','','','','','','','','',''};
%significance = {'','','','',''};
%group_name = {'','','','',''};
%group_name = {'U138','U138','U138'};
text_font  = 24;
n = 5;
n = 6;
n = 9;
n = 14;
n = 8;
gaps = [];
y_axis_label = 'Divisions';
Leg = Name;
plot_inv(data_in, standard_in, Name, Title, significance, group_name, text_font, n, gaps, y_axis_label,Leg)
% save bar plot:
saveas(gcf, sprintf('Cell Divisions %s.png',savename));
savefig(sprintf('Cell Divisions %s.fig',savename))


% Scaling coefficients from MSD - whole time window:
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau,pMeanAll(1:length(VelFieldData(groups(1)).tau),1),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(VelFieldData(i).tau,pMeanAll(1:length(VelFieldData(i).tau),i),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthEast','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count + 1;
    SEM = pStdAll(1:length(VelFieldData(i).tau),i)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    %errorbar(VelFieldData(i).tau((WSize-1)/2:end),nanmean(VelFieldData(i).pAll((WSize-1)/2:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
    shadedErrorBar(VelFieldData(i).tau,pMeanAll(1:length(VelFieldData(i).tau),i),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('\Deltat in min','FontSize',text_font)
ylabel('\alpha','FontSize',text_font)
%ylim([0.90 1.47]);
xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('Scaling Coefficient \alpha','FontSize',text_font+6);
saveas(gcf, sprintf('ScalingCoeffLog %s.png',savename));
savefig(sprintf('ScalingCoeffLog %s.fig',savename))
close all

% Scaling coefficients from MSD - individual time windows:
TimeWin = 6;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau(1:length(pMeanAllTemp)),pMeanAllTemp(1:length(pMeanAllTemp),1,TimeWin),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(VelFieldData(i).tau(1:length(pMeanAllTemp)),pMeanAllTemp(1:length(pMeanAllTemp),i,TimeWin),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthEast','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count + 1;
    SEM = pStdAllTemp(1:length(VelFieldData(i).tau(1:length(pMeanAllTemp))),i,TimeWin)./sqrt(size(VelFieldData(i).MSDAll(1:length(pMeanAllTemp),:),2));
    %errorbar(VelFieldData(i).tau((WSize-1)/2:end),nanmean(VelFieldData(i).pAll((WSize-1)/2:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
    shadedErrorBar(VelFieldData(i).tau(1:length(pMeanAllTemp)),pMeanAllTemp(1:length(pMeanAllTemp),i,TimeWin),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('\Deltat in min','FontSize',text_font)
ylabel('\alpha','FontSize',text_font)
%ylim([0.90 1.47]);
xlim([0 max(VelFieldData(groups(1)).tau(1:length(pMeanAllTemp)))]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('Scaling Coefficient \alpha - 50-70h','FontSize',text_font+6);
saveas(gcf, sprintf('ScalingCoeffLog %dh %s.png',(TimeWin-1)*dt*Sampling/60,savename));
savefig(sprintf('ScalingCoeffLog %dh Windowstart %s.fig',(TimeWin-1)*dt*Sampling/60, savename))
close all

% MSD
figure('units','normalized','outerposition',[0 0 1 1]);
loglog(VelFieldData(groups(1)).tau,PxSize.^2*nanmean(VelFieldData(groups(1)).MSDAll(1:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
hold all
% p = polyfit(log10(VelFieldData(groups(1)).tau),log10(PxSize.^2*nanmean(VelFieldData(groups(1)).MSDAll(1:length(VelFieldData(groups(1)).tau),:),2)),1);
% yfit = polyval(p,log10(VelFieldData(groups(1)).tau));
% loglog(VelFieldData(groups(1)).tau,10.^yfit,'--r','Linewidth',2,'HandleVisibility','off')
%Slope1 = polyval([1.05,0],log10(VelFieldData(groups(1)).tau))-1.95;
%Slope2 = polyval([2,0],log10(VelFieldData(groups(1)).tau))-1.4;
%loglog(VelFieldData(groups(1)).tau,10.^Slope1,'--r','Linewidth',2,'HandleVisibility','off')
%loglog(VelFieldData(groups(1)).tau,10.^Slope2,'--r','Linewidth',2,'HandleVisibility','off')
count = 1;
for i = groups(2:end)
    count = count + 1;
    loglog(VelFieldData(i).tau,PxSize.^2*nanmean(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2, 'Color',color_map(count,:))
    %     p = polyfit(log10(VelFieldData(i).tau),log10(PxSize.^2*nanmean(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2)),1);
    %     yfit = polyval(p,log10(VelFieldData(i).tau));
    %     loglog(VelFieldData(i).tau,10.^yfit,'--r','Linewidth',2,'HandleVisibility','off')
end
legend(LegendNames, 'Location','SouthEast','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count + 1;
    SEM = PxSize.^2*nanstd(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    errorbar(VelFieldData(i).tau,PxSize.^2*nanmean(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
end
xlabel('\Deltat in min','FontSize',text_font)
ylabel('MSD in �m�','FontSize',text_font)
%ylim([0 123]);
xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('MSD','FontSize',text_font+6);
saveas(gcf, sprintf('MSD %s.png',savename));
savefig(sprintf('MSD %s.fig',savename))
close all

% Q - order parameter
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldData(groups(1)).tau,nanmean(VelFieldData(groups(1)).QAll(1:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    semilogx(VelFieldData(i).tau,nanmean(VelFieldData(i).QAll(1:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','SouthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count + 1;
    SEM = nanstd(VelFieldData(i).QAll(1:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    %errorbar(VelFieldData(i).tau,nanmean(VelFieldData(i).QAll(1:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
    shadedErrorBar(VelFieldData(i).tau,nanmean(VelFieldData(i).QAll(1:length(VelFieldData(i).tau),:),2),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('Time in min','FontSize',text_font)
ylabel('Orderparameter','FontSize',text_font)
ylim([0 1]);
xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('Orderparameter','FontSize',text_font+6);
saveas(gcf, sprintf('OrderParameter %s.png',savename));
savefig(sprintf('OrderParameter %s.fig',savename))
close all

% Save Data as new xls sheet:
QMean=NaN(length(nanmean(VelFieldData(1).QAll(1:length(VelFieldData(i).tau),:),2)),length(VelFieldData));
Qsem=NaN(length(nanmean(VelFieldData(1).QAll(1:length(VelFieldData(i).tau),:),2)),length(VelFieldData));
for i = 1:length(VelFieldData)    
    QMean(1:length(nanmean(VelFieldData(i).QAll(1:length(VelFieldData(i).tau),:),2)),i) = ...
        nanmean(VelFieldData(i).QAll(1:length(VelFieldData(i).tau),:),2);
    Qsem(1:length(nanmean(VelFieldData(1).QAll(1:length(VelFieldData(i).tau),:),2)),i) = ...
        nanstd(VelFieldData(i).QAll(1:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
end
writetable(array2table(QMean),'OrderParameterCenter.xlsx')
writetable(array2table(Qsem),'OrderParameterSEMCenter.xlsx')

% Chi - 4 point susceptibility
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldData(groups(1)).tau,nanmean(VelFieldData(groups(1)).ChiAll(1:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    semilogx(VelFieldData(i).tau,nanmean(VelFieldData(i).ChiAll(1:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count + 1;
    SEM = nanstd(VelFieldData(i).ChiAll(1:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    %errorbar(VelFieldData(i).tau,nanmean(VelFieldData(i).ChiAll(1:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
    shadedErrorBar(VelFieldData(i).tau,nanmean(VelFieldData(i).ChiAll(1:length(VelFieldData(i).tau),:),2),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('Time in min','FontSize',text_font)
ylabel('4 Point Suceptibility','FontSize',text_font)
%ylim([0 9.5]);
%ylim([0 4.5]);
xlim([0 max(VelFieldData(groups(1)).tau)])
xlim([0 1300])
xlim([0 t_max]);
set(gca, 'Fontsize',text_font)
title('4 Point Suceptibility','FontSize',text_font+6)
saveas(gcf, sprintf('Chi %s.png',savename));
savefig(sprintf('Chi %s.fig',savename))
close all

% Save Data as new xls sheet:
ChiMean=NaN(length(nanmean(VelFieldData(1).ChiAll(1:length(VelFieldData(i).tau),:),2)),length(VelFieldData));
Chisem=NaN(length(nanmean(VelFieldData(1).ChiAll(1:length(VelFieldData(i).tau),:),2)),length(VelFieldData));
for i = 1:length(VelFieldData)    
    ChiMean(1:length(nanmean(VelFieldData(i).ChiAll(1:length(VelFieldData(i).tau),:),2)),i) = ...
        nanmean(VelFieldData(i).ChiAll(1:length(VelFieldData(i).tau),:),2);
    Chisem(1:length(nanmean(VelFieldData(1).ChiAll(1:length(VelFieldData(i).tau),:),2)),i) = ...
        nanstd(VelFieldData(i).ChiAll(1:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
end
writetable(array2table(ChiMean),'ChiCenter.xlsx')
writetable(array2table(Chisem),'ChiSEMCenter.xlsx')


% Speed in �m/h
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau,PxSize/(dt/60)*nanmean(VelFieldData(groups(1)).Speed(1:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(VelFieldData(i).tau,PxSize/(dt/60)*nanmean(VelFieldData(i).Speed(1:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','SouthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count + 1;
    SEM = PxSize/(dt/60)*nanstd(VelFieldData(i).Speed(1:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    %errorbar(VelFieldData(i).tau,PxSize/(dt/60)*nanmean(VelFieldData(i).Speed(1:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
    shadedErrorBar(VelFieldData(i).tau,PxSize/(dt/60)*nanmean(VelFieldData(i).Speed(1:length(VelFieldData(i).tau),:),2),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('Time in min','FontSize',text_font)
ylabel('Speed in �m/h','FontSize',text_font)
%ylim([0 25.37]);
%ylim([4 17]);
%ylim([1 13]);
%ylim([0 9.37]);
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
xlim([0 t_max])
set(gca, 'Fontsize',text_font)
title('Speed','FontSize',text_font+6);
saveas(gcf, sprintf('Speed %s.png',savename));
savefig(sprintf('Speed %s.fig',savename))
close all

% Root mean squared speed in �m/h
% 3d measurements of astros u138 ln229-> ctl only:
offset = [3,6,0];
tempOffset = [300,0,0]; % in minutes
% else:
offset = zeros(size(groups));
tempOffset = zeros(size(groups));

figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau(tempOffset(1)/dt+1:end)-tempOffset(1),offset(1) + PxSize/(dt/60)*nanmean(VelFieldData(groups(1)).RMSVelAll(tempOffset(1)/dt+1:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(VelFieldData(i).tau(tempOffset(i)/dt+1:end)-tempOffset(i),offset(i) + PxSize/(dt/60)*nanmean(VelFieldData(i).RMSVelAll(tempOffset(i)/dt+1:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count + 1;
    SEM = PxSize/(dt/60)*nanstd(VelFieldData(i).RMSVelAll(tempOffset(i)/dt+1:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(groups(1)).tau),:),2));
    %errorbar(VelFieldData(i).tau,PxSize/(dt/60)*nanmean(VelFieldData(i).RMSVelAll(1:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
    shadedErrorBar(VelFieldData(i).tau(tempOffset(i)/dt+1:end)-tempOffset(i),offset(i) + PxSize/(dt/60)*nanmean(VelFieldData(i).RMSVelAll(tempOffset(i)/dt+1:length(VelFieldData(i).tau),:),2),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)

end
xlabel('Time in min','FontSize',text_font)
ylabel('Speed in �m/h','FontSize',text_font)
%ylim([0 22.37]);
%ylim([12 30]);
%ylim([6 24.37]);
%ylim([3 7]);
%ylim([4 10]);
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
xlim([0 t_max])
set(gca, 'Fontsize',text_font)
title('Speed','FontSize',text_font+6);
saveas(gcf, sprintf('RootMeanSpeed %s.png',savename));
savefig(sprintf('RootMeanSpeed %s.fig',savename))
close all

% Save Data as new xls sheet:
VMean=NaN(length(nanmean(VelFieldData(1).RMSVelAll(1:length(VelFieldData(1).tau),:),2)),length(VelFieldData));
Vsem=NaN(length(nanmean(VelFieldData(1).RMSVelAll(1:length(VelFieldData(1).tau),:),2)),length(VelFieldData));
for i = 1:length(VelFieldData)    
    VMean(1:length(nanmean(VelFieldData(i).RMSVelAll(1:length(VelFieldData(i).tau),:),2)),i) = PxSize/(dt/60)*...
        nanmean(VelFieldData(i).RMSVelAll(1:length(VelFieldData(i).tau),:),2);
    Vsem(1:length(nanmean(VelFieldData(1).RMSVelAll(1:length(VelFieldData(i).tau),:),2)),i) = PxSize/(dt/60)*...
        nanstd(VelFieldData(i).RMSVelAll(1:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
end
writetable(array2table(VMean),'RMS_Center.xlsx')
writetable(array2table(Vsem),'RMS_SEM_Center.xlsx')


% Speed as bar plot:

% F�r die Statistik:
StatData = NaN(length(VelFieldData),100);
for i = 1:length(VelFieldData)    
StatData(i,1:size(VelFieldData(i).RMSVelAll(1,:),2)) = squeeze(PxSize/(dt/60)*nanmedian(VelFieldData(i).RMSVelAll(1:length(VelFieldData(i).tau),:),1));
end
% Save Data as new xls sheet:
writetable(array2table(StatData),'MeanSpeedCenter.xlsx')

[p,tbl,stats] = anova1(StatData');
[c,m,h,gnames] = multcompare(stats,'CType','hsd');


for i = 1:length(VelFieldData)
data_in(i) = nanmean(PxSize/(dt/60)*nanmean(VelFieldData(i).RMSVelAll(1:length(VelFieldData(i).tau),:),1),2);
end
data_in = data_in(groups);
% SEM calculation:
count = 0;
SEM = [];
standard_in = [];
for i = groups
    count = count+1;
    SEM(count) = sqrt(size(VelFieldData(i).MSDAll(1,:),2));
    standard_in(count) = nanstd(PxSize/(dt/60)*nanmean(VelFieldData(i).RMSVelAll(1:length(VelFieldData(i).tau),:),1),[],2)./SEM(count);
end

Name = LegendNames;
%Name = {'A375', '5% MG','10% MG','15% MG','30% MG'};
Title = 'Mean Speed';
%Title = 'Spheroid Size A375 70h';
significance = {'','',''};
group_name = {'','',''};
significance = {'','','','','','','',''};
group_name =   {'','','','','','','',''};
significance = {'','','','','','','','',''};
group_name =   {'','','','','','','','',''};
significance = {'','','','',''};
group_name = {'','','','',''};
significance = {'','','','','','','','','','','','','',''};
group_name =   {'','','','','','','','','','','','','',''};
%significance = {'','','','',''};
%group_name = {'','','','',''};
%group_name = {'U138','U138','U138'};
text_font  = 24;
n = 5;
n = 6;
n = 9;
n = 14;
n = 8;
n = 3;
gaps = [];
y_axis_label = 'Speed in �m/h';
Leg = Name;
plot_inv(data_in, standard_in, Name, Title, significance, group_name, text_font, n, gaps, y_axis_label,Leg)
% save bar plot:
saveas(gcf, sprintf('Mean Speed Bar %s.png',savename));
savefig(sprintf('MeanSpeed Bar %s.fig',savename))


% Correlation length in �m
TimeWin =  [1,t_max/dt];
% Size scaling in units of the velocity velociy correlation function:
LengthScale = (0:size(VelFieldData(2).AutoCorrAll,1)-1)*VelFieldPxSize;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau,smooth(CorrLength(1:length(VelFieldData(groups(1)).tau),groups(1)),15),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(VelFieldData(i).tau,smooth(CorrLength(1:length(VelFieldData(i).tau),i),15),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthEast','FontSize',text_font,'AutoUpdate','off')
count = 0;

for i = groups(1:end)
    count = count+1;
    SEM = [];
    %pos = [];
    for j = 1:length(VelFieldData(i).tau)
    [~,Idx] = min(abs(LengthScale-smooth(CorrLength(j,i),15)));
    %pos(j) = Idx;   
    SEM(j) = AutoCorrStd(j,Idx,i)./sqrt(size(VelFieldData(i).MSDAll(TimeWin(1):TimeWin(2),:),2));
    end
    SEM = 5*SEM.*smooth(CorrLength(1:length(VelFieldData(i).tau),i),15)';
    % errorbar(VelFieldData(i).tau,AutoCorrMean(1:length(VelFieldData(i).tau),Idx,i),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(i,:))
    shadedErrorBar(VelFieldData(i).tau(1:length(VelFieldData(i).tau)),smooth(CorrLength(1:length(VelFieldData(i).tau),i),15)...
        ,SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('Time in min','Fontsize',text_font)
ylabel('Correlation Length in �m','Fontsize',text_font)
%ylim([0 92]);
%ylim([5 23]);
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
xlim([0 t_max])
set(gca, 'Fontsize',text_font)
title('Correlation Length','FontSize',text_font+6);
saveas(gcf, sprintf('CorrelationLength %s.png',savename));
savefig(sprintf('CorrelationLength %s.fig',savename))
close all

% Correlation values at 20�m distance - specified time window:
Dist = 20; % autorcorralation function evaluated at specified distance in �m
TimeWin =  [1,1*t_max/dt];
% Size scaling in units of the velocity velociy correlation function:
LengthScale = (0:size(VelFieldData(2).AutoCorrAll,1)-1)*VelFieldPxSize;
[~,Idx] = min(abs(LengthScale-Dist));
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau(TimeWin(1):TimeWin(2)),AutoCorrMean(TimeWin(1):TimeWin(2),Idx,groups(1)),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(VelFieldData(i).tau(TimeWin(1):TimeWin(2)),AutoCorrMean(TimeWin(1):TimeWin(2),Idx,i),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthEast','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count+1;
    SEM = AutoCorrStd(TimeWin(1):TimeWin(2),Idx,groups(1))./sqrt(size(VelFieldData(i).MSDAll(TimeWin(1):TimeWin(2),:),2));
    % errorbar(VelFieldData(i).tau,AutoCorrMean(1:length(VelFieldData(i).tau),Idx,i),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(i,:))
    shadedErrorBar(VelFieldData(i).tau(TimeWin(1):TimeWin(2)),AutoCorrMean(TimeWin(1):TimeWin(2),Idx,i),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('Time in min','Fontsize',text_font)
ylabel(sprintf('C(\\Deltar = %d �m)',Dist),'Fontsize',text_font)
%ylim([0 92]);
%ylim([5 23]);
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
xlim(TimeWin*dt)
set(gca, 'Fontsize',text_font)
title('Velocity Correlation','FontSize',text_font+6);
saveas(gcf, sprintf('Velocity Correlation at Distance %s.png',savename));
savefig(sprintf('Velocity Correlation at Distance %s.fig',savename))
close all

% V velocity order parameter
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau,nanmean(nanmean(VelFieldData(groups(1)).VTempAll(:,1:length(VelFieldData(groups(1)).tau),:),1),3),...
    'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for j = groups(2:end)
    count = count + 1;
    plot(VelFieldData(j).tau,nanmean(nanmean(VelFieldData(j).VTempAll(:,1:length(VelFieldData(j).tau),:),1),3),...
        'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count + 1;
    SEM = nanstd(nanmean(VelFieldData(i).VTempAll(:,1:length(VelFieldData(i).tau),:),1),[],3)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(groups(1)).tau),:),2));
    errorbar(VelFieldData(i).tau,nanmean(nanmean(VelFieldData(i).VTempAll(:,1:length(VelFieldData(i).tau),:),1),3),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
end
xlabel('Time in min','Fontsize',text_font)
ylabel('Velocity Order Parameter','Fontsize',text_font)
%ylim([0.150 0.45]);
%ylim([0.0 0.60]);
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
xlim([0 t_max])
set(gca,'Fontsize',text_font)
title('Velocity Order Parameter','FontSize',text_font+6);
saveas(gcf, sprintf('VelocityOrderParameter %s.png',savename));
savefig(sprintf('VelocityOrderParameter %s.fig',savename))
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time dependent plots:

% Get different times:
Times = [21:160:600];
Times = 241;
Times = (Sampling * 0)+1;
% Take temporal offset:
offset = 0;%20;

% Take care of constraints imposed by window size:
% WSize = VelFieldData(1).WSize;
% Times(Times <(WSize-1)/2) = [];
% Times(Times>size(VelFieldData(1).AutoCorrAll-2,2)-WSize) = [];

LegendNames = [];
Legends = [];
count = 0;
for i = groups
    count = count + 1;
    LegendNames{count} = VelFieldData(i).Name;
end
for j = 1:length(groups)
    for i = 1:length(Times)
        idx = (j-1)*length(Times) + i;
        Legends{idx} = sprintf('%s t = %0.5g h',LegendNames{j},(Times(i)-1)*3/60 );
    end
end

% Q - order parameter
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldData(groups(1)).tau(1:WSize-1),nanmean(nanmean(VelFieldData(groups(1)).QTempAll(1:end-1,Times(1)-offset:Times(1)+offset,:),2),3),...
    'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
if length(Times>1)
    for i = 2:length(Times)
        semilogx(VelFieldData(groups(1)).tau(1:WSize-1),nanmean(nanmean(VelFieldData(groups(1)).QTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
    end
end
count = 1;
for j = groups(2:end)
    count = count+1;
    for i = 1:length(Times)
        semilogx(VelFieldData(j).tau(1:WSize-1),nanmean(nanmean(VelFieldData(j).QTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
    end
end

count = 0;
for i = groups(1:end)
    count = count + 1;
    for j = 1:length(Times)
        SEM = nanstd(nanmean(VelFieldData(i).QTempAll(1:end-1,Times(j)-offset:Times(j)+offset,:),2),[],3)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
        %errorbar(VelFieldData(i).tau(1:WSize-1),nanmean(nanmean(VelFieldData(i).QTempAll(1:end-1,Times(j)-offset:Times(j)+offset,:),2),3),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
        shadedErrorBar(VelFieldData(i).tau(1:WSize-1),nanmean(nanmean(VelFieldData(i).QTempAll(1:end-1,Times(j)-offset:Times(j)+offset,:),2),3),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)

    end
end
legend(Legends,'Location','SouthWest','Fontsize',text_font)
xlabel('Time in min','Fontsize',text_font)
ylabel('Orderparameter','Fontsize',text_font)
ylim([0 1]);
%xlim([0 dt*2*offset]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca,'Fontsize',text_font)
title('Order Parameter','FontSize',text_font+6);
saveas(gcf, sprintf('OrderParameterTime %s.png',savename));
savefig(sprintf('OrderParameterTime %s.fig',savename))
close all

% Chi - 4 point susceptibility
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldData(groups(1)).tau(1:WSize-1),nanmean(nanmean(VelFieldData(groups(1)).ChiTempAll(1:end-1,Times(1)-offset:Times(1)+offset,:),2),3),...
    'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
if length(Times>1)
    for i = 2:length(Times)
        semilogx(VelFieldData(groups(1)).tau(1:WSize-1),nanmean(nanmean(VelFieldData(groups(1)).ChiTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
    end
end
count = 1;
for j = groups(2:end)
    count = count+1;
    for i = 1:length(Times)
        semilogx(VelFieldData(j).tau(1:WSize-1),nanmean(nanmean(VelFieldData(j).ChiTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
    end
end

count = 0;
for i = groups(1:end)
    count = count + 1;
    for j = 1:length(Times)
        SEM = nanstd(nanmean(VelFieldData(i).ChiTempAll(1:end-1,Times(j)-offset:Times(j)+offset,:),2),[],3)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
        %errorbar(VelFieldData(i).tau(1:WSize-1),nanmean(nanmean(VelFieldData(i).ChiTempAll(1:end-1,Times(j)-offset:Times(j)+offset,:),2),3),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
        shadedErrorBar(VelFieldData(i).tau(1:WSize-1),nanmean(nanmean(VelFieldData(i).ChiTempAll(1:end-1,Times(j)-offset:Times(j)+offset,:),2),3),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)

    end
end
legend(Legends,'Location','NorthWest','Fontsize',text_font)
xlabel('Time in min','Fontsize',text_font)
ylabel('4 Point Suceptibility','Fontsize',text_font)
%ylim([0 1]);
%xlim([0 dt*2*offset]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca,'Fontsize',text_font)
title('4 Point Suceptibility','FontSize',text_font+6);
saveas(gcf, sprintf('ChiTime %s.png',savename));
savefig(sprintf('ChiTime %s.fig',savename))
close all

% MSD in �m� per minute
figure('units','normalized','outerposition',[0 0 1 1]);
loglog(VelFieldData(groups(1)).tau(1:WSize-1),PxSize.^2*nanmean(nanmean(VelFieldData(groups(1)).MSDTempAll(1:end-1,Times(1)-offset:Times(1)+offset,:),2),3),...
    'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
if length(Times>1)
    for i = 2:length(Times)
        loglog(VelFieldData(groups(1)).tau(1:WSize-1),PxSize.^2*nanmean(nanmean(VelFieldData(groups(1)).MSDTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
    end
end
count = 1;
for j = groups(2:end)
    count = count + 1;
    for i = 1:length(Times)
        loglog(VelFieldData(j).tau(1:WSize-1),PxSize.^2*nanmean(nanmean(VelFieldData(j).MSDTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
    end
end
count = 0;
for i = groups(1:end)
    count = count + 1;
    for j = 1:length(Times)
        SEM = PxSize.^2*nanstd(nanmean(VelFieldData(i).MSDTempAll(1:end-1,Times(j)-offset:Times(j)+offset,:),2),[],3)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
        errorbar(VelFieldData(i).tau(1:WSize-1),PxSize.^2*nanmean(nanmean(VelFieldData(i).MSDTempAll(1:end-1,Times(j)-offset:Times(j)+offset,:),2),3),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
    end
end

legend(Legends,'Location','NorthWest','Fontsize',text_font)
xlabel('\Deltat in min','Fontsize',text_font)
ylabel('MSD in �m�','Fontsize',text_font)
%ylim([0 1]);
%xlim([0 dt*2*offset]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca,'Fontsize',text_font)
title('MSD','FontSize',text_font+6);
saveas(gcf, sprintf('MSDTime %s.png',savename));
savefig(sprintf('MSDTime %s.fig',savename))
close all

% Scaling coefficients - temporal
for i = groups
    FitData = PxSize.^2*nanmean(nanmean(VelFieldData(i).MSDTempAll(1:end-1,Times(1)-offset:Times(1)+offset,:),2),3);
    %p = polyfit(log(FitData),log(),1)
    % extract scaling coefficient
    FType = fittype('a*x+b');
    f = fit(log((1:length(FitData)-1)'*dt),log(FitData(1:end-1)),FType)
    uncertainty = confint(f,0.90);
end

% Scaling coefficients - whole time window:
for i = groups
    FitData = PxSize.^2*nanmean(nanmean(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2),3);
    %p = polyfit(log(FitData),log(),1)
    % extract scaling coefficient
    FType = fittype('a*x+b');
    f = fit(log((1:length(FitData)-1)'*dt),log(FitData(1:end-1)),FType)
    uncertainty = confint(f,0.90);
end


%% Implement orientation analysis:

% Create plots of the average vectors for each time track mat:
% Reduce analysis to vectors with average speed of more than maybe 2 �m/h ?
% Same with the quiver plots
% - should adjust for the immobile ln229 and astrocytes.

% Set up probability density functions:
% Von Mises distrubition -> flock:
%   Input:
%   alpha     angles to evaluate pdf at
%   thetahat preferred direction
%   kappa    concentration parameter
pdFlock = @(alpha,kappa,thetahat) ...
    (1/(2*pi*besseli(0,kappa))).* exp(kappa*cos(alpha-thetahat));
% Sum of two Von Mises distributions with their respective prefered
% direction 180degree apart -> stream:
pdStream = @(alpha,kappa,thetahat) ...
    (1/(2*pi*besseli(0,kappa))).* cosh(kappa*cos(alpha-thetahat));
% uniform distribution:
pdUniform = @(data,in) (0*in + 1/(2*pi));

% Distribution type:
DistType = zeros(length(TrackMat),50,4);
% Differences in velocity angles as function of distance:
MeanAngle = NaN(100,length(TrackMat),50,10);
% Images used for analysis of velocity distribution:
im_max = 200;
% Speed constraints -> min 1�m/h
MinSpeed =  0.5;
% Calculate the distance associated with the respective speed -> in pixels:
MinDist = (MinSpeed/PxSize) * im_max * dt/60;
% Number of directional clusters:
NumClusters = 4;
% Domain size:
Area = NaN(length(TrackMat),50000,10);
for i = 1:length(TrackMat)
    i
    for j = 1:length(TrackMat(i).TrackMat)
        for k = 1:length(TrackMat(i).TrackMat{j})
            x = TrackMat(i).TrackMat{j}{1,k}(:,:,2,1);
            % flip y coordinates to adjust the matlab y-axis direction to
            % the more common direction:
            y = flip(TrackMat(i).TrackMat{j}{1,k}(:,:,1,1),2);
            vecX = TrackMat(i).TrackMat{j}{1,k}(:,:,2,im_max ) - TrackMat(i).TrackMat{j}{1,k}(:,:,2,1);
            % invert direction of differences as the matlab counting for
            % the y axis goes into -y direction
            vecY = -TrackMat(i).TrackMat{j}{1,k}(:,:,1,im_max ) + TrackMat(i).TrackMat{j}{1,k}(:,:,1,1);
            
            [xq,yq] = meshgrid(1:ImSize(1),1:ImSize(2));
            VecInX = interp2(x',y',vecX',xq,yq);
            VecInY = interp2(x',y',vecY',xq,yq);
           
             
            
            % Speed:
            sp = sqrt(vecX.^2 + vecY.^2);
            sp = max(sp(:))*PxSize * (60/dt) /im_max;            
%             tic
%             LIC = plotvfield(VecInX,VecInY,1,1,jet,[1 1 size(VecInX,2) size(VecInX,1)],'r');
%             toc
%             figure('units','normalized','outerposition',[0 0 1 1]);
%             imshow(LIC);
%             hold on
%             colormap(jet);
%             c = colorbar('FontSize',24);
%             %c.Label.String = 'Speed in �m/h';
%             %set(c,'YLim',[0,1])
%             ylabel(c,'Speed in �m/h','FontSize',24);
%             caxis([0,sp])
%             %c.TickLabels = linspace(0,1,numel(c.Ticks));
%             quiver(xq(1:30:end,1:30:end),yq(1:30:end,1:30:end),VecInY(1:30:end,1:30:end),VecInX(1:30:end,1:30:end),2,'w','lineWidth',1.2)
%             saveas(gcf,'LIC Astrocytes V3 05 - 60-70h.png');
%             savefig('LIC Astrocytes V3 05 - 60-70h.fig')
            
            vecY(vecY == 0) = NaN;
            vecX(vecX == 0) = NaN;
%             % Create quiver plot:
%             figure('units','normalized','outerposition',[0 0 1 1]);
%             quiver(x,y,vecX,vecY,2,'r','lineWidth',1.2)

            vecX = vecX(:);
            vecY = vecY(:);            
            % Clustering of the vectorfield in 4 different directions:
            [ClustIdx] = kmeans([vecX,vecY], NumClusters,'Distance','cosine');
            
            % Remove NaNs:
            Idx = isnan(vecX) | isnan(vecY);
            vecX(Idx) = [];
            vecY(Idx) = [];
            ClustIdx(Idx) = -1;
            x = x(:);
            y = y(:);
            x(Idx) = [];
            y(Idx) = [];
            % eliminate the non-moving and NaN parts in the clusters:
            % Remove non-moving parts:
            Idx = find(vecnorm([vecX,vecY]')<MinDist);
            vecX(Idx) = [];
            vecY(Idx) = [];
            ClustIdx(Idx) = -1;
            x(Idx) = [];
            y(Idx) = [];
            
            % Get size of domains:
            for m = 1:NumClusters
                Binary = false(size(TrackMat(i).TrackMat{j}{1,k},1),size(TrackMat(i).TrackMat{j}{1,k},2));
                Binary(ClustIdx == m) = true;
                AreaTemp = regionprops(Binary,'Area');
                pos = find(isnan(Area(i,:,k)) == 1,1);
                Area(i,pos:pos+length(AreaTemp)-1,k) = [AreaTemp.Area];
            end
            
            %             % Plot clusters:
            %             ClustSize = [];
            %             for m = 1:max(ClustIdx)
            %                 ClustSize(m) = sum(ClustIdx == m);
            %             end
            %             [ClustSize,MaxIdx] = sort(ClustSize,'descend');
            %             % Plot quiver plot with cluster colorcoding for the largest
            %             % five clusters:
            %             colors = [1,0,0; 0,1,0; 0,0,1; 0.5,0.5,0.5; 0,0,0];
            %             figure;
            %             hold all
            %             for m = 1:max(ClustIdx)
            %                 xTemp = x(ClustIdx == m);
            %                 yTemp = y(ClustIdx == m);
            %                 vecXTemp = vecX(ClustIdx == m);
            %                 vecYTemp = vecY(ClustIdx == m);
            %                 quiver(xTemp,yTemp,vecXTemp,vecYTemp,2,'color',colors(m,:),'lineWidth',1.2)
            %             end
            
            % create angle of velocity field from 0 to 360 degree:
            % Reference (along x-direction) and normal vector:
            RefVec = [1,0,0];
            NormalVec = [0,0,1];
            xCoord = cross(RefVec.*ones(length(vecX(:)),1),[vecX(:),vecY(:),zeros(length(vecX(:)),1)],2);
            c = sign(dot(xCoord,NormalVec.*ones(length(vecX(:)),1),2)).*vecnorm(xCoord,2,2);
            AngleTotal = atan2d(c,dot(RefVec.*ones(length(vecX(:)),1),[vecX(:),vecY(:),zeros(length(vecX(:)),1)],2));
            % transform angle into radians:
            AngleTotal = 2*pi*AngleTotal./360;
            % Shift angles to 0->2pi:
            AngleTotal = AngleTotal+pi;
            % remove NaNs:
            AngleTotal(isnan(AngleTotal)) = [];
            % Get binning for calculating probability distribution:
            %figure;
            [Num,Alpha] = hist(AngleTotal,72);
            % Censoring Interval:
            censInt = Alpha(2)-Alpha(1);
            % Round to multiples of the censoring intervall:
            censData =  1+round((AngleTotal-Alpha(1))/censInt);
            % make tails of distributions to fit:
            censData(censData == 0) = 1;
            censData(censData > length(Num)) = length(Num);
            censData = Alpha(censData)';
            %censData = [censData-censInt/2 censData+censInt/2];
            
            if ~isempty(AngleTotal) && length(Idx)/(length(Idx) + length(vecX))< 1/2 % do not do calculations if immobile fraction is too large.
                % Fit distribution to data:
                [Num,Alpha] = histcounts(AngleTotal,72,'Normalization','pdf');
                % flocking motion:
                ft = fittype(@(kappa,thetahat,alpha) (1/(2*pi*besseli(0,kappa))).* exp(kappa*cos(alpha-thetahat)), ...
                    'independent', 'alpha');
                [curve, goodness] = fit(Alpha(1:length(Num))', Num', ft, 'StartPoint',[min(Alpha(Num == max(Num))),1]);
                RFlock = goodness.rsquare;
                KFlock = curve.kappa;
                % streaming motion:
                ft = fittype(@(kappa,thetahat,alpha) (1/(2*pi*besseli(0,kappa))).*  cosh(kappa*cos(alpha-thetahat)), ...
                    'independent', 'alpha');
                [curve, goodness] = fit(Alpha(1:length(Num))', Num', ft, 'StartPoint',[min(Alpha(Num == max(Num))),1]);
                RStream = goodness.rsquare;
                KStream = curve.kappa;
                % random motion (swarm):
                [curve, goodness] = polyfit(Alpha(1:length(Num))', Num', 0);
                RSwarm = goodness.R;
            end
            
            % Behavior:
            if ~isempty(Idx) && length(Idx)/(length(Idx) + length(vecX))> 1/2 % immobile fraction
                DistType(i,k,4) = DistType(i,k,4) + 1;
            elseif RFlock > RStream && KFlock(1) >0.5
                DistType(i,k,1) = DistType(i,k,1) + 1;
            elseif RFlock <= RStream  && KStream(1) >0.5
                DistType(i,k,2) = DistType(i,k,2) + 1;
            else
                DistType(i,k,3) = DistType(i,k,3) + 1;
            end
            
            %             figure;
            %             h = histogram(censData,72,'Normalization','pdf');
            %             hold on
            %             plot(Alpha(1:length(Num))', Num','ro')
            %             plot(curve,Alpha(1:length(Num))', Num')
            %             hold all
            %             plot(Alpha,pdFlock(Alpha,phatFlock(1),phatFlock(2)),'r-')
            %             plot(Alpha,pdStream(Alpha,phatStream(1),phatStream(2)),'g-')
            % Get maximum likelihood estimates for flocking or streaming:
            %             [phatFlock,pciFlock] = mle(censData,'pdf',pdFlock,'Start',[min(Alpha(Num == max(Num))),1]);
            %             [phatStream,pciStream] = mle(censData,'pdf',pdStream,'Start',[min(Alpha(Num == max(Num))),1]);
            %[phatUniform,pciUniform] = mle(censData,'pdf',pdUniform,'Start',[1]);
            %             figure;
            %             histogram(censData,360,'Normalization','pdf')
            %             hold all
            %             plot(Alpha,pdFlock(Alpha,phatFlock(1),phatFlock(2)),'r-')
            %             plot(Alpha,pdStream(Alpha,phatStream(1),phatStream(2)),'g-')
            %plot(Alpha,pdUniform(Alpha,phatUniform(1)),'m-')
            % Get log-likelihood values of both distributions using the optimal
            % values derived above:
            %             LHoodFlock = -sum(log(pdFlock(Alpha,phatFlock(1),phatFlock(2))));
            %             LHoodStream = -sum(log(pdStream(Alpha,phatStream(1),phatStream(2))));
            %             % Behavior:
            %             if LHoodFlock > LHoodStream && phatFlock(1) >0.5
            %                 DistType(i,j,k) = 1;
            %             elseif LHoodFlock <= LHoodStream  && phatStream(1) >0.5
            %                 DistType(i,j,k) = 2;
            %             else
            %                 DistType(i,j,k) = 3;
            %             end
            
            % Plot similarity of velocity direction as a function of distance:
            % Distance in angles:
            dSpace = pdist2([x,y],[x,y]);
            dAngle = pdist2(AngleTotal,AngleTotal,@circ_dist2);
            dSpace = dSpace(:);
            dAngle = dAngle(:);
            % bin spatial data:
            [binIdx,bins] = discretize(dSpace,50);
            % bin associated angular data:           
            if ~(length(Idx)/(length(Idx) + length(vecX))> 1/2) % calculate angles only if a significant motile fraction is around
                for m = 1:max(binIdx)
                    MeanAngle(m,i,j,k) = nanmean(abs(dAngle(binIdx == m)));
                end
            end

            %figure;
            %plot(1./MeanAngle(2:end,i,j,k))
        end
    end
end

t = 2;
CType = 3;
for t = 1:7
    figure;
    plot(nanmean(MeanAngle(2:end,CType,:,t),3))
end

CType = 3;
figure;
hold on
plot(nanmean(360/(2*pi) * MeanAngle(2:end,CType,:,1),3),'r')
plot(nanmean(360/(2*pi) * MeanAngle(2:end,CType,:,7),3),'b')


% Get mean domain size:
MinDomainSize = 20; % entspricht rund 2000�m�
% For statistics:
MArea = Area;
MArea(Area<MinDomainSize) = NaN;
% time evolution:
for i = 1:3
    i
    [h,p] = ttest2(MArea(i,:,1),MArea(i,:,2))
    [h,p] = ttest2(MArea(i,:,1),MArea(i,:,3))
    [h,p] = ttest2(MArea(i,:,1),MArea(i,:,4))
    [h,p] = ttest2(MArea(i,:,1),MArea(i,:,5))
    [h,p] = ttest2(MArea(i,:,1),MArea(i,:,6))
    [h,p] = ttest2(MArea(i,:,1),MArea(i,:,7))
end
% groups at specific time:
t = 1;
[h,p] = ttest2(MArea(1,:,t),MArea(2,:,t))
[h,p] = ttest2(MArea(1,:,t),MArea(3,:,t))
[h,p] = ttest2(MArea(2,:,t),MArea(3,:,t))

MArea = Area;
MArea(Area<MinDomainSize) = NaN;
% Get sample size:
SampleSize = squeeze(sum(~isnan(MArea),2));
% Get mean domain size:
% teste 20,10,1 -> am besten 10   <- entspricht rund 2000/1000/100�m�
StdArea = squeeze(nanstd(MArea,[],2));
MArea = squeeze(nanmean(MArea,2));


% Pixelsize in �m/"pixel (entry)" in the tracking matrix
TrackMatPxSize = mean((ImSize./size(TrackMat(1).TrackMat{1}{1,1}(:,:,1,1),1:2)))*PxSize;
MArea = MArea * TrackMatPxSize.^2;
StdArea = StdArea * TrackMatPxSize.^2;

Times = [5,15,25,35,45,55,65]; % times analyzed in hours -> centered
Times = [5,15,25,35,45,55]; % times analyzed in hours -> centered
% Plot Domain Size evolution:
figure('units','normalized','outerposition',[0 0 1 1]);
plot(Times,MArea(1,1:length(Times)),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
hold all
% plot(Times,MArea(1,1:length(Times)),'o','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(Times,MArea(i,1:length(Times)),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(count,:))
    % plot(Times,MArea(i,1:length(Times)),'o','MarkerSize',5,'LineWidth',2, 'Color',color_map(count,:))
    %     p = polyfit(log10(VelFieldData(i).tau),log10(PxSize.^2*nanmean(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2)),1);
    %     yfit = polyval(p,log10(VelFieldData(i).tau));
    %     loglog(VelFieldData(i).tau,10.^yfit,'--r','Linewidth',2,'HandleVisibility','off')
end
legend(LegendNames, 'Location','SouthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
SEM = StdArea./sqrt(SampleSize);
for i = groups(1:end)
    count = count + 1;
    %errorbar(Times,MArea(i,1:length(Times)),SEM(i,1:length(Times)),'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
    shadedErrorBar(Times,MArea(i,1:length(Times)),SEM(i,1:length(Times)),'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
    
end
xlabel('Time in h','FontSize',text_font)
ylabel('Domain Size in �m�','FontSize',text_font)
%ylim([0 123]);
%xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('Collective Domain Size','FontSize',text_font+6);
saveas(gcf, sprintf('Domain Size %s.png',savename));
savefig(sprintf('Domain Size %s.fig',savename))
close all



% Normalize distribution type:
DistType = DistType./sum(DistType,3);

Times = [5,15,25,35,45,55,65]; % times analyzed in hours -> centered
Times = [5,15,25,35,45,55]; % times analyzed in hours -> centered
% Plot flocking motion:
figure('units','normalized','outerposition',[0 0 1 1]);
plot(Times,DistType(1,1:length(Times),1),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(Times,DistType(i,1:length(Times),1),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(count,:))
    %     p = polyfit(log10(VelFieldData(i).tau),log10(PxSize.^2*nanmean(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2)),1);
    %     yfit = polyval(p,log10(VelFieldData(i).tau));
    %     loglog(VelFieldData(i).tau,10.^yfit,'--r','Linewidth',2,'HandleVisibility','off')
end
legend(LegendNames, 'Location','SouthEast','FontSize',text_font,'AutoUpdate','off')
% count = 0;
% for i = groups(1:end)
%     count = count + 1;
%     SEM = nanstd(DistType(i,1:length(Times),1),[],2)./sqrt(size(TrackMat(i).TrackMat,1));
%     errorbar(Times,DistType(i,1:length(Times),1),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
% end
xlabel('Time in h','FontSize',text_font)
ylabel('% of Flocking Motion','FontSize',text_font)
%ylim([0 123]);
%xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('Flocking Motion','FontSize',text_font+6);
saveas(gcf, sprintf('Flocking Motion %s.png',savename));
savefig(sprintf('Flocking Motion %s.fig',savename))
close all

% Plot Streaming motion:
figure('units','normalized','outerposition',[0 0 1 1]);
plot(Times,DistType(1,1:length(Times),2),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(Times,DistType(i,1:length(Times),2),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(count,:))
    %     p = polyfit(log10(VelFieldData(i).tau),log10(PxSize.^2*nanmean(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2)),1);
    %     yfit = polyval(p,log10(VelFieldData(i).tau));
    %     loglog(VelFieldData(i).tau,10.^yfit,'--r','Linewidth',2,'HandleVisibility','off')
end
legend(LegendNames, 'Location','SouthEast','FontSize',text_font,'AutoUpdate','off')
% count = 0;
% for i = groups(1:end)
%     count = count + 1;
%     SEM = nanstd(DistType(i,1:length(Times),2),[],2)./sqrt(size(TrackMat(i).TrackMat,1));
%     errorbar(Times,DistType(i,1:length(Times),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
% end
xlabel('Time in h','FontSize',text_font)
ylabel('% of Streaming Motion','FontSize',text_font)
%ylim([0 123]);
%xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('Streaming Motion','FontSize',text_font+6);
saveas(gcf, sprintf('Streaming Motion %s.png',savename));
savefig(sprintf('Streaming Motion %s.fig',savename))
close all

% Plot Swarming motion:
figure('units','normalized','outerposition',[0 0 1 1]);
plot(Times,DistType(1,1:length(Times),3),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(Times,DistType(i,1:length(Times),3),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(count,:))
    %     p = polyfit(log10(VelFieldData(i).tau),log10(PxSize.^2*nanmean(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2)),1);
    %     yfit = polyval(p,log10(VelFieldData(i).tau));
    %     loglog(VelFieldData(i).tau,10.^yfit,'--r','Linewidth',2,'HandleVisibility','off')
end
legend(LegendNames, 'Location','SouthEast','FontSize',text_font,'AutoUpdate','off')
% count = 0;
% for i = groups(1:end)
%     count = count + 1;
%     SEM = nanstd(DistType(i,1:length(Times),3),[],2)./sqrt(size(TrackMat(i).TrackMat,1));
%     errorbar(Times,DistType(i,1:length(Times),3),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
% end
xlabel('Time in h','FontSize',text_font)
ylabel('% of Swarming Motion','FontSize',text_font)
%ylim([0 123]);
%xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('Swarming Motion','FontSize',text_font+6);
saveas(gcf, sprintf('Swarming Motion %s.png',savename));
savefig(sprintf('Swarming Motion %s.fig',savename))
close all

% Plot immobile cells:
figure('units','normalized','outerposition',[0 0 1 1]);
plot(Times,DistType(1,1:length(Times),4),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(Times,DistType(i,1:length(Times),4),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(count,:))
    %     p = polyfit(log10(VelFieldData(i).tau),log10(PxSize.^2*nanmean(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2)),1);
    %     yfit = polyval(p,log10(VelFieldData(i).tau));
    %     loglog(VelFieldData(i).tau,10.^yfit,'--r','Linewidth',2,'HandleVisibility','off')
end
legend(LegendNames, 'Location','SouthEast','FontSize',text_font,'AutoUpdate','off')
% count = 0;
% for i = groups(1:end)
%     count = count + 1;
%     SEM = nanstd(DistType(i,1:length(Times),4),[],2)./sqrt(size(TrackMat(i).TrackMat,1));
%     errorbar(Times,DistType(i,1:length(Times),4),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
% end
xlabel('Time in h','FontSize',text_font)
ylabel('% of Immobile Cells','FontSize',text_font)
%ylim([0 123]);
%xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('Immobile Cells','FontSize',text_font+6);
saveas(gcf, sprintf('Immobile Cells %s.png',savename));
savefig(sprintf('Immobile Cells %s.fig',savename))
close all


% Make stacked bar plots:
CellType = 'U138';
PlotVar = squeeze(DistType(5,1:7,:));

% shape factor over density:
h=figure('units','normalized','outerposition',[0 0 1 1]);
set(h,'DefaultTextFontSize',24)
hold all
%h = boxplot(flip(Young,2),{'10�M Blebbistatin', '40�M Y-27632','U138','5�M Blebbistatin','20�M Y-27632','LN229'},'BoxStyle','outline','Width',0.5,'LabelOrientation','inline','LabelVerbosity','majorminor','Orientation','horizontal');
h = bar(PlotVar,'stacked','FaceColor','flat');
% bar colors:
h(1).CData = [1,1,1];
h(2).CData = 3/4*[1,1,1];
h(3).CData = 2/4*[1,1,1];
h(4).CData = 1/4*[1,1,1];
set(h,{'linewidth'},{2})
xlabel('Time in h','FontSize',24)
xTickLab = cell(1,length(Times)+1);
xTickLab(2:length(Times)+1) = num2cell(Times);
xticklabels(xTickLab)
set(gca,'fontsize', 24);
legend({'Flock','Stream','Swarm','Immobile'},'FontSize',text_font,'Location','NorthEast')
title(sprintf('Migration Type - %s', CellType),'FontSize',30);
saveas(gcf, sprintf('Migration Type %s.png',CellType));
savefig(sprintf('Migration Type %s.fig',CellType))
close all





% compare vector fields:

% generate comparison value for random uniform distribution of angles from
% 0 to 90 degree (as many orientations are bound to that angle):
MaxVar = [];
for i = 1:10000
    MaxVar(i) = circ_std(0.5*pi*rand(8374,1));
end
MaxVar = max(MaxVar)
MaxVar = 0.4563;


% Angular variance and mean between cell orientation field and velocity
% field:
% Calculate Angular variance and mean:
InStdMean = NaN(t_max/dt,size(VecFieldData,2));
InStdStd = NaN(t_max/dt,size(VecFieldData,2));
InMeanMean = NaN(t_max/dt,size(VecFieldData,2));
InMeanStd = NaN(t_max/dt,size(VecFieldData,2));
SampleSize = [];
for j = 1:size(VecFieldData,2)
    j
    SampleSize(j) = size(VecFieldData(j).AngleVecFields,1);
    for i = 1:t_max/dt
        tmpMean = [];
        tmpStd = [];
        for k = 1:size(VecFieldData(j).AngleVecFields,1)
          tmpMean(k) = circ_mean(2*pi*VecFieldData(j).AngleVecFields{k}(:,i)/360);
          tmpStd(k) = circ_std(2*pi*VecFieldData(j).AngleVecFields{k}(:,i)/360);          
        end
        %tmpMean(tmpMean==0) = [];
        %tmpStd(tmpStd==0) = [];        
        InStdMean(i,j) = mean(tmpStd);
        InStdStd(i,j) = std(tmpStd);
        InMeanMean(i,j) = (360/(2*pi))*circ_mean(tmpMean);
        InMeanStd(i,j) = (360/(2*pi))*std(tmpMean);
    end
end

% plot angular standard deviation of relative cell orientation vs cell speed:
figure('units','normalized','outerposition',[0 0 1 1]);
plot(0:dt:t_max -dt,InStdMean(:,1),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
hold all
% plot(Times,MArea(1,1:length(Times)),'o','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(0:dt:t_max -dt,InStdMean(:,i),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(count,:)) 
end
legend(LegendNames, 'Location','SouthWest','FontSize',text_font,'AutoUpdate','off')
plot([0,t_max-dt],[MaxVar,MaxVar],'--','MarkerSize',5,'LineWidth',2, 'Color',[0.5,0.5,0.5])
count = 0;
SEM = InStdStd./sqrt(SampleSize);
for i = groups(1:end)
    count = count + 1;
    %errorbar(Times,MArea(i,1:length(Times)),SEM(i,1:length(Times)),'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
    shadedErrorBar(0:dt:t_max -dt,InStdMean(:,i),SEM(:,i),'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
    
end
xlabel('Time in min','FontSize',text_font)
ylabel('Angular Variance','FontSize',text_font)
%ylim([0 123]);
%xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('Angular Variance of Angle Between Cell Orientation and Velocity','FontSize',text_font+6);
saveas(gcf, sprintf('Angular Variance of Angle Between Cell Orientation and Velocity %s.png',savename));
savefig(sprintf('Angular Variance of Angle Between Cell Orientation and Velocity %s.fig',savename))
close all

% plot angular mean of relative cell orientation vs cell speed:
figure('units','normalized','outerposition',[0 0 1 1]);
plot(0:dt:t_max -dt,InMeanMean(:,1),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
hold all
% plot(Times,MArea(1,1:length(Times)),'o','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(0:dt:t_max -dt,InMeanMean(:,i),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(count,:)) 
end
legend(LegendNames, 'Location','SouthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
SEM = InMeanStd./sqrt(SampleSize);
for i = groups(1:end)
    count = count + 1;
    %errorbar(Times,MArea(i,1:length(Times)),SEM(i,1:length(Times)),'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
    shadedErrorBar(0:dt:t_max -dt,InMeanMean(:,i),SEM(:,i),'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
    
end
xlabel('Time in min','FontSize',text_font)
ylabel('Angular Mean','FontSize',text_font)
%ylim([0 123]);
%xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('Angular Mean of Angle Between Cell Orientation and Velocity','FontSize',text_font+6);
saveas(gcf, sprintf('Angular Mean of Angle Between Cell Orientation and Velocity %s.png',savename));
savefig(sprintf('Angular Mean of Angle Between Cell Orientation and Velocity %s.fig',savename))
close all


% Angular variance of cell orientation field
% field:
% Calculate Angular variance and mean:
InStdMean = NaN(t_max/dt,size(VecFieldData,2));
InStdStd = NaN(t_max/dt,size(VecFieldData,2));
SampleSize = [];
for j = 1:size(VecFieldData,2)
    j
    SampleSize(j) = size(VecFieldData(j).AngleVecFields,1);
    for i = 1:1350%t_max/dt
        tmpStd = [];
        for k = 1:size(VecFieldData(j).AngleOrientation,1)          
          tmpStd(k) = circ_std(2*pi*VecFieldData(j).AngleOrientation{k}(:,i)/360);          
        end    
               
        tmpStd(tmpStd==0) = [];  
        InStdMean(i,j) = mean(tmpStd);
        InStdStd(i,j) = std(tmpStd);
        
        
    end
end

% plot angular standard deviation of relative cell orientation vs cell speed:
figure('units','normalized','outerposition',[0 0 1 1]);
plot(0:dt:t_max -dt,InStdMean(:,1),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
hold all
% plot(Times,MArea(1,1:length(Times)),'o','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(0:dt:t_max -dt,InStdMean(:,i),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(count,:)) 
end
legend(LegendNames, 'Location','SouthWest','FontSize',text_font,'AutoUpdate','off')
plot([0,t_max-dt],[MaxVar,MaxVar],'--','MarkerSize',5,'LineWidth',2, 'Color',[0.5,0.5,0.5])
count = 0;
SEM = InStdStd./sqrt(SampleSize);
for i = groups(1:end)
    count = count + 1;
    %errorbar(Times,MArea(i,1:length(Times)),SEM(i,1:length(Times)),'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
    shadedErrorBar(0:dt:t_max -dt,InStdMean(:,i),SEM(:,i),'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
    
end
xlabel('Time in min','FontSize',text_font)
ylabel('Angular Variance','FontSize',text_font)
%ylim([0 123]);
%xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('Angular Variance of Cell Orientation','FontSize',text_font+6);
saveas(gcf, sprintf('Angular Variance of Cell Orientation %s.png',savename));
savefig(sprintf('Angular Variance of Cell Orientation %s.fig',savename))
close all



% Angular variance of velocity field
% field:
% Calculate Angular variance and mean:
InStdMean = NaN(t_max/dt,size(VecFieldData,2));
InStdStd = NaN(t_max/dt,size(VecFieldData,2));
SampleSize = [];
for j = 1:size(VecFieldData,2)
    j
    SampleSize(j) = size(VecFieldData(j).AngleVelField,1);
    for i = 1:t_max/dt
        tmpStd = [];
        for k = 1:size(VecFieldData(j).AngleOrientation,1)          
          tmpStd(k) = circ_std(2*pi*VecFieldData(j).AngleVelField{k}(:,i)/360);          
        end    
        tmpStd(tmpStd==0) = [];  
        InStdMean(i,j) = mean(tmpStd);
        InStdStd(i,j) = std(tmpStd);
    end
end

% plot angular standard deviation of relative cell orientation vs cell speed:
figure('units','normalized','outerposition',[0 0 1 1]);
plot(0:dt:t_max -dt,InStdMean(:,1),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
hold all
% plot(Times,MArea(1,1:length(Times)),'o','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(0:dt:t_max -dt,InStdMean(:,i),'-','MarkerSize',5,'LineWidth',2, 'Color',color_map(count,:)) 
end
legend(LegendNames, 'Location','SouthWest','FontSize',text_font,'AutoUpdate','off')
plot([0,t_max-dt],[MaxVar,MaxVar],'--','MarkerSize',5,'LineWidth',2, 'Color',[0.5,0.5,0.5])
count = 0;
SEM = InStdStd./sqrt(SampleSize);
for i = groups(1:end)
    count = count + 1;
    %errorbar(Times,MArea(i,1:length(Times)),SEM(i,1:length(Times)),'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
    shadedErrorBar(0:dt:t_max -dt,InStdMean(:,i),SEM(:,i),'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
    
end
xlabel('Time in min','FontSize',text_font)
ylabel('Angular Variance','FontSize',text_font)
%ylim([0 123]);
%xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('Angular Variance of Velocity Field','FontSize',text_font+6);
saveas(gcf, sprintf('Angular Variance of Velocity Field %s.png',savename));
savefig(sprintf('Angular Variance of Velocity Field %s.fig',savename))
close all




















% make the following plots:
% Best made using the AnalyseVectorField script

% plots of speed over cellular orientation (2D plot)

% below are some raw templates listed:
% Plot speed over angle as 2D heatmap
Speed = squeeze((VelField(:,:,1,1:end-1).^2 + VelField(:,:,2,1:end-1).^2).^(0.5));
% Take normalized speed:
%Speed = squeeze(( (VelField(:,:,1,1:end-1)-nanmean(nanmean(VelField(:,:,1,1:end-1),2),1)).^2 + (VelField(:,:,2,1:end-1)-nanmean(nanmean(VelField(:,:,1,1:end-1),2),1)).^2).^(0.5));
Speed = reshape(Speed, size(Speed,1)*size(Speed,2),size(Speed,3));


TempSpeed = Speed(:,1025:1275);
TempAngle = AngleTemp(:,1025:1275);
% Bin data:
%ptsSpeed = linspace(min(TempSpeed(:)),max(Speed(:)),50);
ptsSpeed = linspace(min(TempSpeed(:)),5,50);
%ptsSpeed = linspace(8,16,50);
ptsAngle = linspace(0,90,90);
N = histcounts2(TempSpeed(:),TempAngle(:),ptsSpeed,ptsAngle);
% plot fixed speed as function of angle:
figure; plot(N(1:5:end,:)'./N(1:5:end,1)')

% plot fixed angle as function of speed:
figure; plot(N(:,1:30:end)./N(1,1:30:end))

figure;
imagesc(ptsAngle,ptsSpeed,N)
hold on
colormap(jet)
colorbar
set(gca,'Xlim',ptsAngle([1 end]),'YLim',ptsSpeed([1 end]),'YDir','normal')


ImNum = 1184;
figure;
polarhistogram(2*pi*(AngleTemp(:,ImNum))/360,10,'Normalization', 'probability')
























