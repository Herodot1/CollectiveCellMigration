clearvars -except VelFieldData VecFieldData TrackMat CellDivData;
mPath = mfilename('fullpath');
Idx = max(strfind(mPath,filesep));
mPath = mPath(1:Idx);
addpath(mPath)

% add path for helper functions:
addpath(strcat(mPath,filesep,'HelperFunctions'))
% load parameters
[dt,t_max,Sampling,ImSize,ImPhysSize,VelFieldSize,text_font,groups,...
    savename,color_map] = ParameterFunction;
% Calculate some derived quantities
% pixel size of image in µm/px:
PxSize = ImPhysSize/ImSize;
% Area of image in mm^2:
ImPhysArea = ImPhysSize(1)*ImPhysSize(2)./10^6;
% Pixelsize of velocity field in µm/px:
VelFieldPxSize = mean((ImSize./VelFieldSize))*PxSize;
% Windowsize:
WSize = VelFieldData(1).WSize;
% get legend names:
LegendNames = [];
count = 0;
for i = groups
    count = count + 1;
    LegendNames{count} = VelFieldData(i).Name;
    % Make sure the time frame is set right:
    VelFieldData(i).tau = dt*[0:t_max./dt]';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time evolution of some critical parameters:
% Create scaling coefficient as function of used lag-time for the whole
% time window:
[pMeanAll,pStdAll] = MSDScalingCoeff(VelFieldData);

% Create scaling coefficient as function of used lag-time for different
% time windows only:
[pMeanAllTemp,pStdAllTemp] = MSDScalingCoeffTemp(VelFieldData,Sampling,dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make actual plots
% Plot nearest neighbor stuff:
% Agglomerate Variable:
PlotVar = NaN(length(groups),1700);
count = 0;
for i = groups
    count = count + 1;    
    PlotVar(count,1:size(VelFieldData(i).NumNewNeighborsAll,1)) = 1-VelFieldData(i).NumNewNeighborsAll(:,1);    
end
%Do the plotting: 
h=figure('units','normalized','outerposition',[0 0 1 1]);
set(h,'DefaultTextFontSize',24)
hold all
h = boxplot(PlotVar',LegendNames,'BoxStyle','outline','Width',0.5,'LabelOrientation','inline','LabelVerbosity','majorminor','Orientation','horizontal');
set(h,{'linewidth'},{2})
xlabel('Relative Number of Cells with no Changed Neighbors','FontSize',24)
set(gca,'fontsize', 24);
title('Preserved Neighborhood','FontSize',30);
saveas(gcf, sprintf('%s Relative Number of Cells with no Changed Neighbors.png',savename));
savefig(sprintf('%s Relative Number of Cells with no Changed Neighbors.fig',savename))
close all


% Calculate cumulative cell division number:
[MeanDivs,StdDivs,DivTimes] = CumCellDivs(CellDivData,t_max,dt,groups,ImPhysArea);
% Plot cell divisions
figure('units','normalized','outerposition',[0 0 1 1]);
plot(DivTimes,MeanDivs(:,1),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(DivTimes,MeanDivs(:,i),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count + 1;
    SEM = StdDivs(:,i)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    shadedErrorBar(DivTimes,MeanDivs(:,i),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('Time in min','FontSize',text_font)
ylabel('Number of Divisions per mm²','FontSize',text_font)
xlim([0 t_max]);
set(gca, 'Fontsize',text_font)
title('Number of Cell Divisions','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Divisions %s.png',savename));
savefig(sprintf('Cell Divisions %s.fig',savename))
close all

% Plot Speed around divisions vs. rest of monolayer:
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau,movmean(nanmean(CellDivData(groups(1)).MeanSpeedDiv(1:length(VelFieldData(groups(1)).tau),:)./CellDivData(groups(1)).MeanSpeedNonDiv(1:length(VelFieldData(groups(1)).tau),:),2),21,'omitnan'),...
    'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    plot(VelFieldData(i).tau,movmean(nanmean(CellDivData(i).MeanSpeedDiv(1:length(VelFieldData(i).tau),:)./CellDivData(i).MeanSpeedNonDiv(1:length(VelFieldData(i).tau),:),2),21,'omitnan'),...
        'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','SouthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count + 1;
    SEM = nanstd(CellDivData(i).MeanSpeedDiv(1:length(VelFieldData(i).tau),:)./CellDivData(i).MeanSpeedNonDiv(1:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    shadedErrorBar(VelFieldData(i).tau,movmean(nanmean(CellDivData(i).MeanSpeedDiv(1:length(VelFieldData(i).tau),:)./CellDivData(i).MeanSpeedNonDiv(1:length(VelFieldData(i).tau),:),2),21,'omitnan'),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('Time in min','FontSize',text_font)
ylabel('Speed Division/Speed Layer','FontSize',text_font)
xlim([0 t_max]); 
set(gca, 'Fontsize',text_font)
title('Relative Speed around Divisions and Monolayer','FontSize',text_font+6);
saveas(gcf, sprintf('Ratio Speed Division vs Layer %s.png',savename));
savefig(sprintf('Ratio Speed Division vs Layer %s.fig',savename))
close all


% Plot Speed around divisions vs. rest of monolayer as box plots:
PlotVar = NaN(length(groups),50);
count = 0;
for i = groups
    count = count + 1;    
    PlotVar(count,1:size(CellDivData(i).MeanSpeedDiv,2)) = ...
        nanmean(CellDivData(i).MeanSpeedDiv(1:length(VelFieldData(i).tau),:)./CellDivData(i).MeanSpeedNonDiv(1:length(VelFieldData(i).tau),:),1);
end
% shape factor over density:
h=figure('units','normalized','outerposition',[0 0 1 1]);
set(h,'DefaultTextFontSize',24)
hold all
%h = boxplot(flip(Young,2),{'10µM Blebbistatin', '40µM Y-27632','U138','5µM Blebbistatin','20µM Y-27632','LN229'},'BoxStyle','outline','Width',0.5,'LabelOrientation','inline','LabelVerbosity','majorminor','Orientation','horizontal');
h = boxplot(PlotVar',LegendNames,'BoxStyle','outline','Width',0.5,'LabelOrientation','inline','LabelVerbosity','majorminor','Orientation','horizontal');
set(h,{'linewidth'},{2})
xlabel('Speed Division/Speed Layer','FontSize',24)
set(gca,'fontsize', 24);
xlim([1 1.7])
title('Relative Speed around Divisions and Monolayer','FontSize',text_font+6);
saveas(gcf, sprintf('Ratio Speed Division vs Layer Box Plot %s.png',savename));
savefig(sprintf('Ratio Speed Division vs Layer Box Plot %s.fig',savename))
close all


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
    shadedErrorBar(VelFieldData(i).tau,pMeanAll(1:length(VelFieldData(i).tau),i),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('\Deltat in min','FontSize',text_font)
ylabel('\alpha','FontSize',text_font)
xlim([0 t_max]);
set(gca, 'Fontsize',text_font)
title('Scaling Coefficient \alpha','FontSize',text_font+6);
saveas(gcf, sprintf('ScalingCoeffLog %s.png',savename));
savefig(sprintf('ScalingCoeffLog %s.fig',savename))
close all


% MSD
figure('units','normalized','outerposition',[0 0 1 1]);
loglog(VelFieldData(groups(1)).tau,PxSize.^2*nanmean(VelFieldData(groups(1)).MSDAll(1:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2, 'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    loglog(VelFieldData(i).tau,PxSize.^2*nanmean(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2, 'Color',color_map(count,:))
end
legend(LegendNames, 'Location','SouthEast','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count + 1;
    SEM = PxSize.^2*nanstd(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    errorbar(VelFieldData(i).tau,PxSize.^2*nanmean(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
end
xlabel('\Deltat in min','FontSize',text_font)
ylabel('MSD in µm²','FontSize',text_font)
xlim([0 t_max]);
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
    shadedErrorBar(VelFieldData(i).tau,nanmean(VelFieldData(i).QAll(1:length(VelFieldData(i).tau),:),2),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('Time in min','FontSize',text_font)
ylabel('Orderparameter','FontSize',text_font)
ylim([0 1]);
xlim([0 t_max]); 
set(gca, 'Fontsize',text_font)
title('Orderparameter','FontSize',text_font+6);
saveas(gcf, sprintf('OrderParameter %s.png',savename));
savefig(sprintf('OrderParameter %s.fig',savename))
close all

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
    shadedErrorBar(VelFieldData(i).tau,nanmean(VelFieldData(i).ChiAll(1:length(VelFieldData(i).tau),:),2),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('Time in min','FontSize',text_font)
ylabel('4 Point Suceptibility','FontSize',text_font)
xlim([0 t_max]);
set(gca, 'Fontsize',text_font)
title('4 Point Suceptibility','FontSize',text_font+6)
saveas(gcf, sprintf('Chi %s.png',savename));
savefig(sprintf('Chi %s.fig',savename))
close all

% Speed in µm/h
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
    shadedErrorBar(VelFieldData(i).tau,PxSize/(dt/60)*nanmean(VelFieldData(i).Speed(1:length(VelFieldData(i).tau),:),2),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)
end
xlabel('Time in min','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
xlim([0 t_max])
set(gca, 'Fontsize',text_font)
title('Speed','FontSize',text_font+6);
saveas(gcf, sprintf('Speed %s.png',savename));
savefig(sprintf('Speed %s.fig',savename))
close all

% Root mean squared speed in µm/h
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
    shadedErrorBar(VelFieldData(i).tau(tempOffset(i)/dt+1:end)-tempOffset(i),offset(i) + PxSize/(dt/60)*nanmean(VelFieldData(i).RMSVelAll(tempOffset(i)/dt+1:length(VelFieldData(i).tau),:),2),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)

end
xlabel('Time in min','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
xlim([0 t_max])
set(gca, 'Fontsize',text_font)
title('Speed','FontSize',text_font+6);
saveas(gcf, sprintf('RootMeanSpeed %s.png',savename));
savefig(sprintf('RootMeanSpeed %s.fig',savename))
close all

% compare vector fields:
% generate comparison value for random uniform distribution of angles from
% 0 to 90 degree (as many orientations are bound to that angle):
% MaxVar = [];
% for i = 1:10000
%     MaxVar(i) = circ_std(0.5*pi*rand(8374,1));
% end
% MaxVar = max(MaxVar);
MaxVar = 0.4563;

% Angular variance of cell orientation field
% field:
% Calculate Angular variance and mean:
InStdMean = NaN(t_max/dt,size(VecFieldData,2));
InStdStd = NaN(t_max/dt,size(VecFieldData,2));
SampleSize = [];
for j = 1:size(VecFieldData,2)
    SampleSize(j) = size(VecFieldData(j).AngleOrientation,1);
    for i = 1:t_max/dt-1
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
    shadedErrorBar(0:dt:t_max -dt,InStdMean(:,i),SEM(:,i),'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:)},'patchSaturation',0.2)   
end
xlabel('Time in min','FontSize',text_font)
ylabel('Angular Variance','FontSize',text_font)
set(gca, 'Fontsize',text_font)
title('Angular Variance of Cell Orientation','FontSize',text_font+6);
saveas(gcf, sprintf('Angular Variance of Cell Orientation %s.png',savename));
savefig(sprintf('Angular Variance of Cell Orientation %s.fig',savename))
close all



