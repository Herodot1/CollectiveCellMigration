% Plot comparison plots with shifted time (scratch closure = 0):
a
clearvars -except AreaAll VelFieldData CorrLength;
%% Set constants:
% Groups to compare:
groups = [2,5,10];
% time difference between cell density estimates in minutes:
DensDt = 240;
% time difference between sucessive images in minutes:
dt = 3;
% Julian:
dt = 5;
% Set offset corresponding to the number of images before and after the set
% timepoint that are taken as reference:
offset = 10;
% Area proportion where cells have been counted:
dA = 0.3^2;
% Image size in pixels:
ImSize = [1280,960];
% julian
ImSize = [1040,1392];
% Physical dimension of image in µm:
ImPhysSize = [621.23,465.80];
% Julian
ImPhysSize = [478.18,640.18];
% pixel size of image in µm/px:
PxSize = ImPhysSize/ImSize;
% Area of image in mm^2:
ImPhysArea = ImPhysSize(1)*ImPhysSize(2)./10^6;
% Factor for absolute density (cells/mm^2);
CellDensFac = 1/(ImPhysArea*dA);
% Size Velocity Field in px:
VelFieldSize = [106,79];
% Julian:
VelFieldSize = [86,115];
% Pixelsize of velocity field in µm/px:
VelFieldPxSize = mean((ImSize./VelFieldSize))*PxSize;
% Text font size:
text_font = 24;
% Windowsize:
WSize = VelFieldData(1).WSize;
%% Set up storage structure and some parameters:
VelFieldDataShifted = struct('Name',[],'pAll',[],'VTempAll',[],'AreaAll',[],...
    'RMSVelAll',[],'AngVelAll',[],'MSDAll',[],'MSDTempAll',[],...
    'QAll',[],'QTempAll',[],'ChiAll',[],'ChiTempAll',[],'tau',[],...
    'Speed',[],'AutoCorrAll',[],'CellDensity',[],'WSize',[],...
    'ScratchStartPos',[],'ScratchEndWidth',[], 'ScratchEndTime',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create struct with shifted time:
% start from i = 2, because i = 1 corresponds to positions in the center of
% the monolayer.
for i = 2:length(VelFieldData)
    
    VelFieldDataShifted(i).Name                 = VelFieldData(i).Name;
    VelFieldDataShifted(i).WSize                = VelFieldData(i).WSize;
    VelFieldDataShifted(i).ScratchStartPos      = VelFieldData(i).ScratchStartPos ;
    VelFieldDataShifted(i).ScratchEndWidth      = VelFieldData(i).ScratchEndWidth;
    VelFieldDataShifted(i).ScratchEndTime       = VelFieldData(i).ScratchEndTime;
    
    for j = 1:size(VelFieldData(i).pAll,2)
        pos = VelFieldData(i).ScratchEndTime(j)-0;
        num = size(VelFieldData(i).pAll(pos:end,j),1);
        VelFieldDataShifted(i).pAll(1:num,j)        = VelFieldData(i).pAll(pos:end,j);
        VelFieldDataShifted(i).tau(1:num,1)           = VelFieldData(i).tau(1:num);
        VelFieldDataShifted(i).Speed(1:num,j)       = VelFieldData(i).Speed(pos:end,j);
        VelFieldDataShifted(i).AreaAll(1:num,j)     = VelFieldData(i).AreaAll(pos:end,j);
        VelFieldDataShifted(i).RMSVelAll(1:num,j)   = VelFieldData(i).RMSVelAll(pos:end,j);
        VelFieldDataShifted(i).AngVelAll(1:num,j)   = VelFieldData(i).AngVelAll(pos:end,j);
        VelFieldDataShifted(i).MSDAll(1:num,j)      = VelFieldData(i).MSDAll(pos:end,j);       
        VelFieldDataShifted(i).QAll(1:num,j)        = VelFieldData(i).QAll(pos:end,j);        
        VelFieldDataShifted(i).ChiAll(1:num,j)      = VelFieldData(i).ChiAll(pos:end,j);
        VelFieldDataShifted(i).MSDTempAll(:,1:num,j)    = VelFieldData(i).MSDTempAll(:,pos:end,j);
        VelFieldDataShifted(i).QTempAll(:,1:num,j)      = VelFieldData(i).QTempAll(:,pos:end,j);
        VelFieldDataShifted(i).ChiTempAll(:,1:num,j)    = VelFieldData(i).ChiTempAll(:,pos:end,j);       
        VelFieldDataShifted(i).VTempAll(:,1:num,j)      = VelFieldData(i).VTempAll(:,pos:end,j);
        
        num = size(VelFieldData(i).AutoCorrAll(:,pos:end,j),2);
        VelFieldDataShifted(i).AutoCorrAll   =  VelFieldData(i).AutoCorrAll(:,pos:end,j);
    end
    
    % Replace zeros:
    VelFieldDataShifted(i).VTempAll(find(VelFieldDataShifted(i).VTempAll(:,:,:)==0)) = NaN;
    VelFieldDataShifted(i).pAll(find(VelFieldDataShifted(i).pAll==0)) = NaN;
    VelFieldDataShifted(i).AreaAll(find(VelFieldDataShifted(i).AreaAll==0)) = NaN;
    VelFieldDataShifted(i).RMSVelAll(find(VelFieldDataShifted(i).RMSVelAll==0)) = NaN;
    VelFieldDataShifted(i).AngVelAll(find(VelFieldDataShifted(i).AngVelAll==0)) = NaN;
    VelFieldDataShifted(i).MSDAll(find(VelFieldDataShifted(i).MSDAll==0)) = NaN;
    VelFieldDataShifted(i).MSDTempAll(find(VelFieldDataShifted(i).MSDTempAll==0)) = NaN;
    VelFieldDataShifted(i).QAll(find(VelFieldDataShifted(i).QAll==0)) = NaN;
    VelFieldDataShifted(i).QTempAll(find(VelFieldDataShifted(i).QTempAll==0)) = NaN;
    VelFieldDataShifted(i).ChiAll(find(VelFieldDataShifted(i).ChiAll==0)) = NaN;
    VelFieldDataShifted(i).ChiTempAll(find(VelFieldDataShifted(i).ChiTempAll==0)) = NaN;
    VelFieldDataShifted(i).Speed(find(VelFieldDataShifted(i).Speed==0)) = NaN;
    VelFieldDataShifted(i).AutoCorrAll(find(VelFieldDataShifted(i).AutoCorrAll==0)) = NaN;
end

    % Go back to top level directory:

    
save('VelFieldDataShiftedTime.mat','VelFieldDataShifted')   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time evolution of some critical parameters:
ChiPeakHeight       = NaN(900,length(VelFieldDataShifted));
ChiPeakStd          = NaN(900,length(VelFieldDataShifted));
ChiPeakPos          = NaN(900,length(VelFieldDataShifted));
QDropPos            = NaN(900,length(VelFieldDataShifted));
AutoCorrDropPos     = NaN(900,length(VelFieldDataShifted));
VTemp               = NaN(900,length(VelFieldDataShifted));
for i = 1:length(VelFieldDataShifted)
    % Peak position and value of 4 point susceptibility over time:
    [PeakHeight,PeakPos]                    = max(nanmean(VelFieldDataShifted(i).ChiTempAll,3),[],1);
    PeakStd                                 = nanstd(max(VelFieldDataShifted(i).ChiTempAll,[],1),[],3);
    ChiPeakStd(1:length(PeakStd),i)         = PeakStd;
    ChiPeakHeight(1:length(PeakHeight),i)   = PeakHeight;
    ChiPeakPos(1:length(PeakPos),i)         = PeakPos;
    % Find Point where Q reaches 0.9 first:
    [~,PeakPos]                             = min(abs(nanmean(VelFieldDataShifted(i).QTempAll,3)-0.9),[],1);
    QDropPos(1:length(PeakPos),i)           = PeakPos;
     % Find Point where autocorrelation reaches 1/e  ~ 0.3679 first:
    [~,PeakPos]             = min(abs(nanmean(VelFieldDataShifted(i).AutoCorrAll,3)-1/exp(1)),[],1);
    AutoCorrDropPos(1:length(PeakPos),i)    = PeakPos;
    % Velocity order parameter:
    VTemp(1:size(VelFieldDataShifted(i).VTempAll,2),i) = nanmean(nanmean(VelFieldDataShifted(i).VTempAll,3),1);
    
end

% Calculate Autocorrelation length:
CorrLength = FitAutoCorrelation(VelFieldDataShifted,VelFieldPxSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot means of VelFieldData for different cell densities. Global parameters.
% Groups to compare:
groups = [3,7,11];
groups = [1,2,3];
groups = [1,2,3,4,5,6];
% Savename:
savename = 'Astros+LN+U138';
savename = 'Astros Scratch';


LegendNames = [];
count = 0;
for i = groups
    count = count + 1;
    LegendNames{count} = VelFieldDataShifted(i).Name;
end

% MSD
figure('units','normalized','outerposition',[0 0 1 1]);
loglog(VelFieldDataShifted(groups(1)).tau,PxSize.^2*nanmean(VelFieldDataShifted(groups(1)).MSDAll(1:length(VelFieldDataShifted(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
hold all
p = polyfit(log10(VelFieldDataShifted(groups(1)).tau),log10(PxSize.^2*nanmean(VelFieldDataShifted(groups(1)).MSDAll(1:length(VelFieldDataShifted(groups(1)).tau),:),2)),1);
yfit = polyval(p,log10(VelFieldDataShifted(groups(1)).tau));
loglog(VelFieldDataShifted(groups(1)).tau,10.^yfit,'--r','Linewidth',2,'HandleVisibility','off')

for i = groups(2:end)
    loglog(VelFieldDataShifted(i).tau,PxSize.^2*nanmean(VelFieldDataShifted(i).MSDAll(1:length(VelFieldDataShifted(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)    
    p = polyfit(log10(VelFieldDataShifted(i).tau),log10(PxSize.^2*nanmean(VelFieldDataShifted(i).MSDAll(1:length(VelFieldDataShifted(i).tau),:),2)),1);
    yfit = polyval(p,log10(VelFieldDataShifted(i).tau));
    loglog(VelFieldDataShifted(i).tau,10.^yfit,'--r','Linewidth',2,'HandleVisibility','off')
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font)
xlabel('\Deltat in min','FontSize',text_font)
ylabel('MSD in µm²','FontSize',text_font)
%ylim([0 123]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('MSD','FontSize',text_font+6);
saveas(gcf, sprintf('MSD %s.png',savename));
savefig(sprintf('MSD %s.fig',savename))
close all

% Q - order parameter
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldDataShifted(groups(1)).tau,nanmean(VelFieldDataShifted(groups(1)).QAll(1:length(VelFieldDataShifted(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
hold all
for i = groups(2:end)
    semilogx(VelFieldDataShifted(i).tau,nanmean(VelFieldDataShifted(i).QAll(1:length(VelFieldDataShifted(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font)
xlabel('Time in min','FontSize',text_font)
ylabel('Orderparameter','FontSize',text_font)
ylim([0 1]); 
xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Orderparameter','FontSize',text_font+6);
saveas(gcf, sprintf('OrderParameter %s.png',savename));
savefig(sprintf('OrderParameter %s.fig',savename))
close all

% Chi - 4 point susceptibility
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldDataShifted(groups(1)).tau,nanmean(VelFieldDataShifted(groups(1)).ChiAll(1:length(VelFieldDataShifted(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
hold all
for i = groups(2:end)
    semilogx(VelFieldDataShifted(i).tau,nanmean(VelFieldDataShifted(i).ChiAll(1:length(VelFieldDataShifted(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font)
xlabel('Time in min','FontSize',text_font)
ylabel('4 Point Suceptibility','FontSize',text_font)
%ylim([0 1]); 
xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('4 Point Suceptibility','FontSize',text_font+6)
saveas(gcf, sprintf('Chi %s.png',savename));
savefig(sprintf('Chi %s.fig',savename))
close all

% Scaling coefficients from MSD
WSize = VelFieldDataShifted(groups(1)).WSize;
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldDataShifted(groups(1)).tau((WSize-1)/2:end),nanmean(VelFieldDataShifted(groups(1)).pAll((WSize-1)/2:length(VelFieldDataShifted(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
hold all
for i = groups(2:end)
    semilogx(VelFieldDataShifted(i).tau((WSize-1)/2:end),nanmean(VelFieldDataShifted(i).pAll((WSize-1)/2:length(VelFieldDataShifted(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Time in min','FontSize',text_font)
ylabel('\alpha','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('\alpha','FontSize',text_font+6);
saveas(gcf, sprintf('ScalingCoeffLog %s.png',savename));
savefig(sprintf('ScalingCoeffLog %s.fig',savename))
close all

% Scaling coefficients from MSD
WSize = VelFieldDataShifted(groups(1)).WSize;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldDataShifted(groups(1)).tau((WSize-1)/2:end),nanmean(VelFieldDataShifted(groups(1)).pAll((WSize-1)/2:length(VelFieldDataShifted(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
hold all
for i = groups(2:end)
    plot(VelFieldDataShifted(i).tau((WSize-1)/2:end),nanmean(VelFieldDataShifted(i).pAll((WSize-1)/2:length(VelFieldDataShifted(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Time in min','FontSize',text_font)
ylabel('\alpha','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('\alpha','FontSize',text_font+6);
saveas(gcf, sprintf('ScalingCoeff %s.png',savename));
savefig(sprintf('ScalingCoeff %s.fig',savename))
close all

% Speed in µm/h
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldDataShifted(groups(1)).tau,PxSize/(dt/60)*nanmean(VelFieldDataShifted(groups(1)).Speed(1:length(VelFieldDataShifted(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
hold all
for i = groups(2:end)
    semilogx(VelFieldDataShifted(i).tau,PxSize/(dt/60)*nanmean(VelFieldDataShifted(i).Speed(1:length(VelFieldDataShifted(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','NorthWest','Fontsize',text_font)
xlabel('Time in min','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Speed','FontSize',text_font+6);
saveas(gcf, sprintf('SpeedLog %s.png',savename));
savefig(sprintf('SpeedLog %s.fig',savename))
close all

% Speed in µm/h
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldDataShifted(groups(1)).tau,PxSize/(dt/60)*nanmean(VelFieldDataShifted(groups(1)).Speed(1:length(VelFieldDataShifted(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
hold all
for i = groups(2:end)
    plot(VelFieldDataShifted(i).tau,PxSize/(dt/60)*nanmean(VelFieldDataShifted(i).Speed(1:length(VelFieldDataShifted(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','NorthWest','Fontsize',text_font)
xlabel('Time in min','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Speed','FontSize',text_font+6);
saveas(gcf, sprintf('Speed %s.png',savename));
savefig(sprintf('Speed %s.fig',savename))
close all

% Root mean squared speed in µm/h
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldDataShifted(groups(1)).tau,PxSize/(dt/60)*nanmean(VelFieldDataShifted(groups(1)).RMSVelAll(1:length(VelFieldDataShifted(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
hold all
for i = groups(2:end)
    semilogx(VelFieldDataShifted(i).tau,PxSize/(dt/60)*nanmean(VelFieldDataShifted(i).RMSVelAll(1:length(VelFieldDataShifted(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','NorthWest','Fontsize',text_font)
xlabel('Time in min','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
ylim([0 43]);  
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Speed','FontSize',text_font+6);
saveas(gcf, sprintf('RootMeanSpeedLog %s.png',savename));
savefig(sprintf('RootMeanSpeedLog %s.fig',savename))
close all

% Root mean squared speed in µm/h
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldDataShifted(groups(1)).tau,PxSize/(dt/60)*nanmean(VelFieldDataShifted(groups(1)).RMSVelAll(1:length(VelFieldDataShifted(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
hold all
for i = groups(2:end)
    plot(VelFieldDataShifted(i).tau,PxSize/(dt/60)*nanmean(VelFieldDataShifted(i).RMSVelAll(1:length(VelFieldDataShifted(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','NorthWest','Fontsize',text_font)
xlabel('Time in min','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
ylim([0 13]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Speed','FontSize',text_font+6);
saveas(gcf, sprintf('RootMeanSpeed %s.png',savename));
savefig(sprintf('RootMeanSpeed %s.fig',savename))
close all

% % Angular speed
% figure('units','normalized','outerposition',[0 0 1 1]);
% semilogx(VelFieldDataShifted(groups(1)).tau,nanmean(VelFieldDataShifted(groups(1)).AngVelAll(1:length(VelFieldDataShifted(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
% hold all
% for i = groups(2:end)
%     semilogx(VelFieldDataShifted(i).tau,nanmean(VelFieldDataShifted(i).AngVelAll(1:length(VelFieldDataShifted(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
% end
% legend(LegendNames, 'Location','NorthWest')
% 
% % Angular speed
% figure('units','normalized','outerposition',[0 0 1 1]);
% plot(VelFieldDataShifted(groups(1)).tau,nanmean(VelFieldDataShifted(groups(1)).AngVelAll(1:length(VelFieldDataShifted(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
% hold all
% for i = groups(2:end)
%     plot(VelFieldDataShifted(i).tau,nanmean(VelFieldDataShifted(i).AngVelAll(1:length(VelFieldDataShifted(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
% end
% legend(LegendNames, 'Location','NorthWest')

% Area of quick movers in µm²
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldDataShifted(groups(1)).tau,VelFieldPxSize.^2*medfilt1(nanmean(VelFieldDataShifted(groups(1)).AreaAll(1:length(VelFieldDataShifted(groups(1)).tau),:),2),20),'o','MarkerSize',5,'LineWidth',2)
hold all
for i = groups(2:end)
    plot(VelFieldDataShifted(i).tau,VelFieldPxSize.^2*medfilt1(nanmean(VelFieldDataShifted(i).AreaAll(1:length(VelFieldDataShifted(i).tau),:),2),20),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','SouthWest','FontSize',text_font)
xlabel('Time in min','FontSize',text_font)
ylabel('\xi in µm²','FontSize',text_font)
%ylim([0 43]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('\xi','FontSize',text_font+6);
saveas(gcf, sprintf('AreaQuickMovers %s.png',savename));
savefig(sprintf('AreaQuickMovers %s.fig',savename))
close all

% % Area of quick movers - chose only part with increasing
% % parameter in µm²
% offset = 220;
% figure('units','normalized','outerposition',[0 0 1 1]);
% plot(VelFieldDataShifted(groups(1)).tau(offset:end),VelFieldPxSize.^2*nanmean(VelFieldDataShifted(groups(1)).AreaAll(offset:length(VelFieldDataShifted(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
% hold all
% for i = groups(2:end)
%     plot(VelFieldDataShifted(i).tau(offset:end),VelFieldPxSize.^2*nanmean(VelFieldDataShifted(i).AreaAll(offset:length(VelFieldDataShifted(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
% end
% legend(LegendNames, 'Location','NorthWest')
% xlabel('Time in min','FontSize',22)
% ylabel('\xi in µm²','FontSize',22)

% Correlation length in µm
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldDataShifted(groups(1)).tau,smooth(CorrLength(1:length(VelFieldDataShifted(groups(1)).tau),groups(1)),15),'o','MarkerSize',5,'LineWidth',2)
hold all
for i = groups(2:end)
    plot(VelFieldDataShifted(i).tau,smooth(CorrLength(1:length(VelFieldDataShifted(i).tau),i),15),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','NorthEast','Fontsize',text_font)
xlabel('Time in min','Fontsize',text_font)
ylabel('Correlation Length in µm','Fontsize',text_font)
%ylim([12 22]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Correlation Length','FontSize',text_font+6);
saveas(gcf, sprintf('CorrelationLength %s.png',savename));
savefig(sprintf('CorrelationLength %s.fig',savename))
close all

groups = [9:12];
groups = [1,2,3];
% Correlation length in µm vs root mean squared velocity in µm/h;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(PxSize/(dt/60)*nanmean(VelFieldDataShifted(groups(1)).RMSVelAll(1:length(VelFieldDataShifted(groups(1)).tau),:),2),(CorrLength(1:length(VelFieldDataShifted(groups(1)).tau),groups(1))),'o','MarkerSize',5,'LineWidth',2)
hold all
for i = groups(2:end)
    plot(PxSize/(dt/60)*nanmean(VelFieldDataShifted(i).RMSVelAll(1:length(VelFieldDataShifted(i).tau),:),2),(CorrLength(1:length(VelFieldDataShifted(i).tau),i)),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','NorthEast','Fontsize',text_font)
xlabel('Speed in µm/h','Fontsize',text_font)
ylabel('Correlation Length in µm','Fontsize',text_font)
%ylim([8 22]); 
%xlim([0.5 8]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Correlation Length vs Speed','FontSize',text_font+6);
saveas(gcf, sprintf('CorrelationLength vs Speed %s.png',savename));
savefig(sprintf('CorrelationLength vs Speed %s.fig',savename))
close all
% Works, but gives no helpful information
% % Extract scaling coefficients:
% TempVarX = [];
% TempVarY = [];
% for i = groups
%     TempVarX = [TempVarX,PxSize/(dt/60)*nanmean(VelFieldDataShifted(i).RMSVelAll(1:length(VelFieldDataShifted(groups(1)).tau),:),2)];
%     TempVarY = [TempVarY,CorrLength(1:length(VelFieldDataShifted(groups(1)).tau),i)];
% end
% % Remove NaNs:
% TempVarX(1,:) = [];
% TempVarY(1,:) = [];
% R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
% p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% % extract global scaling coefficient
% FType = fittype('a*x+b');
% f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
% uncertainty = confint(f,0.90);
% 
% % Sort values according to speed:
% [TempVarX,I] = sort(TempVarX);
% TempVarY = TempVarY(I);
% WSize = 41;
% % Take WSize values, starting from the lowest one to calculate speed
% % dependent scaling coefficient:
% for j = 1:length(groups)
%     for i = (WSize-1)/2+1: length(TempVarX)-(WSize-1)/2 
%         f = fit(log(TempVarX(i-(WSize-1)/2:i+(WSize-1)/2,j)),log(TempVarY(i-(WSize-1)/2:i+(WSize-1)/2,j)),FType);
%         ScalingCoeff(i-(WSize-1)/2,j) = f.a;
%         % Corresponding mean Speed:
%         CoeffSpeed(i-(WSize-1)/2,j) = nanmean(TempVarX(i-(WSize-1)/2:i+(WSize-1)/2,j));
%     end
% end
% 
% figure;
% plot(nanmean(CoeffSpeed,2),nanmean(ScalingCoeff,2),'o')
% xlim([27,35])
% ylim([-1,1])


% Drop of Chi Peak Height
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldDataShifted(groups(1)).tau,smooth(ChiPeakHeight(1:length(VelFieldDataShifted(groups(1)).tau),groups(1)),1),'o','MarkerSize',5,'LineWidth',2)
hold all
for i = groups(2:end)
    plot(VelFieldDataShifted(i).tau,smooth(ChiPeakHeight(1:length(VelFieldDataShifted(i).tau),i),1),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','NorthWest','Fontsize',text_font)
xlabel('Time in min','Fontsize',text_font)
ylabel('4 Point Suceptibility Peak Height','Fontsize',text_font)
%ylim([0 43]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('4 Point Suceptibility Peak Height','FontSize',text_font+6);
saveas(gcf, sprintf('ChiDrop %s.png',savename));
savefig(sprintf('ChiDrop %s.fig',savename))
close all

% Drop of Q Drop Pos
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldDataShifted(groups(1)).tau,smooth(QDropPos(1:length(VelFieldDataShifted(groups(1)).tau),groups(1)),1),'o','MarkerSize',5,'LineWidth',2)
hold all
for i = groups(2:end)
    plot(VelFieldDataShifted(i).tau,smooth(QDropPos(1:length(VelFieldDataShifted(i).tau),i),1),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','NorthEast','Fontsize',text_font)
xlabel('Time in min', 'Fontsize',text_font)
ylabel('Order Parameter Half Width', 'Fontsize',text_font)
%ylim([0 43]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Order Parameter Half Width','FontSize',text_font+6);
saveas(gcf, sprintf('OrderParameterDrop %s.png',savename));
savefig(sprintf('OrderParameterDrop %s.fig',savename))
close all

% V velocity order parameter
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldDataShifted(groups(1)).tau,nanmean(nanmean(VelFieldDataShifted(groups(1)).VTempAll(:,1:length(VelFieldDataShifted(groups(1)).tau),:),1),3),...
    'o','MarkerSize',5,'LineWidth',2)
hold all

for j = groups(2:end)    
    semilogx(VelFieldDataShifted(j).tau,nanmean(nanmean(VelFieldDataShifted(j).VTempAll(:,1:length(VelFieldDataShifted(j).tau),:),1),3),...
        'o','MarkerSize',5,'LineWidth',2)    
end
legend(LegendNames,'Location','NorthWest','Fontsize',text_font)
xlabel('Time in min','Fontsize',text_font)
ylabel('Velocity Order Parameter','Fontsize',text_font)
%ylim([0 43]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca,'Fontsize',text_font)
title('Velocity Order Parameter','FontSize',text_font+6);

% V velocity order parameter
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldDataShifted(groups(1)).tau,nanmean(nanmean(VelFieldDataShifted(groups(1)).VTempAll(:,1:length(VelFieldDataShifted(groups(1)).tau),:),1),3),...
    'o','MarkerSize',5,'LineWidth',2)
hold all

for j = groups(2:end)    
    plot(VelFieldDataShifted(j).tau,nanmean(nanmean(VelFieldDataShifted(j).VTempAll(:,1:length(VelFieldDataShifted(j).tau),:),1),3),...
        'o','MarkerSize',5,'LineWidth',2)    
end
legend(LegendNames,'Location','NorthWest','Fontsize',text_font)
xlabel('Time in min','Fontsize',text_font)
ylabel('Velocity Order Parameter','Fontsize',text_font)
%ylim([0 43]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca,'Fontsize',text_font)
title('Velocity Order Parameter','FontSize',text_font+6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time dependent plots:

% Get different times:
Times = [21:160:300];
Times = 21;
% Take temporal offset:
offset = 20;

% Take care of constraints imposed by window size:
WSize = VelFieldDataShifted(1).WSize;
Times(Times <(WSize-1)/2) = []; 
Times(Times>size(VelFieldDataShifted(1).AutoCorrAll-2,2)-WSize) = [];

LegendNames = [];
count = 0;
for i = groups
    count = count + 1;
    LegendNames{count} = VelFieldDataShifted(i).Name;
end
for j = 1:length(groups)
    for i = 1:length(Times)
        idx = (j-1)*length(Times) + i;
        Legends{idx} = sprintf('%s t = %0.5g h',LegendNames{j},(Times(i)-1-offset)*3/60 );
    end
end

% Q - order parameter
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldDataShifted(groups(1)).tau(1:WSize-1),nanmean(nanmean(VelFieldDataShifted(groups(1)).QTempAll(1:end-1,Times(1)-offset:Times(1)+offset,:),2),3),...
    'o','MarkerSize',5,'LineWidth',2)
hold all
if length(Times>1)
    for i = 2:length(Times)
        semilogx(VelFieldDataShifted(groups(1)).tau(1:WSize-1),nanmean(nanmean(VelFieldDataShifted(groups(1)).QTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2)
    end
end
for j = groups(2:end)
    for i = 1:length(Times)
        semilogx(VelFieldDataShifted(j).tau(1:WSize-1),nanmean(nanmean(VelFieldDataShifted(j).QTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2) 
    end
end
legend(Legends,'Location','SouthWest','Fontsize',text_font)
xlabel('Time in min','Fontsize',text_font)
ylabel('Orderparameter','Fontsize',text_font)
ylim([0 1]); 
xlim([0 dt*2*offset]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca,'Fontsize',text_font)
title('Order Parameter','FontSize',text_font+6);
saveas(gcf, sprintf('OrderParameterTime %s.png',savename));
savefig(sprintf('OrderParameterTime %s.fig',savename))
close all

% Chi - 4 point susceptibility
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldDataShifted(groups(1)).tau(1:WSize-1),nanmean(nanmean(VelFieldDataShifted(groups(1)).ChiTempAll(1:end-1,Times(1)-offset:Times(1)+offset,:),2),3),...
    'o','MarkerSize',5,'LineWidth',2)
hold all
if length(Times>1)
    for i = 2:length(Times)
        semilogx(VelFieldDataShifted(groups(1)).tau(1:WSize-1),nanmean(nanmean(VelFieldDataShifted(groups(1)).ChiTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2)
    end
end
for j = groups(2:end)
    for i = 1:length(Times)
        semilogx(VelFieldDataShifted(j).tau(1:WSize-1),nanmean(nanmean(VelFieldDataShifted(j).ChiTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2) 
    end
end
legend(Legends,'Location','NorthWest','Fontsize',text_font)
xlabel('Time in min','Fontsize',text_font)
ylabel('4 Point Suceptibility','Fontsize',text_font)
%ylim([0 1]); 
xlim([0 dt*2*offset]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca,'Fontsize',text_font)
title('4 Point Suceptibility','FontSize',text_font+6);
saveas(gcf, sprintf('ChiTime %s.png',savename));
savefig(sprintf('ChiTime %s.fig',savename))
close all

% MSD in µm² per minute
figure('units','normalized','outerposition',[0 0 1 1]);
loglog(VelFieldDataShifted(groups(1)).tau(1:WSize-1),PxSize.^2*nanmean(nanmean(VelFieldDataShifted(groups(1)).MSDTempAll(1:end-1,Times(1)-offset:Times(1)+offset,:),2),3),...
    'o','MarkerSize',5,'LineWidth',2)
hold all
if length(Times>1)
    for i = 2:length(Times)
        loglog(VelFieldDataShifted(groups(1)).tau(1:WSize-1),PxSize.^2*nanmean(nanmean(VelFieldDataShifted(groups(1)).MSDTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2)
    end
end
for j = groups(2:end)
    for i = 1:length(Times)
        loglog(VelFieldDataShifted(j).tau(1:WSize-1),PxSize.^2*nanmean(nanmean(VelFieldDataShifted(j).MSDTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2) 
    end
end
legend(Legends,'Location','NorthWest','Fontsize',text_font)
xlabel('\Deltat in min','Fontsize',text_font)
ylabel('MSD in µm²','Fontsize',text_font)
%ylim([0 1]); 
xlim([0 dt*2*offset]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca,'Fontsize',text_font)
title('MSD','FontSize',text_font+6);
saveas(gcf, sprintf('MSDTime %s.png',savename));
savefig(sprintf('MSDTime %s.fig',savename))
close all

% % V velocity order parameter
% figure('units','normalized','outerposition',[0 0 1 1]);
% semilogx(VelFieldDataShifted(groups(1)).tau(1:WSize-1),nanmean(nanmean(VelFieldDataShifted(groups(1)).VTempAll(1:end-1,Times(1)-offset:Times(1)+offset,:),2),3),...
%     'o','MarkerSize',5,'LineWidth',2)
% hold all
% if length(Times>1)
%     for i = 2:length(Times)
%         semilogx(VelFieldDataShifted(groups(1)).tau(1:WSize-1),nanmean(nanmean(VelFieldDataShifted(groups(1)).VTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
%             'o','MarkerSize',5,'LineWidth',2)
%     end
% end
% for j = groups(2:end)
%     for i = 1:length(Times)
%         semilogx(VelFieldDataShifted(j).tau(1:WSize-1),nanmean(nanmean(VelFieldDataShifted(j).VTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
%             'o','MarkerSize',5,'LineWidth',2) 
%     end
% end
% legend(Legends,'Location','NorthWest','Fontsize',text_font)
% xlabel('Time in min','Fontsize',text_font)
% ylabel('Velocity Order Parameter','Fontsize',text_font)
% %ylim([0 1]); 
% xlim([0 dt*2*offset]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
% set(gca,'Fontsize',text_font)
% title('Velocity Order Parameter','FontSize',text_font+6);

% velocity auto-correlation
t_end = 66;
% calculate linear offset:
offset_length = 20;
figure('units','normalized','outerposition',[0 0 1 1]);
loglog(VelFieldPxSize*[0:t_end-1],nanmean(nanmean(VelFieldDataShifted(groups(1)).AutoCorrAll(1:t_end,Times(1)-offset:Times(1)+offset,:),2),3),...
    'o','MarkerSize',5,'LineWidth',2)
hold all
if length(Times>1)
    for i = 2:length(Times)
        loglog(VelFieldPxSize*[0:t_end-1],nanmean(nanmean(VelFieldDataShifted(1).AutoCorrAll(1:t_end,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2)
    end
end
for j = groups(2:end)
    for i = 1:length(Times)
        loglog(VelFieldPxSize*[0:t_end-1] ,nanmean(nanmean(VelFieldDataShifted(j).AutoCorrAll(1:t_end,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2) 
    end
end
legend(Legends,'Location','SouthWest','Fontsize',text_font)
xlabel('Distance in µm','Fontsize',text_font)
ylabel('Auto Correlation Function','Fontsize',text_font)
ylim([0 1]); 
xlim([0 dt*2*offset]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca,'Fontsize',text_font)
title('Auto Correlation Function','FontSize',text_font+6);
saveas(gcf, sprintf('AutoCorrelationTime Log %s.png',savename));
savefig(sprintf('AutoCorrelationTime Log %s.fig',savename))
close all


% velocity auto-correlation
t_end = 66;
% calculate linear offset:
offset_length = 20;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldPxSize*[0:t_end-1],nanmean(nanmean(VelFieldDataShifted(groups(1)).AutoCorrAll(1:t_end,Times(1)-offset:Times(1)+offset,:),2),3),...
    'o','MarkerSize',5,'LineWidth',2)
hold all
if length(Times>1)
    for i = 2:length(Times)
        plot(VelFieldPxSize*[0:t_end-1],nanmean(nanmean(VelFieldDataShifted(1).AutoCorrAll(1:t_end,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2)
    end
end
for j = groups(2:end)
    for i = 1:length(Times)
        plot(VelFieldPxSize*[0:t_end-1] ,nanmean(nanmean(VelFieldDataShifted(j).AutoCorrAll(1:t_end,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2) 
    end
end
legend(Legends,'Location','NorthEast','Fontsize',text_font)
xlabel('Distance in µm','Fontsize',text_font)
ylabel('Auto Correlation Function','Fontsize',text_font)
ylim([0 1]); 
xlim([0 dt*2*offset]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca,'Fontsize',text_font)
title('Auto Correlation Function','FontSize',text_font+6);
saveas(gcf, sprintf('AutoCorrelationTime %s.png',savename));
savefig(sprintf('AutoCorrelationTime %s.fig',savename))
close all





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots relative to cell density:
% needs the AreaAll variable
savename = [];
% temporal offset in "images":
offset = 11;
% Get mean and deviation of cell density in 1/mm^2:
% Get values in cells per px:
MeanDens = 1./squeeze(AreaAll(1,:,:,:));
SEMDens = 1./squeeze(AreaAll(2,:,:,:));
% Transform to cells per mm^2:
MeanDens =  MeanDens*(ImSize(1)*ImSize(2))/ImPhysArea;
SEMDens =  SEMDens*(ImSize(1)*ImSize(2))/ImPhysArea .* (1./squeeze(sqrt(AreaAll(3,:,:,:))));



% Get mean values of parameters corresponding to cell density:
count = 0;
for j = Idx
    count = count + 1;
    for i = 1:length(VelFieldDataShifted)        
        Interval = [max(j-offset,1):min(j+offset,length(VelFieldDataShifted(1).tau)+1)];
        VelDens(count,i) = PxSize/(dt/60)*nanmean(nanmean(VelFieldDataShifted(1,i).RMSVelAll(Interval,:),1),2);
        AreaDens(count,i) = VelFieldPxSize.^2*nanmean(nanmean(VelFieldDataShifted(1,i).AreaAll(Interval,:),1),2);
        %CorrLengthDens(count,i) = nanmean(CorrLength(Interval,i));
        VelOrderParamDens(count,i) = nanmean(nanmean(nanmean(VelFieldDataShifted(1,i).VTempAll(:,j,:),1),3));
        ChiPeakHeightDens(count,i) = nanmean(ChiPeakHeight(Interval,i));
    end
end

% cell density versus root mean squared velocity:
% astrocytes:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:,1);
TempVarY = VelDens(:,1:4);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
uncertainty = confint(f,0.90);
hold all
%legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Speed - Astrocytes','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Speed - Astrocytes %s.png',savename));
savefig(sprintf('Cell Density vs Speed - Astrocytes %s.fig',savename))
close all

% LN229:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:,2);
TempVarY = VelDens(:,5:8);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
uncertainty = confint(f,0.90);
plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
hold all
%legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Speed - LN229','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Speed - LN229 %s.png',savename));
savefig(sprintf('Cell Density vs Speed - LN229 %s.fig',savename))
close all

% U138:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:,3);
TempVarY = VelDens(:,9:12);
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
uncertainty = confint(f,0.90);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
hold all
%legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Speed - U138','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Speed - U138 %s.png',savename));
savefig(sprintf('Cell Density vs Speed - U138 %s.fig',savename))
close all

% cell density versus area of quick movers:
% astrocytes:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:,1);
TempVarY = AreaDens(:,1:4);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1);
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
uncertainty = confint(f,0.90);

plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
hold all
%legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('\xi in µm²','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs \xi - Astrocytes','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs AreaQuickMovers - Astrocytes %s.png',savename));
savefig(sprintf('Cell Density vs AreaQuickMovers - Astrocytes %s.fig',savename))
close all

% LN229:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:,2);
TempVarY = AreaDens(:,5:8);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
uncertainty = confint(f,0.90);
plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
hold all
%legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('\xi in µm²','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs \xi - LN229','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs AreaQuickMovers - LN229 %s.png',savename));
savefig(sprintf('Cell Density vs AreaQuickMovers - LN229 %s.fig',savename))
close all

% U138:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:,3);
TempVarY = AreaDens(:,9:12);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
uncertainty = confint(f,0.90);
plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
hold all
%legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('\xi in µm²','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs \xi - U138','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs AreaQuickMovers - U138 %s.png',savename));
savefig(sprintf('Cell Density vs AreaQuickMovers - U138 %s.fig',savename))
close all

% cell density versus correlation length:
% include linear fits:
% astrocytes:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:,1);
TempVarY = CorrLengthDens(:,1:4);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','pearson')
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
uncertainty = confint(f,0.90);
plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
hold all
%legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('Correlation length in µm','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Correlation Length - Astrocytes','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Correlation Length - Astrocytes %s.png',savename));
savefig(sprintf('Cell Density vs Correlation Length - Astrocytes %s.fig',savename))
close all

% LN229:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:,2);
TempVarY = CorrLengthDens(:,5:8);
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
uncertainty = confint(f,0.90);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','pearson')
plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
hold all
%legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('Correlation Length in µm','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Correlation Length - LN229','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Correlation Length - LN229 %s.png',savename));
savefig(sprintf('Cell Density vs Correlation Length - LN229 %s.fig',savename))
close all

% U138:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:,3);
TempVarY = CorrLengthDens(:,9:12);
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
uncertainty = confint(f,0.90);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','pearson')
plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
hold all
%legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('Correlation Length in µm','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Correlation Length - U138','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Correlation Length - U138 %s.png',savename));
savefig(sprintf('Cell Density vs Correlation Length - U138 %s.fig',savename))
close all



R = corr(MeanDens(:),AreaDens(:),'rows', 'pairwise','type','spearman')

R = corr(VelDens(:),AreaDens(:),'rows', 'pairwise','type','spearman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison of heterogeneity length and speed:
% temporal offset in "images":
offset = 11;
savename = [];

% Get mean values of parameters corresponding to cell density:
count = 0;
for j = 1:length(VelFieldDataShifted(1).tau)
    count = count + 1;
    for i = 1:length(VelFieldDataShifted)        
        Interval = [max(j-offset,1):min(j+offset,length(VelFieldDataShifted(1).tau)+1)];
        VelDensAll(count,i) = PxSize/(dt/60)*nanmean(nanmean(VelFieldDataShifted(1,i).RMSVelAll(Interval,:),1),2);
        AreaDensAll(count,i) = VelFieldPxSize.^2*nanmean(nanmean(VelFieldDataShifted(1,i).AreaAll(Interval,:),1),2);        
    end
end

%AreaDensAll = CorrLength(2:end,:);

% cell speed versus area of quick movers:
% astrocytes:
% Groups to compare:
groups = [1:4];
LegendNames = [];
count = 0;
for i = groups
    count = count + 1;
    LegendNames{count} = VelFieldDataShifted(i).Name;
end
figure('units','normalized','outerposition',[0 0 1 1]);
hold all
for i = groups
    TempVarX = VelDensAll(:,i);
    TempVarY = AreaDensAll(:,i);
    plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','NorthWest','Fontsize',text_font)
xlabel('Speed in µm/h²','FontSize',text_font)
ylabel('\xi in µm²','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('\xi vs Speed - Astrocytes','FontSize',text_font+6);
saveas(gcf, sprintf('Speed vs AreaQuickMovers - Astrocytes %s.png',savename));
savefig(sprintf('Speed vs AreaQuickMovers - Astrocytes %s.fig',savename))
close all

% LN229:
% Groups to compare:
groups = [5:8];
LegendNames = [];
count = 0;
for i = groups
    count = count + 1;
    LegendNames{count} = VelFieldDataShifted(i).Name;
end
figure('units','normalized','outerposition',[0 0 1 1]);
hold all
for i = groups
    TempVarX = VelDensAll(:,i);
    TempVarY = AreaDensAll(:,i);
    plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Speed in µm/h²','FontSize',text_font)
ylabel('\xi in µm²','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('\xi vs Speed - LN229','FontSize',text_font+6);
saveas(gcf, sprintf('Speed vs AreaQuickMovers - LN229 %s.png',savename));
savefig(sprintf('Speed vs AreaQuickMovers - LN229 %s.fig',savename))
close all

% U138:
% Groups to compare:
groups = [9:12];
LegendNames = [];
count = 0;
for i = groups
    count = count + 1;
    LegendNames{count} = VelFieldDataShifted(i).Name;
end
figure('units','normalized','outerposition',[0 0 1 1]);
hold all
for i = groups
    TempVarX = VelDensAll(:,i);
    TempVarY = AreaDensAll(:,i);
    plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',2)
end
legend(LegendNames, 'Location','NorthWest','Fontsize',text_font)
xlabel('Speed in µm/h²','FontSize',text_font)
ylabel('\xi in µm²','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldDataShifted(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('\xi vs Speed - U138','FontSize',text_font+6);
saveas(gcf, sprintf('Speed vs AreaQuickMovers - U138 %s.png',savename));
savefig(sprintf('Speed vs AreaQuickMovers - U138 %s.fig',savename))
close all

R = corr(VelDens(:),AreaDens(:),'rows', 'pairwise','type','spearman')




