% make bar plot for speed:
FiltSize = 21;
TimePoint = 0;
TimePoint = 19.0;

Input = {};
count = 0;
for i = groups
    count = count+1;
    Input{count} = VelFieldData(i).Speed;
end

TimePoint = TimePoint*60/dt+1;
[stats,Output] = SpheroidStat(Input,TimePoint,FiltSize);
[c,m,h,gnames] = multcompare(stats,'CType','hsd');

% Generate Bar Plot
%savename = 'A375 + ACM + CBs';
data_in =PxSize/(dt/60)*nanmean(Output,1);
standard_in = PxSize/(dt/60)*nanstd(Output,[],1)./sqrt(size(Output,1) - sum(isnan(Output),1));
Name = LegendNames;
%Name = {'A375', '5% MG','10% MG','15% MG','30% MG'};
Title = 'Speed - LN229 19h';
%Title = 'Spheroid Size A375 70h';
significance = {'','','',''};
group_name = {'','','',''};
significance = {'','','','',''};
group_name = {'','','','',''};
%group_name = {'U138','U138','U138'};
text_font  = 24;
n = 5;
gaps = [];
y_axis_label = 'Speed in µm/h';
Leg = Name;
plot_inv(data_in, standard_in, Name, Title, significance, group_name, text_font, n, gaps, y_axis_label,Leg)

% save bar plot:
saveas(gcf, sprintf('Collective Speed %s.png',savename));
savefig(sprintf('Collective Speed %s.fig',savename))





% Scaling coefficients from MSD
WSize = VelFieldData(groups(1)).WSize;
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldData(groups(1)).tau((WSize-1)/2:end),nanmean(VelFieldData(groups(1)).pAll((WSize-1)/2:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)    
    count = count + 1;
    semilogx(VelFieldData(i).tau((WSize-1)/2:end),nanmean(VelFieldData(i).pAll((WSize-1)/2:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)    
    count = count + 1;
    SEM = nanstd(VelFieldData(i).pAll((WSize-1)/2:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    errorbar(VelFieldData(i).tau((WSize-1)/2:end),nanmean(VelFieldData(i).pAll((WSize-1)/2:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
end
xlabel('Time in min','FontSize',text_font)
ylabel('\alpha','FontSize',text_font)
%ylim([0.90 1.47]); 
xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('\alpha','FontSize',text_font+6);
saveas(gcf, sprintf('ScalingCoeffLog %s.png',savename));
savefig(sprintf('ScalingCoeffLog %s.fig',savename))
close all

% Scaling coefficients from MSD
WSize = VelFieldData(groups(1)).WSize;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau((WSize-1)/2:end),nanmean(VelFieldData(groups(1)).pAll((WSize-1)/2:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)    
    count = count + 1;
    plot(VelFieldData(i).tau((WSize-1)/2:end),nanmean(VelFieldData(i).pAll((WSize-1)/2:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)    
    count = count + 1;
    SEM = nanstd(VelFieldData(i).pAll((WSize-1)/2:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    errorbar(VelFieldData(i).tau((WSize-1)/2:end),nanmean(VelFieldData(i).pAll((WSize-1)/2:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
end
xlabel('Time in min','FontSize',text_font)
ylabel('\alpha','FontSize',text_font)
%ylim([0.9 1.85]); 
%ylim([0.90 1.47]); 
xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('\alpha','FontSize',text_font+6);
saveas(gcf, sprintf('ScalingCoeff %s.png',savename));
savefig(sprintf('ScalingCoeff %s.fig',savename))
close all

% Scaling coefficients from MSD -smoothed
WSize = VelFieldData(groups(1)).WSize;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau((WSize-1)/2:end),smooth(nanmean(VelFieldData(groups(1)).pAll((WSize-1)/2:length(VelFieldData(groups(1)).tau),:),2),15,'sgolay'),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)    
    count = count + 1;
    plot(VelFieldData(i).tau((WSize-1)/2:end),smooth(nanmean(VelFieldData(i).pAll((WSize-1)/2:length(VelFieldData(i).tau),:),2),15,'sgolay'),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthEast','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)    
    count = count + 1;
    SEM = smooth(nanstd(VelFieldData(i).pAll((WSize-1)/2:length(VelFieldData(i).tau),:),[],2),15,'sgolay')./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    errorbar(VelFieldData(i).tau((WSize-1)/2:end),smooth(nanmean(VelFieldData(i).pAll((WSize-1)/2:length(VelFieldData(i).tau),:),2),15,'sgolay'),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
end
xlabel('Time in min','FontSize',text_font)
ylabel('\alpha','FontSize',text_font)
%ylim([0.9 1.85]); 
%ylim([0.90 1.47]); 
xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',text_font)
title('\alpha','FontSize',text_font+6);
saveas(gcf, sprintf('ScalingCoeffSmoothed %s.png',savename));
savefig(sprintf('ScalingCoeffSmoothed %s.fig',savename))
close all

% Speed in µm/h
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldData(groups(1)).tau,PxSize/(dt/60)*nanmean(VelFieldData(groups(1)).Speed(1:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)    
    count = count + 1;
    semilogx(VelFieldData(i).tau,PxSize/(dt/60)*nanmean(VelFieldData(i).Speed(1:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','SouthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)    
    count = count + 1;
    SEM = PxSize/(dt/60)*nanstd(VelFieldData(i).Speed(1:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    errorbar(VelFieldData(i).tau,PxSize/(dt/60)*nanmean(VelFieldData(i).Speed(1:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
end
xlabel('Time in min','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
%ylim([1 9.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
xlim([0 t_max])
set(gca, 'Fontsize',text_font)
title('Speed','FontSize',text_font+6);
saveas(gcf, sprintf('SpeedLog %s.png',savename));
savefig(sprintf('SpeedLog %s.fig',savename))
close all


% Root mean squared speed in µm/h
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldData(groups(1)).tau,PxSize/(dt/60)*nanmean(VelFieldData(groups(1)).RMSVelAll(1:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)    
    count = count + 1;
    semilogx(VelFieldData(i).tau,PxSize/(dt/60)*nanmean(VelFieldData(i).RMSVelAll(1:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','SouthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)    
    count = count + 1;
    SEM = PxSize/(dt/60)*nanstd(VelFieldData(i).RMSVelAll(1:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    errorbar(VelFieldData(i).tau,PxSize/(dt/60)*nanmean(VelFieldData(i).RMSVelAll(1:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
end
xlabel('Time in min','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
%ylim([0 33]); 
%ylim([0 11.37]);
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
xlim([0 t_max])
set(gca, 'Fontsize',text_font)
title('Speed','FontSize',text_font+6);
saveas(gcf, sprintf('RootMeanSpeedLog %s.png',savename));
savefig(sprintf('RootMeanSpeedLog %s.fig',savename))
close all

% % Angular speed
% figure('units','normalized','outerposition',[0 0 1 1]);
% semilogx(VelFieldData(groups(1)).tau,nanmean(VelFieldData(groups(1)).AngVelAll(1:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
% hold all
% for i = groups(2:end)
%     semilogx(VelFieldData(i).tau,nanmean(VelFieldData(i).AngVelAll(1:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
% end
% legend(LegendNames, 'Location','NorthWest')
% 
% % Angular speed
% figure('units','normalized','outerposition',[0 0 1 1]);
% plot(VelFieldData(groups(1)).tau,nanmean(VelFieldData(groups(1)).AngVelAll(1:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
% hold all
% for i = groups(2:end)
%     plot(VelFieldData(i).tau,nanmean(VelFieldData(i).AngVelAll(1:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
% end
% legend(LegendNames, 'Location','NorthWest')

% Area of quick movers in µm²
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau,VelFieldPxSize.^2*medfilt1(nanmean(VelFieldData(groups(1)).AreaAll(1:length(VelFieldData(groups(1)).tau),:),2),20),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)    
    count = count + 1;
    plot(VelFieldData(i).tau,VelFieldPxSize.^2*medfilt1(nanmean(VelFieldData(i).AreaAll(1:length(VelFieldData(i).tau),:),2),20),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)    
    count = count + 1;
    SEM = VelFieldPxSize.^2*medfilt1(nanstd(VelFieldData(i).AreaAll(1:length(VelFieldData(i).tau),:),[],2),20)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(groups(1)).tau),:),2));
    errorbar(VelFieldData(i).tau,VelFieldPxSize.^2*medfilt1(nanmean(VelFieldData(i).AreaAll(1:length(VelFieldData(i).tau),:),2),20),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
end
xlabel('Time in min','FontSize',text_font)
ylabel('\xi in µm²','FontSize',text_font)
%ylim([100 1000]); 
%ylim([0 1200]); 
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
xlim([0 t_max])
set(gca, 'Fontsize',text_font)
title('\xi','FontSize',text_font+6);
saveas(gcf, sprintf('AreaQuickMovers %s.png',savename));
savefig(sprintf('AreaQuickMovers %s.fig',savename))
close all
% 
% % Area of quick movers normalized to scratch area in µm²
% 
% if size(Scratch.StartPosition,2) == 2
%     AbsVal = abs(Scratch.StartPosition(folder_number-2,2)-Scratch.StartPosition(folder_number-2,1));
%     Slope = (Scratch.EndWidth(folder_number-2)- AbsVal)/Scratch.EndTime(folder_number-2);
%     p_lin = [Slope,AbsVal];
%     ScratchCenterT = mean(Scratch.StartPosition(folder_number-2,:));
% else
%     AbsVal = abs(Scratch.StartPosition(folder_number-2,1));
%     Slope = (Scratch.EndWidth(folder_number-2)- AbsVal)/Scratch.EndTime(folder_number-2);
%     p_lin = [Slope,AbsVal];
%     % Denote this is a "virtual" scratch center only (corresponding to
%     % half the empty space from scratch to image boarder):
%     ScratchCenterT = mean([Scratch.StartPosition(folder_number-2,:),0]);
% end
% WidthCalculated(i) = p_lin(1)*i+p_lin(2);

figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau,VelFieldPxSize.^2*medfilt1(nanmean(VelFieldData(groups(1)).AreaAll(1:length(VelFieldData(groups(1)).tau),:),2),20),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)    
    count = count + 1;
    plot(VelFieldData(i).tau,VelFieldPxSize.^2*medfilt1(nanmean(VelFieldData(i).AreaAll(1:length(VelFieldData(i).tau),:),2),20),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)    
    count = count + 1;
    SEM = VelFieldPxSize.^2*medfilt1(nanstd(VelFieldData(i).AreaAll(1:length(VelFieldData(i).tau),:),[],2),20)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    errorbar(VelFieldData(i).tau,VelFieldPxSize.^2*medfilt1(nanmean(VelFieldData(i).AreaAll(1:length(VelFieldData(i).tau),:),2),20),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
end
xlabel('Time in min','FontSize',text_font)
ylabel('\xi in µm²','FontSize',text_font)
%ylim([0 43]); 
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
xlim([0 t_max])
set(gca, 'Fontsize',text_font)
title('\xi','FontSize',text_font+6);
saveas(gcf, sprintf('AreaQuickMovers %s.png',savename));
savefig(sprintf('AreaQuickMovers %s.fig',savename))
close all

% % Area of quick movers - chose only part with increasing
% % parameter in µm²
% offset = 220;
% figure('units','normalized','outerposition',[0 0 1 1]);
% plot(VelFieldData(groups(1)).tau(offset:end),VelFieldPxSize.^2*nanmean(VelFieldData(groups(1)).AreaAll(offset:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
% hold all
% for i = groups(2:end)
%     plot(VelFieldData(i).tau(offset:end),VelFieldPxSize.^2*nanmean(VelFieldData(i).AreaAll(offset:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2)
% end
% legend(LegendNames, 'Location','NorthWest')
% xlabel('Time in min','FontSize',22)
% ylabel('\xi in µm²','FontSize',22)



% groups = [9:12];
% groups = [1,2,3];
% Correlation length in µm vs root mean squared velocity in µm/h;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(PxSize/(dt/60)*nanmean(VelFieldData(groups(1)).RMSVelAll(1:length(VelFieldData(groups(1)).tau),:),2),(CorrLength(1:length(VelFieldData(groups(1)).tau),groups(1))),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)    
    count = count + 1;
    plot(PxSize/(dt/60)*nanmean(VelFieldData(i).RMSVelAll(1:length(VelFieldData(i).tau),:),2),(CorrLength(1:length(VelFieldData(i).tau),i)),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font,'AutoUpdate','off')
xlabel('Speed in µm/h','Fontsize',text_font)
ylabel('Correlation Length in µm','Fontsize',text_font)
%ylim([8 22]); 
%xlim([0.5 8]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Correlation Length vs Speed','FontSize',text_font+6);
saveas(gcf, sprintf('CorrelationLength vs Speed Astros+Inh %s.png',savename));
savefig(sprintf('CorrelationLength vs Speed Astros+Inh %s.fig',savename))
close all
% Works, but gives no helpful information
% Extract scaling coefficients:
for i = groups
TempVarX = [];
TempVarY = [];
%groups = 3
%for i = groups
    TempVarX = [TempVarX,PxSize/(dt/60)*nanmean(VelFieldData(i).RMSVelAll(1:length(VelFieldData(groups(1)).tau),:),2)];
    TempVarY = [TempVarY,CorrLength(1:length(VelFieldData(groups(1)).tau),i)];
%end
% Remove NaNs:
TempVarX(1,:) = [];
TempVarY(1,:) = [];
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% extract global scaling coefficient
FType = fittype('a*x+b')
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType)
uncertainty = confint(f,0.95)
end
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


% Correlation length in µm vs area of quick movers in µm²;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldPxSize.^2*medfilt1(nanmean(VelFieldData(groups(1)).AreaAll(1:length(VelFieldData(groups(1)).tau),:),2),20),(CorrLength(1:length(VelFieldData(groups(1)).tau),groups(1))),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)    
    count = count + 1;
    plot(VelFieldPxSize.^2*medfilt1(nanmean(VelFieldData(i).AreaAll(1:length(VelFieldData(i).tau),:),2),20),(CorrLength(1:length(VelFieldData(i).tau),i)),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font,'AutoUpdate','off')
xlabel('\xi in µm²','Fontsize',text_font)
ylabel('Correlation Length in µm','Fontsize',text_font)
%ylim([8 22]); 
%xlim([0.5 8]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Correlation Length vs \xi','FontSize',text_font+6);
saveas(gcf, sprintf('CorrelationLength vs Speed %s.png',savename));
savefig(sprintf('CorrelationLength vs Speed %s.fig',savename))
close all


% Drop of Chi Peak Height
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau,smooth(ChiPeakHeight(1:length(VelFieldData(groups(1)).tau),groups(1)),1)/5,'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)    
    count = count + 1;
    plot(VelFieldData(i).tau,smooth(ChiPeakHeight(1:length(VelFieldData(i).tau),i),1)/5,'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthEast','FontSize',text_font,'AutoUpdate','off')
count = 0;
for i = groups(1:end)    
    count = count + 1;
    SEM = smooth(ChiPeakStd(1:length(VelFieldData(i).tau),i),1)/5./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(groups(1)).tau),:),2));
    errorbar(VelFieldData(i).tau,smooth(ChiPeakHeight(1:length(VelFieldData(i).tau),i),1)/5,SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
end
xlabel('Time in min','Fontsize',text_font)
ylabel('4 Point Suceptibility Peak Height','Fontsize',text_font)
%ylim([0 32]); 
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
xlim([0 t_max])
set(gca, 'Fontsize',text_font)
title('4 Point Suceptibility Peak Height','FontSize',text_font+6);
saveas(gcf, sprintf('ChiDrop %s.png',savename));
savefig(sprintf('ChiDrop %s.fig',savename))
close all

% Change of Chi Peak Position
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau,dt*smooth(ChiPeakPos(1:length(VelFieldData(groups(1)).tau),groups(1)),5),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)    
    count = count + 1;
    plot(VelFieldData(i).tau,dt*smooth(ChiPeakPos(1:length(VelFieldData(i).tau),i),5),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthEast','FontSize',text_font,'AutoUpdate','off')
% count = 0;
% for i = groups(1:end)    
%     count = count + 1;
%     SEM = smooth(ChiPosStd(1:length(VelFieldData(i).tau),i),5)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(groups(1)).tau),:),2));
%     errorbar(VelFieldData(i).tau,dt*smooth(ChiPeakPos(1:length(VelFieldData(i).tau),i),1),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',color_map(count,:))
% end
xlabel('Time in min','Fontsize',text_font)
ylabel('4 Point Suceptibility Peak Position in min','Fontsize',text_font)
ylim([0 120]); 
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
xlim([0 t_max])
set(gca, 'Fontsize',text_font)
title('4 Point Suceptibility Peak Position','FontSize',text_font+6);
saveas(gcf, sprintf('ChiPos %s.png',savename));
savefig(sprintf('ChiPos %s.fig',savename))
close all

% Drop of Q Drop Pos
figure('units','normalized','outerposition',[0 0 1 1]);
plot(VelFieldData(groups(1)).tau,smooth(QDropPos(1:length(VelFieldData(groups(1)).tau),groups(1)),1),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for i = groups(2:end)    
    count = count + 1;
    plot(VelFieldData(i).tau,smooth(QDropPos(1:length(VelFieldData(i).tau),i),1),'o','MarkerSize',5,'LineWidth',2,'Color',color_map(count,:))
end
legend(LegendNames, 'Location','NorthWest','FontSize',text_font,'AutoUpdate','off')
xlabel('Time in min', 'Fontsize',text_font)
ylabel('Order Parameter Half Width', 'Fontsize',text_font)
%ylim([0 43]); 
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
xlim([0 t_max])
set(gca, 'Fontsize',text_font)
title('Order Parameter Half Width','FontSize',text_font+6);
saveas(gcf, sprintf('OrderParameterDrop %s.png',savename));
savefig(sprintf('OrderParameterDrop %s.fig',savename))
close all

% V velocity order parameter
figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldData(groups(1)).tau,nanmean(nanmean(VelFieldData(groups(1)).VTempAll(:,1:length(VelFieldData(groups(1)).tau),:),1),3),...
    'o','MarkerSize',5,'LineWidth',2,'Color',color_map(1,:))
hold all
count = 1;
for j = groups(2:end)    
    count = count + 1;  
    semilogx(VelFieldData(j).tau,nanmean(nanmean(VelFieldData(j).VTempAll(:,1:length(VelFieldData(j).tau),:),1),3),...
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
%ylim([0 43]); 
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
xlim([0 t_max])
set(gca,'Fontsize',text_font)
title('Velocity Order Parameter','FontSize',text_font+6);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time dependent plots:

% % V velocity order parameter
% figure('units','normalized','outerposition',[0 0 1 1]);
% semilogx(VelFieldData(groups(1)).tau(1:WSize-1),nanmean(nanmean(VelFieldData(groups(1)).VTempAll(1:end-1,Times(1)-offset:Times(1)+offset,:),2),3),...
%     'o','MarkerSize',5,'LineWidth',2)
% hold all
% if length(Times>1)
%     for i = 2:length(Times)
%         semilogx(VelFieldData(groups(1)).tau(1:WSize-1),nanmean(nanmean(VelFieldData(groups(1)).VTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
%             'o','MarkerSize',5,'LineWidth',2)
%     end
% end
% for j = groups(2:end)
%     for i = 1:length(Times)
%         semilogx(VelFieldData(j).tau(1:WSize-1),nanmean(nanmean(VelFieldData(j).VTempAll(1:end-1,Times(i)-offset:Times(i)+offset,:),2),3),...
%             'o','MarkerSize',5,'LineWidth',2) 
%     end
% end
% legend(Legends,'Location','NorthWest','Fontsize',text_font)
% xlabel('Time in min','Fontsize',text_font)
% ylabel('Velocity Order Parameter','Fontsize',text_font)
% %ylim([0 1]); 
% xlim([0 dt*2*offset]); % xlim([0 max(VelFieldData(groups(1)).tau)])
% set(gca,'Fontsize',text_font)
% title('Velocity Order Parameter','FontSize',text_font+6);

% velocity auto-correlation
t_end = 66;
% calculate linear offset:
offset_length = 20;
figure('units','normalized','outerposition',[0 0 1 1]);
loglog(VelFieldPxSize*[0:t_end-1],nanmean(nanmean(VelFieldData(groups(1)).AutoCorrAll(1:t_end,Times(1)-offset:Times(1)+offset,:),2),3),...
    'o','MarkerSize',5,'LineWidth',2)
hold all
if length(Times>1)
    for i = 2:length(Times)
        loglog(VelFieldPxSize*[0:t_end-1],nanmean(nanmean(VelFieldData(1).AutoCorrAll(1:t_end,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2)
    end
end
for j = groups(2:end)
    for i = 1:length(Times)
        loglog(VelFieldPxSize*[0:t_end-1] ,nanmean(nanmean(VelFieldData(j).AutoCorrAll(1:t_end,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2) 
    end
end
legend(Legends,'Location','SouthWest','Fontsize',text_font)
xlabel('Distance in µm','Fontsize',text_font)
ylabel('Auto Correlation Function','Fontsize',text_font)
ylim([0 1]); 
xlim([0 dt*2*offset]); % xlim([0 max(VelFieldData(groups(1)).tau)])
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
plot(VelFieldPxSize*[0:t_end-1],nanmean(nanmean(VelFieldData(groups(1)).AutoCorrAll(1:t_end,Times(1)-offset:Times(1)+offset,:),2),3),...
    'o','MarkerSize',5,'LineWidth',2)
hold all
if length(Times>1)
    for i = 2:length(Times)
        plot(VelFieldPxSize*[0:t_end-1],nanmean(nanmean(VelFieldData(1).AutoCorrAll(1:t_end,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2)
    end
end
for j = groups(2:end)
    for i = 1:length(Times)
        plot(VelFieldPxSize*[0:t_end-1] ,nanmean(nanmean(VelFieldData(j).AutoCorrAll(1:t_end,Times(i)-offset:Times(i)+offset,:),2),3),...
            'o','MarkerSize',5,'LineWidth',2) 
    end
end
legend(Legends,'Location','NorthEast','Fontsize',text_font)
xlabel('Distance in µm','Fontsize',text_font)
ylabel('Auto Correlation Function','Fontsize',text_font)
ylim([0 1]); 
xlim([0 dt*2*offset]); % xlim([0 max(VelFieldData(groups(1)).tau)])
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
%offset = 21;
% number of experiments;
num = 5;
% Get mean and deviation of cell density in 1/mm^2:
% Get values in cells per px for ROCK/Myosin II inhibition experiments:
for i = 1:size(Area,3)
    for j = 1:size(Area,4)
        var1 = Area(:,:,i,j);
        AreaAll(i,j) = nanmean(var1(:));
        STDAll(i,j) = nanstd(var1(:));
        Number(i,j) = size(var1,1)*size(var1,2) - sum(isnan(var1(:)));
    end
end
MeanDens = (ImSize(1)*ImSize(2))./squeeze(AreaAll);
% Proportion of std relative to mean area:
SEMDens = squeeze(STDAll)./(squeeze(AreaAll)*sqrt(num));
MeanDens =  MeanDens./ImPhysArea;
SEMDens =  SEMDens.*MeanDens;

% Get values in cells per px for control experiments:
MeanDens = 1./squeeze(AreaAll(1,:,:,:));
SEMDens = squeeze(AreaAll(2,:,:,:))./(squeeze(AreaAll(1,:,:,:)).*sqrt(num));
% Transform to cells per mm^2:
MeanDens =  MeanDens*(ImSize(1)*ImSize(2))/ImPhysArea;
SEMDens =  SEMDens.*MeanDens;

count = 0;
for j = Idx
    count = count + 1;
    for i = 1:length(VelFieldData)        
        Interval = [max(j-offset,1):min(j+offset,length(VelFieldData(1).tau)+1)];
        VelDens(count,i) = PxSize/(dt/60)*nanmean(nanmean(VelFieldData(1,i).RMSVelAll(Interval,:),1),2);
        VelDensStd(count,i) = PxSize/(dt/60)*nanstd(nanmean(VelFieldData(1,i).RMSVelAll(Interval,:),1),[],2);
        AreaDens(count,i) = VelFieldPxSize.^2*nanmean(nanmean(VelFieldData(1,i).AreaAll(Interval,:),1),2);
        CorrLengthDens(count,i) = nanmean(CorrLength(Interval,i));
        VelOrderParamDens(count,i) = nanmean(nanmean(nanmean(VelFieldData(1,i).VTempAll(:,j,:),1),3));
        ChiPeakHeightDens(count,i) = nanmean(ChiPeakHeight(Interval,i));
    end
end

% For Myosin/ROCK Inhibition Experiments:
% cell density versus root mean squared velocity:
% astrocytes:

figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens;
TempVarY = VelDens;
TempVarStdX = SEMDens;
TempVarStdY = VelDensStd/sqrt(num);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
hold all
for i = 1:size(TempVarX,2)
%plot(TempVarX(:,i),TempVarY(:,i),'o','MarkerSize',5,'LineWidth',5)
errorbar(TempVarX(:,i),TempVarY(:,i),TempVarStdY(:,i),TempVarStdY(:,i),TempVarStdX(:,i),TempVarStdX(:,i),'o','MarkerSize',5,'LineWidth',2)
end
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
uncertainty = confint(f,0.90);
hold all
legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Speed - Astrocytes + Inhibitors','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Speed - Astrocytes  + Inhibitors Errorbars1 %s.png',savename));
savefig(sprintf('Cell Density vs Speed - Astrocytes  + Inhibitors Errorbars1 %s.fig',savename))
close all

figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens;
TempVarY = VelDens;
TempVarStdX = SEMDens;
TempVarStdY = VelDensStd/sqrt(num);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
errorbar(TempVarX(:),TempVarY(:),TempVarStdY(:),TempVarStdY(:),TempVarStdX(:),TempVarStdX(:),'o','MarkerSize',5,'LineWidth',2)

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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Speed - Astrocytes  + Inhibitors','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Speed - Astrocytes  + Inhibitors Errorbars2 %s.png',savename));
savefig(sprintf('Cell Density vs Speed - Astrocytes  + Inhibitors Errorbars2 %s.fig',savename))
close all


% cell density versus area of quick movers:
% LN229:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:);
TempVarY = AreaDens(:,:);
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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs \xi - Astrocytes  + Inhibitors','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs AreaQuickMovers - Astrocytes  + Inhibitors %s.png',savename));
savefig(sprintf('Cell Density vs AreaQuickMovers - Astrocytes  + Inhibitors %s.fig',savename))
close all


% cell density versus correlation length:
% include linear fits:
% LN229:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:);
TempVarY = CorrLengthDens(:,:);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','pearson')
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
uncertainty = confint(f,0.95)
plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
hold all
%legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('Correlation length in µm','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Correlation Length - U138','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Correlation Length - U138  %s.png',savename));
savefig(sprintf('Cell Density vs Correlation Length - U138  %s.fig',savename))
close all




% For Control Experiments:
% cell density versus root mean squared velocity:
% astrocytes:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:,1);
TempVarY = VelDens(:,1:4);
TempVarStdX = SEMDens(:,:,1);
TempVarStdY = VelDensStd(:,1:4)/sqrt(num);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
%plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
errorbar(TempVarX(:),TempVarY(:),TempVarStdY(:),TempVarStdY(:),TempVarStdX(:),TempVarStdX(:),'o','MarkerSize',5,'LineWidth',2)
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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Speed - Astrocytes','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Speed - Astrocytes Errorbar %s.png',savename));
savefig(sprintf('Cell Density vs Speed - Astrocytes Errorbar %s.fig',savename))
close all

% LN229:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:,2);
TempVarY = VelDens(:,5:8);
TempVarStdX = SEMDens(:,:,2);
TempVarStdY = VelDensStd(:,5:8)/sqrt(num);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
uncertainty = confint(f,0.90);
%plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
errorbar(TempVarX(:),TempVarY(:),TempVarStdY(:),TempVarStdY(:),TempVarStdX(:),TempVarStdX(:),'o','MarkerSize',5,'LineWidth',2)
hold all
%legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Speed - LN229','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Speed - LN229 Errorbar %s.png',savename));
savefig(sprintf('Cell Density vs Speed - LN229 Errorbar %s.fig',savename))
close all

% U138:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:,3);
TempVarY = VelDens(:,9:12);
TempVarStdX = SEMDens(:,:,3);
TempVarStdY = VelDensStd(:,9:12)/sqrt(num);
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1)
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType);
uncertainty = confint(f,0.90);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
errorbar(TempVarX(:),TempVarY(:),TempVarStdY(:),TempVarStdY(:),TempVarStdX(:),TempVarStdX(:),'o','MarkerSize',5,'LineWidth',2)
hold all
% for i = 1:size(TempVarX,2)
% %plot(TempVarX(:,i),TempVarY(:,i),'o','MarkerSize',5,'LineWidth',5)
% errorbar(TempVarX(:,i),TempVarY(:,i),TempVarStdY(:,i),TempVarStdY(:,i),TempVarStdX(:,i),TempVarStdX(:,i),'o','MarkerSize',5,'LineWidth',2)
% end
%legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('Speed in µm/h','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Speed - U138','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Speed - U138 Errorbar Colored%s.png',savename));
savefig(sprintf('Cell Density vs Speed - U138 Errorbar Colored%s.fig',savename))
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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Correlation Length - U138','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Correlation Length - U138 %s.png',savename));
savefig(sprintf('Cell Density vs Correlation Length - U138 %s.fig',savename))
close all



R = corr(MeanDens(:),AreaDens(:),'rows', 'pairwise','type','spearman')

R = corr(VelDens(:),AreaDens(:),'rows', 'pairwise','type','spearman')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For Scratch Assay - Plots relative to cell density:
% needs the AreaAll variable 
savename = [];
% temporal offset in "images":
offset = 11;

% For Scratch assays:
% needs CellDensAll variable
MeanDens = CellDensAll' /ImPhysArea;
% Time points of measurement:
Idx = [80,620,950,1200]/dt;
% Get mean values of parameters corresponding to cell density:
count = 0;
for j = Idx
    count = count + 1;
    for i = 6%1:length(VelFieldData)        
        Interval = [max(j-offset,1):min(j+offset,length(VelFieldData(1).tau)+1)];
        VelDens(count,:) = PxSize/(dt/60)*nanmean(VelFieldData(1,i).RMSVelAll(Interval,:),1);
        AreaDens(count,:) = VelFieldPxSize.^2*nanmean(VelFieldData(1,i).AreaAll(Interval,:),1);
        CorrLengthDens(count,i) = nanmean(CorrLength(Interval,i));
        VelOrderParamDens(count,:) = nanmean(nanmean(VelFieldData(1,i).VTempAll(:,j,:),1));
        ChiPeakHeightDens(count,:) = nanmean(ChiPeakHeight(Interval,i));
    end
end

% cell density versus root mean squared velocity:
% astrocytes:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(1:4,:);
TempVarY = VelDens(1:4,:);
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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Speed - Astrocytes Scratch ','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Speed - Astrocytes Scratch %s.png',savename));
savefig(sprintf('Cell Density vs Speed - Astrocytes Scratch %s.fig',savename))
close all


% cell density versus area of quick movers:
% astrocytes:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = MeanDens(:,:);
TempVarY = AreaDens(:,:);
R = corr(TempVarX(:),TempVarY(:),'rows', 'pairwise','type','spearman')
p = polyfit(log(TempVarX(:)),log(TempVarY(:)),1);
% extract scaling coefficient
FType = fittype('a*x+b');
f = fit(log(TempVarX(:)),log(TempVarY(:)),FType)
uncertainty = confint(f,0.90);
plot(TempVarX(:),TempVarY(:),'o','MarkerSize',5,'LineWidth',5)
hold all
%legend(LegendNames, 'Location','SouthWest','Fontsize',text_font)
xlabel('Cell Density per mm²','FontSize',text_font)
ylabel('\xi in µm²','FontSize',text_font)
%ylim([1 1.37]); 
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs \xi - Astrocytes Scratch','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs AreaQuickMovers - Astrocytes Scratch %s.png',savename));
savefig(sprintf('Cell Density vs AreaQuickMovers - Astrocytes Scratch %s.fig',savename))
close all



% cell density versus correlation length:
% include linear fits:
% astrocytes:
figure('units','normalized','outerposition',[0 0 1 1]);
TempVarX = nanmean(MeanDens,2);
TempVarY = CorrLengthDens(:,:);
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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('Cell Density vs Correlation Length - Astrocytes Scratch','FontSize',text_font+6);
saveas(gcf, sprintf('Cell Density vs Correlation Length - Astrocytes Scratch %s.png',savename));
savefig(sprintf('Cell Density vs Correlation Length - Astrocytes Scratch %s.fig',savename))
close all


R = corr(MeanDens(:),AreaDens(:),'rows', 'pairwise','type','spearman')
R = corr(MeanDens(:),VelDens(:),'rows', 'pairwise','type','spearman')
R = corr(VelDens(:),AreaDens(:),'rows', 'pairwise','type','spearman')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison of heterogeneity length and speed:
% temporal offset in "images":
offset = 11;
savename = [];

% Get mean values of parameters corresponding to cell density:
count = 0;
for j = 1:length(VelFieldData(1).tau)
    count = count + 1;
    for i = 1:length(VelFieldData)        
        Interval = [max(j-offset,1):min(j+offset,length(VelFieldData(1).tau)+1)];
        VelDensAll(count,i) = PxSize/(dt/60)*nanmean(nanmean(VelFieldData(1,i).RMSVelAll(Interval,:),1),2);
        AreaDensAll(count,i) = VelFieldPxSize.^2*nanmean(nanmean(VelFieldData(1,i).AreaAll(Interval,:),1),2);        
    end
end

%AreaDensAll = CorrLength(2:end,:);

% cell speed versus area of quick movers:
% astrocytes:
% Groups to compare:
groups = [1:4];
%groups = [1:6];
LegendNames = [];
count = 0;
for i = groups
    count = count + 1;
    LegendNames{count} = VelFieldData(i).Name;
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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
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
    LegendNames{count} = VelFieldData(i).Name;
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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
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
    LegendNames{count} = VelFieldData(i).Name;
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
%xlim([0 1236]); % xlim([0 max(VelFieldData(groups(1)).tau)])
set(gca, 'Fontsize',text_font)
title('\xi vs Speed - U138','FontSize',text_font+6);
saveas(gcf, sprintf('Speed vs AreaQuickMovers - U138 %s.png',savename));
savefig(sprintf('Speed vs AreaQuickMovers - U138 %s.fig',savename))
close all

R = corr(VelDensAll(:),AreaDensAll(:),'rows', 'pairwise','type','spearman')










