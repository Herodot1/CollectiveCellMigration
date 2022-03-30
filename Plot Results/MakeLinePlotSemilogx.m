function MakeLinePlotSemilogx(xData,yData,SampleSize,MarkerSize,LineWidth,ColorMap,LegendLocation,FontSize,xLabel,yLabel,yLim,xLim,Title,Savename)

figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(VelFieldData(groups(1)).tau,nanmean(VelFieldData(groups(1)).QAll(1:length(VelFieldData(groups(1)).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',ColorMap(1,:))
hold all
count = 1;
for i = groups(2:end)
    count = count + 1;
    semilogx(VelFieldData(i).tau,nanmean(VelFieldData(i).QAll(1:length(VelFieldData(i).tau),:),2),'o','MarkerSize',5,'LineWidth',2,'Color',ColorMap(count,:))
end
legend(LegendNames, 'Location','SouthWest','FontSize',FontSize,'AutoUpdate','off')
count = 0;
for i = groups(1:end)
    count = count + 1;
    SEM = nanstd(VelFieldData(i).QAll(1:length(VelFieldData(i).tau),:),[],2)./sqrt(size(VelFieldData(i).MSDAll(1:length(VelFieldData(i).tau),:),2));
    %errorbar(VelFieldData(i).tau,nanmean(VelFieldData(i).QAll(1:length(VelFieldData(i).tau),:),2),SEM,'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',ColorMap(count,:))
    shadedErrorBar(VelFieldData(i).tau,nanmean(VelFieldData(i).QAll(1:length(VelFieldData(i).tau),:),2),SEM,'lineProps',{'.','MarkerSize',0.5,'LineWidth',0.2, 'Color',ColorMap(count,:)},'patchSaturation',0.2)
end
xlabel('Time in min','FontSize',FontSize)
ylabel('Orderparameter','FontSize',FontSize)
ylim([0 1]);
xlim([0 t_max]); % xlim([0 max(VelFieldData(groups(1)).tau)])
%xlim([0 1300])
set(gca, 'Fontsize',FontSize)
title('Orderparameter','FontSize',FontSize+6);
saveas(gcf, sprintf('OrderParameter %s.png',savename));
savefig(sprintf('OrderParameter %s.fig',savename))
close all