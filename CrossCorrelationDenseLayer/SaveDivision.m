function SaveDivision(CellProps,CType,NetName, MedFiltSize,FiltType,...
    BackThresh,WSize,MaxTimeDiff)

save('CellDivisions.mat','CellProps','-v7.3')
save('Settings.mat','CType','NetName', 'MedFiltSize','FiltType',...
    'BackThresh','WSize','MaxTimeDiff')