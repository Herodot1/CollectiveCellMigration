function [DelaunayNeighborsTimeReversed,DelaunayNeighbors,NumNewNeighbors,...
    NewNeighborsPerCell,NewNeigborsPerDiffusiveCell] = NeighborsAndTracks(TrackMat,DataPath,SaveName,PxSize)
% Chose 9 nearest neighbors to include the query point itself and its 8
% nearest neighbors:
NumNeighborsKNN = 9;
% Text font size:
text_font = 24;


% Check if several folders actually exist:
folder = fullfile(strcat(DataPath,'\',SaveName,'\Results'),'Results Path Analysis');
% Checks the existence of the folder with the name "Results". If it
% does not exist it is created.
if (exist(folder) == 0)
    mkdir('Results Path Analysis');
end


% Do Delaunay triangulation to find neighboring cells:
% Ignore first entry, as the squares lead to "edge sharing", giving a
% false number of neighbors. Sleight variations by cell movements are
% enough to cope with that. Thus start with image 5:
% Test the following strategy: start at the end and look at time
% evaluation, as the problem should be time-reversible.
TempVarX = squeeze(TrackMat(2:end-2,2:end-2,2,end));
TempVarY = squeeze(TrackMat(2:end-2,2:end-2,1,end));
dt = delaunayTriangulation(double([TempVarX(:), TempVarY(:)]));
% Get initial configuration for comparison:
InitialNeighbors = NaN(length(dt.Points),10);
for i = 1:length(dt.Points)
    [r,c] = find(dt.ConnectivityList == i);
    Temp = unique(dt.ConnectivityList(r,:));
    InitialNeighbors(i,1:length(Temp)) = Temp;
    NumNeighbors(i) = numel(Temp) - 1;
end
InitialNeighbors(InitialNeighbors==0) = NaN;
% figure
% voronoi(dt)

% Compare each image with the initial configuration:
ConservedNeighbors = NaN(size(TrackMat, 4),1);
for j = 1:size(TrackMat, 4)
    TempVarX = squeeze(TrackMat(2:end-2,2:end-2,2,j));
    TempVarY = squeeze(TrackMat(2:end-2,2:end-2,1,j));
    dt = delaunayTriangulation(double([TempVarX(:), TempVarY(:)]));
    CurrentNeighbors = NaN(length(dt.Points),10);
    ConservedNeighborsTemp = NaN(length(dt.Points),1);
    for i = 1:length(NumNeighbors)%1:length(dt.Points)
        [r,c] = find(dt.ConnectivityList == i);
        Temp = unique(dt.ConnectivityList(r,:));
        CurrentNeighbors(i,1:length(Temp)) = Temp;
        [Val,~,ind] = intersect(InitialNeighbors(i,:),CurrentNeighbors(i,1:length(Temp)),'stable');
        ConservedNeighborsTemp(i) =  (numel(Val) - 1)./NumNeighbors(i);
    end
    ConservedNeighbors(j) = nanmean(ConservedNeighborsTemp);
end
DelaunayNeighborsTimeReversed(1:size(TrackMat,4)) = ConservedNeighbors;

% No time reversion:
% Do Delaunay triangulation to find neighboring cells:
% Ignore first entry, as the squares lead to "edge sharing", giving a
% false number of neighbors. Sleight variations by cell movements are
% enough to cope with that. Thus start with image 5:
% Test the following strategy: start at the end and look at time
% evaluation, as the problem should be time-reversible.
TempVarX = squeeze(TrackMat(2:end-2,2:end-2,2,20));
TempVarY = squeeze(TrackMat(2:end-2,2:end-2,1,20));
dt = delaunayTriangulation(double([TempVarX(:), TempVarY(:)]));
% Get initial configuration for comparison:
InitialNeighbors = NaN(length(dt.Points),10);
for i = 1:length(dt.Points)
    [r,c] = find(dt.ConnectivityList == i);
    Temp = unique(dt.ConnectivityList(r,:));
    InitialNeighbors(i,1:length(Temp)) = Temp;
    NumNeighbors(i) = numel(Temp) - 1;
end

% figure
% voronoi(dt)

% Compare each image with the initial configuration:
ConservedNeighbors = NaN(size(TrackMat, 4),1);
for j = 1:size(TrackMat, 4)
    TempVarX = squeeze(TrackMat(2:end-2,2:end-2,2,j));
    TempVarY = squeeze(TrackMat(2:end-2,2:end-2,1,j));
    dt = delaunayTriangulation(double([TempVarX(:), TempVarY(:)]));
    CurrentNeighbors = NaN(length(dt.Points),10);
    ConservedNeighborsTemp = NaN(length(dt.Points),1);
    for i = 1:length(NumNeighbors)%1:length(dt.Points)
        [r,c] = find(dt.ConnectivityList == i);
        Temp = unique(dt.ConnectivityList(r,:));
        CurrentNeighbors(i,1:length(Temp)) = Temp;
        [Val,~,ind] = intersect(InitialNeighbors(i,:),CurrentNeighbors(i,1:length(Temp)),'stable');
        ConservedNeighborsTemp(i) =  (numel(Val) - 1)./NumNeighbors(i);
    end
    ConservedNeighbors(j) = nanmean(ConservedNeighborsTemp);
end
DelaunayNeighbors(1:size(TrackMat,4)) = ConservedNeighbors;

%% Alterntive nearest neighbors estimate using a simple euclidean metric for estimation of new neighbors
% Use a simple euclidean metric to determine nearest neighbors
% find nearest neighbors in the beginning and end of TrackMat:
% Start with the beginning:
TempVarX = squeeze(TrackMat(2:end-2,2:end-2,2,1));
TempVarY = squeeze(TrackMat(2:end-2,2:end-2,1,1));
XStart = [TempVarX(:),TempVarY(:)];
[IdxStart,DStart] = knnsearch(XStart,XStart,'K',NumNeighborsKNN);

% Continue with track ends:
TempVarX = squeeze(TrackMat(2:end-2,2:end-2,2,end));
TempVarY = squeeze(TrackMat(2:end-2,2:end-2,1,end));
XEnd = [TempVarX(:),TempVarY(:)];
[IdxEnd,DEnd] = knnsearch(XEnd,XEnd,'K',NumNeighborsKNN);

% Get new neighbors:
for i = 1:size(IdxStart,1)
    TempVar = [IdxStart(i,:)',IdxEnd(i,:)'];
    TempVar = TempVar(:);
    [tmp,Idx] = sort(TempVar);
    idp = diff(tmp)>0;
    TempVar = TempVar(Idx([true;idp]&[idp;true]));
    if isempty(TempVar)
        NewNeighbors(i) = 0;
    else
        % Division by 2 is necessary as this approach assignes each changed
        % neigbor 2 new ones. The "old" object is found and the "new" one
        % substituting the "old" one, effectively doubling the number of
        % new neighbors. Thus, divide by two.
        NewNeighbors(i) = numel(TempVar)/2;
    end
end
NewNeighbors = NewNeighbors./(NumNeighborsKNN-1);

% Number of cells with new neighbors;
NumNewNeighbors = numel(find(NewNeighbors>0))./length(NewNeighbors);
% average number of new neighbors per cell:
NewNeighborsPerCell = mean(NewNeighbors).*(NumNeighborsKNN-1);
% average number of new neighbors per cell only of those cells that made
% new neighbors:
NewNeigborsPerDiffusiveCell = mean(NewNeighbors(NewNeighbors>0)).*(NumNeighborsKNN-1);


% Got to save folder and save necessary data:
%cd(folder)
cd('Results Path Analysis')

% Make a plot of each Track:
figure('units','normalized','outerposition',[0 0 1 1]);
hold all
for i = 1:size(TrackMat,1)
    for j = 1:size(TrackMat,2)
        TempVar = PxSize*squeeze(TrackMat(i,j,:,:))';
        DiffX(i,j) = TempVar(1,2) - TempVar(end,2);
        DiffY(i,j) = TempVar(1,1) - TempVar(end,1);
        PathlengthX(i,j) = sum(([0; diff(TempVar(:,2))])); %sum(abs(TempVar(:,2)-TempVar(1,2)));
        PathlengthY(i,j) = sum(([0; diff(TempVar(:,1))]));
        plot(TempVar(:,2),TempVar(:,1),'LineWidth',2)
        plot(TempVar(1,2),TempVar(1,1),'.r','Markersize',20)
        plot(TempVar(end,2),TempVar(end,1),'.b','Markersize',20)
        
    end
end

%     % Set average PathLength in the image to zero:
%     for i = 1:size(TrackMat,1)
%         for j = 1:size(TrackMat,2)
%             TempVar = PxSize*squeeze(TrackMat(i,j,:,:))';
%             PathlengthX(i,j) = sum(([0; diff(TempVar(:,2))])); %sum(abs(TempVar(:,2)-TempVar(1,2)));
%             PathlengthY(i,j) = sum(([0; diff(TempVar(:,1))]));
%         end
%     end
%     figure('units','normalized','outerposition',[0 0 1 1]);
%     hold all
%     for i = 1:size(TrackMat,1)
%         for j = 1:size(TrackMat,2)
%             TempVar = PxSize*squeeze(TrackMat(i,j,:,:))';
%             % Set average PathLength in the image to zero:
%             TempVar(2:end,2) = TempVar(2:end,2) - nanmean(PathlengthX(:)).*[2:size(TempVar,1)]'/size(TempVar,1);
%             TempVar(2:end,1) = TempVar(2:end,1) - nanmean(PathlengthY(:)).*[2:size(TempVar,1)]'/size(TempVar,1);
%             plot(TempVar(:,2),TempVar(:,1),'LineWidth',2)
%             plot(TempVar(1,2),TempVar(1,1),'.r','Markersize',20)
%             plot(TempVar(end,2),TempVar(end,1),'.b','Markersize',20)
%
%         end
%     end
xlabel('Position in µm','FontSize',text_font)
ylabel('Position in µm','FontSize',text_font)
set(gca, 'Fontsize',text_font)
title('Path','FontSize',text_font+6);
saveas(gcf, sprintf('Plot of Path %s.png',SaveName));
savefig(sprintf('Plot of Path %s.fig',SaveName))
close all
% Go back to initial directory:
%cd(strcat(DataPath,'\',SaveName,'\Results'))
cd('..')