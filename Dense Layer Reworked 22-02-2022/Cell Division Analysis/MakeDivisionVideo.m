function MakeDivisionVideo(CellProps,StartPos,jpg_path2,folderList,folder_number,m)

% Creates a video from the images including each cell division event:

% Older but significantly slower (10x) version - get frame is really slow!
% vidobj = VideoWriter('Divisions.avi');
% vidobj.FrameRate = 10;
% open(vidobj);
% % colormap:
% CMap = jet(length(CellProps));
% CMap = CMap(randperm(length(CellProps)),:);
% % visualize tracks:
% for i = 1:m
%     [images]=Reader(strcat(jpg_path2,'\',folderList(folder_number).name),i);
%     Idx = find(StartPos<=i);
%     Divisions = figure('name','Divisions','NumberTitle','off','OuterPosition',[0 0 855 640]);
%     imshow(images)
%     hold on
%     for k = 1:length(Idx)
%         plot(CellProps(Idx(k)).centroid(i,1),CellProps(Idx(k)).centroid(i,2),'r.','MarkerSize',30)
%     end
%     
%     for k = 1:length(Idx)
%         if i>=2 && ~isnan(CellProps(Idx(k)).centroid(i,1))
%             %k
%             for j = 2:i
%                 % Plots the travelled path.
%                 line([CellProps(Idx(k)).centroid(j-1,1),CellProps(Idx(k)).centroid(j,1)], ...
%                     [CellProps(Idx(k)).centroid(j-1,2),CellProps(Idx(k)).centroid(j,2)], ...
%                     'LineWidth',3,'Color',CMap(Idx(k),:));
%             end
%         end
%     end
%     % Gives the images names and save them
%     % saveas(gcf,sprintf('img%04d.tif',i));
%     curr_frame = getframe;
%     curr_frame.cdata = imresize(curr_frame.cdata, [640 855]);
%     writeVideo(vidobj,curr_frame);
%     %close curr_frame;
%     close Divisions;
% end
% close(vidobj);

% newer and 10x faster version:
vidobj = VideoWriter('DivisionsTest.avi');
vidobj.FrameRate = 10;
open(vidobj);
% colormap:
CMap = jet(length(CellProps));
CMap = CMap(randperm(length(CellProps)),:);
% visualize tracks:
for i = 1:m
    [images]=Reader(strcat(jpg_path2,'\',folderList(folder_number).name),i);
    MaxVal = double(max(images(:)));
    Idx = find(StartPos<=i);
    
    
    % get center coordinates:
    x = NaN(length(Idx),1);
    y = x;
    for k = 1:length(Idx)
        x(k) = CellProps(Idx(k)).centroid(i,1);
        y(k) = CellProps(Idx(k)).centroid(i,2);
    end
    
    
    x(isnan(x)) = [];
    y(isnan(y)) = [];
    images = insertShape(images,'FilledCircle',[x,y,12*ones(size(x))],'Color','red','Opacity',1);
    %figure; imshow(images)
    
    % get line coordinates:
    Coords = cell(length(Idx),1);
    NotEmpty = zeros(length(Idx),i);    
    for k = 1:length(Idx)
        if i>=2 && ~isnan(CellProps(Idx(k)).centroid(i,1))
            %k
            for j = 1:i
                % Plots the travelled path.
                Coords{k}(j,1) = CellProps(Idx(k)).centroid(j,1);
                Coords{k}(j,2) = CellProps(Idx(k)).centroid(j,2);
                NotEmpty(k,j) = 1;
            end
        end
    end
    NotEmpty = sum(NotEmpty,2);
    NotEmpty = NotEmpty>0;
    NotEmptyIDX = find(NotEmpty == 1);
    Coords(NotEmpty == 0) = [];
    
    % as insert shape cannot handle NaNs: 
    % remove all (leading??) NaNs
    for k = 1:length(Coords)
        tmp = Coords{k};
        tmp(isnan(tmp(:,1)),:) = []; 
        if size(tmp,1) > 2
            images = insertShape(images,'Line',tmp,'LineWidth',8,'Color',MaxVal.*CMap(Idx(NotEmptyIDX(k)),:),'Opacity',1);   
        end
    end
    % rescale image:
    images = imresize(images, [640 855]);
    curr_frame = im2frame(images);
    
    writeVideo(vidobj,curr_frame);
    %close curr_frame;
    %close Divisions;
end
toc
close(vidobj);