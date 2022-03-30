function [stack_new, m,m_max]=Reader(jpg_path2,varargin)
% the second argument needs to be a vector [x y] that defines the start and
% end image for reading.

%     % Get the full path: Initial folder(jpg_path2)  +
%     % Subfolder(folder_name)
%     jpg_path2 = strcat(jpg_path2,'\',folder_name);
%     % Get the extension of the image Files.
%     extension = dir(jpg_path2);
%     extension = extension(size(extension,1)).name;
%     extension = extension(end-3:end);
    if length(varargin) > 2
        disp('to many input arguments')
        return
    end
     
    jpg_path2 = strcat(jpg_path2);
    files = dir(jpg_path2);
    i=3;
    while i <= size(files,1)
    
         full_file_name = files(i).name;
         extension = full_file_name(end-3:end);
         extension_matrix(i-2) = cellstr(extension);
         i = i + 1;
         
    end
    
    
%     extension_numbers(1) = numel(strmatch('.jpg', extension_matrix,'exact'));
%     extension_numbers(2) = numel(strmatch('.pbm', extension_matrix,'exact'));
%     extension_numbers(3) = numel(strmatch('.pgm', extension_matrix,'exact'));
%     extension_numbers(4) = numel(strmatch('.png', extension_matrix,'exact'));
%     extension_numbers(5) = numel(strmatch('.ppm', extension_matrix,'exact'));
%     extension_numbers(6) = numel(strmatch('.tif', extension_matrix,'exact'));
%     extension_numbers(7) = numel(strmatch('.bmp', extension_matrix,'exact'));
    extension_numbers(1) = sum(strcmpi('.jpg', extension_matrix));
    extension_numbers(2) = sum(strcmpi('.pbm', extension_matrix));
    extension_numbers(3) = sum(strcmpi('.pgm', extension_matrix));
    extension_numbers(4) = sum(strcmpi('.png', extension_matrix));
    extension_numbers(5) = sum(strcmpi('.ppm', extension_matrix));
    extension_numbers(6) = sum(strcmpi('.tif', extension_matrix));
    extension_numbers(7) = sum(strcmpi('.bmp', extension_matrix));

    position = find(extension_numbers == max(extension_numbers));
    
    if position == 1
        
        extension = '.jpg';
        
    elseif position == 2
        
    	extension = '.pbm';   
        
    elseif position == 3
        
        extension = '.pgm';
        
    elseif position == 4
        
        extension = '.png';
        
    elseif position == 5
        
        extension = '.ppm';
        
    elseif position == 6
        
        extension = '.tif';
        
    elseif position == 7
        
        extension = '.bmp';
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Am Ende wieder entfernen!!!
    % extension = '.png';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Accounts for a "Results" folder and the "thumbs.db" file.
    if  sum(extension(1:end) == 's.db') == 4 || strcmp(extension,'ults') == 1 
        
        extension = dir(jpg_path2);
        extension = extension(3).name;
        extension = extension(end-3:end);
    
    end
    
      % Get the filenames
    filenames  = dir(fullfile(jpg_path2, sprintf('*%s',extension)));
    %filenames = dir(fullfile(jpg_path2, '*.png'));
    filenames = {filenames.name};
    %filenames = strrep(filenames,'.png','');   %entfernt die Endung .png aus filename und ersetzt sie mit nichts
    m = numel(filenames);
    m_max = m;
    bilder=[];
    vec = [1:m];
    
    if nargin == 2 
        vec = varargin{1};
        if length(vec) == 1
            vec =  [vec vec];
        end
        if vec(1) < 1
            vec(1) = 1;
        end
        if vec(2)>m
            vec(2) = m;
        end
        
         m = max(vec)- min(vec) + 1;
    end
    % Read all files.
    
    
    for k=vec
        
        % Filename
        d = filenames{k};
        % Get the file
        f = fullfile(jpg_path2 , d);
        % Generate a variable name
        dynamische_variable = genvarname(d(1:end-4)); 
        % Read the file
        bilder.(dynamische_variable) = imread(f);
        eval(sprintf('bild%s = imread(f);',dynamische_variable));
        % Check if image is rgb or grayscale
        
        if size(bilder.(dynamische_variable),3) == 3
            
            stack_new2=bilder.(dynamische_variable);
            pos = k - min(vec) + 1;
            stack_new(:,:,pos)=rgb2gray(stack_new2);
            
        else
            pos = k - min(vec) + 1;
            stack_new(:,:,pos)=bilder.(dynamische_variable);
            
        end
        
        clear bilder.(dynamische_variable);
        clearvars -EXCEPT k jpg_path2 filenames stack_new m vec m_max %stack_new2
        
    end 
    
    
    clear jpg_path;

end