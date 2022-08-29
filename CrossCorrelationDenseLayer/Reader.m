function [Stack, m,m_max]=Reader(FilePath,varargin)
% the second argument needs to be a vector [x y] that defines the start and
% end image for reading.

    if length(varargin) > 2
        disp('Too many input arguments.')
        return
    end
     
    % Get files and extensions in the path:
    FilePath = strcat(FilePath);
    Files = dir(FilePath);
    % Remove directories:
    Files([Files.isdir]==1) = [];
    i=1;
    while i <= size(Files,1)    
         FullFileName = Files(i).name;
         Extension = FullFileName(end-3:end);
         ExtensionMatrix(i) = cellstr(Extension);
         i = i + 1;         
    end
      
    % Get the most frequent extension of jpg, pbm, pgm, png, ppm, tif and
    % bmp and use it for later reading:
    ExtensionNumbers(1) = sum(strcmpi('.jpg', ExtensionMatrix));
    ExtensionNumbers(2) = sum(strcmpi('.pbm', ExtensionMatrix));
    ExtensionNumbers(3) = sum(strcmpi('.pgm', ExtensionMatrix));
    ExtensionNumbers(4) = sum(strcmpi('.png', ExtensionMatrix));
    ExtensionNumbers(5) = sum(strcmpi('.ppm', ExtensionMatrix));
    ExtensionNumbers(6) = sum(strcmpi('.tif', ExtensionMatrix));
    ExtensionNumbers(7) = sum(strcmpi('.bmp', ExtensionMatrix));
    position = find(ExtensionNumbers == max(ExtensionNumbers));    
    if position == 1        
        Extension = '.jpg';        
    elseif position == 2        
    	Extension = '.pbm';          
    elseif position == 3        
        Extension = '.pgm';        
    elseif position == 4        
        Extension = '.png';        
    elseif position == 5        
        Extension = '.ppm';        
    elseif position == 6        
        Extension = '.tif';        
    elseif position == 7        
        Extension = '.bmp';        
    end
    
    % Get the filenames
    FileNames  = dir(fullfile(FilePath, sprintf('*%s',Extension)));
    FileNames = {FileNames.name};
    m = numel(FileNames);
    m_max = m;
    Image=[];
    vec = [1:m];
    % If additional argument is given restrict reading to the image numbers
    % given
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
        d = FileNames{k};
        % Get the file
        f = fullfile(FilePath , d);
        % Generate a variable name
        DynamischeVariable = genvarname(d(1:end-4)); 
        % Read the file
        Image.(DynamischeVariable) = imread(f);
        eval(sprintf('bild%s = imread(f);',DynamischeVariable));
        % Check if image is rgb or grayscale. If rgb transform to grayscale        
        if size(Image.(DynamischeVariable),3) == 3            
            StackTemp=Image.(DynamischeVariable);
            pos = k - min(vec) + 1;
            Stack(:,:,pos)=rgb2gray(StackTemp);            
        else
            pos = k - min(vec) + 1;
            Stack(:,:,pos)=Image.(DynamischeVariable);            
        end        
        clear bilder.(dynamische_variable);
        clearvars -EXCEPT k FilePath FileNames Stack m vec m_max %stack_new2
    end       
end