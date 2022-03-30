function plot_inv(varargin)
%   plot_inv is a function to generate plots x and y axes 
%   with specified errors modified from codes written by Nils Sjöberg 
%   (http://www.mathworks.com/matlabcentral/fileexchange/5444-xyerrorbar)
%  
%   plot_inv(data) plots "data" as bar diagram with default settings.
%   plot_inv(data, standard) plots the deviation "standard" in y direction
%   as well.
%   plot_inv(data, standard, name) plots the data with given xlabel in
%   namestring of the form ['asd','',...,''] with length(data)elements.
%   plot_inv(data, standard, name, significance) plots the data with a
%   given significance string of the form ['*','',...,''] with length(data)
%   elements with the significance string on top of the error bars.
%   plot_inv(data, standard, name, significance,group_name) plots the names
%   of experimental groups in the bars of the bar plot. Length and style of
%   the string is identical to the "signifcance" variable.
%   plot_inv(data, standard, name, significance,group_name, text_font) adds
%   a specified text font.
%   plot_inv(data, standard, name, significance,group_name, text_font,n,
%   gaps) plots data introducing gaps between data points given by "gaps"
%   and "n" gives the maximal number of elements per group. The gaps
%   between data points are of the same size as one bar. 
  
if length(varargin)==1
    error('Insufficient number of inputs');
    return;
end
data = varargin{1};


if length(varargin) == 1
   
    standard = zeros(1,length(data));
    n = length(data);
    name = [];
    title2 = [];
    significance = cell(length(data),1);
    group_name = [];
    text_font = 24;
    gaps = [];
    y_axis_label = 'Normalized Covered Area';
    
elseif length(varargin) == 2
    
    standard = varargin{2};
    n = length(data);
    name = [];
    title2 = [];
    significance = cell(length(data),1);
    group_name = [];
    text_font = 24;
    gaps = [];  
    y_axis_label = 'Normalized Covered Area';
    
elseif length(varargin) == 3

    standard = varargin{2};
    n = length(data);
    name = varargin{3};
    title2 = [];
    significance = cell(length(data),1);
    group_name = [];
    text_font = 24;
    gaps = [];
    y_axis_label = 'Normalized Covered Area';
    
elseif length(varargin) == 4
    
    standard = varargin{2};
    n = length(data);
    name = varargin{3};
    title2 = varargin{4};
    significance = cell(length(data),1);
    group_name = [];
    text_font = 24;
    gaps = [];
    y_axis_label = 'Normalized Covered Area';
        
elseif length(varargin) == 5
    
    standard = varargin{2};
    n = length(data);
    name = varargin{3};
    title2 = varargin{4};
    significance = varargin{5};
    group_name = [];
    text_font = 24;
    gaps = [];
    y_axis_label = 'Normalized Covered Area';

elseif length(varargin) == 6
        
    standard = varargin{2};
    n = length(data);
    name = varargin{3};
    title2 = varargin{4};
    significance = varargin{5};
    group_name = varargin{6};
    text_font = 24;
    gaps = [];
    y_axis_label = 'Normalized Covered Area';
    
elseif length(varargin) == 7
        
    standard = varargin{2};
    n = length(data);
    name = varargin{3};
    title2 = varargin{4};
    significance = varargin{6};
    group_name = varargin{6};
    text_font = varargin{7};
    gaps = [];
    y_axis_label = 'Normalized Covered Area';
        
elseif length(varargin) == 8   
    
    error('Inout must contain a gap string and a group number')
    return;
    
elseif length(varargin) == 9
        
    standard = varargin{2};
    name = varargin{3};
    title2 = varargin{4};
    significance = varargin{5};
    group_name = varargin{6};
    text_font = varargin{7};
    n = varargin{8};
    gaps = varargin{9};  
    gaps = gaps-1;
    y_axis_label = 'Normalized Covered Area';
    
elseif length(varargin) == 10
    
    standard = varargin{2};
    name = varargin{3};
    title2 = varargin{4};
    significance = varargin{5};
    group_name = varargin{6};
    text_font = varargin{7};
    n = varargin{8};
    gaps = varargin{9};  
    gaps = gaps-1;
    y_axis_label = varargin{10};
    
elseif length(varargin) == 11
    
    standard = varargin{2};
    name = varargin{3};
    title2 = varargin{4};
    significance = varargin{5};
    group_name = varargin{6};
    text_font = varargin{7};
    n = varargin{8};
    gaps = varargin{9};  
    gaps = gaps-1;
    y_axis_label = varargin{10};
    leg = varargin{11};
    
elseif length(varargin) > 11
    
    error('Too many input arguments!')
    return;
    
end

if length(data) ~= length(standard)
    error('Data and standard deviation must have the same number of elements!')
    return;
end
if length(data) ~= length(significance)
    error('Data and significance must have the same number of elements!')
    return;
end

colorbar = ones(n,3);
for i = 1:3
colorbar(:,i) =  transpose([1:n]/n)  .*colorbar(:,i);
end
colorbar = colorbar - min(colorbar(:));

% c1 = [0 0 0]; % Black Color for segment 1.
% c2 = [0.2 0.2 0.2];
% c3 = [0.4 0.4 0.4];
% c4 = [0.6 0.6 0.6];

var = [];
var = data(:)+standard(:);
var = var(:);
data_num = length(data);
figure('units','normalized','outerposition',[0 0 1 1]);
hold all
for i=1:data_num
    
  
  
  if isempty(gaps) == 1 || i<=min(gaps)  
    % Make all errorbars of equal formating (width of errorbars dependent
    % on length of x-axis)
    col = colorbar(i,:);
    a = 1:data_num+length(gaps);
    b = zeros(1,data_num+length(gaps));
    b(i) = data(i);
    c = zeros(1,data_num+length(gaps));
    c(i) = standard(i);
    errorbar(a,b,c, 'linestyle','none', 'LineWidth', 2.005,'Color',col,'HandleVisibility','off');
    %errorbar(i,data(i,2),standard(i,1), 'linestyle','none', 'LineWidth', 0.005,'Color',col);
    h = bar(i,data(i));
    var(i) = data(i)+standard(i);
    text(i,(var(i)+0.03*max(var)), ...
    significance(i),'horizontalalignment','center','verticalalignment',...
    'middle','Rotation',0,'FontSize',text_font,'Color',[0 0 0])
  elseif i>min(gaps)    
    
    % Make all errorbars of equal formating (width of errorbars dependent
    % on length of x-achses)   
    gap_num = find(gaps<i,1,'last');
    col = colorbar(i-gaps(gap_num),:);
    a = 1:data_num+length(gaps);
    b = zeros(1,data_num+length(gaps));
    b(i+gap_num) = data(i);
    c = zeros(1,data_num+length(gaps));
    c(i+gap_num) = standard(i);
    errorbar(a,b,c, 'linestyle','none', 'LineWidth', 2.005,'Color',col,'HandleVisibility','off'); 
    %errorbar(i+1,data(i-3,3),standard(i-3,2), 'linestyle','none', 'LineWidth', 0.005,'Color',col);
    h = bar(i+gap_num,data(i));
    var(i) = data(i)+standard(i);    
    text(i+gap_num,(var(i)+0.03*max(var)), ...
    significance(i),'horizontalalignment','center','verticalalignment',...
    'middle','Rotation',0,'FontSize',text_font,'Color',[0 0 0])
  end
  
  set(h, 'FaceColor', col,'LineWidth',0.5,'EdgeColor','black')

end
set(gca, 'XTickLabel', '','Fontsize',text_font)
testvar = max(data(:)+standard(:));
ymax = max([0,1.1*testvar]);
%ymax = 1.1*max(var);
testvar= min(data(:)-standard(:));
ymin = min([0,1.1*testvar]);
xlabetxt{1,length(name)+length(gaps)} = [];
group_name2{1,length(group_name)+length(gaps)} = [];
%group_name2 = cell(1,length(group_name)+length(gaps));

for i = 1:length(name)

    if isempty(name) == 1
    
    elseif isempty(gaps) == 1 || i<=min(gaps)  
        xlabetxt(i) = name(i); 
        
    elseif i>min(gaps)
        gap_num = find(gaps<i,1,'last');
        xlabetxt(i+gap_num) = name(i); 
        
    end
    
end
        
for i = 1:length(group_name)

    if isempty(group_name) == 1
    
    elseif isempty(gaps) == 1 || i<=min(gaps)  
        group_name2(i) = group_name(i); 
        
    elseif i>min(gaps)
        gap_num = find(gaps<i,1,'last');
        group_name2(i+gap_num) = group_name(i); 
        
    end
    
end              

ylim([ymin ymax*1.30]); 
%ylim([0 25]); 
ypos = -max(ylim)/50;
text(1:data_num+length(gaps),repmat(ypos,data_num+length(gaps),1), ...
     xlabetxt','horizontalalignment','right','verticalalignment','middle','Rotation',15,'FontSize',text_font-2)
 
text(1:data_num+length(gaps),repmat(-ypos,data_num+length(gaps),1), ...
     group_name2,'horizontalalignment','center','verticalalignment','middle','Rotation',0,'FontSize',text_font-2,'Color',[1 0 0])
 
ylabel(y_axis_label,'FontSize',text_font)
title(title2,'FontSize',text_font+6);
 
if length(varargin) == 11
legend(leg,'FontSize',text_font,'Location','NorthWest')
end

% Plot a line inside the bar plot (this works, because a "hold all" is
% still set):
% plot([1,2],[3 3], 'r')


