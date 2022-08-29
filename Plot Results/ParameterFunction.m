function [dt,t_max,Sampling,ImSize,ImPhysSize,VelFieldSize,text_font,groups,...
    savename,color_map] = ParameterFunction
%% Set constants for plotting
% Groups to compare:
% time difference between sucessive images in minutes:
dt = 3;
% max time to plot in min:
t_max = 3*13;%3*1200; 
% Sampling of MSD, Chi and Q;
Sampling = 200;
% Image size in pixels:
ImSize = [1280,960];
% Physical dimension of image in µm:
ImPhysSize = [621.23,465.80];
% Size Velocity Field in px:
VelFieldSize = [97,72];

% Text font size for plots:
text_font = 24;
% Groups to compare:
groups = [1:5];
% Savename:
savename = 'MultiCellStreaming CTLs';
% Colormap. Matrix of size nx3, with n corresponding to the number of
% groups.
color_map = [1,0,0; 0,0,1; 0,0,0; 0.5,0,0; 0.5,0.5,0.5;];


