function [s_size,im_size,ImPhysSize,dt,CSize,OvThresh,WSize,...
    CenterSpeed,SubPxResolution,Sampling,MaxVisibleTime] = ParameterFunctionMain()

%% Set Parameters
% Define points of interest. Used to create a grid with the given spacing.
% Points on the grid will be used for particle tracking, calculation of the
% order parameter, etc. Low values yield a high resolution but may increase
% computation time strongly. 
s_size = 80;
% image size:
im_size = [960,1280];
% Physical dimension of image in µm:
ImPhysSize = [465.80,621.23];
% Time between successive images in minutes:
dt = 3; 
% Maximal number of images a single cell division is considered to be a
% true positive and not a cell extrusion, apoptosis event, etc.
MaxVisibleTime = 150/dt; 
% Estimate of cell diameter (in px):
CSize = 80;
% Overlap threshold used for calculating the order parameter and 4-point
% susceptibility:
OvThresh = 0.1;
% Center velocity field if value is 1, do not center otherwise:
% Centering is done as follows: Vx = Vx-mean(Vx(:)), etc
CenterSpeed = 0; 
% Sub-Pixel resolution for the tracking matrix. Values in range [0,1]. 1 =
% no sub-pixel resolution. 0.1 = 10x resolution. Attention, low values
% might drastically slow down calculations.
SubPxResolution = 1; 
% Calculate time dependent variables only at every x-th time point/image.
% See below for an example of how the calculation works.
% Low values may lead to large resulting data files and may slow down the
% calculations.
Sampling = 100; 
% Window Size -> number of images used to determine time sensitive
% parameters, such as MSD, order parameter Q and 4-point susceptibility
% chi. Denote that all these parameters will also be calculated for the
% full measurement time as well.
WSize = 41; 

% example how WSize and Sampling interact:
% Sampling = 100; WSize = 41; 
% Using these values all relevant parameters will be calculated for
% image 1,101,201,etc. As MSD, Q, Chi, etc. do need "temporal averaging",
% WSize is used to determine the time frame. Thus, parameters are
% calculated for images 1,2,...42, than for 101,...142, etc


