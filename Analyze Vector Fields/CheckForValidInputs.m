function CheckForValidInputs(s_size,im_size,ImPhysSize,dt,CSize,OvThresh,WSize,...
    CenterSpeed,SubPxResolution,Sampling,MaxVisibleTime)

if ~(isa(s_size,'double') || isa(s_size,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''s_size'' must be a double or integer, not a %s.',class(s_size))
end

if sum(size(s_size) ~= [1,1]) ~= 0 
    error('Error in ParameterFunctionMain.m. Input of ''s_size'' must be a scalar of size 1x1.')
end

if ~(isa(im_size,'double') || isa(im_size,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''im_size'' must be a double or integer, not a %s.',class(im_size))
end

if sum(size(im_size) ~= [1,2]) ~= 0 
    error('Error in ParameterFunctionMain.m. Input of ''im_size'' must be a vector of size 1x2.')
end

if s_size>min(im_size)
    error('Error in ParameterFunctionMain.m. ''s_size'' must be smaller than the minimum of im_size.')
end

if ~(isa(ImPhysSize,'double') || isa(ImPhysSize,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''ImPhysSize'' must be a double or integer, not a %s.',class(ImPhysSize))
end

if sum(size(ImPhysSize) ~= [1,2]) ~= 0 
    error('Error in ParameterFunctionMain.m. Input of ''ImPhysSize'' must be a vector of size 1x2.')
end

if ~(isa(dt,'double') || isa(dt,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''dt'' must be a double or integer, not a %s.',class(dt))
end

if sum(size(dt) ~= [1,1]) ~= 0 
    error('Error in ParameterFunctionMain.m. Input of ''dt'' must be a scalar of size 1x1.')
end

if ~(isa(MaxVisibleTime,'double') || isa(MaxVisibleTime,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''MaxVisibleTime'' must be a double or integer, not a %s.',class(MaxVisibleTime))
end

if sum(size(MaxVisibleTime) ~= [1,1]) ~= 0 
    error('Error in ParameterFunctionMain.m. Input of ''MaxVisibleTime'' must be a scalar of size 1x1.')
end

if ~(isa(CSize,'double') || isa(CSize,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''CSize'' must be a double or integer, not a %s.',class(CSize))
end

if sum(size(CSize) ~= [1,1]) ~= 0 
    error('Error in ParameterFunctionMain.m. Input of ''CSize'' must be a scalar of size 1x1.')
end

if ~(isa(OvThresh,'double') || isa(OvThresh,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''OvThresh'' must be a double or integer, not a %s.',class(OvThresh))
end

if sum(size(OvThresh) ~= [1,1]) ~= 0 
    error('Error in ParameterFunctionMain.m. Input of ''OvThresh'' must be a scalar of size 1x1.')
end

if OvThresh<=0
    error('Error in ParameterFunctionMain.m. Input of ''OvThresh'' must be larger than zero.')
end

if ~(isa(WSize,'double') || isa(WSize,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''WSize'' must be a double or integer, not a %s.',class(WSize))
end

if sum(size(WSize) ~= [1,1]) ~= 0 
    error('Error in ParameterFunctionMain.m. Input of ''OvThresh'' must be a scalar of size 1x1.')
end

if WSize<=1
    error('Error in ParameterFunctionMain.m. Input of ''WSize'' must be larger than one.')
end

if ~(isa(CenterSpeed,'double') || isa(CenterSpeed,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''CenterSpeed'' must be a double or integer, not a %s.',class(CenterSpeed))
end

if sum(size(CenterSpeed) ~= [1,1]) ~= 0 
    error('Error in ParameterFunctionMain.m. Input of ''CenterSpeed'' must be a scalar of size 1x1.')
end

if ~(isa(SubPxResolution,'double') || isa(SubPxResolution,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''SubPxResolution'' must be a double or integer, not a %s.',class(SubPxResolution))
end

if sum(size(SubPxResolution) ~= [1,1]) ~= 0 
    error('Error in ParameterFunctionMain.m. Input of ''SubPxResolution'' must be a scalar of size 1x1.')
end

if (SubPxResolution>1 ||SubPxResolution<=0)
    error('Error in ParameterFunctionMain.m. Input of ''SubPxResolution'' must be in range [0,1].')
end

if ~(isa(Sampling,'double') || isa(Sampling,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''Sampling'' must be a double or integer, not a %s.',class(Sampling))
end

if sum(size(Sampling) ~= [1,1]) ~= 0 
    error('Error in ParameterFunctionMain.m. Input of ''Sampling'' must be a scalar of size 1x1.')
end
