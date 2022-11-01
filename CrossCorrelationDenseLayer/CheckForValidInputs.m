function CheckForValidInputsPIVlab(UsePIV,s,p,r,UseCellDivDetect,...
    CType,NetName,NetNameSeg,MinDivArea,BlockSize,MakeVideo,MedFiltSize,...
    FiltType,BackThresh,WSize,MaxTimeDiff,UseOrientationAnalysis,MSize,...
    rho,UseBM3D,UseDriftCorrect,WinSizeDrift)

% Checks inputs in the parameter function on generell validity:
% Errors for "UseModulXYZ":
if ~(isa(UsePIV,'double') || isa(UsePIV,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''UsePIV'' must be a double or integer, not a %s.',class(UsePIV))
end

if sum(size(UsePIV) ~= [1,1]) ~= 0
    error('Error in ParameterFunctionMain.m. Input of ''UsePIV'' must be a scalar of size 1x1.')
end

if ~(isa(UseCellDivDetect,'double') || isa(UseCellDivDetect,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''UseCellDivDetect'' must be a double or integer, not a %s.',class(UseCellDivDetect))
end

if sum(size(UseCellDivDetect) ~= [1,1]) ~= 0
    error('Error in ParameterFunctionMain.m. Input of ''UseCellDivDetect'' must be a scalar of size 1x1.')
end

if ~(isa(UseOrientationAnalysis,'double') || isa(UseOrientationAnalysis,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''UseOrientationAnalysis'' must be a double or integer, not a %s.',class(UseOrientationAnalysis))
end

if sum(size(UseOrientationAnalysis) ~= [1,1]) ~= 0
    error('Error in ParameterFunctionMain.m. Input of ''UseOrientationAnalysis'' must be a scalar of size 1x1.')
end

if ~(isa(UseBM3D,'double') || isa(UseBM3D,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''UseBM3D'' must be a double or integer, not a %s.',class(UseBM3D))
end

if sum(size(UseBM3D) ~= [1,1]) ~= 0
    error('Error in ParameterFunctionMain.m. Input of ''UseBM3D'' must be a scalar of size 1x1.')
end

if ~(isa(UseDriftCorrect,'double') || isa(UseDriftCorrect,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''UseDriftCorrect'' must be a double or integer, not a %s.',class(UseDriftCorrect))
end

if sum(size(UseDriftCorrect) ~= [1,1]) ~= 0
    error('Error in ParameterFunctionMain.m. Input of ''UseDriftCorrect'' must be a scalar of size 1x1.')
end

if ~(isa(MakeVideo,'double') || isa(MakeVideo,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''MakeVideo'' must be a double or integer, not a %s.',class(MakeVideo))
end

if sum(size(MakeVideo) ~= [1,1]) ~= 0
    error('Error in ParameterFunctionMain.m. Input of ''MakeVideo'' must be a scalar of size 1x1.')
end

% Error for chosen ANN and cell type:
if ~isa(NetName,'char')
    error('Error in ParameterFunctionMain.m. \nInput of ''NetName'' must be a char, not a %s.',class(NetName))
end

if ~isa(NetNameSeg,'char')
    error('Error in ParameterFunctionMain.m. \nInput of ''NetNameSeg'' must be a char, not a %s.',class(NetNameSeg))
end

if ~isa(CType,'char')
    error('Error in ParameterFunctionMain.m. \nInput of ''CType'' must be a char, not a %s.',class(CType))
end
% Drift Correction Parameters:
if isempty(WinSizeDrift) || ~(isa(WinSizeDrift,'double') || isa(WinSizeDrift,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''WinSizeDrift'' must not be empty and a double or integer, not a %s.',class(WinSizeDrift))
end

% PIV Parameters:
if isempty(s{1,2}) || ~(isa(s{1,2},'double') || isa(s{1,2},'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''s{1,2}'' must not be empty and a double or integer, not a %s.',class(s{1,2}))
end
if ~(isa(s{2,2},'double') || isa(s{2,2},'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''s{2,2}'' must be a double or integer, not a %s.',class(s{2,2}))
end

if  ~( (s{3,2}==2) || (s{3,2}==1) )
    error('Error in ParameterFunctionMain.m. Input of ''s{3,2}'' must be equal to 1 or 2.')
end

if  (~isempty(s{4,2}) && ~isa(s{4,2},'logical'))
    error('Error in ParameterFunctionMain.m. \nInput of ''s{4,2}'' must be a logical or empty, not a %s.',class(s{4,2}))
end

if  (~isempty(s{5,2}) && (~isa(s{5,2},'integer') && ~isa(s{5,2},'double')) )
    error('Error in ParameterFunctionMain.m. \nInput of ''s{5,2}'' must be a double or integer, not a %s.',class(s{5,2}))
end

if  ~( (s{6,2}==4) || (s{6,2}==3) || (s{6,2}==2) || (s{6,2}==1) )
    error('Error in ParameterFunctionMain.m. Input of ''s{6,2}'' must be equal to 1,2,3 or 4.')
end

if  (~isa(s{7,2},'integer') && ~isa(s{7,2},'double'))
    error('Error in ParameterFunctionMain.m. \nInput of ''s{7,2}'' must not be a double or integer, not a %s.',class(s{7,2}))
end

if  (~isa(s{8,2},'integer') && ~isa(s{8,2},'double'))
    error('Error in ParameterFunctionMain.m. \nInput of ''s{8,2}'' must not be a double or integer, not a %s.',class(s{8,2}))
end

if  (~isa(s{9,2},'integer') && ~isa(s{9,2},'double'))
    error('Error in ParameterFunctionMain.m. \nInput of ''s{9,2}'' must not be a double or integer, not a %s.',class(s{9,2}))
end

if ~isa(s{10,2},'char')
    error('Error in ParameterFunctionMain.m. \nInput of ''s{10,2}'' must be a char, not a %s.',class(NetName))
end

if ~(strcmp(s{10,2},'*linear') || (strcmp(s{10,2},'*spline')))
    error('Error in ParameterFunctionMain.m. \nInput of ''s{10,2}'' must be either ''*linear'' or ''*spline''.')
end

if isempty(s{11,2}) || (~isa(s{11,2},'integer') && ~isa(s{11,2},'double')) || ~( (s{11,2}==1) || (s{11,2}==0)  )
    error('Error in ParameterFunctionMain.m. Input of ''s{11,2}'' must be a double or integer of value 0 or 1 and not empty.')
end

if isempty(s{12,2}) || (~isa(s{12,2},'integer') && ~isa(s{12,2},'double')) || ~( (s{12,2}==1) || (s{12,2}==0)  )
    error('Error in ParameterFunctionMain.m. Input of ''s{12,2}'' must be a double or integer of value 0 or 1 and not empty.')
end

if isempty(s{13,2}) || (~isa(s{13,2},'integer') && ~isa(s{13,2},'double')) || ~( (s{13,2}==1) || (s{13,2}==0)  )
    error('Error in ParameterFunctionMain.m. Input of ''s{13,2}'' must be a double or integer of value 0 or 1 and not empty.')
end

if isempty(s{14,2}) || (~isa(s{14,2},'integer') && ~isa(s{14,2},'double')) || ~( (s{14,2}==1) || (s{14,2}==0)  )
    error('Error in ParameterFunctionMain.m. Input of ''s{14,2}'' must be a double or integer of value 0 or 1 and not empty.')
end

if isempty(s{15,2}) || (~isa(s{15,2},'integer') && ~isa(s{15,2},'double'))
    error('Error in ParameterFunctionMain.m. Input of ''s{15,2}'' must be a double or integer and not empty.')
end

% Pre processing parameters
if  ~((p{2,2}==0) || (p{2,2}==1) )
    error('Error in ParameterFunctionMain.m. Input of ''p{2,2}'' must be equal to 0 or 1.')
end

if (isempty(p{3,2}) && p{2,2}==1) || (~isa(p{3,2},'integer') && ~isa(p{3,2},'double'))  || p{3,2} <= 0
    error('Error in ParameterFunctionMain.m. Input of ''p{3,2}'' must be a double or integer larger than 0.')
end

if  ~((p{4,2}==0) || (p{4,2}==1) )
    error('Error in ParameterFunctionMain.m. Input of ''p{4,2}'' must be equal to 0 or 1.')
end

if (isempty(p{5,2}) && p{4,2}==1) || (~isa(p{5,2},'integer') && ~isa(p{5,2},'double'))  || p{5,2} <= 0
    error('Error in ParameterFunctionMain.m. Input of ''p{5,2}'' must be a double or integer larger than 0.')
end

if  ~((p{6,2}==0) || (p{6,2}==1) )
    error('Error in ParameterFunctionMain.m. Input of ''p{6,2}'' must be equal to 0 or 1.')
end

if  ~((p{7,2}==0) || (p{7,2}==1) )
    error('Error in ParameterFunctionMain.m. Input of ''p{7,2}'' must be equal to 0 or 1.')
end

if (isempty(p{8,2}) && p{7,2}==1) || (~isa(p{8,2},'integer') && ~isa(p{8,2},'double'))  || p{8,2} <= 0
    error('Error in ParameterFunctionMain.m. Input of ''p{8,2}'' must be a double or integer larger than 0.')
end

if isempty(p{9,2}) || ~isa(p{9,2},'double')  || p{9,2} < 0 || p{9,2} > 1
    error('Error in ParameterFunctionMain.m. Input of ''p{9,2}'' must be a double in between 0 and 1 and not empty.')
end

if isempty(p{10,2}) || ~isa(p{10,2},'double')  || p{10,2} < 0 || p{10,2} > 1
    error('Error in ParameterFunctionMain.m. Input of ''p{10,2}'' must be a double in between 0 and 1 and not empty.')
end

if p{9,2} > p{10,2}
    error('Error in ParameterFunctionMain.m. Input of ''p{9,2}'' must be smaller than ''p{10,2}''.')
end


%Post processing parameters

if ~(isa(r{1,2},'double') || isa(r{1,2},'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''r{1,2}'' must be a double or integer, not a %s.',class(r{1,2}))
end

if ~(isa(r{2,2},'double') || isa(r{2,2},'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''r{2,2}'' must be a double or integer, not a %s.',class(r{2,2}))
end

if  ~((r{4,2}==0) || (r{4,2}==1) )
    error('Error in ParameterFunctionMain.m. Input of ''r{4,2}'' must be equal to 0 or 1.')
end

if (isempty(r{5,2}) && r{4,2}==1) || (~isa(r{5,2},'integer') && ~isa(r{5,2},'double'))  || r{5,2} <= 0
    error('Error in ParameterFunctionMain.m. Input of ''r{5,2}'' must be a double or integer larger than 0.')
end

if  ~((r{6,2}==0) || (r{6,2}==1) )
    error('Error in ParameterFunctionMain.m. Input of ''r{6,2}'' must be equal to 0 or 1.')
end

if (isempty(r{7,2}) && r{6,2}==1) || (~isa(r{7,2},'integer') && ~isa(r{7,2},'double'))  || r{7,2} <= 0
    error('Error in ParameterFunctionMain.m. Input of ''r{7,2}'' must be a double or integer larger than 0.')
end

if ~(isa(r{3,2},'double') || isa(r{3,2},'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''r{3,2}'' must be a double or integer, not a %s.',class(r{3,2}))
end

if ~xor((size(r{3,2},1) == 4 && size(r{3,2},2) == 1) , (size(r{3,2},1) == 1 && size(r{3,2},2) == 4))
    error('Error in ParameterFunctionMain.m. Input of ''r{3,2}'' must be a a vector of size 4x1 or 1x4.')
end

% Cell division parameters:
if ~(isa(MinDivArea,'double') || isa(MinDivArea,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''MinDivArea'' must be a double, not a %s.',class(MinDivArea))
end

if ~(isa(BlockSize,'double') || isa(BlockSize,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''BlockSize'' must be a double, not a %s.',class(BlockSize))
end

if sum(size(BlockSize) ~= [1,2]) ~= 0
    error('Error in ParameterFunctionMain.m. Input of ''BlockSize'' must be a vector of size 1x2.')
end

if ~(isa(MedFiltSize,'double') || isa(MedFiltSize,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''MedFiltSize'' must be a double, not a %s.',class(MedFiltSize))
end

if sum(size(MedFiltSize) ~= [1,2]) ~= 0
    error('Error in ParameterFunctionMain.m. Input of ''MedFiltSize'' must be a vector of size 1x2.')
end

if ~(isa(FiltType,'double') || isa(FiltType,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''FiltType'' must be a NxN matrix of type double, not a %s.',class(FiltType))
end

if ~(isa(BackThresh,'double') || isa(BackThresh,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''BackThresh'' must be a double, not a %s.',class(BackThresh))
end

if sum(size(BackThresh) ~= [1,1]) ~= 0
    error('Error in ParameterFunctionMain.m. Input of ''BackThresh'' must be a scalar of size 1x1.')
end

if ~(isa(WSize,'double') || isa(WSize,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''WSize'' must be a double corresponding to half the diameter of the images used for training of the ANN called in ''NetName'', not a %s.',class(WSize))
end

if sum(size(WSize) ~= [1,1]) ~= 0
    error('Error in ParameterFunctionMain.m. Input of ''WSize'' must be a scalar of size 1x1.')
end

if ~(isa(MaxTimeDiff,'double') || isa(MaxTimeDiff,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''MaxTimeDiff'' must be an integer >0, not a %s.',class(MaxTimeDiff))
end

if sum(size(MaxTimeDiff) ~= [1,1]) ~= 0
    error('Error in ParameterFunctionMain.m. Input of ''MaxTimeDiff'' must be a scalar of size 1x1.')
end

if ~(isa(MSize,'double') || isa(MSize,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''MSize'' must be a double, not a %s.',class(MSize))
end

if sum(size(MSize) ~= [1,1]) ~= 0
    error('Error in ParameterFunctionMain.m. Input of ''MSize'' must be a scalar of size 1x1.')
end

if ~(isa(rho,'double') || isa(rho,'integer'))
    error('Error in ParameterFunctionMain.m. \nInput of ''rho'' must be a double, not a %s.',class(rho))
end

if sum(size(rho) ~= [1,1]) ~= 0
    error('Error in ParameterFunctionMain.m. Input of ''rho'' must be a scalar of size 1x1.')
end
