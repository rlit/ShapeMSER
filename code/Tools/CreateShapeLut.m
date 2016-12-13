function lut = CreateShapeLut(shape_full_path,shape,null_shape)

if isempty(shape_full_path),error('shape file not found'),end

[shape_path,shape_name] = fileparts(shape_full_path);

corr_path = fullfile(shape_path ,'corr');
assert(isdir(corr_path),'ground-truth directory does not exist. contact me to get it')

lutPath = fullfile(corr_path,[shape_name '.corr.lut.mat']);

if ~exist(lutPath,'file')

    corr_path = fullfile(corr_path,[shape_name '.corr']);
    
    corr_data = load(corr_path);
    lut = CreateLutFromCorr(null_shape,shape,corr_data);
    
    save(lutPath,'lut')
    
else
    loadRes = load(lutPath);
    lut = loadRes.lut;
end
