function CreateShrecVolumeData(shrecVersion,nPoints)
if nargin < 1
    shrecVersion = 'SHREC11';

end
shapes2dPath = GetShrecFolderLocation(shrecVersion);
if nargin < 2
    nPoints = [150 135 110];
end

corrPath2d  = fullfile(shapes2dPath,'corr');
assert(isdir(corrPath2d),'ground-truth directory does not exist. contact me to get it')

volPath = fullfile(shapes2dPath,'volume');
corrPathVol = fullfile(volPath,'corr');
MakeDir(corrPathVol)

% voxelize nulls & then all transformed shapes
for iNull = 1:numel(dir(fullfile(shapes2dPath,'*null.0.off')))
    numberStr = sprintf('%04d',iNull);
    nullName  = sprintf('%04d.null.0.off',iNull);
    nullShape = load_shape(fullfile(shapes2dPath,nullName));
    
    VerboseDisp('voxelizing %s...',nullName)
    nullVol = struct;
    [nullVol.G nullVol.d nullVol.g nullVol.points nullVol.index] = mesh2volume(nullShape,nPoints(iNull));
    nullD = nullVol.d;
    SaveVolShape(volPath,nullName,nullVol)
    nullVol = AppendFields(nullVol,nullShape);
    
    VerboseDisp('creating symmetry ground-truth for %s...',nullName)
    sym2dFile = which(sprintf('%04d.null.0.sym.mat',iNull));
    sym2d  = load(sym2dFile);
    symVol = GenerateVolCorr(nullVol,nullVol,sym2d.isym);
    SaveVolCorr(strrep(sym2dFile,'sym','vol.sym'),symVol)
    
    dirShapes = dir(fullfile(shapes2dPath,[numberStr '*.off']));
    dirShapes = setdiff({dirShapes.name},nullName);
    for iShape = 1:numel(dirShapes)
        shapeName = dirShapes{iShape};
        currShape = load_shape(fullfile(shapes2dPath,shapeName));
        VerboseDisp('voxelizing %s...',shapeName)
        
        shapeVol = struct;
        if isempty(strfind(shapeName,'.scale.'))
            [shapeVol.G shapeVol.d shapeVol.g shapeVol.points shapeVol.index] = mesh2volume_D(currShape,nullD);
        else
            iso1 = load_shape(fullfile(shapes2dPath,[numberStr '.isometry.1.off']));            
            [shapeVol.G shapeVol.d shapeVol.g shapeVol.points shapeVol.index] = mesh2volume_D(currShape,nullD * ScaleRatio(iso1,currShape));
        end
        SaveVolShape(volPath,shapeName,shapeVol)
        shapeVol = AppendFields(shapeVol,currShape);
        
        VerboseDisp('creating ground-truth for %s...',shapeName)
        corr2d      = load(fullfile(corrPath2d,strrep(shapeName,'off','corr.lut.mat')));
        corrVol     = GenerateVolCorr(shapeVol,nullVol,corr2d.lut);
        corrVolName = strrep(shapeName,'off','vol.corr.mat');
        SaveVolCorr(fullfile(corrPathVol,corrVolName),corrVol)
    end
    
end


%% calc laplacian
% CalcVolLaplacian(volPath)


%%



end
function s3d = AppendFields(s3d,s2d)
s3d.X    = s2d.X;
s3d.Y    = s2d.Y;
s3d.Z    = s2d.Z;
s3d.TRIV = s2d.TRIV;
end
function SaveVolShape(volPath,sName,s)
save(fullfile(volPath,strrep(sName,'off','vol.mat')),'-struct','s','G','d','g','points','index')
end
function SaveVolCorr(n,volCorr)
save(n,'volCorr')
end
function ratio = ScaleRatio(original,scaled)
mX = min(original.X);  MX = max(original.X);
mY = min(original.Y);  MY = max(original.Y);
mZ = min(original.Z);  MZ = max(original.Z);
d1 = max([(MX-mX) (MY-mY) (MZ-mZ)]);

mX = min(scaled.X);  MX = max(scaled.X);
mY = min(scaled.Y);  MY = max(scaled.Y);
mZ = min(scaled.Z);  MZ = max(scaled.Z);
d2 = max([(MX-mX) (MY-mY) (MZ-mZ)]);

ratio = d2/d1;
end