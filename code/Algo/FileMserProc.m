function [mserRes shape_loaded] = FileMserProc(shape_file_name,params)
%
%
%

%% init
tStart = now;

shape_folder = params.data_path;
shape_full_path = fullfile(shape_folder,shape_file_name); 
assert(exist(shape_full_path,'file')>0,'shape file "%s" not found',shape_file_name)
[~,shape_name] = fileparts(shape_full_path);

%% try to load MSERs from mat-file
p             = fullfile(shape_folder ,'mser_res',params.benchName);
mser_res_file = fullfile(p,[shape_name '.mser.mat']);
if params.load_msers && exist(mser_res_file,'file')
    VerboseDisp('loading MSER result for shape %s',shape_name)
    mserRes = load(mser_res_file);
    if nargout > 1
        shape_loaded = load_shape(shape_full_path);
    end
    return
elseif ~isdir(p),mkdir(p)
end
clear p

%% load shape
VerboseDisp('loading shape %s',shape_name)
shape_loaded = load_shape(shape_full_path);
% TODO - save a ".mat" copy of the shape

%% remesh if necessary
nVertices = length(shape_loaded.X);
if nVertices > params.maxVertices
    shape = remesh(shape_loaded, struct(...
        'verbose',false,...
        'vertices', params.maxVertices));
    
    % remember what was the original vertex
    tree = ann('init',[shape_loaded.X shape_loaded.Y shape_loaded.Z]');
    shape.origIdx = ann('search',tree,[shape.X shape.Y shape.Z]',1);
    ann('deinit', tree); clear tree;
    
else
    shape = shape_loaded;
    shape.origIdx = 1:nVertices;

end


%%
mserRes = ShapeMserProc(shape,params,shape_name);
nMsers  = numel(mserRes.msers);

%%
VerboseDisp('Finalizing...')
save(mser_res_file, '-struct','mserRes','shape','msers')


VerboseDisp('finished shape %s, found %d MSERs, took %s \n',....
    shape_name,nMsers,datestr(now-tStart,13) )

