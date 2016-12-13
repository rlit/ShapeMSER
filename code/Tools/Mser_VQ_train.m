function Mser_VQ_train
%% should run on 64-bit
if ~strcmp(computer,'PCWIN64')
    %error('!!! You better run this in 64 bit !!!!! ' )
end
%% get data - assuming all needed files are in ..\data
code_path = fileparts(mfilename('fullpath'));
root_path = fileparts(code_path);
data_path = fullfile(root_path,'data');
shapesPath = fileparts(which('0001.null.0.off'));
if isempty(shapesPath)
    addpath(genpath(data_path))
    shapesPath = fileparts(which('0001.null.0.off'));
end
InitPath;


%%
timeSamples = 2.^(4:0.5:7);

shapeFiles = dir(fullfile(shapesPath,'*.off'));
shapeFiles = {shapeFiles.name};
desc       = cell(numel(shapeFiles),1);
siDesc     = cell(numel(shapeFiles),1);
for iShape = 1:numel(shapeFiles)
    shapeName = shapeFiles{iShape};
    
    laplacian_path = strrep(shapeName,'off','laplacian.mat');
    laplacianRes = load(laplacian_path);
    
    evals = laplacianRes.evals;
    evecs = laplacianRes.evecs;
    
    siHksDesc = Eigen2SIHKS(evals,evecs);
    %siHksDesc = normalize(siHksDesc, 'l2', 2);

    siDesc{iShape} = siHksDesc;
    
    %     % % save...?    
    %     area = laplacianRes.area;
    %     save(laplacian_path,'evals','evecs','area','siHksDesc');
    
    %     isZeroVal = evals < eps(max(evals));
    %     evals = evals(~isZeroVal);
    %     evecs = evecs(:,~isZeroVal);
    
    desc{iShape} = Eigen2HKS(evals,evecs,timeSamples);
    
    
    
    
end

%%
desc   = cell2mat(desc);
siDesc = cell2mat(siDesc);
vs = 6:2:12;% vocabulary size
for ii=1:numel(vs)
    disp(['starting ' int2str(vs(ii)) ])
    %% HKS
    [idx, vocab,sum_dist] = kmeans(desc, vs(ii), ...
        'maxiter', 1000, ...
        'display', 1, ...
        'replicates', 5, ...
        'randstate', 0, ...
        'outlierfrac', 1e-3); %#ok<ASGLU>
    
    % Compute sigma
    tree = ann('init', vocab');
    [sig, min_dist] = ann('search', tree, desc', 1);
    sigma = median(min_dist);
    ann('deinit', tree);
    
    vocab_file_name = ['vocab_vs_' int2str(vs(ii)) '.mat'];
    save(fullfile(data_path,vocab_file_name), ...
        'vocab', ...
        'sigma', ...
        'sum_dist', ...
        'timeSamples')
    
    %% SI-HKS
    [idx, vocab,sum_dist] = kmeans(siDesc, vs(ii), ...
        'maxiter', 1000, ...
        'display', 1, ...
        'replicates', 5, ...
        'randstate', 0, ...
        'outlierfrac', 1e-3); %#ok<ASGLU>
    
    % Compute sigma
    tree = ann('init', vocab');
    [sig, min_dist] = ann('search', tree, desc', 1);
    sigma = median(min_dist);
    ann('deinit', tree);
    
    vocab_file_name = ['vocab_SI_vs_' int2str(vs(ii)) '.mat'];
    save(fullfile(data_path,vocab_file_name), ...
        'vocab', ...
        'sigma', ...
        'sum_dist')
        
    
end
ann('close');

