function InitPath


code_root = fileparts(mfilename('fullpath'));

addpath(code_root)

addpath(fullfile(code_root,'common'))
addpath(fullfile(code_root,'Tools'))
% addpath(fullfile(code_root,'volume'))
addpath(fullfile(code_root,'Bench'))
addpath(fullfile(code_root,'Algo'))

addpath(fullfile(code_root,'io'))
% addpath(fullfile(code_root,'io','obj_io'))
% addpath(fullfile(code_root,'io','ply_io'))

% addpath(fullfile(code_root,'3rdParty','FastMarching_version3b'))
% addpath(fullfile(code_root,'3rdParty','FastMarching_version3b','shortestpath'))
% addpath(fullfile(code_root,'3rdParty','FastMarching_version3b','functions'))

data_root = GetDataFolderLocation();
addpath(data_root)

% 
% if isempty(which('0001.null.0.off'))
%     warning('adding HUGE path....')
%     addpath(genpath(data_root))
% end