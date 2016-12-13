function  data_root = GetDataFolderLocation


code_root = fileparts(fileparts(mfilename('fullpath')));
proj_root = fileparts(code_root);
data_root = fullfile(proj_root,'data');


