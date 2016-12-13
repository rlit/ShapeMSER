function shrec_path = GetShrecFolderLocation(version)

data_root = GetDataFolderLocation();

switch version
    
    case {'2010','SHREC10'}
        dataLocation = 'http://tosca.cs.technion.ac.il/data/shrec_correspondence.zip';
        dataSize = '(~300MB)';
        name = 'SHREC 2010';
        shrec_path = fullfile(data_root,'SHREC10','shrec_correspondence');
        nFiles = (5 * 9 + 1) * 3;

    case {'2011','SHREC11'}
        dataLocation = 'http://dl.dropbox.com/u/8686064/shrec2011_hires.zip';
        dataSize = '(~185MB)';
        name = 'SHREC 2011';
        shrec_path = fullfile(data_root,'SHREC11','shrec2011_hires');
        nFiles = (5 * 11 + 1) * 1;        
        
    otherwise
        error;
end


if ~isdir(shrec_path) || numel(dir(fullfile(shrec_path,'000*.off'))) < nFiles

    button = questdlg(...
        [name ' dataset not found. Do you wish to download it? ' dataSize],...
        'Download SHREC dataset','Yes','No','No');
    
    if strcmp(button,'Yes')
        MakeDir(shrec_path)
        unzip(dataLocation,shrec_path)
        
    else
        error('Shrec data not found')
    end
    
    
end
