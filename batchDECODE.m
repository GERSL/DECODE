function batchDECODE(task, tasks, tilename)
    addpath(fullfile(pwd));
    addpath(fullfile(pwd, 'Detection'));

    if ~exist('task', 'var')
        task = 1;
    end
    if ~exist('tasks', 'var')
        tasks = 1;
    end
    if ~exist('tilename', 'var')
        tilename = [];
        %% can also be defined here
        ARDTiles = globalsets.getARDTiles('decode');
%         ARDTiles = ARDTiles(1:5,:);
    else
        ARDTiles = {tilename}; % get Landsat ARD tile from .sh
    end
    
    
%    ARDTiles = {'h028v009'}; % A case tile
    %% Indicate the regions
    region = globalsets.processedRegion;
    
    %% Locate to the working folder, in which all outputing data can be found
    path_working = fullfile(globalsets.PathDECODE,region);
    
    %% Path to mask image as a reference source
    pathMask = fullfile(globalsets.PathMask,region);
    %% run COLD one by one Landsat ARD
    % The COLD function will process the stacking line data.
    for iARD = 1: length(ARDTiles)
        hv_name = ARDTiles{iARD};
        
        fprintf('Start to detect change for %s\n', hv_name);
        folderpath_tiledeocde = fullfile(path_working, hv_name);
        wetlandMask = fullfile(pathMask, [hv_name, '_mask.tif']);
        DECODE(folderpath_tiledeocde, wetlandMask,...
            'onceread', true, 'msg', true, 'mask', true,'ephemeralDetect',false,...
            'task', task ,'ntasks', tasks);
    end
end
