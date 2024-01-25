function batchDECODE_Phase3_S1_IsolatedRemoval_SingleTile(task,ntasks)
%% Remove the isolated patch and misclassification of the tidal wetlands by using Minimum Mapping Unit and contextual analysis
    
    if ~exist('task', 'var')
        task = 1;
    end
    if ~exist('ntasks', 'var')
        ntasks = 1;
    end
    restoredefaultpath;
    folderpath_mfile = fileparts(mfilename('fullpath'));
    addpath(folderpath_mfile);
    addpath(fullfile(folderpath_mfile,'Packages', 'GRIDobj'));
    addpath(fullfile(folderpath_mfile,'Export'));

    ARDTiles = globalsets.getARDTiles('decode');
    
    path_working = globalsets.PathDECODE;

    mapType = 'cover_MMU';

    years = 1986:2020;

    %% Path to mask image as a reference source
    pathMask = fullfile(globalsets.PathMask);
    
    % HPC process    
    patchSize = 100; % 5 Pixels
     % HPC process
    num_t = length(ARDTiles);
    tasks_per = ceil(num_t/ntasks);
    start_i = (task-1)*tasks_per + 1;
    end_i = min(task*tasks_per, num_t);
    tic
    for i_task = start_i: end_i
        hv_name = ARDTiles{i_task};
        fprintf('Begin to process %s in task %d/%d (Totally %d)\n',hv_name,i_task,end_i,end_i-start_i+1);
        pathMap = fullfile(path_working, hv_name,[globalsets.FolderMap]);

        maskImage = fullfile(pathMask, [hv_name, '_mask.tif']); 
        for yr = years
            coverImage = fullfile(pathMap,sprintf('%s_%d.tif',mapType,yr));
            patchImage = replace(coverImage,'MMU','Patch');
            if isfile(patchImage)
                fprintf('Skip %s\n',patchImage);        
                continue
            else              
                fprintf('Remove isolated patches... \n');
                DECODE_Phase3S_removeIsolatedTidalWetlandOverestimation(coverImage,patchImage,maskImage,'patchSize', patchSize);   
            end    
        end
    end
    fprintf('--------------------------------------\n');
end