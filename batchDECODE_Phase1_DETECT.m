function batchDECODE_Phase1_DETECT(task, tasks, tilename)
    restoredefaultpath;
    addpath(fullfile(pwd));
    addpath(fullfile(pwd, 'Detection'));    
    if ~exist('task', 'var')
        task = 1;
    end
    if ~exist('tasks', 'var')
        tasks = 1;
    end
    
    if ~exist('tilename', 'var')
        %% can also be defined here
        ARDTiles = globalsets.getARDTiles('decode');        
    else
        ARDTiles = {tilename}; % get Landsat ARD tile from .sh
    end
    
    beginYear = 1984;
    %% Locate to the working folder, in which all outputing data can be found
    path_working = fullfile(globalsets.PathDECODE);
    
    %% Path to mask image as a reference source
    pathMask = fullfile(globalsets.PathMask);
    %% run COLD one by one Landsat ARD
    % The COLD function will process the stacked data.
%     parfor iARD = 1: length(ARDTiles)
    for iARD = 1: length(ARDTiles)    
        hv_name = ARDTiles{iARD};
        folderpath_tiledecoder = fullfile(path_working, hv_name);
        % Load the metadata
        metadata = load(fullfile(folderpath_tiledecoder,sprintf('metadata.mat')));
        nblocks = metadata.metadata.nblocks;
        fprintf('Start to detect change for %s with %d blocks (%d/%d)\n', hv_name, nblocks, iARD, length(ARDTiles));
     
        % Check existence of TSFits
        folderpath_tsf = fullfile(folderpath_tiledecoder, [globalsets.FolderDetection]);
   
        if exist(folderpath_tsf,'dir')
            numTSFits = dir(fullfile(folderpath_tsf,'rec_cg*.mat'));
            if length(numTSFits)==nblocks
                fprintf('------SKIP %s (%d) - All of the lines have been calculated\n',hv_name,nblocks);
                continue
            else
                fprintf('------Remain %d/%d of lines\n',nblocks-length(numTSFits));
            end
        else
            mkdir(folderpath_tsf);
            % make TSFitLine folder for storing coefficients of Time Series Fitting
            fprintf('------First to generate the TSFit\n'); 
        end
        
        wetlandMask = fullfile(pathMask, sprintf('%s_mask.tif',hv_name));
        
        %% First step is to detect the spectral breaks and refine the model        
        DECODE_Phase1_DETECT(folderpath_tiledecoder,folderpath_tsf,wetlandMask,...
            'onceread', true, 'orbitpath','single','msg', true, 'mask', true, 'beginYear',beginYear, ...            
            'modelRefinement',false,...
            'task', task ,'ntasks', tasks);    
        
    end
end
