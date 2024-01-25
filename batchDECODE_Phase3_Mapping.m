function batchDECODE_Phase3_Mapping(task,ntasks)
    if ~exist('task', 'var')
        task = 1;
    end
    if ~exist('ntasks', 'var')
        ntasks = 1;
    end

    folderpath_mfile = fileparts(mfilename('fullpath'));
    addpath(folderpath_mfile);
    addpath(fullfile(folderpath_mfile, 'Export'));
    addpath(fullfile(folderpath_mfile, 'Characterization'));
    addpath(fullfile(folderpath_mfile, 'Detection'));
    addpath(fullfile(folderpath_mfile,'Packages','GRIDobj'));
    
    mapTypes = {'disturbance','cover'};
    MMU = 4; % Minimum mapping unit of landsat pixels
    ARDTiles = globalsets.getARDTiles('decode');

    isCheckExist = 1; % If check, then dont rewrite    

    % Locate to the DECODE results
    path_working = globalsets.PathDECODE;    

    folderTSFit = [globalsets.FolderDetection,'_R'];
    
    %% Indicate the regionssqueu 
    years = 1986:2020;

    %% Path to mask image as a reference source
    pathMask = fullfile(globalsets.PathMask);
    
    % HPC process
    num_t = length(ARDTiles);
    tasks_per = ceil(num_t/ntasks);
    start_i = (task-1)*tasks_per + 1;
    end_i = min(task*tasks_per, num_t);
    tic
    for i_task = start_i: end_i
        hv_name = ARDTiles{i_task};
        fprintf('Begin to process %s in task %d/%d (Totally %d)\n',hv_name,i_task,end_i,end_i-start_i+1);
        folderpath_tileDecode = fullfile(path_working, hv_name);           
        pathMap = fullfile(path_working, hv_name,[globalsets.FolderMap]);
        maskImage = fullfile(pathMask, [hv_name, '_mask.tif']);
 
        %% 3-1 export annual map
        fprintf(' -> 1 Export Annual Maps (%d min)... \n',round(toc/60));
        % Check the existence of the maps
        mapCount = 0;
        for i_yr = 1: length(years)
            yr = years(i_yr);
            for i_mt = 1:length(mapTypes)
                if isfile(fullfile(pathMap, sprintf('%s_%d.tif',mapTypes{i_mt},yr)))
                    mapCount = mapCount + 1;
                end
            end
        end
        if isCheckExist && mapCount == length(mapTypes)*length(years)
            fprintf('\nHaving generated all the maps in this tile\n Skip the tile\n');
        else
            %% Export four kinds of maps for mutiple years
            DECODE_Phase3_1_exportAnnualMap(folderpath_tileDecode,folderTSFit,maskImage, years);
        end
        
        %% 3-2 export MMU cover map
        fprintf(' -> 2 MMU Land cover (%d min)... \n',round(toc/60));
        for yr = years
            fprintf('%d ...    Land Cover',yr);
            DECODE_Phase3_2_MMUCoverMap(pathMap,maskImage,yr,'mmu', MMU); 
            fprintf('    Disturbance  \n');
            DECODE_Phase3_3_ObjectBasedChangeMap(pathMap,maskImage,yr,'mmu', MMU);
        end
        
    end
    fprintf('------------ COMPLETE ---------------\n');
end