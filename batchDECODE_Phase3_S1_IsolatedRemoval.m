function batchDECODE_Phase3_S1_IsolatedRemoval(task,ntasks)
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
    folderMosaic = 'Mosaic';

    folderMap = globalsets.FolderMap;
    mapType = 'cover_MMU';

    years = 1986:2020;
    regions = {'Northwest','West','Southwest','Texas','WestLouisiana','LouisianaDelta','South','NorthFlorida','CentralFlorida','SouthFlorida','Southeast','East','Northeast'};

    dataType = 'uint8';
    %% Path to mask image as a reference source
    pathMask = fullfile(globalsets.PathMask);
    
    %% Conda env and functions
    %% Python modules and functions
    pathCODE = fullfile(fileparts(mfilename('fullpath')),'PyTool');
    
    condaPath = '. $HOME/miniconda3/etc/profile.d/conda.sh; conda activate wetland; python ';
    condaDeactivate = 'conda deactivate';    
    mosaicFunction = [pathCODE,'/mosaicTiles.py'];    
    compressFunction = [pathCODE,'/compressStack.py'];    
    clipFunction = [pathCODE,'/cust_warp.py'];
    % HPC process    
    patchSize = 100; % 5 Pixels
    tasks = [];
    itask = 1;
    
    for yr = years   
        for i_r = 1:length(regions) 
            tasks(itask).yr = yr;
            tasks(itask).region = regions{i_r};
            itask = itask + 1;
        end
    end
    
    totalTasks = length(tasks);
    tasks_per = ceil(totalTasks/ntasks);
    start_i = (task-1)*tasks_per + 1;
    end_i = min(task*tasks_per, totalTasks);
    fprintf('In total %d images needed to be processed\n',end_i-start_i+1);
    tic
    
    count = 1;
    for i_task = start_i: end_i
        tasknow = tasks(i_task);
        yr = tasknow.yr;
        region = tasknow.region;  
        yearlyMapType = sprintf('%s_%d',mapType,yr);
        pathMap = fullfile(fileparts(fileparts(path_working)),folderMosaic,'TempPatchProcess',region);
        if ~exist(pathMap,'dir')
            mkdir(pathMap);
        end
        ARDTilesInRegion = divideZones(ARDTiles,region); 
        maskImage = fullfile(pathMap,sprintf('mask_%s.tif',region));
        coverImage = fullfile(pathMap,sprintf('%s_%s_%d.tif',region,mapType,yr));
        patchImage = replace(coverImage,'MMU','Patch');
        if isfile(patchImage)
            fprintf('Skip %s\n',patchImage);        
            continue
        else
            fprintf('Process %02d -> %s (%d tiles)\n',count,patchImage,length(ARDTilesInRegion));
            count = count +1;
%           continue
        
            %% Create the mask layer for the region if not existing
            if ~isfile(maskImage)
                for i_t = 1:length(ARDTilesInRegion)
                    hv_name = ARDTilesInRegion{i_t};
                    mapFiles = dir(fullfile(pathMask,[hv_name,'_mask.tif']));
                    tempPathMap = fullfile(pathMap,'Mask');
                    if ~exist(tempPathMap,'dir')
                        mkdir(tempPathMap)
                    end
                    for i_f = 1:length(mapFiles)
                        copyfile (fullfile(mapFiles(i_f).folder,mapFiles(i_f).name),tempPathMap);
                    end
                    fprintf('%s...',hv_name);
                end
    
                outputMosaic = replace(maskImage,'mask','mask_temp');
                outputCompress = maskImage;
    
                % By using Conda environment
                commandStack = sprintf("%s%s --input '%s' --output '%s' --strName %s --dataType %s;%s",condaPath,mosaicFunction,...
                    pathMap,outputMosaic,'h',dataType,condaDeactivate);
    
                system(commandStack);
    
                fprintf('  Compress...\n');
                commandCompress = sprintf("%s%s --input '%s' --output '%s;%s'",condaPath,compressFunction,outputMosaic,outputCompress,condaDeactivate);
                system(commandCompress);
            end
    %         continue
            %% Create the mosaic image if not existing
            if ~isfile(coverImage)
                %% remove the maps to a temporal folder 
                for i_t = 1:length(ARDTilesInRegion)
                    hv_name = ARDTilesInRegion{i_t};
                    mapFiles = dir(fullfile(path_working,hv_name,folderMap,[yearlyMapType,'*']));
                    tempPathMap = fullfile(pathMap,hv_name);
                    if ~exist(tempPathMap,'dir')
                        mkdir(tempPathMap)
                    end
                    for i_f = 1:length(mapFiles)
                        copyfile (fullfile(mapFiles(i_f).folder,mapFiles(i_f).name),tempPathMap);
                    end
                    fprintf('%s...',hv_name);
                end
    
                %% Mosaic the cover maps
                fprintf('\nMosaic %s - %s\n',yearlyMapType,region);
    
                outputMosaic = fullfile(pathMap,sprintf('%s_Temp_%s.tif',region,yearlyMapType));
                outputCompress = coverImage;
    
                fprintf('Image (%s) composite -> %s\n',yearlyMapType,pathMap);
                commandStack = sprintf("%s%s --input '%s' --output '%s' --strName %s --dataType %s;%s",condaPath,mosaicFunction,...
                    pathMap,outputMosaic,yearlyMapType,dataType,condaDeactivate);            
                system(commandStack);
    
                fprintf('  Compress...\n');
                commandCompress = sprintf("%s%s --input '%s' --output '%s;%s'",condaPath,compressFunction,outputMosaic,outputCompress,condaDeactivate);
                
                system(commandCompress);
    
                delete (outputMosaic);
                fprintf('Completed -> %s\n',outputCompress);
            end
            fprintf('Remove isolated patches... \n');
            DECODE_Phase3S_removeIsolatedTidalWetlandOverestimation(coverImage,patchImage,maskImage,'patchSize', patchSize);   
        end
        % Clip the region-based map to tiles
        for i_t = 1:length(ARDTilesInRegion)
            hv_name = ARDTilesInRegion{i_t};            
            outputPatchImagePath = fullfile(path_working,hv_name,folderMap,sprintf('cover_Patch_%d.tif',yr));
            inputPatchImagePath = patchImage;
            imgPathLike = fullfile(pathMask,sprintf('%s_mask.tif',hv_name));
            if isfile(outputPatchImagePath)
                 fprintf('Delete %s\n',outputPatchImagePath);
                 delete (outputPatchImagePath);
            end
            commandStack = sprintf("%s%s %s %s --like %s;%s",condaPath,clipFunction,inputPatchImagePath, outputPatchImagePath, imgPathLike,condaDeactivate);
                
            system(commandStack);
            fprintf('%s...',hv_name);
        end
    end
    fprintf('--------------------------------------\n');
end