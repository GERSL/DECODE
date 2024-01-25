function batchDECODE_Phase3_S2_CoverChangeMap(task,ntasks,interval)
% Generate the cover change maps by comparing the cover maps between
% adjacent years

    if ~exist('task', 'var')
        task = 1;
    end
    if ~exist('ntasks', 'var')
        ntasks = 1;
    end
    if ~exist('interval', 'var')
        interval = 1;
    end
    restoredefaultpath;
    folderpath_mfile = fileparts(mfilename('fullpath'));
    addpath(folderpath_mfile);
    addpath(fullfile(folderpath_mfile, 'Export'));
    addpath(fullfile(folderpath_mfile,'Packages', 'GRIDobj'));
    
    ARDTiles = globalsets.getARDTiles('decode');
    
    % Locate to the COLD results
    folderpath_parentwork = globalsets.PathDECODE;
    
    classIDs = num2cell([1:6 0]);
    [Tidalmarsh,Mangrove,Dieback,Tidalflats,Openwater,Others,Background] = deal(classIDs{:});
    targetClassIDs = 1:4;
    %% Path to mask image as a reference source
    pathMask = fullfile(globalsets.PathMask);
    patchMap = 'cover_Patch'; % map types % Cover maps
    coverMap = 'cover'; % Original Classification results
    % HPC process
    % interval = 5;
    % 
    if interval==1
        years = 1986:interval:2020;
    else
        years = [1986,1990:interval:2020];
    end
    MMU = 5; % 5 Pixels
    conn = 8; % Connection is 4-direction or 8-
    tasks = [];
    itask = 1;
    
    for i_y = 1:length(years)-1
        % No cover change for the first year
        for i_ARD = 1:length(ARDTiles)
            tasks(itask).startYear = years(i_y);
            tasks(itask).endYear = years(i_y+1);
            tasks(itask).hv_name = ARDTiles{i_ARD};
            itask = itask + 1;
        end
    end
    
    totalTasks = length(tasks);
    tasks_per = ceil(totalTasks/ntasks);
    start_i = (task-1)*tasks_per + 1;
    end_i = min(task*tasks_per, totalTasks);
    fprintf('In total %d images needed to be processed\n',end_i-start_i+1);
    for i_task = start_i: end_i
        tasknow = tasks(i_task);
        startYear = tasknow.startYear;
        endYear = tasknow.endYear;
    
        hv_name = tasknow.hv_name;
        fprintf('Land cover change of %d to %d - %s\n',startYear,endYear,hv_name);
        
        pathMap = fullfile(folderpath_parentwork, hv_name,globalsets.FolderMap);
        if isfile(fullfile(pathMap,sprintf('cover_Conversion_%d_%d.tif',startYear,endYear)))
            fprintf('  Skip\n');
            continue
        end
    
        maskImage = fullfile(pathMask, [hv_name, '_mask.tif']);
        maskLayer = GRIDobj(maskImage);
        %% Cover change for the postprocessed maps
        % Read the cover maps
        formerCoverPatchMap = imread(fullfile(pathMap,sprintf('%s_%d.tif',patchMap,startYear)));
        latterCoverPatchMap = imread(fullfile(pathMap,sprintf('%s_%d.tif',patchMap,endYear)));
        
        % Only consider the tidal wetland classes
        formerCoverPatchMap(~ismember(formerCoverPatchMap,targetClassIDs)) = 0;
        latterCoverPatchMap(~ismember(latterCoverPatchMap,targetClassIDs)) = 0;
    
        idxChange = uint8(formerCoverPatchMap~=latterCoverPatchMap);
     
        % Cover change map after post-processing
        coverChangeMap = formerCoverPatchMap.*idxChange*10+latterCoverPatchMap.*idxChange;
    
    
        if MMU>0
        % Remove the isolated noise
            coverValue = coverChangeMap;
            coverValue(coverValue>0)=1;
    
            CC = bwconncomp(coverValue, conn);
            Stat = regionprops(CC, 'area','BoundingBox','SubarrayIdx');
            L = labelmatrix(CC);
            SmallIds = find([Stat.Area] < MMU);
            SmallObj = ismember(L, SmallIds);
            coverChangeMap = coverChangeMap.*(uint8(~SmallObj));
        end
    
        maskLayer.Z = coverChangeMap;
        GRIDobj2geotiff(maskLayer,fullfile(pathMap,sprintf('cover_Conversion_%d_%d.tif',startYear,endYear)));
    end
    fprintf('--------------------------------------\n');
end


