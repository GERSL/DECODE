function batchDECODE_Phase2_Characterize(task,ntasks)
    %% Before characterize the mangrove status, a random forest classifier is needed 
    % Here we have built the random forst models for Eastern and Western US, respectively
    % If the approach would like to be transferred to outside US regions,
    % it is better to re-built the random forest model
    if ~exist('task', 'var')
        task = 1;
    end
    if ~exist('ntasks', 'var')
        ntasks = 1;
    end

    folderpath_mfile = fileparts(mfilename('fullpath'));
    addpath(folderpath_mfile);
    addpath(fullfile(folderpath_mfile, 'Characterization'));
    addpath(fullfile(folderpath_mfile, 'Packages','PackageRandomForest'));

    ARDTiles = globalsets.getARDTiles('decode');
    years = globalsets.Years;
    
    % Add land cover type and disturbance type into the rec_cg
    trainingDataName = 'RandomForest_CoverSample';    

    %% Locate to the working folder, in which all outputing data can be found
    path_working = fullfile(globalsets.PathDECODE);
   
    tic
    for iARD = 1: length(ARDTiles)
        hv_name = ARDTiles{iARD};
        fprintf('Start to characterize change and cover for %s\n', hv_name);
        % Location of the per-tile results
        folderpath_tileDecoder = fullfile(path_working, hv_name);
        
        folderpath_tsf = fullfile(folderpath_tileDecoder, [globalsets.FolderDetection]); %Input Coeff from change detection
        
        % folder to reserve the characterization results
        folderpath_tsf_r = fullfile(folderpath_tileDecoder, [globalsets.FolderDetection,'_R']); % Output Coeff
        if ~isfolder(folderpath_tsf_r)
            mkdir(folderpath_tsf_r);
        end
        
        records = dir(fullfile(folderpath_tsf,'rec_cg*.mat')); % folder names
        num_block = size(records,1);

        %% get metadata
        load(fullfile(path_working, hv_name, sprintf('metadata.mat')));
        
        if num_block~=metadata.nblocks
            fprintf('Some rows are missing for %s\r',folderpath_tsf);
            continue
        end
       
        %% Load RF models based on east/west coasts
        if str2double(hv_name(2:4))<10
            load(sprintf('%s_West_RF',trainingDataName)); % Load modelRF
            fprintf('Loaded the model for west region\n');
        else            
            load(sprintf('%s_East_RF',trainingDataName)); % Load modelRF
            fprintf('Loaded the model for east region\n');
        end
        %% Prior knowledge of the latitude limitation for mangrove distribution
        if str2double(hv_name(end-2:end))<16
            nonMangrove = 1; % cannot have mangrove
        else
            nonMangrove = 0;
        end

        %% Parallel tasks on the row datasats
        tasks_per = ceil(num_block/ntasks);
        start_i = (task-1)*tasks_per + 1;
        end_i = min(task*tasks_per, num_block);
        fprintf('-----Total %d blocks need to be processed in this core between %d and %d ----\n',end_i-start_i+1,start_i,end_i);
        for i_task = start_i:end_i
            block = i_task;        
            blockName = records(i_task).name;
            fprintf('Processing %.2f percent for block %d/%d (%d mins)\n',100*(block/num_block),i_task,end_i,round(toc/60));
            filepath_rcg = fullfile(folderpath_tsf_r, blockName); % r: row    
            if isfile(filepath_rcg)
                fprintf('Having existed block %s\n',filepath_rcg);
                continue
            end
            % load one line of time series models
            load(fullfile(folderpath_tsf,blockName)); %#ok<LOAD>

            rec_cg =DECODE_Phase2_Characterize(rec_cg,modelRF,nonMangrove,years);

            save([filepath_rcg, '.part'] ,'rec_cg'); % save as .part
            clear rec_cg;
            movefile([filepath_rcg, '.part'], filepath_rcg);  % and then rename it as normal format
        end
    end
end