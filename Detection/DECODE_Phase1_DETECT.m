function DECODE_Phase1_DETECT(folderpath_decoder,folderpath_tsf, wetlandMask,varargin)
    %DECODE This function is to detect breaks using the DECODE algorithem based on
    %the inputs created by funtion <stackLandsatARDTide2Line>.
    % Create the record of the rows that have not been processed yet before parallel.
    %
    % INPUT:
    %
    %   folderpath_DECODE:        Locate to DECODE working folder, in which the
    %                           metadata (metadata.mat) and the stacking folder
    %                           <StackData> are necessery, and if single path
    %                           Landsat will used, the single path layer
    %                           (singlepath_landsat.tif) is also requried. All
    %                           those can be generated using the function
    %                           <stackLandsatARDTide2Line>
    %   wetlandMask (optional): This is the mask layer, if provided, only these pixels would be processed.
    %   orbitpath (optional):   Orbit path of Landsat ('single' or 'all').
    %                           'single' means change will be identified based
    %                           on the single Landsat path layer with geotiff
    %                           format. 'all' means not to do that. (default
    %                           value: single)
    %
    %   onceread (optional):    To load all the time series data once time. It
    %                           is recommended to set as true, if enough
    %                           computer memory available. This will benefit
    %                           the efficiency of I/O. For instance, to process
    %                           a stacked row dataset with 10 rows by 5000
    %                           columns by 8 bands by 2000 images, with uint16
    %                           format, more 1.5 G memory (1,600,000,000 bytes)
    %                           is required. Usually, we do not to set this as
    %                           true in UCONN HPC because of the memory
    %                           limitation. (default value is false)
    %
    %   delstack (optional):    To delete the stack row folder once change
    %                           detection done. (default value is false)
    %
    %   cprob (optional):       Change probability threshold (default value is 0.99)
    %
    %   conse (optional):       Number of consecutive observations (default value is 6)
    %
    %   maxc (optional):        Maximum number of coefficients used (default value is 6)
    %
    %   task (optional):        Task ID of parallel computation
    %
    %   ntasks (optional):      Total number of tasks of parallel computation
    %
    %   msg (optional)          [false/true] Display processing status (default
    %                           value: false)
    %
    %
    %
    % RETURN:
    %
    %   null
    %
    % REFERENCE(S):
    %   Yang, Xiucheng et al. "".
    
    %   Zhu, Zhe, et al. "Continuous monitoring of land disturbance based on
    %   Landsat time series." Remote Sensing of Environment 238 (2020): 111116.
    %
    %
    % AUTHOR(s): Zhe Zhu, Shi Qiu & Xiucheng Yang
    % DATE: June. 6, 2023
    % COPYRIGHT @ GERSLab
    
    
    %% Have user's inputs
    % requried
    
    if ~exist('folderpath_decoder', 'var')
        folderpath_decoder = pwd;
    end
    
    % optional
    p = inputParser;
    addParameter(p,'cprob', 0.975); % probability for detecting surface change
    addParameter(p,'conse', 5); % number of consecutive observation
    addParameter(p,'maxc', 9); % number of maximum coefficients
    addParameter(p,'onceread', false); % read the landsat data line by line, if enough memory, please set as true (optional)
    addParameter(p,'orbitpath', 'single'); % based on non-overlap landsat data
    addParameter(p,'task', 1); % 1st task
    addParameter(p,'ntasks', 1); % single task to compute
    addParameter(p,'delstack', false); % delete stack data or not
    addParameter(p,'msg', false); % not to display info
    addParameter(p,'mask',true); % Use mask
    addParameter(p,'beginYear',1984); % Begin year default is 1984
    addParameter(p,'modelRefinement',false); % Whether to conduct the model refinement; designed for further recovery analysis
    % request user's input
    parse(p,varargin{:});
    
    task = p.Results.task;
    ntasks = p.Results.ntasks;
    msg = p.Results.msg;
    T_cg = p.Results.cprob;
    conse = p.Results.conse;
    max_c = p.Results.maxc;
    beginYear = p.Results.beginYear;
    modelRefinement = p.Results.modelRefinement;
    % char will be better understood for users
    switch lower(p.Results.orbitpath)
        case 'single'
            singlepath = true;
            landsatpath = 'Single';
        case 'all'
            singlepath = false;
            landsatpath = 'All';
    end
    
    %% Update the TSFit using single or all observations
    if singlepath
        fprintf('\n Using Sinple path observations \n %s \n',folderpath_tsf);
    else
        fprintf('\n Using all path observations \n %s \n',folderpath_tsf);
    end
    
    %% add the matlab search path of GLMnet
    % addpath(fullfile(fileparts(mfilename('fullpath')), 'GLMnet'));
    
    %% Constants:
    % Bands for detection change
    B_detect = 2:6;
    % Treshold of noisercg_new
    Tmax_cg = 1-1e-5;
    
    %% set paths and folders
    folderpath_stack = fullfile(folderpath_decoder, 'StackData');
    
    %% Parallel tasks on the row datasats
    % if strcmp(stackMode,'Pixel')
    %     stackblocks = dir(fullfile(folderpath_stack, 'P*'));
    % elseif strcmp(stackMode,'Line')
    %     stackblocks = dir(fullfile(folderpath_stack, 'R*'));
    % else
    %     error('Select the correct stack mode!\n');
    % end
    stackblocks = dir(fullfile(folderpath_stack, 'P*'));
    num_stacks = length(stackblocks);
    tasks_per = ceil(num_stacks/ntasks);
    start_i = (task-1)*tasks_per + 1;
    end_i = min(task*tasks_per, num_stacks);
    tic
    %% Locate to a certain task, one task for one row folder
%     for i_task = start_i:end_i
    for i_task = start_i:end_i
        %% according to the name of stacking row dataset, the rows # at start and
        % end can be known well.
        foldername_stackblock = stackblocks(i_task).name;
        folderpath_stackrows = fullfile(folderpath_stack, foldername_stackblock);
    
        %% load metadata.mat for having the basic info of the dataset that is in proccess
        metadata = load(fullfile(folderpath_stackrows, sprintf('metadata.mat')));
        metadata = metadata.metadata;
        %% Check the completeness of the stacking
        n_images = dir(fullfile(folderpath_stackrows,'L*'));
        if length(n_images)~=metadata.nimages
            fprintf('\n----> %s \n ....................! Stacking is not completed !........\n',folderpath_stackrows);            
            return
        else
            fprintf('In total %d images to run change detection\n',length(n_images));
        end
    
        filepath_rcg = fullfile(folderpath_tsf, sprintf('rec_cg_%s.mat', foldername_stackblock));
        %% check exsiting record files, and remain the id of the rows that have no results
        if isfile(filepath_rcg)
            if msg
                fprintf('\nAlready exist -> %s\n', foldername_stackblock);
            end
            continue;
        end
        if msg
            fprintf('\nProcessing block #%s at task# %d/%d (%d mins)\n',foldername_stackblock, i_task, end_i,round(toc/60));
        end
    
        %% report log of CCD only for the first first task
        if task == 1 && i_task == 1
            reportLog(folderpath_decoder, ntasks, folderpath_decoder, metadata.nimages, landsatpath, T_cg, conse, max_c);
        end
    
    
        %% read single path data
        obsInfo = readStackBlockData(folderpath_stackrows,metadata.ncols,beginYear);
        proc_cols = obsInfo.pos;
        sdate = obsInfo.sdate;
        line_t = table2array(obsInfo(:,{'b1','b2','b3','b4','b5','b7','b6','tide','qa'}));
        path_r = [];
        path_t = [];
        if singlepath
            singlepathlayer = imread(fullfile(folderpath_decoder, 'singlepath_landsat.tif'));
            path_t = obsInfo.path;
            % Generate the reference path for each pixel
            imageLocs = sub2ind([metadata.nrows,metadata.ncols],obsInfo.row,obsInfo.col);
    
            path_r = singlepathlayer(imageLocs);
        end
        if msg
            fprintf('  Loaded (%d)...',round(toc/60));
        end
    
        %% based on single path Landsat data, ...
        rec_cg = TrendSeasonalFit_DECODE_Block(sdate, line_t, path_t, path_r, ...
            proc_cols, ... % process each pixel vis columns
            T_cg, Tmax_cg, conse, max_c, metadata.nbands, B_detect);
        if msg
            fprintf('  DECODE (%d)...\n',round(toc/60));
        end
        %% Post-Processing of the segments and append more info
        if modelRefinement
            fprintf('       ... Model Refine for recovery analysis (Dont need for cover change)...\n');
            if ~isempty([rec_cg.t_start])
                %% Refine the model fits by residual judgment
                rec_cg = refineModelFit(rec_cg,conse);
            end
        end
        %% save record of time series segments
        save([filepath_rcg, '.part'] ,'rec_cg'); % save as .part
        clear rec_cg;
        movefile([filepath_rcg, '.part'], filepath_rcg);  % and then rename it as normal format
        close all;
    
        if msg
            fprintf('ProcesingTimeSingleRow = %0.2f mins for %s with %d images\r\n', toc/60, foldername_stackblock, metadata.nimages);
        end
    
    end % end of all tasks

end % end of function