function DECODE(folderpath_cold,wetlandMask, varargin)
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
% DATE: June. 6, 2021
% COPYRIGHT @ GERSLab


%% Have user's inputs
% requried

if ~exist('folderpath_decode', 'var')
    folderpath_decode = pwd;
end

% optional
p = inputParser;
addParameter(p,'cprob', 0.99); % probability for detecting surface change
addParameter(p,'conse', 6); % number of consecutive observation
addParameter(p,'maxc', 8); % number of maximum coefficients
addParameter(p,'onceread', false); % read the landsat data line by line, if enough memory, please set as true (optional)
addParameter(p,'orbitpath', 'single'); % based on non-overlap landsat data
addParameter(p,'task', 1); % 1st task
addParameter(p,'ntasks', 1); % single task to compute
addParameter(p,'delstack', false); % delete stack data or not
addParameter(p,'msg', false); % not to display info
addParameter(p,'mask',true); % Use mask 
addParameter(p,'ephemeralDetect',false); % Detection ephemeral changes or not
% request user's input
parse(p,varargin{:});

onceread = p.Results.onceread;
task = p.Results.task;
ntasks = p.Results.ntasks;
msg = p.Results.msg;
T_cg = p.Results.cprob;
conse = p.Results.conse;
max_c = p.Results.maxc;
delstack =  p.Results.delstack;
mask = p.Results.mask;
ephemeralDetect = p.Results.ephemeralDetect;
% char will be better understood for users
switch lower(p.Results.orbitpath)
    case 'single'
        singlepath = true;
        landsatpath = 'Single';
    case 'all'
        singlepath = false;
        landsatpath = 'All';
end
    
%% add the matlab search path of GLMnet
% addpath(fullfile(fileparts(mfilename('fullpath')), 'GLMnet'));

%% Constants:
% Bands for detection change
B_detect = 2:6;
% Treshold of noisercg_new
Tmax_cg = 1-1e-5;

%% set paths and folders
% folderpath_cold = '/lustre/scratch/qiu25856/TestGERSToolbox/h029v005/';
folderpath_stack = fullfile(folderpath_cold, 'StackData');
folderpath_tsf = fullfile(folderpath_cold, globalsets.FolderDetection);
% make TSFitLine folder for storing coefficients of Time Series Fitting
if ~isfolder(folderpath_tsf)
    mkdir(folderpath_tsf);
end

%% Parallel tasks on the row datasats
stackrows = dir(fullfile(folderpath_stack, 'R*'));
num_stacks = length(stackrows);
tasks_per = ceil(num_stacks/ntasks);
start_i = (task-1)*tasks_per + 1;
end_i = min(task*tasks_per, num_stacks);

%% Locate to a certain task, one task for one row folder
for i_task = start_i:end_i
    %% according to the name of stacking row dataset, the rows # at start and
    % end can be known well.
    foldername_stackrows = stackrows(i_task).name;
    % name format: R xxxxx xxxxx
    row_start = str2num(foldername_stackrows(2:6));
    row_end = str2num(foldername_stackrows(7:11));
    rows = row_start: row_end;
    folderpath_stackrows = fullfile(folderpath_stack, foldername_stackrows);
    
    %% check exsiting record files, and remain the id of the rows that have no results
    ids_more = [];
    for ir = 1: length(rows)
        filepath_rcg = fullfile(folderpath_tsf, sprintf('record_change_r%05d.mat', rows(ir))); % r: row
        if isfile(filepath_rcg) 
            ids_more = [ids_more; ir];
            if msg
                fprintf('\nAlready exist change results for row #%d\n', rows(ir));
            end
            continue;
        end
    end
    rows(ids_more) = [];
    
    % if all the rows at current task are well done, just skip to next
    % process
    if isempty(rows)
        continue;
    end
    
    %% load metadata.mat for having the basic info of the dataset that is in proccess
    load(fullfile(folderpath_stackrows, 'metadata.mat'));
    
    %% report log of CCD only for the first first task
    if task == 1 && i_task == 1
        reportLog(folderpath_cold, ntasks, folderpath_cold, metadata.nimages, landsatpath, T_cg, conse, max_c);
    end
    
    %% read single path data
    if singlepath
        singlepathlayer = imread(fullfile(folderpath_stackrows, 'singlepath_landsat.tif'));
        singlepathlayer = singlepathlayer(rows,:); % only sub rows for lowing memery
    end
    
    %% read all the time series data ahead of time if enough memory
    if onceread
        %     For instance, to process a stacked row dataset with 10 rows by 5000
        %     columns by 8 bands by 2000 images, with uint16 format, more 1.5 G
        %     memory (1,600,000,000 bytes) is required
        %     1) 10 rows by 5000 pixels by 1300 images will need 1 G mememory
        %     2) 10 rows by 5000 pixels by 2600 images will need 2 G mememory for
        %     the places that have very dense landsat data

        % ~8 secs for loading 2000 images; if single row, ~2 secs
        [sdate, line_t_all, path_t] = readStackLineData(folderpath_stackrows, metadata.ncols, metadata.nbands, rows, singlepathlayer);
    end

    %% for each row, CCD
    for ir = 1: length(rows)
        tic % start to count computing time
        if msg
            fprintf('\nProcessing row #%d at task# %d/%d\n', rows(ir), task, ntasks);
        end
        %% based on single path Landsat data, ...
        if singlepath
            % singlepathlayer(ir,:) means the single path value for the
            % pixel in the current row

            path_r = singlepathlayer(ir,:);
            if onceread
                line_t = line_t_all(:,:,1);
                line_t_all(:,:,1) = []; % when loading successfully, empty the row data for saving memeory
            else
                [sdate, line_t, path_t] = readStackLineData(folderpath_stackrows, metadata.ncols, metadata.nbands, rows(ir), path_r);
            end
            % ccd processing
            line_t = double(line_t); % Modify the type from int16 to float
            % Processed cols for each row (default: the whole rows)
            if mask
                sampleMaskImage = geotiffread(wetlandMask);
                proc_cols = find(sampleMaskImage(rows(ir),:)>0);
            else
                proc_cols = 1:metadata.ncols;
            end
            rec_cg = TrendSeasonalFit_DECODE_Line(sdate, line_t, path_t, path_r, ...
            metadata.ncols, rows(ir),proc_cols, ... % process each pixel vis columns
            T_cg, Tmax_cg, conse, max_c, metadata.nbands, B_detect,ephemeralDetect);
        else
            %% based on all path Landsat data, ...
            if onceread
                line_t = line_t_all(:,:,1);
                line_t_all(:,:,1) = []; % when loading successfully, empty the row data for saving memeory
            else
                [sdate, line_t] = readStackLineData(folderpath_stackrows, metadata.ncols, metadata.nbands, rows(ir), []);
            end
            % ccd processing
            line_t = double(line_t); % Modify the type from int16 to float
            % Processed cols for each row (default: the whole rows)
            if mask
                sampleMaskImage = geotiffread(wetlandMask);
                proc_cols = find(sampleMaskImage(rows(ir),:)>0);
            else
                proc_cols = 1:metadata.ncols;
            end
            
            rec_cg = TrendSeasonalFit_DECODE_Line(sdate, line_t, [], [], ...
                metadata.ncols, rows(ir),proc_cols, ... % process each pixel vis columns
                T_cg, Tmax_cg, conse, max_c, metadata.nbands, B_detect,ephemeralDetect);
        end
        
        %% save record of time series segments
        filepath_rcg = fullfile(folderpath_tsf, sprintf('record_change_r%05d.mat', rows(ir))); % r: row
        save([filepath_rcg, '.part'] ,'rec_cg'); % save as .part
        clear rec_cg;
        movefile([filepath_rcg, '.part'], filepath_rcg);  % and then rename it as normal format
        close all;
        
        if msg
            fprintf('ProcesingTimeSingleRow = %0.2f mins for row #%d with %d images\r\n', toc/60, rows(ir), metadata.nimages); 
        end
    end
    if delstack
        rmdir(folderpath_stackrows, 's');
        if msg
            fprintf('Finished deleting the stack row dataset %s\r\n', foldername_stackrows); 
        end
    end
end % end of all tasks

end % end of function