function batchClipDailyTide2ARDExtentLandsatDate(task, tasks)
%% This is an example of clip intepolated daily tide into ARD tiles (only needed when plot the time series)
%
% 
% AUTHOR(s): Xiucheng Yang
% DATE: Jun. 8, 2021
% COPYRIGHT @ GERSLab

    folderpath_mfile = fileparts(mfilename('fullpath'));
    addpath(folderpath_mfile);
    folderpath_mfile = fileparts(folderpath_mfile);
    addpath(folderpath_mfile);
    if ~exist('task', 'var')
        task = 1;
    end
    if ~exist('tasks', 'var')
        tasks = 1;
    end
    %% Indicate the regions
    region = globalsets.processedRegion;
    
    %% Locate to the working folder, in which all outputing data can be found
    path_working = fullfile(globalsets.PathDECODE,region);
    
     %% Location of intepolated tide info
    pathTide = fullfile(globalsets.PathTide,region);
    
    %% Tiles
%     ARDTiles = globalsets.getARDTiles('stack');
    ARDTiles = {'h027v019'};
    fprintf('********* In total, %d of tiles need to process ... *********\n', length(ARDTiles));
    
    %% For each Landsat ARD tile
    for iARD = 1: length(ARDTiles)
        % switch to a certain Landsat ARD tile
        hv_name = ARDTiles{iARD};
        % display
        fprintf('Start to clip tide information for %d - tile - %s\n', iARD, hv_name);

        %% Path to mask image as a reference source
        if isfile(fullfile(path_working, hv_name,[hv_name, '_mask_1km.tif']))
            filepath_likeimage = fullfile(path_working, hv_name,[hv_name, '_mask_1km.tif']);
        else
            %% Path to mask image as a reference source
            pathMask = fullfile(globalsets.PathMask,region);
            ref30m = GRIDobj(fullfile(pathMask, [hv_name, '_mask.tif']));
            ref1km = resample(ref30m,1000);
            GRIDobj2geotiff(ref1km,fullfile(folderpath_decode,[hv_name, '_mask_1km.tif']));
            filepath_likeimage = fullfile(folderpath_decode,[hv_name, '_mask_1km.tif']);
        end
                
        % switch to the outputing folder <DailyTide>
        folderpath_out = fullfile(path_working, hv_name, 'DailyTide');
        
        if ~isfolder(folderpath_out)
            mkdir(folderpath_out);
        end
        % triger the function to clip the tide information into a ARD-like
        % image
        clipDailyTide2ARDExtentLandsatDate(pathTide, folderpath_out,filepath_likeimage,...
           task,tasks);
    end
end

function clipDailyTide2ARDExtentLandsatDate(pathTide, folderpath_out, filepath_likeimage, task,ntasks) 
    imfs = dir(fullfile(pathTide,'wl*tif'));
    imfs = {imfs.name}; 
    if isempty(imfs)
        warning('No daily tide data at %s!', pathTide);
        return;
    end
    % convert to char array
    imfs = vertcat(imfs{:});
  
    num_t = size(imfs,1);
    fprintf('A total of %04d doys at %s\r\n',num_t, pathTide);

    
    %% Assign stacking tasks to each core
    tasks_per = ceil(num_t/ntasks);
    start_i = (task-1)*tasks_per + 1;
    end_i = min(task*tasks_per, num_t);

    fprintf('At task# %d/%d to process %d images\n', task, ntasks, end_i- start_i + 1);

    %% Start to stack the Landsat ARD from start_i th to end_i th
    tic % record all running time
    for i = start_i:end_i
        %% locate to a certain image
        tide = imfs(i,:);
        inputImagePath = fullfile(pathTide, tide);

        outputImagePath = fullfile(folderpath_out, tide);
        if isfile(outputImagePath)                                                                                                                                                                                                                     
            fprintf('(%d/%d) Existing %s\r', i, end_i, outputImagePath);
        else
            clipRasterRIOWarp(inputImagePath, filepath_likeimage, outputImagePath);
            fprintf('(%d/%d) Finished clipping %s with total %0.2f mins\n',  i, end_i, outputImagePath,  toc/60);
        end
    end
end