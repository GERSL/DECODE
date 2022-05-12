function batchClipWLData2ARDExtentLandsatDate(task, tasks)
%% This is an example of clip images to reference of Landsat ARD tiles
%
% 
% AUTHOR(s): Shi Qiu & Xiucheng Yang
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
    %% Locate to the Landsat ARD dataset
    path_ard_tars = globalsets.PathLandsatARD;
    %% Indicate the regions
    region = globalsets.processedRegion;
    
    %% Locate to the working folder, in which all outputing data can be found
    path_working = fullfile(globalsets.PathDECODE,region);
    
     %% Location of intepolated tide info
    pathTide = fullfile(globalsets.PathTide,region);
    
    %% Path to mask image as a reference source
    pathMask = fullfile(globalsets.PathMask,region);
    
    %% Tiles
    ARDTiles = globalsets.getARDTiles('stack');
    
%     ARDTiles = {'h026v018','h026v019','h027v018','h027v019'};
%    ARDTiles =  {'h027v009','h027v010','h028v008','h028v009','h028v010','h029v006','h029v007','h029v008','h029v009','h030v004','h030v005','h030v006','h030v007','h031v003','h031v004','h031v006','h032v003'};
   
    fprintf('********* In total, %d of tiles need to process ... *********\n', length(ARDTiles));
    
    %% For each Landsat ARD tile
    for iARD = 1: length(ARDTiles)
        % switch to a certain Landsat ARD tile
        hv_name = ARDTiles{iARD};
        % display
        fprintf('Start to clip tide information for %d - tile - %s\n', iARD, hv_name);
  
        filepath_likeimage = fullfile(pathMask, [hv_name, '_mask.tif']);
        
        % switch to the final destination
        folderpath_ard =  fullfile(path_ard_tars, hv_name);

        % switch to the outputing folder <TideData>
        folderpath_out = fullfile(path_working, hv_name, 'TideData');
        
        if ~isfolder(folderpath_out)
            mkdir(folderpath_out);
        end
        % triger the function to clip the tide information into a ARD-like
        % image
        clipWLData2ARDExtentLandsatDate(folderpath_ard, folderpath_out,pathTide,filepath_likeimage,...
           task,tasks);
    end
end

function clipWLData2ARDExtentLandsatDate(folderpath_ard, folderpath_out,dir_waterlevel, filepath_likeimage, task,ntasks) 
    imfs = dir(fullfile(folderpath_ard,'L*SR.tar'));
    % filter for Landsat folders
    % espa data
    %imf = regexpi({imf.name}, 'L(T05|T04|E07|C08)(\w*)\-(\w*)', 'match'); 
    % Landsat ARD
    imfs = regexpi({imfs.name}, 'L(T05|T04|E07|C08)(\w*)', 'match'); % no expand name
    % imfs = regexpi({imfs.name}, 'LE07(\w*)', 'match'); % no expand name and only for Landsat 7 data

    imfs = [imfs{:}];
    if isempty(imfs)
        warning('No Landsat ARD data at %s!', folderpath_ard);
        return;
    end
    % convert to char array
    imfs = vertcat(imfs{:});
    % no _SR at the end, to locate other filenames such as xxx_SR.tar,
    % xxx_BT.tar, xxx.xml, xxx_SRBx.tif, xxx_BT.tif
    imfs(:, end-2:end) = []; 

    % sort according to yeardoy, and this ensure each core would have same
    % image list if parallel.
    yyyymmdd = str2num(imfs(:, 16:23)); % should change for different sets
    
    % Water level info for cover mapping
    coverYear = (1984:2021)';
    coverDates = [coverYear*10000+0101;coverYear*10000+0701];
    idx = ismember(coverDates,yyyymmdd);
    coverDates(idx) = [];
    yyyymmdd = [yyyymmdd;coverDates];
    [sortedDates, ~] = sort(yyyymmdd);
    sortedDates = sortedDates(sortedDates>19840000,:);
    clear imfs yyyymmdd;
    % number of folders start with "L"
    num_t = size(sortedDates,1);
    fprintf('A total of %04d Landsat images at %s\r\n',num_t, folderpath_ard);

    
    %% Assign stacking tasks to each core
    tasks_per = ceil(num_t/ntasks);
    start_i = (task-1)*tasks_per + 1;
    end_i = min(task*tasks_per, num_t);

    fprintf('At task# %d/%d to process %d images\n', task, ntasks, end_i- start_i + 1);

    %% Start to stack the Landsat ARD from start_i th to end_i th
    tic % record all running time
    for i = start_i:end_i
        %% locate to a certain image
        selectDate = sortedDates(i);   
        selectDOY = num2str(year(datetime(num2str(selectDate),'InputFormat','yyyyMMdd'))*1000+day(datetime(num2str(selectDate),'InputFormat','yyyyMMdd'),'dayofyear'));
        inputImagePath = fullfile(dir_waterlevel, ['wlDoy',selectDOY,'.tif']);

        outputImagePath = fullfile(folderpath_out, ['tide',num2str(selectDate),'.tif']);
        
        if isfile(outputImagePath) 
            % check file size
            filesize = dir(outputImagePath);
            filesize = filesize.bytes;
            if filesize>100000000
                fprintf('(%d/%d) Existing %s\r', i, end_i, outputImagePath);
            else
                delete (outputImagePath);
                clipRasterRIOWarp(inputImagePath, filepath_likeimage, outputImagePath);
                fprintf('(%d/%d) Finished updating %s with total %0.2f mins\n',  i, end_i, outputImagePath,  toc/60);
            end
        else
            clipRasterRIOWarp(inputImagePath, filepath_likeimage, outputImagePath);
            fprintf('(%d/%d) Finished clipping %s with total %0.2f mins\n',  i, end_i, outputImagePath,  toc/60);
        end
    end
end