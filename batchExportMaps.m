function batchExportMaps(task,ntasks)
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

ARDTiles = globalsets.getARDTiles('decode');
% ARDTiles = {'h027v019','h026v019','h027v018','h026v018'};
% ARDTiles = ARDTiles(46:end,:);
% Locate to the DECODE results
folderpath_parentwork = globalsets.PathDECODE;
%% Indicate the regions
region = globalsets.processedRegion;


%% Path to mask image as a reference source
pathMask = fullfile(globalsets.PathMask,region);

% HPC process
num_t = length(ARDTiles);
tasks_per = ceil(num_t/ntasks);
start_i = (task-1)*tasks_per + 1;
end_i = min(task*tasks_per, num_t);

for i_task = start_i: end_i
    hv_name = ARDTiles{i_task};
    folderpath_tileDecode = fullfile(folderpath_parentwork, region, hv_name);
    
    maskImage = fullfile(pathMask, [hv_name, '_mask.tif']);
    %% To clear the spectral break file and replace them by cover/change type files
    % temporal folder to reserve the updated rec_cg
    
    folderpath_tsfTemp = fullfile(folderpath_tileDecode, [globalsets.FolderDetection,'Type']);
    if exist(folderpath_tsfTemp,'dir')
    % The first time to replace the original DECODE break detection results
    % with cover and change type
        folderpath_tsf = fullfile(folderpath_tileDecode, globalsets.FolderDetection);
        
        records = dir(fullfile(folderpath_tsfTemp,'record_change*.mat')); % folder names
        num_records = size(records,1);
        nrows = 5000;
        if num_records==nrows
            rmdir(folderpath_tsf, 's');
            movefile(folderpath_tsfTemp,folderpath_tsf);
        else
            fprintf('Lack of cover/change type results%s (existing-%d)\r',folderpath_tsfTemp,num_records);
            records_break = dir(fullfile(folderpath_tsf,'record_change*.mat'));
            fprintf('Num of spectral break results%s (existing-%d)\r',folderpath_tsf,size(records_break,1));
            if num_records==records_break
                rmdir(folderpath_tsf, 's');
                movefile(folderpath_tsfTemp,folderpath_tsf);
                fprintf('Equal numbers less than 5000');
            else
                fprintf('Not equal numbers\r');
                continue
            end
        end
    end
    %% Export four kinds of maps for mutiple years
     exportMaps(folderpath_tileDecode,maskImage, [1986: 2020], 'msg', true);
%     exportChangeMagnitude(folderpath_tileDecode,maskImage, [1986: 2020], 'msg', true);
     accumulateChangeMap(folderpath_tileDecode, maskImage, [1986: 2020], 'msg', true,'ctype', 'all');
end
end