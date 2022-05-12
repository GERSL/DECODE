function batchPostProcessCoverImage(task,ntasks)
if ~exist('task', 'var')
    task = 1;
end
if ~exist('ntasks', 'var')
    ntasks = 36;
end


folderpath_mfile = fileparts(mfilename('fullpath'));
addpath(folderpath_mfile);
addpath(fullfile(folderpath_mfile, 'Export'));
addpath(fullfile(folderpath_mfile, 'Characterization'));
addpath(fullfile(folderpath_mfile, 'Detection'));

ARDTiles = globalsets.getARDTiles('decode');
% ARDTiles = {'h027v019','h026v019','h027v018'};
ARDTiles = ARDTiles(11:end,:);

% Locate to the COLD results
folderpath_parentwork = globalsets.PathDECODE;
%% Indicate the regions
region = globalsets.processedRegion;


%% Path to mask image as a reference source
pathMask = fullfile(globalsets.PathMask,region);

% HPC process
years = 1986:2020;

tasks = [];
itask = 1;
    
for year = years
    for i_ARD = 1:length(ARDTiles)
        tasks(itask).year = year;
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
    year = tasknow.year;
    hv_name = tasknow.hv_name;
    fprintf('Process %d - %s\n',year,hv_name);
    
    pathMap = fullfile(folderpath_parentwork, region, hv_name,globalsets.FolderMap);
    maskImage = fullfile(pathMask, [hv_name, '_mask.tif']);
    
    seiveSmallObjects(pathMap,maskImage,year);
end
fprintf('--------------------------------------\n');
end