function batchStackLandatARDTide2Line(task, tasks)
%% This is an example of stacking data for mutilple tiles

    folderpath_mfile = fileparts(mfilename('fullpath'));
    addpath(folderpath_mfile);
    addpath(fullfile(folderpath_mfile, 'Stack'));


    if ~exist('task', 'var')
        task = 1;
    end
    if ~exist('tasks', 'var')
        tasks = 1;
    end
    
    %% Tiles
   % ARDTiles =  {'h027v009','h027v010','h028v008','h028v009','h028v010','h029v006','h029v007','h029v008','h029v009','h030v004','h030v005','h030v006','h030v007','h031v003','h031v004','h031v006','h032v003'};
    %% Tiles
    ARDTiles = globalsets.getARDTiles('decode');
%     ARDTiles = ARDTiles(32:44,:);
    ARDTiles = ARDTiles(18:end,:);
%     ARDTiles = {'h026v018','h026v019','h027v018','h027v019'};

    %% Locate to the Landsat ARD dataset
    path_ard_tars = globalsets.PathLandsatARD;
    
    %% Indicate the regions
    region = globalsets.processedRegion;
    
    %% Locate to the working folder, in which all outputing data can be found
    path_working = fullfile(globalsets.PathDECODE,region);
    
    %% For each Landsat ARD tile
    for iARD = 1: length(ARDTiles)
        % switch to a certain Landsat ARD tile
        hv_name = ARDTiles{iARD};
        % display
        fprintf('Start to stack ARD Tile %s\n', hv_name);
        
        % switch to the final destination
        folderpath_ard =  fullfile(path_ard_tars, hv_name);
        
        % folder of the intepolated water level images
        folderpath_wl = fullfile(path_working,hv_name,'TideData');
        
        % switch to the outputing folder <StackData>
        folderpath_out = fullfile(path_working, hv_name, 'StackData');
        
        % triger the function of stacking Landsat ARD into multiple row
        % dataset with BIP format
        stackLandsatARDTide2Line(folderpath_ard, folderpath_out,folderpath_wl,...
            'orbitpath', 'single', 'check', true, ...
            'task', task ,'ntasks', tasks);
    end
end
