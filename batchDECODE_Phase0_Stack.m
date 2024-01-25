function batchDECODE_Phase0_Stack(task, tasks)
%% This is an example of stacking data for mutilple tiles
%%% 
    % DATE: Jun. 8, 2023
    % COPYRIGHT @ GERSLab
    restoredefaultpath;
    folderpath_mfile = fileparts(mfilename('fullpath'));
    addpath(folderpath_mfile);
    addpath(fullfile(folderpath_mfile, 'Stack'));
    addpath(fullfile(folderpath_mfile, 'Packages/GRIDobj'));

    if ~exist('task', 'var')
        task = 1;
    end
    if ~exist('tasks', 'var')
        tasks = 1;
    end

    Collection = 1;
    %% Tiles
    ARDTiles = globalsets.getARDTiles('decode');    
    %% Locate to the Landsat ARD dataset
    path_ard_tars = globalsets.PathLandsatARD;
    beginYear = globalsets.Years(1); 

    %% Locate to the working folder, in which all outputing data can be found
    path_working = fullfile(globalsets.PathDECODE);
    
    %% For each Landsat ARD tile
    tic % record all running time   
    for iARD = 1: length(ARDTiles)
        % switch to a certain Landsat ARD tile
        hv_name = ARDTiles{iARD};        
        fprintf('Stack %s (%d/%d) %d mins\n',hv_name,iARD,length(ARDTiles),round(toc/60));
        % switch to the final destination
        folderpath_ard =  fullfile(path_ard_tars, hv_name);
        
        % folder of the intepolated water level imag
            % Mask image es
        folderpath_wl = fullfile(path_working, hv_name,'TideData');
        
        % switch to the outputing folder <StackData>
        folderpath_out = fullfile(path_working, hv_name, 'StackData');
   
        pathMaskImage = fullfile(globalsets.PathMask,sprintf('%s_mask.tif',hv_name));
        %% Check the completeness of the stack
        if isfile(fullfile(fullfile(path_working, hv_name,sprintf('metadata.mat'))))
            load(fullfile(fullfile(path_working, hv_name,sprintf('metadata.mat'))));
            nblocks = metadata.nblocks;
            nimages = metadata.nimages;
            blockFolders = dir(fullfile(folderpath_out,'P*'));
            flag = 1;
            if isempty(blockFolders)
                flag = 0;    
            else
                i_b = nblocks;
                try
                    imgs = dir(fullfile(folderpath_out,blockFolders(i_b).name,'L*mat'));
                    if length(imgs)~=nimages
                        fprintf('   Lack %d images to stack\n',nimages-length(imgs))
                        flag = 0;                           
                    end
                catch 
                    flag = 0;
                end
            end 
            if flag == 1
                fprintf('---- Skip --> %s with %d blocks and %d images \n',hv_name,nblocks,nimages);
                continue
            end
        end
        %% Stack the images
        stackLandsatARD2MaskRegion(folderpath_ard, folderpath_out,folderpath_wl,pathMaskImage,...
                'orbitpath', 'single', 'check', true, 'rigorousCheck', true,'beginYear',beginYear, ...
                'Collection',Collection,...
                'task', task ,'ntasks', tasks);
        
    end
end