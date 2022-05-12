function batchCoverChangeMaps(task,ntasks)
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
    ARDTiles = ARDTiles(46:end,:);
    % ARDTiles = {'h026v018','h026v019','h027v018','h027v019'};

    % Locate to the COLD results
    folderpath_parentwork = globalsets.PathDECODE;
    %% Indicate the regions
    region = globalsets.processedRegion;

    %% Path to mask image as a reference source
    pathMask = fullfile(globalsets.PathMask,region);

    for iARD = 1: length(ARDTiles)
        hv_name = ARDTiles{iARD};
        fprintf('Start to characterize change and cover for %s\n', hv_name);
        folderpath_tileDecode = fullfile(folderpath_parentwork, region, hv_name);
        % Locate to the water level folder
        folderpath_Tide = fullfile(folderpath_tileDecode,'TideData'); 
        %% Add land cover type and disturbance type into the rec_cg
        recordCoverChangeToRecCg(task,ntasks,folderpath_tileDecode,folderpath_Tide,[1984:2020],'msg', true);
    end
end