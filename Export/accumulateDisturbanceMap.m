function accumulateDisturbanceMap(folderpath_map, wetlandMask, years,windowSize)
% ACCUMULATECHANGEMAP This is to create the accumulated change map based on
% the yearly change maps in the folder <ChangeMap>
%
% INPUT:
%
%   folderpath_DECODE:        Locate to DECODE working folder, in which the
%                           change folder <TSFitLine> is necessery, and
%                           this folder was created by <DECODE.m>.
%   wetlandMask (optional): If use the mask regions      
%   years:                  [Array] The years of change map. 
%
%   ctype (optional):       {'all', 'cover', 'condition'} which change type?
%
%   minsize (optional):     Mininum size of change object in each year.(default value is 1). 
%
%   msg (optional):         [false/true] Display processing status (default
%                           value: false)
% 
%
% OUTPUT:
%
% accuchangemap_ctype_yyyy_yyyy.tif in folder <outMap> (ctype means the change type and yyyy means the year)
% uint16
% pixel value: xdoy (x indicates the type of change, and doy indicates DOY) 
% x's range is between 1 to 3
% 1 => cover change
% 2 => condition change
% 3 => ephemeral condition change
%
% i.e., 1002 means the regrowth break occured in the 2nd day
%
%
% REFERENCE(S):
%
%   Zhu, Zhe, et al. "Continuous monitoring of land disturbance based on
%   Landsat time series." Remote Sensing of Environment 238 (2020): 111116.
%
% 
% AUTHOR(s): Zhe Zhu, Shi Qiu & Xiucheng yang
% DATE: July. 7, 2021
% COPYRIGHT @ GERSLab
    tic
    [~, foldername_working] = fileparts(fileparts(folderpath_map));
    fprintf('Start to accumulate change types with between %d and %d for %s\r\n', min(years), max(years), foldername_working);
    
    %% typedoy data first
    distmap_name = 'object_disturbance';
    yearlymaps = dir(fullfile(folderpath_map, [distmap_name,'*.tif'])); % this has change types
  
    %% Check files are ready
    if isempty(yearlymaps)
        fprintf('No yearly change maps at %s!\r\n', folderpath_map);6
        return;
    else
        years_nomap = [];
        for imap = 1: length(years)
            yr = years(imap);
            yearlymap_filepath = fullfile(folderpath_map, [distmap_name, '_', num2str(yr),'.tif']);
            if ~isfile(yearlymap_filepath)
                years_nomap = [years_nomap; yr];
            end
        end
        if ~isempty(years_nomap)
            fprintf('No yearly change maps for Years %d!\r', years_nomap);
            return;
        end
    end

    %% Calculate the accumulated map
    for imap = 1: length(years)
        yr = years(imap);
        
        fprintf('Accumulating the change map for year %d\r', yr);
        % Load the annual map
        yearlymap_filepath = fullfile(folderpath_map, [distmap_name, '_', num2str(yr),'.tif']);
        yearlymap = double(imread(yearlymap_filepath)); % Load the change map

        % Change pixels
        idxChangeLocs = find(yearlymap<9 & yearlymap>0);
        
        % Update the accumalated map
        if imap == 1
            lastChangeYear = zeros(size(yearlymap));
            numOfChange = zeros(size(yearlymap));
        end
        lastChangeYear(idxChangeLocs) = yearlymap(idxChangeLocs)*yr; % Updated with newly year 
        numOfChange(idxChangeLocs) = numOfChange(idxChangeLocs) + yearlymap(idxChangeLocs);
    end

    fprintf('Mask is used: %s \n',wetlandMask);
    mask = GRIDobj(wetlandMask);
    
    outputMap = mask;
    lastChangeYear(mask.Z==0) = 9999; % Non-wetland background
    
    outputMap.Z = uint16(lastChangeYear);    
    GRIDobj2geotiff(outputMap, fullfile(folderpath_map,sprintf('%s_%d_%d.tif', 'lastDisturbanceYear', min(years), max(years))));
    
    numOfChange(mask.Z==0) = 255; % Non-wetland background
    outputMap.Z = uint8(numOfChange);
    GRIDobj2geotiff(outputMap, fullfile(folderpath_map,sprintf('%s_%d_%d.tif', 'numberOfDisturbance', min(years), max(years))));
    
    %% Export the resized map
    % Load the original image and conduct the resample
    [Z_DECODE,R] = readgeoraster(wetlandMask); 
    [Z_resize,R_resize] = mapresize(Z_DECODE,R,1/windowSize,'nearest'); % resize to the calculation scale
    info = geotiffinfo(wetlandMask); % Load geotiff info for write

    % Area within the GRID and percentage within the GRID
    lastChangeYearResize = zeros(size(Z_resize));
    numOfChangeResize = zeros(size(Z_resize));
 
    for i_row = 1:size(Z_resize,1)
        for i_col = 1:size(Z_resize,2)
            % Aggregate the last change year to GRID
            lastChangeWindow = lastChangeYear((i_row-1)*windowSize+1:i_row*windowSize,(i_col-1)*windowSize+1:i_col*windowSize);
            
            numChangeWindow = numOfChange((i_row-1)*windowSize+1:i_row*windowSize,(i_col-1)*windowSize+1:i_col*windowSize);
                        
            % Mask regions
            maskWindow = Z_DECODE((i_row-1)*windowSize+1:i_row*windowSize,(i_col-1)*windowSize+1:i_col*windowSize);
            idxMask = find(maskWindow==1);

            if length(unique(lastChangeWindow))==1
                lastChangeYearResize(i_row,i_col) = unique(lastChangeWindow); 
            else
                aggregatedValue = mode(lastChangeWindow(idxMask),'all');
                lastChangeYearResize(i_row,i_col) = aggregatedValue;
            end
            if length(unique(numChangeWindow))==1
                numOfChangeResize(i_row,i_col) = unique(numChangeWindow); 
            else
                aggregatedValue = mode(numChangeWindow(idxMask),'all');
                numOfChangeResize(i_row,i_col) = aggregatedValue;
            end          
        end
    end
    geotiffwrite(fullfile(folderpath_map,sprintf('%s_%d_%d_%d.tif', 'lastDisturbanceYear', min(years), max(years),windowSize)),lastChangeYearResize,R_resize,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(fullfile(folderpath_map,sprintf('%s_%d_%d_%d.tif', 'numberOfDisturbance', min(years), max(years),windowSize)),uint8(numOfChangeResize),R_resize,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
      
    fprintf('Finished accumluting change maps for %s with %0.2f mins\r\n', foldername_working, toc/60); 
end

