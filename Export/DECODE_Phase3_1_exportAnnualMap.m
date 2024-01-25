function DECODE_Phase3_1_exportAnnualMap(folderpath_Decode,folderTSFit,maskImage, years, varargin)
% EXPORTCHANGEMAP This is to export the change map based on the DECODE
% results.
%
%   
% DATE: May. 7, 2023
% COPYRIGHT @ GERSLab
%
    %% cover types
    classes = num2cell([1:10 6]);
    [Tidalmarsh,Mangrove,Dieback,Tidalflats,Openwater,Woody,Herbaceous,Cultivated,Developedbarren,Landwaterinterface,Others] = deal(classes{:});
    
    %% Begin to run

    % optional
    p = inputParser;  
    addParameter(p,'disturbance',true); % Export disturbance map
    addParameter(p,'cover',true); % Export annual     
    addParameter(p,'minsize', 4); % all size change object
    parse(p,varargin{:});
    isDisturbance =  p.Results.disturbance;
    isCover = p.Results.cover;
    minsize = p.Results.minsize;
    [~, foldername_working] = fileparts(folderpath_Decode);

    fprintf('Start to export maps for %s\r\n', foldername_working);

    folderpath_tsf = fullfile(folderpath_Decode, folderTSFit);   
    folderpath_map = fullfile(folderpath_Decode, [globalsets.FolderMap]);
   
    if ~isfolder(folderpath_map)
        mkdir(folderpath_map);
    end

    tic

    %% get metadata
    metadata = load(fullfile(folderpath_Decode, sprintf('metadata.mat')));
    metadata = metadata.metadata;
    nrows = metadata.nrows;
    ncols = metadata.ncols;
    nblocks = metadata.nblocks;

    nbands = metadata.nbands; % 7 Landsat bands + 1 QA band + 1 wl band 
    nIds = 5; % NDVI, MNDWI, TC Transformation
    % dimension and projection of the image
    jiDim = [ncols,nrows];
    % max number of maps
    max_n = length(years);

    % Produce break map
    disturbanceMap = 255*ones(nrows,ncols,max_n,'uint8'); % press/pulse disturbance
      
    % produce cover map
    coverMap = 255*ones(nrows,ncols,max_n,'uint8'); % annual type map

    % cd to the folder for storing recored structure
    % cd(v_input.name_rst);
    records = dir(fullfile(folderpath_tsf,'rec*.mat')); % folder names
    num_blocks = size(records,1);
    if num_blocks~=nblocks
        fprintf('Some rows are missing for %s\r',folderpath_tsf);
        return
    else        
        fprintf('%d blocks to load for %s\r',nblocks,folderpath_tsf);
    end
    for block = 1: num_blocks
        % show processing status
        if block/num_blocks < 1
            fprintf('Processing %.2f percent for block %d/%d (%d mins)\r',100*(block/num_blocks),block,num_blocks,round(toc/60));
        else
            fprintf('Processing %.2f percent for block %d%d (%d mins)\n',100*(block/num_blocks),block,num_blocks,round(toc/60));
        end

        % load one line of time series models
        load(fullfile(folderpath_tsf,records(block).name)); %#ok<LOAD>

        % postions
        pos = [rec_cg.pos];

        % continue if there is no model available
        l_pos = length(pos);
        if l_pos == 0
            continue
        end

        % break time
        t_break = [rec_cg.t_break];
        % change probability
        change_prob = [rec_cg.change_prob];
        changeType = [rec_cg.changeType];       
        
        for i = 1:l_pos % -1: segment of time series will not have change record!
            % get row and col
            [I,J] = ind2sub(jiDim,pos(i));

            %% initialize pixels have at least one model
            if sum(coverMap(J,I,:) == 255) == max_n
                % Mask regions               
                disturbanceMap(J,I,:) = 9;
                coverMap(J,I,:) = 15; % Define as 9 for the background (others)               
            end
            
            %%  Produce the press disturbance related map
            if length(rec_cg(i).severity)>1 && sum(ismember([Tidalmarsh,Mangrove,Dieback,Tidalflats],[rec_cg(i).type]))
                %% Transitional cover change of tidal wetlands                 
                if any(rec_cg(i).transitionYear)
                    % Leading to cover transitional change
                    tranYrs = rec_cg(i).transitionYear;
                    for i_tr = 1:length(tranYrs)
                        tranYr = tranYrs(i_tr);
                        if rec_cg(i).type(i_tr)<Openwater || rec_cg(i).type(i_tr+1)<Openwater
                            % Only record the tidal wetlands-related
                            % gradual trasition
                            [~,locs] = ismember(tranYr,years);
                            if locs % to confirm the transitional year within the time span                                
                                disturbanceMap(J,I,locs) = 3; % cover press Disturbance
                            end
                        end
                    end                    
                else
                    % Not yet lead to cover tranditional
                    [intersectYr,locs] = ismember(rec_cg(i).years,years);                    
                    disturbanceMap(J,I,locs(intersectYr)) = 4; % condition press Disturbance                    
                end
            end
            
            %% Produce the change type and disturbance type maps for pulse
            % disturbance
            if change_prob(i) == 1 % Break point --> abrupt disturbance   
                break_year = datevecmx(t_break(i));
                break_year = break_year(1);
                
                if changeType(i)<4 && changeType(i)>0 
                    if changeType(i) < 3
                        disturbanceMap(J,I,years == break_year) = 1;  % --> cover pulse disturbance
                    elseif changeType(i) == 3
                        disturbanceMap(J,I,years == break_year) = 2;  % --> conditional pulse disturbance
                    end                    
                elseif changeType(i) == 8
                    disturbanceMap(J,I,years == break_year) = 5;   % Non-tidal wetlands abrupt disturbance   
                end
           
            end
            
            %% Produce the cover maps
            % define time span for each time series segment inside the input years
            yr_end = max(year(rec_cg(i).t_end),rec_cg(i).years(end));
            yr_end = min(yr_end,years(end));
            if yr_end(1) < years(1)
                continue
            end
            yr_start = min(year(rec_cg(i).t_start),rec_cg(i).years(1));
            yr_start = max(yr_start,years(1));
            if yr_start(1) > years(end)
                continue
            end
            y_span = find(yr_start==years):find(yr_end==years);

            if length(rec_cg(i).type)==1
                coverMap(J,I,y_span) = rec_cg(i).type; % stable
            else
                tranYr = rec_cg(i).transitionYear;
                startIdx = find(yr_start==years);
                for i_yr = 1:length(tranYr)
                    endIdx = find(tranYr(i_yr)==years)-1;
                    coverMap(J,I,startIdx:endIdx) = rec_cg(i).type(i_yr); % stable
                    startIdx = find(tranYr(i_yr)==years);
                end
                endIdx = find(yr_end==years);
                coverMap(J,I,startIdx:endIdx) = rec_cg(i).type(end); % stable
            end
        end
    end

    geotif_obj = metadata.GRIDobj;    

    for i_yr = 1: length(years)
        yr = years(i_yr);
        if isDisturbance
            geotif_obj.Z = uint8(disturbanceMap(:,:,i_yr));
            GRIDobj2geotiff(geotif_obj, fullfile(folderpath_map, ['disturbance_',num2str(yr),'.tif']));
        end        
        if isCover
            geotif_obj.Z = uint8(coverMap(:,:,i_yr));
            GRIDobj2geotiff(geotif_obj, fullfile(folderpath_map, ['cover_',num2str(yr),'.tif']));
        end       
    end
    fprintf('Finished exporting annual map for %s with %0.2f mins\r\n', foldername_working, toc/60); 
end