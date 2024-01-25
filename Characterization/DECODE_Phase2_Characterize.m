function rec_cg = DECODE_Phase2_Characterize(rec_cg,modelRF,nonMangrove,years)
% This function is used to update the rec_cg to (1)remove the change caused
% by water quality and (2) offer the type of the wetland
%
% INPUT:
%
%   folderpath_Decode:        Locate to DECODE working folder, in which the
%                           change folder <TSFitLine> is necessery, and
%                           this folder was created by <DECODE.m>.
%   folderpath_Tide:          Locate to folder of daily water level images
%
%   years:                  [Array] The years of change map. 
%
%
%   msg (optional):         [false/true] Display processing status (default
%                           value: false)
%
% OUTPUT:
%
% updated rec_cg: add one field of land cover types:
% Values: rec_cg.type
% 1 => Open water
% 2 => vegetated wetlands
% 9 => others, tidal flats, soil, uplands
% updated rec_cg: add one field of disturbance types:
% Values: rec_cg.dist
% 0 => no change
% 1 => cover change
% 2 => condition change

% AUTHOR(s): Xiucheng Yang, Shi Qiu and Zhe Zhu 
% DATE: May. 7, 2021
% COPYRIGHT @ GERSLab
    predictYear = 50;
    %% cover types
    classes = num2cell([1:10 6]);
    [Tidalmarsh,Mangrove,Dieback,Tidalflats,Openwater,Woody,Herbaceous,Cultivated,Developedbarren,Landwaterinterface,Others] = deal(classes{:});
        
%     if isfield(rec_cg,'predictedValue')
%         rec_cg = rmfield(rec_cg,{'predictedValue','type','obs','bands','tide'});        
%     end
%     if isfield(rec_cg,'dist')
%         rec_cg = rmfield(rec_cg,{'dist'});        
%     end
    if isfield(rec_cg,'types')
        rec_cg = rmfield(rec_cg,{'types','changeType','transitionYear'});        
    end
    % postions
    pos = [rec_cg.pos];
    l_pos = length(pos);

    % continue if there is no model available for this row
    if l_pos == 0
        rec_cg.years = [];
        rec_cg.types = [];
        rec_cg.type = [];
        rec_cg.changeType = [];
        rec_cg.transitionYear = [];
        rec_cg.severity = [];
        filepath_rcg = fullfile(folderpath_tsfTemp, sprintf('record_change_r%05d.mat', line)); % r: row
        save([filepath_rcg, '.part'] ,'rec_cg'); % save as .part
        clear rec_cg;
        movefile([filepath_rcg, '.part'], filepath_rcg);  % and then rename it as normal format
        return
    end
      
    % Two fields to reserve the land cover type and disturbance type
    types = num2cell(zeros(l_pos,1));
    changeTypes = num2cell(zeros(l_pos,1));  
    [rec_cg(:).years] = deal(types{:});
    [rec_cg(:).types] = deal(types{:});
    [rec_cg(:).type] = deal(types{:});
    [rec_cg(:).changeType] = deal(changeTypes{:});
    [rec_cg(:).transitionYear] = deal(changeTypes{:});
    [rec_cg(:).severity] = deal(types{:});
    % break time
    t_break = [rec_cg.t_break];
    % change probability
    change_prob = [rec_cg.change_prob];
    
    %% Initialize the annual cover types
    for i = 1:l_pos
        %% Define the years for each temporal segments
        i_start = rec_cg(i).t_start;
        yr_start = ymd(datetime(i_start,'ConvertFrom','datenum'));
        i_end = rec_cg(i).t_end;
        yr_end = ymd(datetime(i_end,'ConvertFrom','datenum'));
        yr_end = min(yr_end(1),years(end));

        yr_plot = (yr_start:yr_end)';
        rec_cg(i).years = yr_plot; % During of the temporal segment
        %% Generate the annual class based on RF model
        coef = rec_cg(i).coefs;
        rmse =  rec_cg(i).rmse';

        inputVariables = [];
        for iYr = 1:length(rec_cg(i).years)            
            coefYear = coef;
            coefYear(1,:) = coef(1,:) + coef(2,:)*datenum(rec_cg(i).years(iYr),7,1);
            inputVariables(iYr,:) = [reshape(coefYear,9*12,[])',rmse];
        end
        [predictedClass,~,votes] = classRF_predict(inputVariables,modelRF); % class
        
        if nonMangrove == 1
            % Update the misclassified mangrove into the most possible
            % catogery
            idx = predictedClass==Mangrove|predictedClass==Dieback;
            if sum(idx)
                votes(:,[Mangrove,Dieback]) = 0;
                [~,updatedClass] = max(votes,[],2); 
                predictedClass(idx) = updatedClass(idx); % Noth regions dont have mangrove
            end
        else %Mangrove dieback only occurs when for mangrove after hurricane
            idx = predictedClass==Dieback;
            votes(:,Dieback) = 0;
            [~,updatedClass] = max(votes,[],2); 
            if i>1 && change_prob(i-1) == 1 && rec_cg(i-1).types(end)==Mangrove   
                % In this case, dieback exists
            elseif sum(idx)
                predictedClass(idx) = updatedClass(idx); % dieback cannot exist for non-mangrove region
            end
        end
        rec_cg(i).types = predictedClass;
       
        %% Check consistent of the cover type        
        if rec_cg(i).types(1)==rec_cg(i).types(end)
            rec_cg(i).types(:)=mode(rec_cg(i).types(:)); % Define the same type for the segment with highest confidence
            rec_cg(i).type = rec_cg(i).types(1);                     
        else
            % Count the number of consecutive identical elements
            d = [1;diff(rec_cg(i).types(:))~=0;1];
            np = diff(find(d));
            np = repelem(np,np);
            % remove the isolated cover type
            stableIdx = find(np>1);
            isolatedIdx = find(np==1);
            if ismember(1,isolatedIdx)
                % The beginning year of the time series
                rec_cg(i).types(1) = rec_cg(i).types(2);
                isolatedIdx(isolatedIdx==1)=[];
                stableIdx = [1; stableIdx];
            end
            if ~isempty(isolatedIdx)
                for i_isolated = 1:length(isolatedIdx)
                    replaceId = max(stableIdx(stableIdx<isolatedIdx(i_isolated)));
                    rec_cg(i).types(isolatedIdx(i_isolated)) = rec_cg(i).types(replaceId);
                end
            end
            % Check if the change is uni-direction and exclude the
            % bi-directional changes
            dirIdx = find([1;diff(rec_cg(i).types(:))]~=0);
            if length(rec_cg(i).types(dirIdx)) ~= length(unique(rec_cg(i).types(dirIdx)))
                rec_cg(i).types(:) = mode(rec_cg(i).types(:));
                rec_cg(i).type = rec_cg(i).types(1);
            elseif length(unique(rec_cg(i).types(dirIdx)))==1 % single directional cover change
                rec_cg(i).type = rec_cg(i).types(1);
            else
                % Confirm the transition change
                startIdx = 1;
                tranTypes = unique(rec_cg(i).types(:),'stable'); % the types according to time
                rec_cg(i).type = tranTypes;
                for i_type = 2:length(tranTypes)
                    transitionLoc = find(rec_cg(i).types(:)==tranTypes(i_type),1,'first');
%                     if rec_cg(i).years(transitionLoc) ==2020
%                         TempBreak = 1;
%                     end
                    rec_cg(i).types(startIdx:transitionLoc-1) = tranTypes(i_type-1);
                    startIdx = transitionLoc;
                    rec_cg(i).transitionYear(i_type-1) = rec_cg(i).years(transitionLoc); % Add the transitional year

                end
                rec_cg(i).types(startIdx:end) = rec_cg(i).types(end);               
            end
        end
%         if ~sum(ismember([1,2,3,4,5,6,7,8,9],rec_cg(i).type))
%             fprintf('Check');
%         end
    end
    %% Cover type refinement for each pixel based on the whole segment
    pixelIDs = unique([rec_cg.pos]);
    for i_p = 1:length(pixelIDs)
        pixelIdx = find([rec_cg.pos]==pixelIDs(i_p));
        
        types = [];
        segmentIDs = []; %To save the segment ID
        for i_s = pixelIdx            
            types = [types;rec_cg(i_s).types];
            segmentIDs = [segmentIDs;i_s*ones(length(rec_cg(i_s).types),1)];
        end
        
        TargetTypes = [Tidalmarsh,Tidalflats];
        for i_ty = 1:length(TargetTypes)
            TargetType = TargetTypes(i_ty);
            if sum(ismember(types,TargetType))
            % possible tidal marsh pixel to refine the classes
                idxTarget = ismember(types,[TargetType,Openwater]); % Temporal change to open water is possible
                dirChange = find([0;diff(idxTarget)]~=0); % Find the conversion time
                if length(dirChange)==2
                % There are two times of cover conversion between tidal
                % marsh/others
                    if dirChange(2)-dirChange(1)>10
                        % The middle status is longer than 10 years then do nothing
                    elseif types(1)==types(end)
                        % The beginning and end of the types are same
                        idxNonWater = ~ismember(types,Openwater);
                        majVoteType = mode(types(idxNonWater));
                        % Majority vote to define the consistent cover type
                        types(idxNonWater) = majVoteType;
                        for i_s = pixelIdx                 
                            rec_cg(i_s).types = types(ismember(segmentIDs,i_s));
                            rec_cg(i_s).type = unique(rec_cg(i_s).types,'stable');% the types according to time
                            if rec_cg(i_s).transitionYear(1)>0 && length(rec_cg(i_s).type)==1
                            % No transitional year
                                rec_cg(i_s).transitionYear = 0;
                            elseif length(rec_cg(i_s).type)>1
                                % Update transitional year
                                rec_cg(i_s).transitionYear = [];
                                startIdx = 1;
                                tranTypes = unique(rec_cg(i_s).types(:),'stable'); % the types according to time                            
                                for i_type = 2:length(tranTypes)
                                    transitionLoc = find(rec_cg(i_s).types(:)==tranTypes(i_type),1,'first');
                                    rec_cg(i_s).types(startIdx:transitionLoc-1) = tranTypes(i_type-1);
                                    startIdx = transitionLoc;
                                    rec_cg(i_s).transitionYear(i_type-1) = rec_cg(i_s).years(transitionLoc); % Add the transitional year
                
                                end
                                rec_cg(i_s).types(startIdx:end) = rec_cg(i_s).types(end);     
    
                            end
                        end
                    end
                elseif length(dirChange)>2
                % Multiple times of conversion occur  
                    idxNonWater = ~ismember(types,Openwater);
                    majVoteType = mode(types(idxNonWater));
                    % Majority vote to define the consistent cover type
                    types(idxNonWater) = majVoteType;
                    for i_s = pixelIdx                 
                        rec_cg(i_s).types = types(ismember(segmentIDs,i_s));
                        rec_cg(i_s).type = unique(rec_cg(i_s).types,'stable');% the types according to time
                        if rec_cg(i_s).transitionYear(1)>0 && length(rec_cg(i_s).type)==1
                        % No transitional year
                            rec_cg(i_s).transitionYear = 0;
                        elseif length(rec_cg(i_s).type)>1
                            % Update transitional year
                            rec_cg(i_s).transitionYear = []; % Initialize
                            startIdx = 1;
                            tranTypes = unique(rec_cg(i_s).types(:),'stable'); % the types according to time                            
                            for i_type = 2:length(tranTypes)
                                transitionLoc = find(rec_cg(i_s).types(:)==tranTypes(i_type),1,'first');
                                rec_cg(i_s).types(startIdx:transitionLoc-1) = tranTypes(i_type-1);
                                startIdx = transitionLoc;
                                rec_cg(i_s).transitionYear(i_type-1) = rec_cg(i_s).years(transitionLoc); % Add the transitional year
            
                            end
                            rec_cg(i_s).types(startIdx:end) = rec_cg(i_s).types(end);     

                        end
                    end
                end
                % If cover change occurs only one time between tidal marsh and
                % others; It is a real change; Do nothing
            end
        end
    end

%     
    %% Change characterization for pulse disturbance and calculate the severity
    for i = 1:l_pos

      %% LCMAP Rule
%             avg_ref_NIR = coef(1,4) + coef(2,4) * datenum(rec_cg(i).years,7,1); 
%             avg_ref_SWIR1 = coef(1,5) + coef(2,5) * datenum(rec_cg(i).years,7,1);
%             BRs = (avg_ref_NIR-avg_ref_SWIR1)./(avg_ref_NIR+avg_ref_SWIR1+eps);
%             if abs(BRs(end) - BRs(1)) > 0.05 
        % Define the gradual change with no cover conversion by using the
        % future prediction 
        % IF grudual change is occurring although cover conversion not
        % yet happened
        if rec_cg(i).types(1)==rec_cg(i).types(end)
            coefYear = coef;
            coefYear(1,:) = coef(1,:) + coef(2,:)*datenum(rec_cg(i).years(1)+predictYear,7,1);
            inputVariables = [];
            inputVariables(1,:) = [reshape(coefYear,9*12,[])',rmse];
            [predictType] = classRF_predict(inputVariables,modelRF); % class
            if predictType ~=  rec_cg(i).type
                % Type conversion occurs in future
                if rec_cg(i).type < Openwater
                    rec_cg(i).changeType = 5*10 + predictType; % Gradual conditional change of tidal wetlands without cover change
                else
                    rec_cg(i).changeType = 5*10; % Gradual conditional change of others
                end
            end
        end
        % Define the gradual cover conversion     
        if length(rec_cg(i).type)>1 && sum(ismember(rec_cg(i).type,[Tidalmarsh,Mangrove,Dieback,Tidalflats]))
            rec_cg(i).changeType = 4; % Gradual cover change of tidal wetlands occurs
        else
            rec_cg(i).changeType = 9; % other gradual cover changes
        end
        % Define the abrupt change
        if change_prob(i) == 1            
            if rec_cg(i).type(end)>Tidalflats  && rec_cg(i+1).type(1)>Tidalflats
                %% Sometimes, the segment includes the press disturbance then remain its change type
                %% Ignore the condition change of Open water and uplands
                if rec_cg(i).type(end) ~= Openwater && rec_cg(i+1).type(1)~=Openwater && ~sum(ismember([Tidalmarsh,Mangrove],rec_cg(i).type))
                    rec_cg(i).changeType = 8; % conversion between non-water-related covers
                end
            else
                %% Coastal wetlands-related change
                if rec_cg(i).type(end) == rec_cg(i+1).type(1)
                    rec_cg(i).changeType = 3; % Condition change within the same cover tye                
                elseif rec_cg(i).type(end) <=Tidalflats && rec_cg(i+1).type(1) <=Tidalflats
                    rec_cg(i).changeType = 2; % coversion within the tidal wetlands
                else
                    rec_cg(i).changeType = 1; % Cover change of tidal wetland and others
                end
            end
        end
    end
    
    %% Calculate the severity
    for i = 1:l_pos
        severity = 0; % Initialize        
        % Sometimes, one segment would exist both pulse and press disturbance
        if any(rec_cg(i).transitionYear) && sum(ismember([Tidalmarsh,Mangrove],rec_cg(i).type))
        % Press disturbance with cover transition
        % Calculate the accumulated severity for vegetated tidal wetlands
            [~,idxTranYr] = ismember(rec_cg(i).transitionYear,rec_cg(i).years);
            idxTranYr = [1,idxTranYr] ; % Add the first year to generate the segments
            severity = zeros(length(rec_cg(i).years),1); 
            for i_tr = 2:length(idxTranYr)
                idxYr = idxTranYr(i_tr-1):idxTranYr(i_tr);
                % For each segment
                if rec_cg(i).type(i_tr-1)< Dieback
                % Cover change of vegetated tidal wetlands 
                    avg_ref_NDVI = rec_cg(i).coefs(1,8) + rec_cg(i).coefs(2,8) * datenum(rec_cg(i).years(idxYr),7,1);
                    severity(idxYr) = abs(round((avg_ref_NDVI(1)-avg_ref_NDVI)./(avg_ref_NDVI(1)+eps)*100));
                    % abs --> NDVI could be either increase or decrease
                    % when converting to other types
                    severity(idxTranYr(i_tr)) = 100; % Cover conversion and severity comes to 100
                elseif rec_cg(i).type(i_tr-1) >= Dieback && rec_cg(i).type(i_tr)<Dieback
                    severity(idxTranYr(i_tr)) = -100; % Recovered to tidal wetland
                end
            end  
%         elseif rec_cg(i).changeType == 5 && rec_cg(i).type<3 % LCMAP-based
%             avg_ref_NDVI = rec_cg(i).coefs(1,8) + rec_cg(i).coefs(2,8) * datenum(rec_cg(i).years,7,1);
%             severity = round((avg_ref_NDVI(1)-avg_ref_NDVI)./(avg_ref_NDVI(1)+eps)*100);
        elseif rec_cg(i).changeType >50 && rec_cg(i).type<Dieback % Prediction-based    
        %    No transition years
        % Press disturbance without cover transition for marsh and mangrove
        % Calculate the accumulated severity (negative means recovery)
            avg_ref_NDVI = rec_cg(i).coefs(1,8) + rec_cg(i).coefs(2,8) * datenum(rec_cg(i).years,7,1);
            severity = abs(round((avg_ref_NDVI(1)-avg_ref_NDVI)./(avg_ref_NDVI(1)+eps)*100));
        end
        
        %% Process the abrupt disturbance
        if change_prob(i) == 1 
        % Pulse disturbance
            if rec_cg(i).type(end)<Dieback
                % Tidal wetland before the change
                if rec_cg(i).changeType==1 || rec_cg(i).changeType==2
                % Cover change 
                    severity(end) = 100;
                else
                % Condition change (positive means severity; negative means recovery)
                % Use magnitude
                    relativeSeverity = -rec_cg(i).magnitude(8)/(rec_cg(i).coefs(1,8) + rec_cg(i).coefs(2,8)*datenum(rec_cg(i).years(end),7,1)+eps);
                % Use July1st
                    greenPrevious = rec_cg(i).coefs(1,8) + rec_cg(i).coefs(2,8) * datenum(rec_cg(i).years(end),7,1);
                    greenDisturbed = rec_cg(i+1).coefs(1,8) + rec_cg(i+1).coefs(2,8) * datenum(rec_cg(i).years(end),7,1);
                    relativeSeverity = (greenPrevious-greenDisturbed)/greenPrevious;
                    severity(end) = round(relativeSeverity*100);                
                end
            elseif rec_cg(i+1).type(1)<Dieback && rec_cg(i).changeType==1
                % Change to tidal wetlands from non-tidal wetlands
                severity(end) = -100; % Totally recovery
            end
        end
        % Define the boundary of the severity
        severity(severity>100) = 100;
        severity(severity<-100) = -100;
        rec_cg(i).severity = severity;
    end
     %% Check whether gradual conditional change occurs under press disturbance (LCMAP-like 0.05) and severity
  
    %% Fill in the gaps 
    pixelLocs = unique(pos);
    for i_loc = 1:length(pixelLocs)
        pixelLoc = pixelLocs(i_loc);
        pixelIdx = find(pos==pixelLoc);
        
        % Exist change point
        for i_p = 1:length(pixelIdx)
            pixelId = pixelIdx(i_p);
            if i_p == 1  
                % The initialization begins later than the interested
                % years  
                if years(1) ~= rec_cg(pixelId).years(1)
                    gapYear = find(years==rec_cg(pixelId).years(1));
                    rec_cg(pixelId).years = [years(1:gapYear-1)';rec_cg(pixelId).years];
                    rec_cg(pixelId).types = [ones(gapYear-1,1)*rec_cg(pixelId).types(1);rec_cg(pixelId).types];
                    if length(rec_cg(pixelId).severity)>1
                        rec_cg(pixelId).severity = [zeros(gapYear-1,1);rec_cg(pixelId).severity];
                    end
                end
            else
                % No time series model for some years and gap exists
                % between two segments
                gapYearNum = rec_cg(pixelId).years(1) - rec_cg(pixelId-1).years(end) - 1;
                if  gapYearNum > 0                       
                    gapYear = find(years==rec_cg(pixelId).years(1));
                    rec_cg(pixelId).years = [years(gapYear-gapYearNum:gapYear-1)';rec_cg(pixelId).years];
                    rec_cg(pixelId).types = [ones(gapYearNum,1)*rec_cg(pixelId).types(1);rec_cg(pixelId).types];
                    if length(rec_cg(pixelId).severity)>1
                        rec_cg(pixelId).severity = [zeros(gapYearNum,1);rec_cg(pixelId).severity];
                    end
                end                    
            end
        end
        
    end
    
end
 