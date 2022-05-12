function recordCoverChangeToRecCg(task,ntasks,folderpath_Decode,folderpath_Tide,years,varargin)
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

% optional
p = inputParser;

addParameter(p,'msg', true); % not to display info
parse(p,varargin{:});
msg = p.Results.msg;

[~, foldername_working] = fileparts(folderpath_Decode);
if msg
    fprintf('Start to calculate the category and change type for %s\r\n', foldername_working);
end
% 
folderpath_tsf = fullfile(folderpath_Decode, globalsets.FolderDetection);

% temporal folder to reserve the updated rec_cg
folderpath_tsfTemp = fullfile(folderpath_Decode, [globalsets.FolderDetection,'Type']);

if ~isfolder(folderpath_tsfTemp)
    mkdir(folderpath_tsfTemp);
end

tic

%% get metadata
load(fullfile(folderpath_Decode, 'metadata.mat'));
nrows = metadata.nrows;
ncols = metadata.ncols;
nIds = 5;
nbands = metadata.nbands; %7 Landsat bands + 1 QA band + 1WL band + 1 SPEI band
% dimension and projection of the image
jiDim = [ncols,nrows];

% water level bands in January and July (used for classification)
dir_files=[dir(fullfile(folderpath_Tide,'*0101.tif'));dir(fullfile(folderpath_Tide,'*0701.tif'))];
fileNames = {dir_files.name};
fileNames= char(fileNames);

% % sort according to yeardoy
ymdDate = str2num(fileNames(:, 5:12)); 
[~, sort_order] = sort(ymdDate);
fileNames = fileNames(sort_order, :);

% Filter according to the interested years
yr = str2num(fileNames(:, 5:8)); %#ok<*ST2NM>
idxYr = ismember(yr,years);
month = str2num(fileNames(:, 9:10));
dy = str2num(fileNames(:, 11:12));
yr = yr(idxYr);
mt = month(idxYr);
dy = dy(idxYr);
fileNames = fileNames(idxYr,:);

idxJan = mt==1;
idxJul = mt==7;
sdate = datenum(yr, mt, dy);
yr = unique(yr);
sJan = sdate(idxJan);
sJul = sdate(idxJul);

num_y = size(sJul,1);
wlJan = zeros(num_y,ncols);
wlJul = zeros(num_y,ncols);


tempWLJan = zeros(nrows,ncols,num_y,'uint32');
tempWLJul = zeros(nrows,ncols,num_y,'uint32');

for i = 1:num_y
    imJan = fullfile(folderpath_Tide,fileNames(i*2-1, :));
    imJul = fullfile(folderpath_Tide,fileNames(i*2, :));

    tempWLJan(:,:,i) = imread(imJan); 
    tempWLJul(:,:,i) = imread(imJul); 
end

records = dir(fullfile(folderpath_tsf,'record_change*.mat')); % folder names
num_line = size(records,1);

if num_line~=nrows
    fprintf('Some rows are missing for %s\r',folderpath_tsf);
    return
end

%% Parallel tasks on the row datasats
tasks_per = ceil(num_line/ntasks);
start_i = (task-1)*tasks_per + 1;
end_i = min(task*tasks_per, num_line);
fprintf('-----Totall %d lines need to be processed in this core between %d and %d ----\n',end_i-start_i+1,start_i,end_i);
for i_task = start_i:end_i
% for line = 1:num_line
    line = i_task;
    % show processing status
    if line/num_line < 1
        fprintf('Processing %.2f percent for Line %d\r',100*(line/num_line),line);
    else
        fprintf('Processing %.2f percent for Line %d\n',100*(line/num_line),line);
    end
    if isfile(fullfile(folderpath_tsfTemp, sprintf('record_change_r%05d.mat.part', line)))
        delete (fullfile(folderpath_tsfTemp, sprintf('record_change_r%05d.mat.part', line)));
    end
    
    if isfile(fullfile(folderpath_tsfTemp, sprintf('record_change_r%05d.mat', line)))
        fprintf('Having existed line %d\n',line);
        continue
    end
    % load one line of time series models
    load(fullfile(folderpath_tsf,records(line).name)); %#ok<LOAD>

    % postions
    pos = [rec_cg.pos];
    l_pos = length(pos);

    % continue if there is no model available for this row
    if l_pos == 0
        rec_cg.type = [];
        rec_cg.dist = [];
        filepath_rcg = fullfile(folderpath_tsfTemp, sprintf('record_change_r%05d.mat', line)); % r: row
        save([filepath_rcg, '.part'] ,'rec_cg'); % save as .part
        clear rec_cg;
        movefile([filepath_rcg, '.part'], filepath_rcg);  % and then rename it as normal format
        continue
    end
    
    % remove the old results
    if isfield(rec_cg,'transitionYear')
        rec_cg = rmfield(rec_cg,{'years','pred_jan','pred_jul','types','transitionYear'});
    end
    % Two fields to reserve the land cover type and disturbance type
    types = num2cell(zeros(l_pos,1));
    dists = num2cell(zeros(l_pos,1));
    transYr = num2cell(zeros(l_pos,1));
    [rec_cg(:).type] = deal(types{:});
    [rec_cg(:).dist] = deal(dists{:});
    [rec_cg(:).transitionYear] = deal(transYr{:});
    
    [~,j] = ind2sub(jiDim,pos(1));
    for i = 1:num_y
        wlJan(i,:) = tempWLJan(j,:,i); 
        wlJul(i,:) = tempWLJul(j,:,i); 
    end
    % break time
    t_break = [rec_cg.t_break];
    % change probability
    change_prob = [rec_cg.change_prob];
    % change vector magnitude
    mag = [rec_cg.magnitude];
    % reshape magnitude
    mag = reshape(mag,nbands+nIds-2,[]);
    % coefficients
    coefs = [rec_cg.coefs];
    coefs = reshape(coefs,9,nbands+nIds-2,[]);
    
    for i = 1:l_pos
    % sub class type classification based on rules
        
        [I,J] = ind2sub(jiDim,pos(i));

        i_start = rec_cg(i).t_start;
        yr_start = ymd(datetime(i_start,'ConvertFrom','datenum'));
        i_end = rec_cg(i).t_end;
        yr_end = ymd(datetime(i_end,'ConvertFrom','datenum'));
        
        yr_end = min(yr_end(1),years(end));
        
        x_plot= rec_cg(i).t_start:rec_cg(i).t_end;
        yr_plot = (yr_start:yr_end)';
        rec_cg(i).years = yr_plot;
        
        yr_idx = yr_plot - yr_start + 1;
        [~,idxYr] = ismember(yr_plot,yr);
%         idxYr = find(yr==yr_plot);
        date_jan = sJan(idxYr);
        date_jul = sJul(idxYr);
        wl_jan = wlJan(idxYr,I);
        wl_jul = wlJul(idxYr,I);
        pred_jan = autoTSPredWL(date_jan,rec_cg(i).coefs(:,:,:),wl_jan);
        pred_jul = autoTSPredWL(date_jul,rec_cg(i).coefs(:,:,:),wl_jul);
        rec_cg(i).pred_jan = pred_jan;
        rec_cg(i).pred_jul = pred_jul;
        seasonalCoefs = rec_cg(i).coefs([3,4],:,:);
        seasonalAmplitude = sqrt(seasonalCoefs(1,:).*seasonalCoefs(1,:)+seasonalCoefs(2,:).*seasonalCoefs(2,:));
        for y = 1:length(yr_plot)
            if pred_jul(y,9)>0 && pred_jan(y,9)>0 && pred_jul(y,9)-pred_jul(y,8)>0
                rec_cg(i).types(y) = 1; %Yearlong water body
            elseif pred_jul(y,8)> 2000 && seasonalAmplitude(8) > 500
                rec_cg(i).types(y) = 2; % Vegetation
            elseif pred_jul(y,8)> 2000 && pred_jan(y,8)> 2000 && seasonalAmplitude(8) > 250
                rec_cg(i).types(y) = 2; % wetland affected by   
            elseif pred_jul(y,8)> 3500 && pred_jan(y,8)> 3500   
                rec_cg(i).types(y) = 2; % dense vegetation in the whole year in the rainforest
            else
                rec_cg(i).types(y) = 9; % others
            end
        end
        
        %% To check the transitional cover change
        if length(yr_plot)<5 % Transactional change cover occurs for segment longer than 4 years
            rec_cg(i).types(:)=mode(rec_cg(i).types(:));
            rec_cg(i).type = rec_cg(i).types(1);
        elseif ~isempty(intersect(rec_cg(i).types(1:2),rec_cg(i).types(end-1:end)))
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
                rec_cg(i).types(1) = mode(rec_cg(i).types(2:4));
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
                    if yr_plot(transitionLoc) ==2020
                        TempBreak = 1;
                    end
                    rec_cg(i).types(startIdx:transitionLoc-1) = tranTypes(i_type-1);
                    startIdx = transitionLoc;
                    rec_cg(i).transitionYear(i_type-1) = yr_plot(transitionLoc); % Add the transitional year
                end
                rec_cg(i).types(startIdx:end) = rec_cg(i).types(end);
            end
        end
    end

    for i = 1:l_pos-1
        if change_prob(i) == 1
            if rec_cg(i).type(end) ~= rec_cg(i+1).type(1) % if the types are different
                rec_cg(i).dist = 1; % cover change for every class
                
            elseif rec_cg(i).type(end)~= 1 % Water body has no condition change
                rec_cg(i).dist = 2; % condition change for vegetation and others
            end
        elseif change_prob(i) == 3 && rec_cg(i).type(end)==2
            rec_cg(i).dist = 3; % ephemeral change for vegetation
        end
    end

    filepath_rcg = fullfile(folderpath_tsfTemp, sprintf('record_change_r%05d.mat', line)); % r: row
    save([filepath_rcg, '.part'] ,'rec_cg'); % save as .part
    clear rec_cg;
    movefile([filepath_rcg, '.part'], filepath_rcg);  % and then rename it as normal format
end

%% if the whole rows are calculated, remove the temporal file and rename the folder name
% records = dir(fullfile(folderpath_tsfTemp,'record_change*.mat')); % folder names
% num_records = size(records,1);
% if num_records==nrows
%     rmdir(folderpath_tsf, 's');
%     movefile(folderpath_tsfTemp,folderpath_tsf);
% else
%     fprintf('Wrong occurs when processing %s\r',folderpath_tsfTemp);
% end
end