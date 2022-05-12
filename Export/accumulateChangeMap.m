function accumulateChangeMap(folderpath_Decode, wetlandMask, years, varargin)
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
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'GRIDobj'));

% optional
p = inputParser;
addParameter(p,'minsize', 1); % all size change object
addParameter(p,'ctype', 'all'); % include all change types
addParameter(p,'msg', false); % not to display info
parse(p,varargin{:});
minsize = p.Results.minsize;
msg = p.Results.msg;
ctype = p.Results.ctype;

tic
[~, foldername_working] = fileparts(folderpath_Decode);
if msg
    fprintf('Start to accumulate change maps with %s change between %d and %d for %s\r\n', ctype, min(years), max(years), foldername_working);

end

mapfolder = fullfile(folderpath_Decode, globalsets.FolderMap);
        
accmap_outfilepath = fullfile(mapfolder, ...
    sprintf('%s_%s_%d_%d.tif', 'accChange', ctype, min(years), max(years)));

accMapGridobj = [];

%% typedoy data first
distmap_name = 'changeType';
yearlymaps = dir(fullfile(mapfolder, 'changeType_*.tif')); % this has change types
% if isempty(yearlymaps)
%     distmap_name = 'break';
%     yearlymaps = dir(fullfile(mapfolder, 'break_*.tif')); % this has break doy
% end

%% Check files are ready
if isempty(yearlymaps)
    fprintf('No yearly change maps at %s!\r\n', mapfolder);
    return;
else
    years_nomap = [];
    for imap = 1: length(years)
        yr = years(imap);
        yearlymap_filepath = fullfile(mapfolder, [distmap_name, '_', num2str(yr),'.tif']);
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
    if msg
        fprintf('Accumulating the change map for year %d\r', yr);
    end
    yearlymap_filepath = fullfile(mapfolder, [distmap_name, '_', num2str(yr),'.tif']);
    if isempty(accMapGridobj)
        accMapGridobj = GRIDobj(yearlymap_filepath); % first map gives to trget map
        accMapGridobj.Z = int16(accMapGridobj.Z);
        % select which change type
        switch ctype
            case 'cover' % type = 1
                accMapGridobj.Z(accMapGridobj.Z~=1 | accMapGridobj.Z~=4) = 0; % remain cover change
            case 'condition' % type = 2 || 3
                accMapGridobj.Z(accMapGridobj.Z~=2 | accMapGridobj.Z~=3) = 0; % exclude cover change
            case 'transition' 
                accMapGridobj.Z(accMapGridobj.Z~=4) = 0; % remain cover change
            case 'abrupt'
                accMapGridobj.Z(accMapGridobj.Z>2) = 0; % Considering the abrupt changes
        end
    
        mask_changmap = accMapGridobj.Z > 0 & accMapGridobj.Z < 9; % Stable = 9 ;
        
        if minsize > 1
            mask_changmap = mask_changmap & bwareaopen(mask_changmap, minsize, 8); % remove small objects with 8-conn
        end
        accMapGridobj.Z(mask_changmap) = yr;
        accMapGridobj.Z(~mask_changmap) = 9999; % label as 9999 for the removed pixels. 
        clear yearlymap_filepath mask_changmap;
        continue;
    end
    yearlymap = GRIDobj(yearlymap_filepath);
    
    % select which change type
    switch ctype
%         case 'all' % do not need to select
        case 'cover' % type = 1
            yearlymap.Z(yearlymap.Z~=1 | yearlymap.Z~=4) = 0; % remain cover change
        case 'condition' % type = 2 || 3
            yearlymap.Z(yearlymap.Z~=2 | yearlymap.Z~=3) = 0; % exclude cover change
        case 'transition' 
            yearlymap.Z(yearlymap.Z~=4) = 0; % remain cover change
        case 'abrupt'
            yearlymap.Z(yearlymap.Z>2) = 0; % Considering the abrupt changes
    end
    
    mask_changmap = yearlymap.Z >0& yearlymap.Z < 9;
            
    if minsize > 1
        mask_changmap = mask_changmap & bwareaopen(mask_changmap, minsize, 8);  % remove small objects with 8-conn
    end
    accMapGridobj.Z(mask_changmap) = yr;
    clear yearlymap yearlymap_filepath mask_changmap;
end

if isfile(wetlandMask)
    fprintf('Mask is used: %s \n',wetlandMask);
    mask = GRIDobj(wetlandMask);
    accMapGridobj.Z(accMapGridobj.Z==9999 & mask.Z==1) = 0;
end

accMapGridobj.Z = uint16(accMapGridobj.Z);
% save out
GRIDobj2geotiff(accMapGridobj, accmap_outfilepath);


% % % save out as a folder
% % if ~isempty(outputfolder)
% %     copyfile(accmap_outfilepath, outputfolder);
% % end

if msg
    fprintf('Finished accumluting change maps for %s with %0.2f mins\r\n', foldername_working, toc/60); 
end



