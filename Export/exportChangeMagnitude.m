function exportChangeMagnitude(folderpath_Decode,maskImage, years, varargin)
% EXPORTCHANGEMAP This is to export the change map based on the DECODE
% results.
%
% INPUT:
%
%   folderpath_decode:        Locate to DECODE working folder, in which the
%                           change folder <TSFitLine> is necessery, and
%                           this folder was created by <DECODE.m>.
%
%   years:                  [Array] The years of change map. 

%   Outputs:
%   break (optional):            [false/true] Export the break doy or not.
%   break_yyyy.tif in folder <outMap>
%   pixel value: doy (doy indicates DOY) 
%   uint16
%   changeType (optional):       [false/true] Export the change type or not.
%   changeType_yyyy.tif in folder <outMap> (yyyy means the year)
%   uint8
%   x's range is between 0 to 2
%   0 => stable
%   1 => cover change
%   2 => condition change


%   cover (optional):       [false/true] Export the annual cover map or not.  
%   cover_yyyy.tif in folder <outMap>
%   uint 8
%   pixel value: indicates the type of the pixel
%   1 => permanent water body
%   2 => vegetated wetland
%   9 => others: tidal flat upland soil beach
%   
%   tideInfluence (optional)   Export the tidal influence map or not
%   tide_yyyy.tif in folder <outMap>
%   unit 16
%   pixel value: indicates the degree of the tidal influence
%   coefficient of c2 * 10e3

%   msg (optional):         [false/true] Display processing status (default
%                           value: false)


%
% REFERENCE(S):
%   
%   Yang, Xiucheng et al. 
%
%   Zhu, Zhe, et al. "Continuous monitoring of land disturbance based on
%   Landsat time series." Remote Sensing of Environment 238 (2020): 111116.

 
% AUTHOR(s): Xiucheng Yang, Zhe Zhu, & Shi Qiu
% DATE: May. 7, 2021
% COPYRIGHT @ GERSLab
%

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'GRIDobj'));

% optional
p = inputParser;
addParameter(p,'changeMag', true); % Export change type map
addParameter(p,'msg', true); % not to display info
parse(p,varargin{:});
msg = p.Results.msg;
isChangeMag = p.Results.changeMag;

[~, foldername_working] = fileparts(folderpath_Decode);
if msg  
        fprintf('Start to export maps for %s\r\n', foldername_working);
end



folderpath_tsf = fullfile(folderpath_Decode, globalsets.FolderDetection);
% folderpath_tsf = fullfile(folderpath_Decode, [globalsets.FolderDetection,'Temp']);
% folderpath_tsf = strrep(folderpath_tsf,'DECODEResults','DECODETest');

folderpath_map = fullfile(folderpath_Decode, [globalsets.FolderMap]);
% folderpath_map = strrep(folderpath_map,'DECODEResults','DECODETest');

if ~isfolder(folderpath_map)
    mkdir(folderpath_map);
end

tic

%% get metadata
load(fullfile(folderpath_Decode, 'metadata.mat'));
nrows = metadata.nrows;
ncols = metadata.ncols;
nbands = metadata.nbands; % 7 Landsat bands + 1 QA band + 1 wl band 
nIds = 5; % NDVI, MNDWI, TC Transformation
% dimension and projection of the image
jiDim = [ncols,nrows];
% max number of maps
max_n = length(years);

% produce tidal influence map
changeMagMap = -32768*ones(nrows,ncols,max_n,'int16'); % water level coefficient

% cd to the folder for storing recored structure
% cd(v_input.name_rst);
records = dir(fullfile(folderpath_tsf,'record_change*.mat')); % folder names
num_line = size(records,1);
if num_line~=nrows
    fprintf('Some rows are missing for %s\r',folderpath_tsf);
    return
end
for line = 1: num_line
    % show processing status
    if line/num_line < 1
        fprintf('Processing %.2f percent for Line %d\r',100*(line/num_line),line);
    else
        fprintf('Processing %.2f percent for Line %d\n',100*(line/num_line),line);
    end
   
    % load one line of time series models
    load(fullfile(folderpath_tsf,records(line).name)); %#ok<LOAD>
    
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
    % change type
    dist = [rec_cg.dist];
    
%     % cover type
%     type = [rec_cg.type];
    % change vector magnitude
    mag = [rec_cg.magnitude];
    % reshape magnitude
    mag = reshape(mag,nbands+nIds-2,[]);
    % coefficients
    coefs = [rec_cg.coefs];
    coefs = reshape(coefs,9,nbands+nIds-2,[]);
    
    
    for i = 1:l_pos % -1: segment of time series will not have change record!
        % get row and col
        [I,J] = ind2sub(jiDim,pos(i));
        
        % initialize pixels have at least one model
        if sum(changeMagMap(J,I,:) == -32768) == max_n
            % Mask regions
            changeMagMap (J,I,:) = 0;
        end
        
        % Produce the change and break maps
        if change_prob(i) == 1 % ephemeral change
            break_type = dist(i);
            break_year = datevecmx(t_break(i));
            break_year = break_year(1);
                        
            changeMagMap(J,I,years == break_year) = mag(8,i); %doy
    
        end
    end
end

geotif_obj = metadata.GRIDobj;
% Read the mask layer
sampleMaskImage = geotiffread(maskImage);

for i_yr = 1: length(years)
    yr = years(i_yr);
    if isChangeMag
        geotif_obj.Z = int16(changeMagMap(:,:,i_yr));
        geotif_obj.Z(sampleMaskImage==0) = -32768;
        GRIDobj2geotiff(geotif_obj, fullfile(folderpath_map, ['changeMag_',num2str(yr),'.tif']));
    end  
end
if msg
    fprintf('Finished exporting change map for %s with %0.2f mins\r\n', foldername_working, toc/60); 
end
end