function createSyntheticImage(folderpath_decode, folderpath_wl, yy, mm, dd, varargin)
%createSyntheticImage This function is to create synthetic Landsat images
%based on the results of continous change detection (CCD).
%
% INPUT:
%
%   folderpath_cold:        Locate to COLD working folder, in which the
%                           change folder <TSFitLine> is necessery, and
%                           this folder was created by <COLD.m>.
%
%   year:                   [Number] The year of synthtic image. 
%
%   mm:                     [Number] The month or DOY of synthtic image. 
%
%   dd (optional):          [Number] The day of synthtic image. When dd =
%                           [], this function will consider mm as DOY.
%
%   compress (optional):    [true/false] Compress all the files into a tar.
%                           (default value: false)
%
%   msg (optional):         [false/true] Display processing status (default
%                           value: false)
%
% REFERENCE(S):
%
%   Zhu, Zhe, et al. "Generating synthetic Landsat images based on all
%   available Landsat data: Predicting Landsat surface reflectance at any
%   given time." Remote Sensing of Environment 162 (2015): 67-83.
%
% EXAMPLE OF USE:
%
%   > To create a synthetic image for 07/01/2000
%
%     createSyntheticImage('/lustre/scratch/qiu25856/TestGERSToolbox/h029v005/',
%     2000, 7, 1)
%
%   > To create a synthetic image for DOY 171, 2000
%
%     createSyntheticImage('/lustre/scratch/qiu25856/TestGERSToolbox/h029v005/',
%     2000, 171)
% 
% AUTHOR(s): Zhe Zhu and Shi Qiu
% DATE: Feb. 7, 2021
% COPYRIGHT @ GERSLab

    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Detection'));
    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'GRIDobj'));
    
    % optional
    p = inputParser;
    addParameter(p,'compress', false); % not to generate the qa band of synthetic data
    addParameter(p,'msg', true); % not to display info
    parse(p,varargin{:});
    msg = p.Results.msg;
    compress = p.Results.compress;
    
    % define date
    if yy>0
        if mm>0&&~isempty(dd)&&dd>0 % using yyyy/mm/dd
            % converted to julian date
            j_date=datenum(yy,mm,dd);
%             j0_date=datenum(yy,1,0);
%             doy = j_date-j0_date;
            % date for show i.e. 19990815
            s_date=yy*10000+mm*100+dd;
            s_doy = yy*1000 + day(datetime(yy,mm,dd),'dayofyear');
%             s_date = yy*1000 + doy;
        elseif mm>0 % using yyyy DOY
            j_date = datenum(yy,1,0) + mm;
            s_date = yy*1000 + mm;
            s_doy = s_date;
        end
    end
    
    
    [~, foldername_working] = fileparts(folderpath_decode);
    if msg
        fprintf('Start to synthetic image on %d for %s\r\n', s_date, foldername_working);
    end

    folderpath_tsf = fullfile(folderpath_decode, 'TSFitLine');
%     if ~isfolder(folderpath_tsf) % also support the old name <TSFitMap>
%         folderpath_tsf = fullfile(folderpath_decode, 'TSFitMap');
%     end
    %% Image folder
    folderpath_syn = fullfile(folderpath_decode, 'SyntheticImage');

    
    %% get metadata
    load(fullfile(folderpath_decode, 'metadata.mat'));
    
    filename_syn = sprintf('SL_%s_%s_%d', metadata.region, [metadata.tile(2:4), metadata.tile(6:8)], s_date);  % S: Synthetic % L: Landsat
    
    %% check exist: .tar
    if isfile(fullfile(folderpath_syn, [filename_syn, '.tar']))
        if msg
            fprintf('Already exist synthetic image %s\r\n', [filename_syn, '.tar']);
        end
        return;
    end
    
    %% check exist: folder with uncompress
    if isfolder(fullfile(folderpath_syn, filename_syn))
        files_have = dir(fullfile(folderpath_syn, filename_syn, 'SL*.tif'));
        bytes = [files_have.bytes];
        bytes = uint16(floor(bytes/10000));
        if length(files_have) == 6 && length(unique(bytes))==1% metadata.nbands
            if msg
                fprintf('Already exist synthetic image %s with bytes %d *10e4\r\n', filename_syn,unique(bytes));
            end
            return;
        end
    end
    
    %% create folder if no
    if ~isfolder(folderpath_syn)
        mkdir(folderpath_syn);
    end
    % synthetic image's name
    dir_out_syn = fullfile(folderpath_syn, filename_syn);
    if ~isfolder(dir_out_syn)
        mkdir(dir_out_syn);
    end
    
    %% Sync starts ...
    geotif_obj = metadata.GRIDobj;
    geotif_obj.Z = uint16(geotif_obj.Z);
    nrows = metadata.nrows;
    ncols = metadata.ncols;
    % dimension and projection of the image
    jiDim = [ncols,nrows];
    nbands = metadata.nbands-1; % 7 Landsat bands + 1 QA band + tide band
    ncoefs = 9; % number of model coeffs 
    
    % Load tide value this day
    path_tide = dir(fullfile(folderpath_wl,['*',num2str(s_doy),'*']));
    tide = imread(fullfile(folderpath_wl,path_tide.name));
    tide = double(tide);
    % produce synthetic data (1, 2, 3, 4, 5, 7, 6, and QA)
    SR = -9999*ones(nrows,ncols,nbands,'int16'); % 1. change here for excluding thermal
    imf=dir(fullfile(folderpath_tsf, 'record_change*.mat')); % folder names
    num_line = size(imf,1);

    %% line by line
    tic
    for line = 1: num_line
        % show processing status
        if msg
            if line/num_line < 1
                fprintf('Processing %.2f percent\r',100*(line/num_line));
            else
                fprintf('Processing %.2f percent\n',100*(line/num_line));
            end
        end
        load(fullfile(folderpath_tsf,imf(line).name));

        % postions & coefficients
        pos = [rec_cg.pos];
        % continue if there is no model available
        l_pos = length(pos);
        if l_pos == 0
            continue;
        end
        % get coefficients matrix
        coefs = [rec_cg.coefs];
        % reshape coefs
%         coefs = reshape(coefs,nbands*ncoefs,l_pos);
        coefs = reshape(coefs,size(rec_cg(1).coefs,1)*size(rec_cg(1).coefs,2),l_pos);
        % get category matrix
        category = [rec_cg.category];
        % take the first digit (number of coefficients)
        category = category - 10*floor(category/10);

        % model start, end, & break time for prediction
        model_start = [rec_cg.t_start];
        model_end = [rec_cg.t_end];
        model_break = [rec_cg.t_break];
        clear rec_cg;
        % model on the right
        ids_right = model_start > j_date;
        % model on the left
        ids_left = (model_end < j_date & model_break == 0)|(model_break <= j_date & model_break > 0);
        % id within model interval
        ids_middle = model_start <= j_date & ((model_end >= j_date & model_break == 0) | (model_break > j_date & model_break > 0));
        clear model_start model_end model_break;

        % position for model in the middle
        pos_middle = pos(ids_middle);
        % coefficients for model in the middle
        coefs_middle = coefs(:,ids_middle);
        % category for model in the middle
        category_middle = category(ids_middle);
        clear ids_middle;

        % positions for the nearest model on the right
        pos_right = pos(ids_right);
        [pos_near_right,ids_near_right] = unique(pos_right,'first');
        % coefficients for the nearest model on the right
        coefs_right = coefs(:,ids_right);
        coefs_near_right = coefs_right(:,ids_near_right);
        % category for the nearest model on the right
        category_right = category(ids_right);
        clear ids_right;
        category_near_right = category_right(ids_near_right);
        clear category_right ids_near_right;

        % postions for the nearest model on the left
        pos_left = pos(ids_left);
        [pos_near_left,ids_near_left] = unique(pos_left,'last');
        % coefficients for the nearest model on the left
        coefs_left = coefs(:,ids_left);
        clear coefs;
        coefs_near_left = coefs_left(:,ids_near_left);
        clear coefs_left;
        % category for the nearest model on the left
        category_left = category(ids_left);
        clear ids_left category;
        category_near_left = category_left(ids_near_left);
        clear category_left ids_near_left;
        % pass if there is no nearest model on the left 
        l_pos=length(pos_near_left);
        if l_pos > 0  
            % providing predictions
            for i=1:l_pos
                [I,J]=ind2sub(jiDim,pos_near_left(i));
                for j_b=1:nbands-1 % excluding QA
                    SR(J,I,j_b)=autoTSPredWL(j_date,coefs_near_left(((j_b-1)*ncoefs+1):j_b*ncoefs,i),tide(J,I));
                end
                % QA band
                SR(J,I,end) = 20 + category_near_left(i); % model forward predicted values         
            end
        end
        clear category_near_left pos_near_left;

        % pass if there is no nearest model on the right 
        l_pos=length(pos_near_right);
        if l_pos > 0  
            % providing predictions
            for i=1:l_pos
                [I,J]=ind2sub(jiDim,pos_near_right(i));
                for j_b=1:nbands-1 % excluding QA
                    SR(J,I,j_b)=autoTSPredWL(j_date,coefs_near_right(((j_b-1)*ncoefs+1):j_b*ncoefs,i),tide(J,I));
                end
                % QA band
                SR(J,I,end) = 10 + category_near_right(i); % model backward predicted values         
            end
        end
        clear category_near_right pos_near_right;

        % pass if there is no nearest model in the middle 
        l_pos=length(pos_middle);
        if l_pos > 0  
            % providing estimations
            for i=1:l_pos
                [I,J]=ind2sub(jiDim,pos_middle(i));
                for j_b=1:nbands-1 % excluding QA
                    SR(J,I,j_b)=autoTSPredWL(j_date,coefs_middle(((j_b-1)*ncoefs+1):j_b*ncoefs,i),tide(J,I));
                end
                % QA band
                SR(J,I,end) = 0 + category_middle(i); % model estimated values         
            end
        end  
        clear category_middle pos_middle l_pos;  
    end
    
    % save disturbance image
    for ib = 1: nbands
        geotif_obj.Z = SR(:,:,ib);
        if ib <=5
            bandname = sprintf('SRB%d', ib);
        end
        if ib == 6
            bandname = sprintf('SRB%d', 7); % swir 2 band of Landsat
        end
        if ib == 7
            bandname = sprintf('BTB%d', 6); % thermal band
        end
        if ib == 8
            bandname = 'QA';
        end
%         if ib<2 ||ib>4
%             continue
%         end
        if ib > 6
            continue
        end
        GRIDobj2geotiff(geotif_obj, fullfile(dir_out_syn, sprintf('%s_%s.tif', filename_syn, bandname)));
    end
    clear geotif_obj SR;
    if compress
        tar( fullfile(folderpath_syn, [filename_syn, '.tar']), {fullfile(dir_out_syn, 'SL*.tif')});
        rmdir(fullfile(folderpath_syn, filename_syn), 's'); % delete the data folder
    end
    if msg
        fprintf('Finished creating synthetic Landsat image for %s with %0.2f mins\r\n', filename_syn, toc/60); 
    end

end
