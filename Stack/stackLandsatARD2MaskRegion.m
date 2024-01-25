function stackLandsatARD2MaskRegion(folderpath_ard, folderpath_stack, folderpath_wl,pathMaskImage,varargin)
%STACKLANDSATARD2LINE This function is to stack Landsat ARD time series
% into few row data with BIP format. Running time is estimated to ~0.25
% mins per Landsat ARD.
%
% INPUT:
%
%   folderpath_ard:         Locate to Landsat ARD folder that have all the
%                           surface reflectance (*_SR.tar) and brightness
%                           temperature (*_BT.tar). Usually we have
%                           hundreds or thousands of Landsat images.
%
%   folderpath_stack:       Locate to output folder where the stack data
%                           will be stored by lines (rows).
%
%   task (optional):        Task ID of parallel computation
%
%   ntasks (optional):      Total number of tasks of parallel computation
%
%   clear (optional):       Percentage of clear pixels (uint: %). The image
%                           with < this value will be not processed.
%                           (default value: 0 for stacking all the images in the ARD folder)
%
%   nsubrow (optional):     Number of seperated lines (rows) (default value
%                           : 10)
%
%   orbitpath (optional):   Orbit path of Landsat ('single' or 'all').
%                           'single' means generate the single Landsat path
%                           layer with geotiff format. 'all' means not to
%                           do that. (default value: single)
%
%   msg (optional)          [false/true] Display processing status (default
%                           value: false)
%
%   check (optional):       [false/true] Check the existing files, and if
%                           exist, this function will skip to process it
%                           (default value: false). Note the default
%                           "false" will let the process not to SCAN all
%                           the already stack files or folders, and this
%                           will be more efficient. However, if one more
%                           time stacking process needed, please set it as
%                           "true" for avoiding to stack the already exist
%                           data.
%
%
% RETURN:
%
%   null
%
% REFERENCE(S):
%
%   null
%
%%% 
% DATE: Feb. 5, 2023
% COPYRIGHT @ GERSLab
    % add the matlab search path of GRIDobj
%     addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Packages/GRIDobj'));

    %% Have user's inputs
    % requried
    if isempty(folderpath_ard)
        warning('No Landsat ARD folder input!');
        return;
    end
    if isempty(folderpath_stack)
        warning('No stack folder input!');
        return;
    end

    % optional
    p = inputParser;
    addParameter(p,'clear', 0); % dedault as 0 , that will stack all the images in the input folder
    addParameter(p,'npixelper', 1000); % default pixel number per file (optional)
    addParameter(p,'orbitpath', 'all'); % to create single path layer
    addParameter(p,'task', 1); % 1st task
    addParameter(p,'ntasks', 1); % single task to compute
    addParameter(p,'msg', true); % not to display info
    addParameter(p,'check', false); % not to scan history data; but once more time, please set as true, and if exist, just skip
    addParameter(p,'rigorousCheck',false); % Rigorous check slow 
    addParameter(p,'beginYear',1984); % The beginning year of the model fitting
    addParameter(p,'Collection',1); % The beginning year of the model fitting
    % request user's input
    parse(p,varargin{:});
    
    clr_pct_min = p.Results.clear;
    npixelper = p.Results.npixelper;
    
    task = p.Results.task;
    ntasks = p.Results.ntasks;
    msg = p.Results.msg;
    checkexist = p.Results.check;
    rigorousCheck = p.Results.rigorousCheck;
    beginYear = p.Results.beginYear;
    Collection = p.Results.Collection;

    % char will be better understood for users
    switch lower(p.Results.orbitpath)
        case 'single'
            singlepath = true;
        case 'all'
            singlepath = false;
    end
    
    %% Load the mask regions
    maskInfo = GRIDobj(pathMaskImage);
    maskImage = imread(pathMaskImage);
    ncols = size(maskImage,2);
%     pos = (nrows-1)*ncols+i_ids;
    [rows,cols] = find(maskImage==1);
    poss = (rows-1)*ncols + cols;
    [poss,idx] = sort(poss);
    rows = rows(idx);
    cols = cols(idx);
    pixelLocs = sub2ind(size(maskImage),rows,cols);
    maskPixelsInfo = [poss,rows,cols,pixelLocs];

    numPixels = length(poss);
    numStackFiles = ceil(numPixels/npixelper);
    fprintf('-------%d pixels and divided into %d data blocks--------\n',numPixels,numStackFiles);
    %% constant parameter    
    nbands = 9; % 7 spectral bands + wl  + Fmask QA band
    
    %% Filter for Landsat folders
    % get num of total folders start with "L"
    imfs = dir(fullfile(folderpath_ard,'L*SR.tar'));
    % filter for Landsat folders
    imfs = regexpi({imfs.name}, 'L(T05|T04|E07|C08|C09)(\w*)', 'match'); % no expand name
    % imfs = regexpi({imfs.name}, 'LE07(\w*)', 'match'); % no expand name and only for Landsat 7 data
    
    imfs = [imfs{:}];
    if isempty(imfs)
        warning('No Landsat ARD data at %s!', folderpath_ard);
        return;
    end
    % convert to char array
    imfs = vertcat(imfs{:});
    % no _SR at the end, to locate other filenames such as xxx_SR.tar,
    % xxx_BT.tar, xxx.xml, xxx_SRBx.tif, xxx_BT.tif
    imfs(:, end-2:end) = []; 
    
    % sort according to yeardoy, and this ensure each core would have same
    % image list if parallel.
    yyyymmdd = str2num(imfs(:, 16:23)); % should change for different sets

    [sortedYear, sort_order] = sort(yyyymmdd);
    imfs = imfs(sort_order, :);
    imfs = imfs(sortedYear>beginYear*10000,:);
    clear yyyymmdd sort_order;
    % number of folders start with "L"
    num_t = size(imfs,1);
    if msg
        fprintf('A total of %04d Landsat images at %s\r\n',num_t, folderpath_ard);
    end
    
    %% Assign stacking tasks to each core
    tasks_per = ceil(num_t/ntasks);
    start_i = (task-1)*tasks_per + 1;
    end_i = min(task*tasks_per, num_t);
    if msg
        fprintf('At task# %d/%d to process %d images\n', task, ntasks, end_i- start_i + 1);
    end

    % create a task folder, that will be uesed to store the divided row
    % data for all the images @ the currenr core
    folderpath_task = fullfile(folderpath_stack, sprintf('tmptaskfolder_%d_%d', task, ntasks));
    if ~isfolder(folderpath_task)
        mkdir(folderpath_task);
    end
    
    %% Start to stack the Landsat ARD from start_i th to end_i th
    
    for i = start_i:end_i %start_i:end_i
        %% locate to a certain image
        imf = imfs(i, :);
        % new filename in format of LXSHHHVVVYYYYDOYLLTTT_PPP
        stackname = subFunStackDataName(imf);
        
        %% check exist or not by DIR the last row folder
        if checkexist
            filepath_meta = fullfile(fileparts(folderpath_stack), sprintf('metadata.mat'));
            if isfile(filepath_meta) % we must have the metadata first
                if ~exist('metadata', 'var')
                    load(filepath_meta); % only load once
                end
                if ~rigorousCheck
                    %% First mode of the check, loose check and would miss some errors
                    % Assume if the last block is existing
                    irow = max(1: metadata.npixelper: metadata.npixels); % the last row folder
                    irow_end = min(irow + metadata.npixelper -1, metadata.npixels);
                    foldername_rowdata = subFunBlockFolderName(irow, irow_end);
                    rowdata_exist = dir(fullfile(folderpath_stack, foldername_rowdata, [stackname, '_*']));
                    isexist = false;
                    for iexist = 1: length(rowdata_exist)
                        try 
                            load(fullfile(folderpath_stack, foldername_rowdata,rowdata_exist(iexist).name));
                            isexist = true;
                        catch
                        end
                    end
                    if isexist
                        if msg
                            fprintf('\nAlready exist the %dth image (%s)\n', i, imf);
                        end
                        continue;
                    end
                else
                    %% Second mode of the check, rigorous check
                    flag = 0;
                    for irow = 1: npixelper: numPixels % total of 5000 rows
                        % @ the current task, only parts of the images will be reconstructed
                        irow_end = min(irow + npixelper -1, numPixels);
                        % row data's foldername R0000100010, which means the row 00001 to 00010
                        foldername_rowdata = subFunBlockFolderName(irow, irow_end);     
                        rowdata_exist = dir(fullfile(folderpath_stack, foldername_rowdata, [stackname, '_*']));
    
                        try 
                            load(fullfile(folderpath_stack, foldername_rowdata,rowdata_exist(1).name));       
                            flag = flag + 1;
                        catch
                            break
                        end
                    end
                    % If all blocks are successsfully stacked then skip it
                    if flag == metadata.nblocks
                        if msg
                            fprintf('\nAlready exist the %dth image (%s)\n', i, imf);
                        end
                        continue
                    end
                end
            end
        end
        
        %% check surface ref. and bright temperature already
        filepath_sr = fullfile(folderpath_ard, [imf, '_SR.tar']);
        filepath_bt = fullfile(folderpath_ard, [imf, '_BT.tar']);
        if isfile(filepath_sr)
            if ~isfile(filepath_bt)
                warning('No BT for %s', imf);
                continue;
            end
        else
            warning('No SR for %s', imf);
            continue;
        end
        
        if msg
            fprintf('\nProcessing the %dth image (%s)\n', i, imf);
        end
 
        %% Uncompress data
        % all .tar files will be uncompressed to a temp folder named by tmp_*
        folderpath_tmp = fullfile(folderpath_task, sprintf('tmpimage_%d', i));
        % unstar *_SR.tar and *_BT.tar to the temp folder
        try
            untar(filepath_sr, folderpath_tmp);
            if ~isfile(filepath_bt)
%                 warning('No BT for %s', imf);
%                 continue;
            else
                untar(filepath_bt, folderpath_tmp);
            end
            
        catch me
            fprintf('Error for stacking %s!\n',imf);
            fprintf([me.message, '\r\n']);
            continue;
        end
        
        
        %% Read CFmask QA band
        % @ the 1st image at 1st task, to create a metadata to save, which
        % is to backup the geo info of the geotiff.
        if task == 1 && i == 1
            % metaset generator            
%             cfmask0 = GRIDobj(fullfile(folderpath_tmp, [imf, '_QA_PIXEL.TIF']));
            metadata = [];
            metadata.GRIDobj = maskInfo;
            metadata.GRIDobj.Z = []; % set as [] for saving storage
            metadata.GRIDobj.name = []; % set as [] for saving storage
            metadata.region = imf(6:7); % CU: CONUS
            metadata.tile = sprintf('h%sv%s', imf(9:11), imf(12:14)); % tile name
            metadata.nrows = maskInfo.size(1); % record of the row size
            metadata.ncols = maskInfo.size(2); % record of the column size.
            metadata.nbands = nbands; % record of the number of bands in stack data
            metadata.npixelper = npixelper; % record of the number of rows per file
            metadata.npixels = numPixels; % record of the number of rows per file
            metadata.nblocks = ceil(numPixels/npixelper); % record of the number of rows per file
            metadata.nimages = num_t; % record of the total number of Landsat images
            metadata.cloudcover = 100 - clr_pct_min; % record of max percentage of cloud cover
            metadata.createtime = datestr(now); % record of the create time
            
            save(fullfile(fileparts(folderpath_stack), sprintf('metadata.mat')), 'metadata'); % saveas metadata
            
            % single path layer generatorni h
            if singlepath
                singlepathlayer = GRIDobj(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'LandsatPathData', 'singlepath_landsat_conus.tif'));  % at parent code's folder
                path_ard = reproject2utm(singlepathlayer, maskInfo,'method', 'nearest');
                path_ard.Z = uint8(path_ard.Z);
                GRIDobj2geotiff(path_ard, fullfile(fileparts(folderpath_stack), 'singlepath_landsat.tif'));
                clear path_ard;
            end          
        end
        if Collection==2
            cfmask = readgeoraster(fullfile(folderpath_tmp, [imf, '_QA_PIXEL.TIF']));
        elseif Collection==1
            cfmask = readgeoraster(fullfile(folderpath_tmp, [imf, '_PIXELQA.tif']));
        else
            error('Check Collection data!');
        end
        cfmask0 = cfmask;
        % convert pixel QA to fmask values
%         cfmask(bitget(cfmask0,1) == 1) = 255;
%         cfmask(bitget(cfmask0,2) == 1) = 0;
%         cfmask(bitget(cfmask0,3) == 1) = 1;
%         cfmask(bitget(cfmask0,4) == 1) = 2;
%         cfmask(bitget(cfmask0,5) == 1) = 3;
%         cfmask(bitget(cfmask0,6) == 1) = 4;
        % see more details from USGS document at https://prd-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/atoms/files/LSDS-1328_Landsat8-9-OLI-TIRS-C2-L2-DFCB-v6.pdf
        % minus one (-1)
        if Collection==2
            cfmask(bitget(cfmask0,1) == 1) = 255; % Filled
            cfmask(bitget(cfmask0,7) == 1) = 0; % Clear Land and Water [No Cloud & No Dialted Cloud]
            cfmask(bitget(cfmask0,8) == 1) = 1; % Water
            cfmask(bitget(cfmask0,5) == 1) = 2; % Cloud Shadow
            cfmask(bitget(cfmask0,6) == 1) = 3; % Snow
            cfmask(bitget(cfmask0,2) == 1) = 4; % Dilated Cloud
            cfmask(bitget(cfmask0,4) == 1) = 4; % Cloud
        elseif Collection==1
            % convert pixel QA to fmask values
            cfmask(bitget(cfmask0,1) == 1) = 255;
            cfmask(bitget(cfmask0,2) == 1) = 0;
            cfmask(bitget(cfmask0,3) == 1) = 1;
            cfmask(bitget(cfmask0,4) == 1) = 2;
            cfmask(bitget(cfmask0,5) == 1) = 3;
            cfmask(bitget(cfmask0,6) == 1) = 4;
        end
        clear cfmask0;
        clr_pct = sum(cfmask(:)<=1)/sum(cfmask(:)<255);
        clr_pct = 100*clr_pct;
        % less than 20% clear observations (default, but optional)
        if clr_pct < clr_pct_min 
            rmdir(folderpath_tmp,'s');
            fprintf('Clear observation less than %.2f percent (%.2f)\r\n', clr_pct_min, clr_pct);
            continue;
        end
        
        %% Find the collection 1's path from .xml, and rename stack name with end of path
        filepath_xml = fullfile(folderpath_tmp, [imf, '.xml']);
        fid_in=fopen(filepath_xml,'r');
        xml_str=strread(fscanf(fid_in,'%c',inf)','%s');

        if Collection == 2
            path_char=char(xml_str(strmatch('<WRS_PATH>',xml_str))); %<WRS_PATH>000
    %         path = str2double(path_char(1,11:13)); % only 1 item
            path = sscanf(path_char(1,:), '<WRS_PATH>%d<WRS_PATH>');% only 1 item
%         time_char= char(xml_str(strmatch('<SCENE_CENTER_TIME>',xml_str))); %<SCENE_CENTER_TIME>HH:MM:SS.
%         time = str2double([time_char(1,20:21),time_char(1,23:24)]); % only 1 item
        elseif Collection ==1 
            path_char=char(xml_str(strmatch('path=',xml_str))); % path=
            path = str2double(path_char(1,7:end-1)); % only 1 item
        end
        stackname = [stackname, sprintf('_%03d.mat', path)]; % landsat path range 1 to 233
       
        if isempty(path)
           fprintf('%s\n',stackname);
        end
        
        %% Load all bands
        l_num = str2double(stackname(3)); % which Landsat satellite?
        if Collection==2
            switch l_num
                case {4, 5, 7} % Landsat 4, 5, and 7
                    b1 = subFunLoadSingleBand(folderpath_tmp, imf, 'SR_B1.TIF');
                    b2 = subFunLoadSingleBand(folderpath_tmp, imf, 'SR_B2.TIF');
                    b3 = subFunLoadSingleBand(folderpath_tmp, imf, 'SR_B3.TIF');
                    b4 = subFunLoadSingleBand(folderpath_tmp, imf, 'SR_B4.TIF');
                    b5 = subFunLoadSingleBand(folderpath_tmp, imf, 'SR_B5.TIF');
                    b7 = subFunLoadSingleBand(folderpath_tmp, imf, 'SR_B7.TIF');
                    if ~isfile(filepath_bt)
                        b6 = zeros(size(b1));
                    else
                        b6 = subFunLoadSingleBand(folderpath_tmp, imf, 'BT_B6.TIF');
                    end
                    
                case {8, 9} % Landsat 8 and 9
                    b1 = subFunLoadSingleBand(folderpath_tmp, imf, 'SR_B2.TIF');
                    b2 = subFunLoadSingleBand(folderpath_tmp, imf, 'SR_B3.TIF');
                    b3 = subFunLoadSingleBand(folderpath_tmp, imf, 'SR_B4.TIF');
                    b4 = subFunLoadSingleBand(folderpath_tmp, imf, 'SR_B5.TIF');
                    b5 = subFunLoadSingleBand(folderpath_tmp, imf, 'SR_B6.TIF');
                    b7 = subFunLoadSingleBand(folderpath_tmp, imf, 'SR_B7.TIF');
                    if ~isfile(filepath_bt)
                        b6 = zeros(size(b1));
                    else
                        b6 = subFunLoadSingleBand(folderpath_tmp, imf, 'BT_B10.TIF');
                    end
            end
            
            b1 = 10000.*(b1.*0.0000275 - 0.2); % rescale to 0 ~10000
            b2 = 10000.*(b2.*0.0000275 - 0.2); % rescale to 0 ~10000
            b3 = 10000.*(b3.*0.0000275 - 0.2); % rescale to 0 ~10000
            b4 = 10000.*(b4.*0.0000275 - 0.2); % rescale to 0 ~10000
            b5 = 10000.*(b5.*0.0000275 - 0.2); % rescale to 0 ~10000
            b7 = 10000.*(b7.*0.0000275 - 0.2); % rescale to 0 ~10000
            b6 =    10.*(b6.*0.00341802 + 149); % rescale same as Landsat Collection 1 ARD (0.1)
        elseif Collection==1
            switch l_num
                case {4, 5, 7} % Landsat 4, 5, and 7
                    b1 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB1.tif');
                    b2 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB2.tif');
                    b3 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB3.tif');
                    b4 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB4.tif');
                    b5 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB5.tif');
                    b7 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB7.tif');
                    b6 = subFunLoadSingleBand(folderpath_tmp, imf, 'BTB6.tif');
                case 8 % Landsat 8
                    b1 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB2.tif');
                    b2 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB3.tif');
                    b3 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB4.tif');
                    b4 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB5.tif');
                    b5 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB6.tif');
                    b7 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB7.tif');
                    b6 = subFunLoadSingleBand(folderpath_tmp, imf, 'BTB10.tif');
            end
        end
%         refImage = GRIDobj(fullfile(folderpath_tmp, [imf, '_SR_B4.TIF']));
%         refImage.Z = int16(b4);
%         GRIDobj2geotiff(refImage,fullfile(folderpath_tmp,['stack_',imf,'_SR_B4.tif']));
%         refImage.Z = int8(cfmask);
%         GRIDobj2geotiff(refImage,fullfile(folderpath_tmp,['stack_',imf,'_QA.tif']));
        
        %% Load tide band into the stacked file
        imageDate = datetime(str2double(imf(16:23)),'ConvertFrom','yyyymmdd');
        filepath_wl = fullfile(folderpath_wl,sprintf('wlDoy%04d%03d.tif',year(imageDate),day(imageDate,'dayofyear')));
        %%---- Temporal to delete  
        %%-----
        if ~isfile(filepath_wl)
            warning('No water level band for %s', imf);
            return;
        end
        [waterLevel,R] = readgeoraster(filepath_wl);
        % Resize the water levels from 1000 m resolution to 30 m
        [waterLevel,~] = mapresize(waterLevel,R,1000/30,'nearest'); % resize to the calculation scale
        waterLevel = int16(waterLevel);

        %% Step 4. Stack all the bands together row by row into the temp foler of current task
%         [nrows, ncols] = size(b1); % have the dimension from b1
        for irow = 1: npixelper: numPixels % total of 5000 rows
            % @ the current task, only parts of the images will be reconstructed
            irow_end = min(irow + npixelper -1, numPixels);
            nPixelInBlock = irow_end-irow+1;
            % row data's foldername R0000100010, which means the row 00001 to 00010
            foldername_rowdata = subFunBlockFolderName(irow, irow_end);
            folderpath_out = fullfile(folderpath_stack, foldername_rowdata);
            if ~exist(folderpath_out,'dir')
                mkdir(folderpath_out);
            end
            
            % @ the 1st image at 1st task, to copt the metadata to each sub
            % row folder, that will be loaded for proceseesing further
            if task == 1 && i == 1
                copyfile(fullfile(fileparts(folderpath_stack), sprintf('metadata.mat')),  fullfile(folderpath_stack, foldername_rowdata)); % saveas metadata
                % single path layer generator
                if singlepath
                    copyfile(fullfile(fileparts(folderpath_stack), 'singlepath_landsat.tif'), fullfile(folderpath_stack, foldername_rowdata));
                end
            end
            
            % if exsit, skip
%             if isfile(fullfile(folderpath_out, stackname)) % ignore to check the file broken or not, since the row file is too small to write incorrectly (at least low possibly)
%                 continue;
%             end
            
            % as for a certain row data, the data will be in format of BIP
            % (data format)    [b1 b2 b3 b4 b5 b7 b6 qa b1 b2 b3 b4 b5 b7 b6 qa b1 b2 b3 b4 b5 b7 b6 qa .... b1 b2 b3 b4 b5 b7 b6 qa]
            % (index in array) [1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 .....] 
            blockdata = zeros(nPixelInBlock, 2+nbands, 'uint16'); % row, col, nrowsper * 5000 cols * 8 bands + WL band info, orbit, sensor
            locInBlock = maskPixelsInfo(irow:irow_end,4);
            blockdata(:,end-1:end) = maskPixelsInfo(irow:irow_end,2:3);
            blockdata(:,1) = b1(locInBlock); % blue
            blockdata(:,2) = b2(locInBlock); % green
            blockdata(:,3) = b3(locInBlock); % red
            blockdata(:,4) = b4(locInBlock); % nir
            blockdata(:,5) = b5(locInBlock); % swir1
            blockdata(:,6) = b7(locInBlock); % swir2
            blockdata(:,7) = b6(locInBlock); % thermal band!
            blockdata(:,8) = waterLevel(locInBlock); % tide water level band
            blockdata(:,9) = cfmask(locInBlock); % qa
            blockdata = array2table(blockdata,'VariableNames',{'b1','b2','b3','b4','b5','b7','b6','tide','qa','row','col'});
            save(fullfile(folderpath_out, stackname),'blockdata');
        end
        
        %% Step 5. Clear temp image folder once done
%         rmdir(folderpath_tmp, 's');
        
        if msg
            fprintf('Finished %0.2f percent at task# %d/%d with total %0.2f mins\n', ...
                100*(1 + i - start_i)/(end_i - start_i + 1), task, ntasks, toc/60);
        end
    end % end of processing all the images at the current task
    fclose('all'); % close the IO for reading .xml files
    rmdir(folderpath_task, 's'); % clear temp task folder once done
end

%% sub function of creating stack data name
function stackname = subFunStackDataName(imf)
    % i.e., LT50290051984251C1V01_012
    % LT5 means Landsat 5 TM
    % 029005 means the Landsat ARD tile h029v005
    % 1984 means the year of 1984
    % 251 means the DOY of 251
    % C1 means Landsat Collection 1
    % V01 means the Landsat ARD version V01
    % PPP means the path of the Landsat orbit.
    yr = str2double(imf(16:19));
    mm = str2double(imf(20:21));
    dd = str2double(imf(22:23));
    doy = datenummx(yr,mm,dd)-datenummx(yr,1,0);
      
    stackname = [imf([1, 2, 4,9:14]),num2str(yr,'%04d'), num2str(doy,'%03d'), 'C', imf(34:35)]; 
    
end

%% sub function of creating row folder name
function foldername_rowdata = subFunBlockFolderName(irow, irow_end)
    % row data's foldername R0000100010, which means the row 00001 to 00010
    foldername_rowdata = sprintf('P%05dK%05dK', floor(irow/1000), ceil(irow_end/1000));
end

%% sub function of loading a certain band from the landsat data folder
function surf_b = subFunLoadSingleBand(imgfoler, imgname, specifyname)
    surf_b = double(readgeoraster(fullfile(imgfoler, [imgname,'_' ,specifyname])));
end