function stackLandsatARDTide2Line(folderpath_ard, folderpath_stack, folderpath_wl, varargin)
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
% EXAMPLE OF USE:
%
%   > To stack Landsat ARD at task # 1/20
%
%   stackLandsatARD2Line('/lustre/scratch/qiu25856/DataLandsatARD/CONUS/h029v005', '/lustre/scratch/qiu25856/TestGERSToolbox/h029v005/StackData', 'task', 1 ,'ntasks', 20)
%
%   > To stack Landsat ARD at task # 1/20 one more time. then, the check as
%   true is recommanded.
%
%   stackLandsatARD2Line('/lustre/scratch/qiu25856/DataLandsatARD/CONUS/h029v005', '/lustre/scratch/qiu25856/TestGERSToolbox/h029v005/StackData', 'task', 1 ,'ntasks', 20 , 'check', true)
%
% 
% AUTHOR(s): Shi Qiu
% DATE: Feb. 5, 2021
% COPYRIGHT @ GERSLab

    % add the matlab search path of GRIDobj
    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'GRIDobj'));

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
    addParameter(p,'nsubrow', 10); % default rows per file (optional)
    addParameter(p,'orbitpath', 'single'); % to create single path layer
    addParameter(p,'task', 1); % 1st task
    addParameter(p,'ntasks', 1); % single task to compute
    addParameter(p,'msg', true); % not to display info
    addParameter(p,'check', false); % not to scan history data; but once more time, please set as true, and if exist, just skip
    % request user's input
    parse(p,varargin{:});
    
    clr_pct_min = p.Results.clear;
    nrowsper = p.Results.nsubrow;
    task = p.Results.task;
    ntasks = p.Results.ntasks;
    msg = p.Results.msg;
    checkexist = p.Results.check;
    % char will be better understood for users
    switch lower(p.Results.orbitpath)
        case 'single'
            singlepath = true;
        case 'all'
            singlepath = false;
    end

    
    %% constant parameter    
    nbands = 9; % 7 spectral bands + wl  + Fmask QA band
    
    %% Filter for Landsat folders
    % get num of total folders start with "L"
    imfs = dir(fullfile(folderpath_ard,'L*SR.tar'));
    % filter for Landsat folders
    % espa data
    %imf = regexpi({imf.name}, 'L(T05|T04|E07|C08)(\w*)\-(\w*)', 'match'); 
    % Landsat ARD
    imfs = regexpi({imfs.name}, 'L(T05|T04|E07|C08)(\w*)', 'match'); % no expand name
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
    imfs = imfs(sortedYear>19840000,:);
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
    tic % record all running time
    for i = start_i:end_i
        %% locate to a certain image
        imf = imfs(i, :);
        % new filename in format of LXSHHHVVVYYYYDOYLLTTT_PPP
        stackname = subFunStackDataName(imf);
        
        %% check exist or not by DIR the last row folder
        if checkexist
            filepath_meta = fullfile(fileparts(folderpath_stack), 'metadata.mat');
            if isfile(filepath_meta) % we must have the metadata first
                if ~exist('metadata', 'var')
                    load(filepath_meta); % only load once
                end
%                 %% normal checking
%                 irow = max(1: metadata.nsubrows: metadata.nrows); % the last row folder
%                 irow_end = min(irow + metadata.nsubrows -1, metadata.nrows);
%                 foldername_rowdata = subFunRowFolderName(irow, irow_end);
%                 rowdata_exist = dir(fullfile(folderpath_stack, foldername_rowdata, [stackname, '_*']));
%                 isexist = false;
%                 for iexist = 1: length(rowdata_exist)
%                     if rowdata_exist(iexist).bytes == metadata.nsubrows*metadata.nrows*metadata.nbands*2 % bit size at uint 16
%                         isexist = true;
%                     end
%                 end
                %% rigorous checking
                isexist = true;
                for irow = 1: metadata.nsubrows: metadata.nrows % total of 5000 rows
                    % @ the current task, only parts of the images will be reconstructed
                    irow_end = min(irow +  metadata.nsubrows -1,  metadata.nrows);

                    % row data's foldername R0000100010, which means the row 00001 to 00010
                    foldername_rowdata = subFunRowFolderName(irow, irow_end);
                    rowdata_exist = dir(fullfile(folderpath_stack, foldername_rowdata, [stackname, '_*']));
                    if ~isempty(rowdata_exist)
                         if rowdata_exist(1).bytes ~= metadata.nsubrows*metadata.nrows*metadata.nbands*2 % bit size at uint 16
                             isexist = false;
                             delete(fullfile(folderpath_stack, foldername_rowdata, rowdata_exist.name))
                             fprintf('Need to be updated for %s\n',rowdata_exist.name);
                             break
                         end
                    else
                        isexist = false;
                        fprintf('Lack %s\n',rowdata_exist.name);
                        break
                    end
                end
                %% if check ok then skip
                if isexist
                    if msg
                        fprintf('\nAlready exist the %dth image (%s)\n', i, imf);
                    end
                    continue;
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
            untar(filepath_bt, folderpath_tmp);
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
            cfmask0 = GRIDobj(fullfile(folderpath_tmp, [imf, '_PIXELQA.tif']));
            metadata = [];
            metadata.GRIDobj = cfmask0;
            metadata.GRIDobj.Z = []; % set as [] for saving storage
            metadata.GRIDobj.name = []; % set as [] for saving storage
            metadata.region = imf(6:7); % CU: CONUS
            metadata.tile = sprintf('h%sv%s', imf(9:11), imf(12:14)); % tile name
            metadata.nrows = cfmask0.size(1); % record of the row size
            metadata.ncols = cfmask0.size(2); % record of the column size.
            metadata.nbands = nbands; % record of the number of bands in stack data
            metadata.nsubrows = nrowsper; % record of the number of rows per file
            metadata.nimages = num_t; % record of the total number of Landsat images
            metadata.cloudcover = 100 - clr_pct_min; % record of max percentage of cloud cover
            metadata.createtime = datestr(now); % record of the create time
            
            save(fullfile(fileparts(folderpath_stack), 'metadata'), 'metadata'); % saveas metadata
            
            % single path layer generatorni h
            if singlepath
                singlepathlayer = GRIDobj(fullfile( fileparts(fileparts(mfilename('fullpath'))), 'ToolData', 'singlepath_landsat_conus.tif'));  % at parent code's folder
                path_ard = reproject2utm(singlepathlayer, cfmask0,'method', 'nearest');
                path_ard.Z = uint8(path_ard.Z);
                GRIDobj2geotiff(path_ard, fullfile(fileparts(folderpath_stack), 'singlepath_landsat.tif'));
                clear path_ard;
            end
            % convert normal image array with uint16
            cfmask0 = uint16(cfmask0.Z); % only array and uint16 same as before
        else
            % normal process
            cfmask0 = geotiffread(fullfile(folderpath_tmp, [imf, '_PIXELQA.tif']));
        end
        cfmask = cfmask0;
        % convert pixel QA to fmask values
        cfmask(bitget(cfmask0,1) == 1) = 255;
        cfmask(bitget(cfmask0,2) == 1) = 0;
        cfmask(bitget(cfmask0,3) == 1) = 1;
        cfmask(bitget(cfmask0,4) == 1) = 2;
        cfmask(bitget(cfmask0,5) == 1) = 3;
        cfmask(bitget(cfmask0,6) == 1) = 4;
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
        path_char=char(xml_str(strmatch('path=',xml_str))); % path=
        path = str2double(path_char(1,7:end-1)); % only 1 item
        stackname = [stackname, sprintf('_%03d', path)]; % landsat path range 1 to 233
        
        
        %% Load all bands
        l_num = str2double(stackname(3)); % which Landsat satellite?
        switch l_num
            case {4, 5, 7} % Landsat 4, 5, and 7
                b1 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB1');
                b2 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB2');
                b3 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB3');
                b4 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB4');
                b5 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB5');
                b7 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB7');
                b6 = subFunLoadSingleBand(folderpath_tmp, imf, 'BTB6');
            case 8 % Landsat 8
                b1 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB2');
                b2 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB3');
                b3 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB4');
                b4 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB5');
                b5 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB6');
                b7 = subFunLoadSingleBand(folderpath_tmp, imf, 'SRB7');
                b6 = subFunLoadSingleBand(folderpath_tmp, imf, 'BTB10');
        end
        
        %% Load tide band into the stacked file
        filepath_wl = dir(fullfile(folderpath_wl,['tide',imf(16:23),'.tif']));
        filepath_wl = fullfile(folderpath_wl,filepath_wl.name);
        if ~isfile(filepath_wl)
            warning('No water level band for %s', imf);
            return;
        end
        wl = geotiffread(filepath_wl);
        
        %% Step 4. Stack all the bands together row by row into the temp foler of current task
        [nrows, ncols] = size(b1); % have the dimension from b1
        for irow = 1: nrowsper: nrows % total of 5000 rows
            % @ the current task, only parts of the images will be reconstructed
            irow_end = min(irow + nrowsper -1, nrows);
 
            % row data's foldername R0000100010, which means the row 00001 to 00010
            foldername_rowdata = subFunRowFolderName(irow, irow_end);
            folderpath_out = fullfile(folderpath_stack, foldername_rowdata);
            if ~isfolder(folderpath_out)
                mkdir(folderpath_out);
            end
            
            % @ the 1st image at 1st task, to copt the metadata to each sub
            % row folder, that will be loaded for proceseesing further
            if task == 1 && i == 1
                copyfile(fullfile(fileparts(folderpath_stack), 'metadata.mat'),  fullfile(folderpath_stack, foldername_rowdata)); % saveas metadata
                % single path layer generator
                if singlepath
                    copyfile(fullfile(fileparts(folderpath_stack), 'singlepath_landsat.tif'), fullfile(folderpath_stack, foldername_rowdata));
                end
            end
            
            % if exsit, skip
            if isfile(fullfile(folderpath_out, stackname)) % ignore to check the file broken or not, since the row file is too small to write incorrectly (at least low possibly)
                continue;
            end
            
            % as for a certain row data, the data will be in format of BIP
            % (data format)    [b1 b2 b3 b4 b5 b7 b6 qa b1 b2 b3 b4 b5 b7 b6 qa b1 b2 b3 b4 b5 b7 b6 qa .... b1 b2 b3 b4 b5 b7 b6 qa]
            % (index in array) [1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 .....] 
            linedata = zeros(length(irow: irow_end), ncols, nbands, 'int16'); % nrowsper * 5000 cols * 8 bands
            linedata(:,:,1) = b1(irow: irow_end, :); % blue
            linedata(:,:,2) = b2(irow: irow_end, :); % green
            linedata(:,:,3) = b3(irow: irow_end, :); % red
            linedata(:,:,4) = b4(irow: irow_end, :); % nir
            linedata(:,:,5) = b5(irow: irow_end, :); % swir1
            linedata(:,:,6) = b7(irow: irow_end, :); % swir2
            linedata(:,:,7) = b6(irow: irow_end, :); % thermal band!
            linedata(:,:,8) = wl(irow: irow_end, :); % tide water level band
            linedata(:,:,end) = cfmask(irow: irow_end, :); % qa
            
%             % Check if there is no data available
%             isNanData = linedata(:,:,end);
%             isNanData(isNanData<10) = 1;
%             isNanData(isNanData>10) = 0;
%             if ~sum(isNanData(:))
%                 continue
%             end
            multibandwrite(linedata, fullfile(folderpath_out, stackname), 'bip');
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
    yr = str2num(imf(16:19));
    mm = str2num(imf(20:21));
    dd = str2num(imf(22:23));
    doy = datenummx(yr,mm,dd)-datenummx(yr,1,0);
    stackname = [imf([1,2,4,9:14,16:19]),num2str(doy,'%03d'),imf([34,16,38:40])]; % Landsat Path will be given from .xml file later
end

%% sub function of creating row folder name
function foldername_rowdata = subFunRowFolderName(irow, irow_end)
    % row data's foldername R0000100010, which means the row 00001 to 00010
    foldername_rowdata = sprintf('R%05d%05d', irow, irow_end);
end

%% sub function of loading a certain band from the landsat data folder
function surf_b = subFunLoadSingleBand(imgfoler, imgname, specifyname)
    surf_b = geotiffread(fullfile(imgfoler, [imgname,'_' ,specifyname, '.tif']));
end