function line_t = readStackBlockData(folderpath_stackrows,ncols,beginYear)
%READSTACKLINEDATA This function is to read the stack Landsat ARD time
%series. ~ 10 seconds to loading all 1500 images without I/O waiting.
%
% INPUT:
%
%   folderpath_row:         Locate to the stacking row dataset folder.
%
%   ncols:                  Number of columns. If Landsat ARD, ncols =
%                           5000.
%
%   nbands:                 Number of bands in the stacking data. nbands =
%                           8 usually, 7 Landsat spectral bands + 1 QA band
%
%   rows:                   [Array] Index of rows. i.e., [1], [2], or [1, 2] [1: nrows]
%
%   path_r (optional):      Single Landsat path labels per row. If set as empty
%                           ([]), all Landsat data including overlap will be used.
%
%
% RETURN:
%
%   sdate:                  [1-dimension array], series dates of all images
%
%   line_t:                 [3-dimension array], 1st dimension reponses the
%                           sdate, 2nd dimension indicates the landsat data
%                           with BIP format, and 3rd dimension depends on
%                           the inputing rows.
%
%   path_t:                 [1-dimension array], the single path value via
%                           row. of which size is same as ncols.
%   
%
% EXAMPLE OF USE:
%
%   See COLD.m
% 
% AUTHOR(s): Shi Qiu
% DATE: Feb. 5, 2021
% COPYRIGHT @ GERSLab

%     folderpath_stackrows = '/scratch/xiy19029/FloridaMangrove/DECODERResults/h026v018/StackData/P00000K00010K/';
%     ncols = 5000;
    if ~exist('beginYear','var')
        beginYear = 1984;
    end
       
    %% get num of total files start with "L"
    imf = dir(fullfile(folderpath_stackrows, 'L*')); % folder names
    % filter for Landsat folders
    imf = regexpi({imf.name}, 'L(\w*)_(\w*)', 'match'); %Lxxx_Path
    imf = [imf{:}];
    imf = vertcat(imf{:});
    
    %% sort according to yeardoy
    yeardoy = str2num(imf(:, 10:16)); 
    [~, sort_order] = sort(yeardoy);
    imf = imf(sort_order, :);

    %% Filter the images with target years
    yr = str2num(imf(:, 10:13)); %#ok<*ST2NM>
    idxTarget = yr>=beginYear;
    imf = imf(idxTarget,:);
    % number of folders start with "L"
    num_t = size(imf,1);

    % Load the first record
    load(fullfile(folderpath_stackrows,imf(1,:)));
    num_p = height(blockdata);
    n_fields = width(blockdata);
    fieldName = blockdata.Properties.VariableNames;

    %% Find dates and paths
    yr = str2num(imf(:, 10:13)); %#ok<*ST2NM>
    doy = str2num(imf(:, 14:16));
    sdate = datenum(yr, 1, 0) + doy;
    path_t = str2num(imf(:,21:23)); % record the path of Landsat swath %C2
    
    %% Read Landsat time series data
    line_t = double(zeros(num_p*num_t,n_fields+3)); 
    % at each Landsat image
    for i = 1:num_t
        load(fullfile(folderpath_stackrows,imf(i,:)));
        blockdata = double(table2array(blockdata));
        poss = (blockdata(:,end-1)-1)*ncols + blockdata(:,end);
        blockdata = [poss,ones(num_p,1)*sdate(i),blockdata,ones(num_p,1)*path_t(i)];
        line_t((i-1)*num_p+1:i*num_p,:) = blockdata;       
    end
    line_t = sortrows(line_t,[1,2]);
    line_t = array2table(line_t,'VariableNames',['pos','sdate',fieldName,'path']); 
    fclose('all'); % close all files
end

