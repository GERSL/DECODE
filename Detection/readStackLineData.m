function [sdate, line_t, path_t]= readStackLineData(folderpath_row, ncols, nbands, rows, path_r)
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

%     folderpath = '/lustre/scratch/qiu25856/COLDResults/TestBRDF_DensityAdjust_StackLine/h014v009/';
%     nrowsper = 10;
%     nrows = 5000;
%     ncols = 5000;
%     nbands = 8;
%     row = 4037;
    
%     folderpath = '/lustre/scratch/qiu25856/TestGERSToolbox/h029v005/StackData/R0329103300/';
%     nbands = 8;
%     ncols = 5000;
%     row = 1;
%     singlepath = [];

    %% get the real row# at start and end
    [~, foldername_stackrows] = fileparts(folderpath_row); % R xxxxx xxxxx
    fprintf('rowstart:  %s\n',foldername_stackrows);
    row_start = str2num(foldername_stackrows(2:6));
    row_end = str2num(foldername_stackrows(7:11));

    %% get num of total files start with "L"
    imf = dir(fullfile(folderpath_row, 'L*')); % folder names
    % filter for Landsat folders
    imf = regexpi({imf.name}, 'L(\w*)_(\w*)', 'match'); %Lxxx_Path
    imf = [imf{:}];
    imf = vertcat(imf{:});
    
    %% Filter out the images out of single path layer per single row or current row cluster
    if ~isempty(path_r)
        % only remain the single path images to load
        WRSPaths = uint8(str2num(imf(:,23:25))); % record the path of Landsat path (1:233)
        % Exclude the images that do not at the current path
        path_t = unique(path_r);
        path_t(path_t<1) = []; % 1 to 233 for Landsat
        % Update the obj
        ids_more = ~ismember(WRSPaths, path_t);
        imf(ids_more,:) = [];
    else
        path_t = []; % return an empty var
    end
    
    %% sort according to yeardoy
    yeardoy = str2num(imf(:, 10:16)); 
    [~, sort_order] = sort(yeardoy);
    imf = imf(sort_order, :);
    % number of folders start with "L"
    num_t = size(imf,1);
    
    %% Find dates and paths
    yr = str2num(imf(:, 10:13)); %#ok<*ST2NM>
    doy = str2num(imf(:, 14:16));
    sdate = datenum(yr, 1, 0) + doy;
    if ~isempty(path_r)
        path_t = str2num(imf(:,23:25)); % record the path of Landsat swath
    end
    
    %% Read Landsat time series data
    line_t = zeros(num_t,nbands*ncols, length(rows), 'uint16'); %Ys
    % at each Landsat image
    for i = 1:num_t
        fid_t = fopen(fullfile(folderpath_row, imf(i, :)),'r'); % get file ids
%         fprintf('%s\n',fullfile(folderpath_row, imf(i, :)))
        %% Load all data once time
        if length(rows) == length(row_start: row_end)
            line_t(i,:,:) = fread(fid_t, [nbands*ncols, length(rows)],'int16','ieee-le'); % get Ys
            continue;
        end
        
        %% Load parts of the data
        irows = find(ismember(row_start: row_end, rows));
        for ir = 1: length(irows)
            % move forward the irows(ir) that is the real index of the row
            % dataset
            fseek(fid_t,2*(irows(ir)-1)*ncols*nbands,'bof'); % 2 means uint16
            % ir is the index of line_t
            line_t(i,:,ir) = fread(fid_t,nbands*ncols,'int16','ieee-le'); % get Ys
        end
    end
    fclose('all'); % close all files

    
end

