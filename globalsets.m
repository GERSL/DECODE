classdef globalsets
    properties (Constant)
        %% Sets for DECODE
        %% Paths
        % [Required] Directory of original Landsat ARD with .tar, including surface reflectance and brightness temperature
        PathLandsatARD = '/shared/cn449/DataLandsatARDCONUS/';
        % Directory of DECODE workspace
        PathDECODE = '/scratch/xiy19029/CONUS/DECODEResults/';
        % Path to intepolated water level (not need to be cliped to Tiles)
        PathTide = '/shared/cn451/Yang/Tide/';
        % Path to mask image of the tidal wetland (scaled to Landsat ARD)
        PathMask = '/scratch/xiy19029/layers/Mask/';
        %% process region
        processedRegion = 'northeast'; % 'all'; 'southeast'; 'northeast';'west'
        %% Input File
        TableARDTiles = 'CCAPGrid.xlsx'; %ARD tiles in different locations of CONUS
        
        %% Folders to save the results
        FolderDetection = 'TSFitLine';
        FolderTSFit = 'TSFitMap';
        FolderMap = 'outMaps';
        FolderMosaic = 'OutMapCONUS'; 
        FolderSample = 'samples';

        %% Filenames
        FilenameRecordRemainRows = 'RecordRemainRows';
        
        %% Parameters
        Years = 1986:2020;
        LandsatCollectionVersion = '01';
        %% End of sets for CCD ###################
    end
    methods (Static)
        function id_status = checkRequirementsToCCD(ARDTileName)
            id_status = 0;
            % Requirements:
            % 1) stack data ready?
            % 2) extent.tif and nooverlap.tif ready?
            if 5000 == length(dir(fullfile(globalsets.PathDECODE, ARDTileName, globalsets.FolderTSFit, 'record_change*.mat')))
                id_return = 9;
                return;
            end
            % check STACK ok or not
            if length(dir(fullfile(globalsets.PathLandsatARD, hv_name, 'L*_SR.tar'))) == ...
                  length(dir(fullfile(globalsets.PathDECODE, hv_name, 'L*')))
    %                                 fprintf('Skip: Finished stacking for ARD Tile %s\n', hv_name);
                idsmore = [idsmore; iARD];
            end
        end
        function [ARDTiles, ARDTilesDone] = getARDTiles(checkStatus)
%             % ARDTiles = textread(globalsets.PathARDTiles, '%s');
            ARDTiles = readtable(fullfile(globalsets.TableARDTiles));
           
            if strcmp(globalsets.processedRegion,'all')
                ARDTiles = ARDTiles.tile;
            else
                % index rows with the correct variables in column
                ARDTiles.region = string(ARDTiles.region);
                regions = unique(ARDTiles.region);

                idx = ARDTiles.region == globalsets.processedRegion;
                ARDTiles = ARDTiles.tile(idx,:);
            end
            if exist('checkStatus', 'var')
                idsmore = [];
                switch lower(checkStatus)
                    case 'stack'
                        for iARD = 1: length(ARDTiles)
                            hv_name = ARDTiles{iARD};
                            % check DECODE ok or not
                            if 5000 == length(dir(fullfile(globalsets.PathDECODE, hv_name, globalsets.FolderTSFit, 'record_change*.mat')))
                                fprintf('Skip: Finished CCD for ARD Tile %s\n', hv_name);
                                idsmore = [idsmore; iARD];
                                continue;
                            end
                            % check STACK ok or not
                            if length(dir(fullfile(globalsets.PathLandsatARD, hv_name, 'L*_SR.tar'))) == ...
                                  length(dir(fullfile(globalsets.PathDECODE, hv_name, 'L*')))
                                fprintf('Skip: Finished stacking for ARD Tile %s\n', hv_name);
                                idsmore = [idsmore; iARD];
                                continue;
                            end
                            % check in the mask region or not
                            if ~isfile(fullfile(globalsets.PathMask,globalsets.processedRegion,[hv_name,'_mask.tif']))
                                idsmore = [idsmore; iARD];
                                fprintf('Skip: No potential region for ARD Tile %s\n', hv_name);
                                continue
                            end
                        end
                    case 'decode'
                        % check DECODE ok or not
                        for iARD = 1: length(ARDTiles)
                            hv_name = ARDTiles{iARD};
                            if 5000 == length(dir(fullfile(globalsets.PathDECODE, hv_name, globalsets.FolderTSFit, 'record_change*.mat')))
                                 fprintf('Skip: Finished CCD for ARD Tile %s\n', hv_name);
                                idsmore = [idsmore; iARD];
                                continue;
                            end
                             % check in the mask region or not
                            if ~isfile(fullfile(globalsets.PathMask,[hv_name,'_mask.tif']))
                                fprintf('Skip: No potential region for ARD Tile %s\n', hv_name);
                                idsmore = [idsmore; iARD];
                                continue
                            end
                        end
                       
                end
                % eleliminate the tile name we have finished already
                if ~isempty(idsmore)
                    ARDTilesDone = ARDTiles(idsmore);
                    ARDTiles(idsmore) = [];
                else
                    ARDTilesDone = [];
                end
            end
        end
    end
end

