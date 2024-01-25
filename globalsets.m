classdef globalsets
    properties (Constant)
        %% Sets for DECODE
        %% Paths
        % -------------------------------Modify the path here-------------------------------
        % Path to mask image of the tidal wetland (scaled to Landsat ARD)
        PathMask = '/scratch/xiy19029/AssociatedData/layers/Mask/';
        % Directory of DECODE workspace
        PathDECODE = '/scratch/xiy19029/AssociatedData/DECODEResults/';    
        % -----------------------------------------------------------------------------------       

        % [Required] Directory of original Landsat ARD with .tar, including surface reflectance and brightness temperature
        PathLandsatARD = '/shared/cn449/DataLandsatARDCONUSTemp/';
        % Path to intepolated water level (not need to be cliped to Tiles)
        PathTide = '/shared/zhulab/Yang/Tide/';
           
        %% Input File
        TableARDTiles = 'CONUS_Coastal_Tiles_CaseStudy.csv'; %ARD tiles in different locations of CONUS     
%       IdxARDTiles = [13,20,45,61,72,74,89];

        %% window size to conduct the analyses 
        windowSize = 10; % pixel width per grid
        
        %% Folders to save the results
%         FolderDetection = 'TSFitLine_Single_composite';
        FolderDetection = 'TSFitBlock';
        FolderTSFit = 'TSFit';
        FolderMap = 'outMaps';
        FolderMosaic = 'OutMapCONUS';
        FolderSample = 'samples';
        folderCoverAnalysis = 'CoverAnalyses';
        %% Filenames
        FilenameRecordRemainRows = 'RecordRemainRows';
        
        %% Colors for different types
        colors = {'#1b7837','#7fbf7b','#cc79a7','#fc8d59','#4575b4','#404040'};

        %% Parameters
        Years = 1986:2021;
        LandsatCollectionVersion = '01';
        %% End of sets for DECODE ###################
    end
    methods (Static)
    %     function id_status = checkRequirementsToCCD(ARDTileName)
    %         id_status = 0;
    %         % Requirements:
    %         % 1) stack data ready?  
    %         % check STACK ok or not
    %         if length(dir(fullfile(globalsets.PathLandsatARD, hv_name, 'L*_SR.tar'))) == ...
    %               length(dir(fullfile(globalsets.PathDECODE, hv_name, 'L*')))
    % %                                 fprintf('Skip: Finished stacking for ARD Tile %s\n', hv_name);
    %             idsmore = [idsmore; iARD];
    %         end
    %     end
        function [ARDTiles, ARDTilesDone] = getARDTiles(checkStatus)
%             % ARDTiles = textread(globalsets.PathARDTiles, '%s');
            ARDTiles = readtable(fullfile(globalsets.TableARDTiles));
           
            
            ARDTiles = ARDTiles.tile;
           
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
                            if ~isfile(fullfile(globalsets.PathMask,[hv_name,'_mask.tif']))
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

