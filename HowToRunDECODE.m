%% This page demonstrates the examples of running the Matlab functions in DECODE version 1.1
% If any questions, please feel free to contact Xiucheng Yang (xiucheng.yang@uconn.edu) in GERS lab (https://gerslab.uconn.edu/).

%% We provide a case study area (15,731 Landsat pixels), in Florida, to run the DECODE approach. 
% The study area experienced gradual shift from tidal marsh to mangrove
% The input data 

%% Global setting of the path and parameters
globalsets.m
%% IMPORTANT!!! 
%% Before to run!!!
% (1) UNZIP Characterization/RandomForest_CoverSample_East_RF.zip in <Characterization>
% (2) To modify the two path variables in globalsets:
% PathMask = '/scratch/xiy19029/CaseStudyOfUSTidalWetland/layers/Mask/';
% % Directory of DECODE workspace
% PathDECODE = '/scratch/xiy19029/CaseStudyOfUSTidalWetland/DECODEResults/';    


%% After the modification of these paths, the proposed DECODE approach can automatically run:
Main_AutoRun_DECODE.m 


%% ------------------------------------------DETAILS-------------------------------------------------------------------------
%% Input data including all available Landsat data and NOAA guage station measurement
% Landsat data can be downloaded from https://earthexplorer.usgs.gov/ by using https://m2m.cr.usgs.gov/ batch downloading
% NOAA gauge station prediction can be downloaded from:
% https://tidesandcurrents.noaa.gov/ by using https://www.tidesandcurrents.noaa.gov/api-helper/url-generator.html
% The sample-based measures are interpolated to raster images 
% --> See Subfolder <prepareWL> for the related codes
%% -------------------------- Step 0.Stack Landsat ARD by blocks ----------------------------------------
% Proposed time series analysis approach is pixel-based approach and needs
% to load all available historic Landsat observations for each pixel to
% complete the calculation. It is difficult to use the original Landsat ARD
% images (5000*5000 pixels) to complete the computing. To enable the
% computing, a block processing is necessary to firstly divide the whole
% pixels into blocks (such as 1000 pixels per block).

% Here, as a case study, we provide a case study area in Florida with stacked
% Landsat data and tide predictions.
batchDECODE_Phase0_Stack.m         % Folder --> <Stack>

%% -------------------------- Step 1. Spectral break detection using DECODE ---------------------------------
batchDECODE_Phase1_DETECT.m        % Folder --> <Detection>
% For the case study area, it spends ~ 1 hour to complete in a local
% computer without parallel compusing; If parfor is active to use, we
% highly recommended it 

%% -------------------------- Step 2. cover classification: Tidal marsh, Mangrove, Dieback, Tidal flats, Open water, and Others -------------
% ---------------------------------- Characterize change - condition/cover changes without/with categorical class conversion
batchDECODE_Phase2_Characterize.m  % Folder --> <Characterization>
% For the case study area, it spends ~ 3 minutes to complete in a local computer without parallel compusing

%% -------------------------- Step 3. Export the pixel-level map && Post-processing to patch-level --------------------------
batchDECODE_Phase3_Mapping.m       % Folder --> <Export>

% - Postprocess - Generate the patch-based map by MMU and contextual algorithm
% The result might be slightly different from the interactive map, because
% the contextual algrithm depends on the surrounding open water bodies. The
% case study area only involves a small region compared to US-wide analysis.
batchDECODE_Phase3_S1_IsolatedRemoval.m % For US-wide
% For single-tile based object process, like the case study attached here,
% Please use the alternative function:
batchDECODE_Phase3_S1_IsolatedRemoval_SingleTile.m

% - Postprocess - Cover change map
batchDECODE_Phase3_S2_CoverChangeMap.m

%% -------------------------- submit parallel jobs from Matlab program --------------------------
% see HPC folder for batch parallel job submission
% The hpc slurm file is prepared for UCONN HPC platform
% If using HPC platform, submit the jobs one by one with automatic processing
% Suggesting 100-200 cores per tile to calculate
Phase0_Stack.sh 
Phase1_DETECT.sh
Phase2_Characterize.sh
Phase3_Mapping.sh
Phase3_S1_GeneratePatch.sh
Phase3_S2_CoverChange.sh

%% -------------------------- Main figures with analysis --------------------------
% Folder <FigureGenerateForCoverAnalysis> --> Analysis based on the results in the US wide.


%% --------------------------------   Input samples ----------------------------------------------
% Folder <Samples> --> calibration change sample; classification training data


%% ------------------------------- Additional supportive functions -------------------------------
% Folder <Packages> --> GRIDobj; RandomForest; ktaub; accuareaadjust(Good Practice);
% Folder <LandsatPathData> --> Landsat Single-path file to ensure the consistent density of Landsat observations in spatial domain



