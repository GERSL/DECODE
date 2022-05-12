%% This page demonstrates the examples of running the Matlab functions in DECODE version 1.1
% If any questions, please feel free to contact Xiucheng Yang (xiucheng.yang@uconn.edu) in GERS lab (https://gerslab.uconn.edu/).

%% Step 0. Prepare tide data based on NOAA CO-OPS or other tide information
% see folder prepareWL to predict the daily water level informations 

%% Step 1.Stack Landsat ARD by lines
% see batchStackLandatARDTide2Line.m
% see 1stack.sh at <HPC> for job submission in HPC platform for parallel
% computing

%% Step 2. Spectral break detection using DECODE
% see batchDECODE.m
% see 2BreakDetection.sh at <HPC> for parallel

%% Step 3. cover classification: vegetated wetlands, open water and others for temporal segments
%%         And characterize change - condition/cover changes without/with cover change
% see batchCoverChangeMaps.m
% see 3coverChangeType.sh at <HPC> for parallel

%% Step 4. Export the change and cover maps
% see batchExportMaps.m
% see 4export.sh at <HPC> for parallel

%% Step 5. Postprocess - seive the isolated pixels in the cover classificaiton map
% see batchPostProcessCoverImage.m
% see 5seive.sh at <HPC> for parallel

%% submit parallel jobs from Matlab program
% see HPC folder for batch parallel job submission
% The hpc slurm file is prepared for UCONN HPC platform