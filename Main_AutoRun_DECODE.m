%% See Readme.m for Details

%% The whole processing would need ~ 1 hour, owing to the dense time series analysis to detect 

%% !IMPORTANT!: Modify the two path variables ["PathMask" "pathDECODE"] in the globalset.m 
  
clc;
clear;
fprintf('+++++++++++++++++++++++ Begin to Run DECODE +++++++++++++++++++++++\n');

fprintf('\n\n-----------Phase 1 --> Dense time series change detection -----------------\n');
batchDECODE_Phase1_DETECT();
% For the case study area, it spends ~ 1 hour to complete in a local computer without parallel computing

fprintf('\n\n-----------Phase 2 --> Characterization of cover and change -----------------\n');
batchDECODE_Phase2_Characterize();  

fprintf('\n\n-----------Phase 3 --> Exporting annual maps -----------------\n');
batchDECODE_Phase3_Mapping();  

fprintf('\n\n-----------Phase 3-S1 --> Object-based postprocessing -----------------\n');
batchDECODE_Phase3_S1_IsolatedRemoval_SingleTile();

fprintf('\n\n-----------Phase 3-S2 --> Export cover change maps -----------------\n');
batchDECODE_Phase3_S2_CoverChangeMap();

pathworking = globalsets.PathDECODE;
fprintf('--------- Completed the DECODE process!\n--------- Check results -> %s\n',pathworking);