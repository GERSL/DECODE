function [rec_cg,clrx,clry] = TrendSeasonalFit_DECODE_Block(sdate_all, line_t, path_t_all, path_r, proc_cols, T_cg,Tmax_cg,conse,num_c,nbands,B_detect,varargin)
%% Spectral Break detection componenet of DEtection and Characterization of cOastal tiDal wEtland change (DECODE) 
% GERS Lab, University of Connecticut, Storrs
%
% INPUT:
%
% nrows:     should be the row ID, not that total number of rows
%% DECODE Version: $ Date: 20/07/2021 $ Copyright: GERS Lab, UCONN
% Version 1.1   Spectral break detection with five indices and additional
% water lever regressor; designed for coastal tidal wetlands
% COLD Revisions: $ Date: 1/1/2021 $ Copyright: GERS Lab, UCONN
%  Version 14.01  Parallel computing based on new stacking row data folders (2/6/2021)
% Version 14.00  COLD Version 2 for producing CONUS disturbance maps (1/1/2021) (Qiu and Zhu et al., 2021)
%                 Specially,
%                 1) Single Landsat swath data with minimum viewing zenith angle (do not consider overlap dara from adjacent orbits).
%                 2) Update frequency of the time series model based on 3% of the number of observations of previous time fit (withdraw version 13.04).
%                 3) Stop scanning files when loading Landsat images.
%  Version 13.06  Record during-change such slope and RMSE for the continous comparions (6/26/2020)
%  Version 13.05  Mat file with .mat.part was used to save locally and then rename it as .mat (6/1/2020)
%                 This can avoid the broken mat file caused by the dirsturbed saving process of mat file (e.g., shut down of compute)
%  Version 13.04  Update model for every three observations (09/03/2020)
%  Version 13.03  Modified model intitializatin test (10/29/2018)
%  Version 13.02  Add included angle to exlcude false positive change (10/05/2018)
%  Version 13.01  Do not need clear observations more than 25 percent (03/28/2018)
% Version 13.00  Optimized parameters for monitoring disturbance - COLD (03/20/2018) (Zhu et al., 2020, RSE)
%  Version 12.36  Reduce 75 percent of memory needed for line_t (11/22/2017)
%  Version 12.35  Adjust T_cg based on conse (10/03/2017)
%  Version 12.34  Add data density requirement (09/20/2017)
%  Version 12.33  Adjust conse based on delta time (09/19/2017)
%  Version 12.32  Fix bud for not enough data (09/17/2017)
%  Version 12.31  Read input varibles in txt file (07/05/2016)
%  Version 12.30  Fixed a bug for pixels without minimum observations (11/20/2015)
%  Version 12.29  Modified fit for perennial snow and Fmask failed pixels (09/20/2015)
%  Version 12.29  Do not fit disturbed time period (09/18/2015)
%  Version 12.28  Fixed a bug for missing values in land cover maps (09/16/2015)
%  Version 12.27  Fixed bugs for persistent snow and falied Fmask pixels (06/17/2015)
%  Version 12.26  Connected time for all models (05/28/2015)
%  Version 12.25  Bug fixed in snow percent (05/19/2015)
%  Version 12.24  Change T_const in Tmask (03/31/2015)
%  Version 12.23  Update iteratively before 24 observations (03/22/2015)
%  Version 12.22  Adjust mini RMSE based on temporal variability (03/22/2015)
%  Version 12.21  Add more categories and update i_start in the end (03/14/2015)
%  Version 12.20  Convert BT to from K to C before analysis (03/12/2015)
%  Version 12.19  Fit for permanent snow if is more than 75% (03/12/2015)
%  Version 12.18  No change detection if clear observation less than 25% (03/12/2015)
%  Version 12.17  Use median value for very simple model & change magnitude (02/24/2015)
%  Version 12.16  Finding changes in all water pixels (02/24/2015)
%  Version 12.15  Use the original multitemporal cloud mask (02/15/2015)
%  Version 12.14  Do not need example_img in images folder (02/09/2015)
%  Version 12.13: More infromation in "category" (11/10/2014)
%  This version (12.13) is used for the third round of the LCMS project.
%  Command: TrendSeasonalFit_v12Plot(N_row,N_col,min=0.5,T_cg=0.99,n_times=3,conse=6,B_detect=2:6)
%  Version 12.12: Fit for pixels where Fmask fails (11/09/2014)
%  Version 12.11: Bug fixed in num_fc (11/09/2014)
%  Version 12.10: Better multietmporal cloud detection at the beginning (11/06/2014)
%  Version 12.9:  Detect change only for land pixel (water/snow speical case) (10/31/2014)
%  Version 12.8:  Speed up by reducing time for RMSE and model computing (10/17/2014)
%  Version 12.7:  mini rmse should be larger than 10% of the mean (10/13/2014)
%  Version 12.6:  Fit model again when there are a 33.3% more data (10/08/2014)
%  Version 12.5:  Use subset of bands (2-6) for detecting surface change (10/01/2014)
%  Version 12.4:  Only apply multitemporal cloud masking during model initialization (09/29/2014)
%  Version 12.3:  Use subset of bands (3-5) to balance change in diferent dimensions (09/01/2014)
%  This version (12.3) is used for the second round of the LCMS project.
%  Command: TrendSeasonalFit_v12Plot(N_row,N_col,min=1,T_cg=0.99,n_times=3,conse=5,B_detect=3:6)
%  Version 12.2:  Bug fixed in model intialization (08/14/2014)
%  Version 12.1:  Use subset of bands (3-6) to avoid atmosphere influences (08/04/2014)
% Version 12.0   Detecting change based on probability (07/19/2014)
%  Version 11.6:  No need to change folder name & faster in speed (by Christ Holden 06/06/2014)
%  Version 11.5:  Improved calculation of temporally adjusted RMSE (04/23/2014)
%  Version 11.4:  Revise "rec_cg.category" to better seperate different fit processes (04/01/2014)
%  This version (11.4) is used for generating synthetic data for ACRE project and
%  detecting change for LCMS project.
%  Command: TrendSeasonalFit_v11Plot(N_row,N_col,min=1,T_cg=2,n_times=3,conse=6,B_detect=1:6)
%  Version 11.3:  Add "rec_cg.magnitude" as indicator of change magnitude (04/01/2014)
%  Version 11.2:  Change very simple fit with mean value for start and end of timeseries (04/01/2014)
%  Version 11.1:  Do not need metadata in the image folder to run CCDC (03/25/2014)
% Version 11.0:  Use change vector magnitude as threshold for detecting change (03/25/2014)
%  Version 10.13: Use subset of bands (1-6) to avoid atmosphere influences (01/31/2014)
%  Version 10.12: More accurate number of days per year "num_yrs" (01/30/2014)
%  Version 10.11: RMSE up% agriculture activity (9)dates with time series fit (01/26/2014)
%  Version 10.10: Update temperature extreme in recent studies (01/16/2014)
%  Version 10.9:  Find break in max value in any of the band (01/08/2014)
%  Version 10.8:  Add very simple fit with median value for start and end of timeseries (10/21/2013)
%  This version (10.8) is used for generating synthetic data for the LCMS project.
%  Command: TrendSeasonalFit_v10Plot('stack',N_row,N_col,mini=0.5,T_cg=3,n_times=3,conse=6,B_detect=2:6)
%  Version 10.7:  Better multitemporal cloud detection (10/19/2013)
%  Version 10.6:  Add "Tmax_cg" for last step noise removal (10/18/2013)
%  Version 10.5:  Use subset of bands (2-6) to avoid atmosphere influences (10/18/2013)
%  Version 10.4:  Let dynamic fitting for pixels at the beginning (09/23/2013)
%  Version 10.3:  Able to detect change at the verying beginning (09/06/2013)
%  Version 10.2:  Add mini years "mini_yrs" in model intialization (09/03/2013)
%  Version 10.1:  Reduce time for calcuating "v_dif" (09/02/2013)
%% Version 10.0:  Fit for beginning and end of the time series (08/31/2013)
%  Version 9.9:   Only fit more than 50% of Landat images overlap area (08/28/2013)
%  Version 9.8:   Force model fit for persistent snow pixels (08/27/2013)
%  Version 9.7:   Add "rec_cg.category" as indicator of fitting procudure (08/20/2013)
%                 Add rec_cg.change_prob as indicator of change probability (08/20/2013)
%                 Add rec_cg.num_obs ad indicator of number of observations (08/20/2013)
%  Version 9.6:   Remove mininum rmse "mini" and minimum years "mini_yrs" (08/16/2013)
%  Version 9.5:   Model gets more coefficients with more observations (08/16/2013)
%  Version 9.4:   Bug fixed in calculating temporally adjusted rmse (08/01/2013)
%  Version 9.3:   Fit curve again after one year (03/28/2013)
%  This version (9.3) is used for mapping land cover for the IDS project.
%  Command: TrendSeasonalFit_v9Plot('stack',N_row,N_col,T_cg=2,n_times=3,conse=4)
%  Version 9.2:   Use "mini = T_const/T_cg" for small rmse cases (03/26/2013)
%  Version 9.1:   Remove out of range pixels before time series analysis (02/09/2013)
% Version 9.0:   Using 8 coefficients and lasso fit (02/01/2013)
%  Version 8.4:   Use "max v_slope" instead of "average v_slope" (01/16/2013)
%  Version 8.3:   Start initialization when "time_span" > 1 year (01/16/2013)
%  Version 8.2:   Bug fix% agriculture activity (9)ed in not fitting models at the begining (01/16/2013)
%  Version 8.1:   Bug fixed in counting "i" and "i_span"(01/13/2013)
% Version 8.0:   Temporally changing RMSE (01/09/2013)
% Version 7.3:   Continuous Change Detection and Classification (CCDC) (07/11/2012)
%  This version (7.3) is explained by Zhu, Z. & Woodcock, C.E., Continuous Change
%  Detection and Classification (CCDC) of land cover using all available
%  Landsat data, Remote Sensing of Environment (2014).
%  Command: TrendSeasonalFit_v7Plot('stack',N_row,N_col,T_cg=3,n_times=3,conse=3)
% Version 1.0:   Continous Monitoring of Forest Disturbance Algorithm (CMFDA) (07/13/2010)
%  This version (1.0) is explained by Zhu, Z., Woodcock, C.E., Olofsson, P.,
%  Continuous monitoring of forest disturbance using all available Landsat
%  data, Remote Sensing of Environment (2012).
%
% Inputs:
% stk_n='stack'; stack image name
% ncols = 8021; % number of pixels processed per line
% nrows=1; % the nrowsth lines
% for example    1 2 3 4 5
%                6 7 8 9 10
%
% Outputs:
%
% rec_cg RECord information about all curves between ChanGes
% rec_cg(i).t_start record the start of the ith curve fitting (julian_date)
% rec_cg(i).t_end record the end of the ith curve fitting (julian_date)
% rec_cg(i).t_break record the first observed break time (julian_date)
% rec_cg(i).coefs record the coefficients of the ith curve
% rec_cg(i).pos record the position of the ith pixel (pixel id)
% rec_cg(i).magnitude record the change vector of all spectral bands
% rec_cg(i).category record what fitting procudure and model is used
% cateogry category 5x: persistent snow    4x: Fmask fails
% cateogry category 3x: modified fit       2x: end fit
% category category 1x: start fit           x: normal procedure
% cateogry category x1: mean value         x4: simple model
% category category x6: advanced model     x8: full model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  defining variables
warning('off', 'all'); % no warnings for linear fit 
%% Constants
p = inputParser;
addParameter(p,'compositeInterval', 24); % probability for detecting surface change
addParameter(p,'recordSR', 1); % number of consec

parse(p,varargin{:});
compositeInterval = p.Results.compositeInterval; % Composite interval 
recordSR = p.Results.recordSR; % Whether to record the observation values

% number of indices used as additional variables
n_indices = 5;

% ADD indices 
extend_B_detect = 8:(8+n_indices-1);
B_detect = [B_detect, extend_B_detect];
% maximum number of coefficient required
% 2 for tri-modal; 2 for bi-modal; 2 for seasonality; 2 for linear; 1 for
% water level 
min_num_c = 5;
mid_num_c = 7;
max_num_c = 9;

% Used to label the pixel location
ncols = 0; % To keep same as the line-based processing
nrows = 1; % To keep same as the line-based processing

% number of clear observation / number of coefficients
n_times = 3;
def_n_times = n_times;
% update frequency
prct_update_model = 0.03; % 3% of the previous model fit
% initialize NUM of Functional Curves
num_fc = 0;
% number of days per year
num_yrs = 365.25;
% Band for multitemporal cloud/snow detection (Green)
num_B1 = 2;
% Band for multitemporal shadow/snow shadow detection (SWIR){[5,6,7,8,9],[10]},[5,6]
num_B2 = 5;
% Threshold for cloud, shadow, and snow detection.
T_const = 4.42;
% minimum year for model intialization
mini_yrs = 1;
% no change detection for permanent snow pixels
t_sn = 0.75;
% threshold (degree) of mean included angle 
nsign = 45;


% consecutive number
def_conse = conse;

% Tmasking of noise
Tmax_cg = chi2inv(Tmax_cg,length(B_detect));
% adjust threshold based on chi-squared distribution
def_pT_cg = T_cg;
def_T_cg = chi2inv(def_pT_cg,length(B_detect));

% initialize the struct data of RECording of ChanGe (rec_cg)
rec_cg = struct('t_start',[],'t_end',[],'t_break',[],'coefs',[],'rmse',[],...
    'pos',[],'change_prob',[],'num_obs',[],'category',[],'magnitude',[],'durchange',[],'data',[]);

%% Some additional process for CSV data
plotIDs = proc_cols;
proc_cols = unique(proc_cols);
sdate_all_samples = sdate_all;
if ~isempty(path_t_all)
    % Single path process at first to remain the target path data
    line_t(path_t_all ~= path_r,end) = 254;
end

%% End of CSV process


%% @ each pixel
% for i_ids = 1:ncols % move pixel via the column demension
% making sure proc_cols's demension OK for FOR function
if size(proc_cols,1) > 1
    proc_cols = proc_cols';
end

flagPixels = 0;

for i_ids = proc_cols % define in the parent funtion, which can be 1:ncols, and one col for locating a certain pixel
    flagPixels = flagPixels+1;
    if mod(flagPixels,500)==0
        fprintf('%d (%d)..', flagPixels,round(toc/60));
%         return % for test
    end
    % get default conse & T_cg
    conse = def_conse;
    T_cg = def_T_cg;
    n_times = def_n_times;
    %% Modification for CSV-based process
    idxSample = plotIDs == i_ids;
    obs = line_t(idxSample,:);
    sdate_all = sdate_all_samples(idxSample,:);    
    obs = [obs,sdate_all];


    % Generate the composite data
    %% Take the composite
%     fprintf('Composite the observations within %d intervals\n',compositeInterval);   
    if compositeInterval>0 
        if mod(compositeInterval,16)==0
            intervals = sdate_all(1):compositeInterval:sdate_all(end)+compositeInterval;
            compositeObs = zeros(length(intervals)-1,size(obs,2));
    
            for i_v = 1:length(intervals)-1
                interval = intervals(i_v):intervals(i_v+1)-1;
                idx = ismember(sdate_all,interval);
    
                tempObs = obs(idx,:);
                clearIdx = tempObs(:,nbands) < 2; % Clear observations based on Fmask
                clearObs = tempObs(clearIdx,:);
    
                % Composite
                compositeInfo = zeros(1,size(obs,2));
                if size(clearObs,1)>1
                    ratios = clearObs(:,4)./(clearObs(:,1)+eps);
                    compositeInfo = clearObs(ratios==max(ratios),:);
                elseif size(clearObs,1)==1
                    compositeInfo = clearObs;
                elseif ~isempty(tempObs)
                    compositeInfo = tempObs(1,:);
                else            
                    compositeInfo(1,nbands) = 254;
                end
                compositeObs(i_v,:) = compositeInfo(1,:);
            end
    
            obs = compositeObs(:,1:nbands);
            sdate = compositeObs(:,end);
        else
            % Combined composite values
            compositeIntervalSparse = compositeInterval*2/3;
            compositeIntervalDense = compositeInterval*4/3;
            maxAnnualClearObs = ceil(365/compositeIntervalDense);
           
            biSatelliteDate = datenum(1999,4,1);
            landsat5End = datenum(2011,12,1);
            landsat8Start = datenum(2013,3,18);
    
            intervals = [sdate_all(1):compositeIntervalSparse:biSatelliteDate,...
                biSatelliteDate:compositeIntervalDense:landsat5End,...
                landsat5End:compositeIntervalSparse:landsat8Start,...
                landsat8Start:compositeIntervalDense:sdate_all(end)+compositeIntervalDense];
            
            compositeObs = zeros(length(intervals)-1,size(obs,2));
    
            for i_v = 1:length(intervals)-1
              
                interval = intervals(i_v):intervals(i_v+1)-1;
                idx = ismember(sdate_all,interval);
    
                tempObs = obs(idx,:);
                clearIdx = tempObs(:,nbands) < 2; % Clear observations based on Fmask
                clearObs = tempObs(clearIdx,:);
    
                % Composite
                compositeInfo = zeros(1,size(obs,2));
                if size(clearObs,1)>1
                    % Use composite approach select the optimal value
                    ratios = clearObs(:,4)./(clearObs(:,1)+eps);
                    compositeInfo = clearObs(ratios==max(ratios),:);
                elseif size(clearObs,1)==1
                    % Only one clear observation available
                    compositeInfo = clearObs;
                elseif ~isempty(tempObs)
                    % No clear observation yet cloud observation available
                    compositeInfo = tempObs(1,:);
                else
                    % No Landsat observation available
                    compositeInfo(1,nbands) = 254; % No composite available in the interval
                end
                compositeObs(i_v,:) = compositeInfo(1,:);
            end
        
            %% To confirm that the there are no more than monthly frequency every year
            obsYears = year(compositeObs(:,end-1));
            idxValidObs = compositeObs(:,nbands)<2;
            for yr = years
                if sum(idxValidObs(obsYears==yr)) > maxAnnualClearObs
                    obsWithInYear = idxValidObs(obsYears ==yr);
                    y = randsample(find(obsWithInYear==1),maxAnnualClearObs);
                    obsWithInYear(:) = 0;
                    obsWithInYear(y) = 1;
                    idxValidObs(obsYears == yr) = obsWithInYear;
                end
            end
            compositeObs(~idxValidObs,nbands) = 255;
            
            %% Update of the observation values, observation dates, and path information
            obs = compositeObs(:,1:nbands);
            sdate = compositeObs(:,end);
        end
        clear compositeObs;
    else
%         fprintf('Use the all observations without composite \n');  
        sdate = obs(:,end);
        obs = obs(:,1:nbands);
    end
    
    
    % mask data
    line_m = obs(:,end);
    % maintain single path data for the current pixel 
    % if we only consider single path Landsat data
    % path_r(i_ids) means the path value of the current pixel
    
    
    % calculate the indices
    blue = obs(:,1);
    green = obs(:,2);
    red = obs(:,3);
    nir = obs(:,4);
    swir1 = obs(:,5);
    swir2 = obs(:,6);
    ndvi = (nir - red)./(red + nir+eps)*10000;
    tcg = -0.2838*blue+0.2435*green-0.5436*red+0.7243*nir+0.0840*swir1-0.1800*swir2;
    mndwi = (green - swir1)./(green + swir1+eps)*10000;    
    tcw = 0.1509*blue+0.1973*green+0.3279*red+0.3406*nir-0.7112*swir1-0.4572*swir2;
    tcb = 0.3037*blue+0.2793*green+0.4743*red+0.5585*nir+0.5082*swir1+0.1863*swir2;
%     NBR = (nir-swir2)./(swir2 + nir+eps)*10000;
%     kNDVI = tanh(ndvi.*ndvi/10000/10000)*10000;
%     NIRv = nir.*ndvi/10000;

    indices = [ndvi,mndwi,tcg,tcw,tcb];
%     indices = [ndvi,tcg,NBR,kNDVI,NIRv];
    % clear pixel should have reflectance between 0 and 1
    % brightness temperature should between -93.2 to 70.7 celsius degree
    % convert Kelvin to Celsius
    obs(:,7) = 10*(obs(:,7) - 2732);
    idrange = obs(:,1)>0&obs(:,1)<10000&...
        obs(:,2)>0&obs(:,2)<10000&...
        obs(:,3)>0&obs(:,3)<10000&...
        obs(:,4)>0&obs(:,4)<10000&...
        obs(:,5)>0&obs(:,5)<10000&...
        obs(:,6)>0&obs(:,6)<10000&...
        obs(:,7)>-9320&obs(:,7)<7070;
    
    % # of clear observatons
    idclr = line_m < 2;
    
    % # of all available observations
    % idall = line_m < 255;
    % clear observation percentage
    % clr_pct = sum(idclr)/sum(idall);
    % snow pixels
    idsn = line_m == 3;
    % percent of snow observations
    sn_pct = sum(idsn)/(sum(idclr)+sum(idsn)+0.01);
    
    % not enough clear observations for change detection
    if sum(idclr) < n_times*max_num_c
        % permanent snow pixels
        if sn_pct > t_sn
            % snow observations are "good" now
            idgood = idsn|idclr;
            % number of snow pixel within range
            n_sn = sum(idgood);
            
            if n_sn < n_times*min_num_c % not enough snow pixels
                continue
            else
                % Xs & Ys for computation
                clrx = sdate(idgood);
                % bands 1-5,7,6
                clry = obs(idgood,1:end-2);
                clry = double(clry);
                
                % combine indices into the reflectance values
                indices = indices(idgood,:);
                clry = [clry,indices];
                
                % tide info 
                wl = obs(idgood,end-1);
                wl = double(wl);
                
                % find repeated ids
                [clrx,uniq_id,~] = unique(clrx);
                % mean of repeated values
                tmp_y = zeros(length(clrx),nbands+n_indices-2);
                tmp_wl = zeros(length(clrx));
                
                % get the mean values
                for i = 1:nbands+n_indices-2
                    tmp_y(:,i) = clry(uniq_id,i);
                end
                clry = tmp_y;
                
                tmp_wl = wl(uniq_id);
                wl = tmp_wl;
                
                % the first observation for TSFit
                i_start = 1;
                % identified and move on for the next curve
                num_fc = num_fc + 1; % NUM of Fitted Curves (num_fc)
                
                % defining computed variables
                fit_cft = zeros(max_num_c,nbands+n_indices-2);
                % rmse for each band
                rmse = zeros(nbands+n_indices-2,1);
                % snow qa = 50
                qa = 50;
                
                for i_B=1:nbands+n_indices-2
                    if i_B ~= nbands-2 % treat saturated and unsaturated pixels differently
                        idgood = clry(:,i_B) < 10000; % saturate if ref > 1 or NBR NDVR > 1
                        i_span = sum(idgood);
                        if i_span < min_num_c*n_times % fill value for frequently saturated snow pixels
                            fit_cft(1,i_B) = 10000; % give constant value
                        else % fit for enough unsaturat snow pixels
                            [fit_cft(:,i_B),rmse(i_B)]=autoTSFitWL(clrx(idgood),clry(idgood,i_B),min_num_c,wl(idgood));
                        end
                    else % fit for temperature band
                        idgood = clry(:,i_B)>-9320&clry(:,i_B)<7070;
                        [fit_cft(:,i_B),rmse(i_B)]=autoTSFitWL(clrx(idgood),clry(idgood,i_B),min_num_c,wl((idgood)));
                    end
                end
                
                % updating information at each iteration
                % record time of curve start
                rec_cg(num_fc).t_start=clrx(i_start);
                % record time of curve end
                rec_cg(num_fc).t_end=clrx(end);
                % record break time
                rec_cg(num_fc).t_break = 0; % no break at the moment
                % record postion of the pixel
                rec_cg(num_fc).pos = (nrows-1)*ncols+i_ids;
                % record fitted coefficients
                rec_cg(num_fc).coefs = fit_cft;
                % record rmse of the pixel
                rec_cg(num_fc).rmse = rmse;
                % record change probability
                rec_cg(num_fc).change_prob = 0;
                % record number of observations
                rec_cg(num_fc).num_obs = n_sn;
                % record fit category
                rec_cg(num_fc).category = qa + min_num_c;
                % record change magnitude
                rec_cg(num_fc).magnitude = zeros(1,nbands+n_indices-2);
                % record durchange
                rec_cg(num_fc).durchange = zeros(2,nbands+n_indices-2);

                if recordSR
%                     % record clear obs dates
%                     rec_cg(num_fc).obs = clrx(idgood);
%                     
%                     % record SR and indices
%                     rec_cg(num_fc).bands = clry(idgood,:);
%                     % record the tide information
%                     rec_cg(num_fc).tide = wl(idgood);                    
                    data = [clrx,clry,wl];
                    rec_cg(num_fc).data = array2table(data(idgood,:),'VariableNames',{'sdate','blue','green','red','nir','swir1','swir2','tem','ndvi','mndvi','tcg','tcw','tcb','tide'});
                end
            end
        else % no change detection for clear observations
            % within physical range pixels
            idgood = idrange;
            
            % Xs & Ys for computation
            clrx = sdate(idgood);
            % bands 1-5,7,6
            clry = obs(idgood,1:end-2);
            clry = double(clry);
            % combine indices into the reflectance values
            indices = indices(idgood,:);
            clry = [clry,indices];
            % tide info
            wl = obs(idgood,end-1);
            wl = double(wl);
            
            % find repeated ids
            [clrx,uniq_id,~] = unique(clrx);
            % mean of repeated values
            tmp_y = zeros(length(clrx),nbands+n_indices-2);
            tmp_wl = zeros(length(clrx));
            % get the mean values
            for i = 1:nbands+n_indices-2
                tmp_y(:,i) = clry(uniq_id,i);
            end
            clry = tmp_y;
            
            tmp_wl = wl(uniq_id);
            wl = tmp_wl;
            idclr = clry(:,num_B1) < median(clry(:,num_B1)) + 400;
            n_clr = sum(idclr);
            
            if n_clr < n_times*min_num_c % not enough clear pixels
                continue
            else
                % Xs & Ys for computation
                clrx = clrx(idclr);
                clry = clry(idclr,:);
                wl = wl(idclr);
                % the first observation for TSFit
                i_start = 1;
                % identified and move on for the next curve
                num_fc = num_fc + 1; % NUM of Fitted Curves (num_fc)
                
                % defining computed variables
                fit_cft = zeros(max_num_c,nbands+n_indices-2);
                % rmse for each band
                rmse = zeros(nbands+n_indices-2,1);
                % Fmask fail qa = 40
                qa = 40;
                
                for i_B = 1:nbands+n_indices-2
                    % fit basic model for all within range snow pixels
                    [fit_cft(:,i_B),rmse(i_B)] = autoTSFitWL(clrx,clry(:,i_B),min_num_c,wl);
                end
                
                % record time of curve start
                rec_cg(num_fc).t_start = clrx(i_start);
                % record time of curve end
                rec_cg(num_fc).t_end = clrx(end);
                % record break time
                rec_cg(num_fc).t_break = 0;
                % record postion of the pixel
                rec_cg(num_fc).pos = (nrows-1)*ncols+i_ids;
                % record fitted coefficients
                rec_cg(num_fc).coefs = fit_cft;
                % record rmse of the pixel
                rec_cg(num_fc).rmse = rmse;
                % record change probability
                rec_cg(num_fc).change_prob = 0;
                % record number of observations
                rec_cg(num_fc).num_obs = length(clrx);
                % record fit category
                rec_cg(num_fc).category = qa + min_num_c;
                % record change magnitude
                rec_cg(num_fc).magnitude = zeros(1,nbands+n_indices-2);
                % record during change
                rec_cg(num_fc).durchange = zeros(2,nbands+n_indices-2);

                if recordSR
                    clear obs dates
%                     rec_cg(num_fc).obs = clrx;
%                     % record SR and indices
%                     rec_cg(num_fc).bands = clry;
%                     % record the tide information
%                     rec_cg(num_fc).tide = wl;

                    data = [clrx,clry,wl];
                    rec_cg(num_fc).data = array2table(data,'VariableNames',{'sdate','blue','green','red','nir','swir1','swir2','tem','ndvi','mndvi','tcg','tcw','tcb','tide'});
                end
            end
        end
    else % normal CCDC procedure
        % clear and within physical range pixels
        idgood = idclr & idrange;

        % Xs & Ys for computation
        clrx = sdate(idgood);
        % bands 1-5,7,6
        clry = obs(idgood,1:end-2);
        clry = double(clry);
        indices = indices(idgood,:);
        clry = [clry,indices];
        
        wl = obs(idgood,end-1);
        wl = double(wl);
        
        % find repeated ids
        [clrx,uniq_id,~] = unique(clrx);
        
        % continue if not enough clear pixels
        if length(clrx) < n_times*min_num_c+conse
            continue
        end
        
        % mean of repeated values
        tmp_y = zeros(length(clrx),nbands+n_indices-2);
        tmp_wl = zeros(length(clrx));
        
        % get the mean values
        for i = 1:nbands+n_indices-2
            tmp_y(:,i) = clry(uniq_id,i);
        end
        clry = tmp_y;
        clear tmp_y;
        
        tmp_wl = wl(uniq_id);
        wl = tmp_wl;
        clear tmp_wl;
        
        % caculate median variogram
        var_clry = clry(2:end,:)-clry(1:end-1,:);
        adj_rmse = median(abs(var_clry),1);
        
        % start with the miminum requirement of clear obs
        i = n_times*min_num_c;
        
        % initializing variables
        % the first observation for TSFit
        i_start = 1;
        % record the start of the model initialization (0=>initial;1=>done)
        BL_train = 0;
        % identified and move on for the next curve
        num_fc = num_fc + 1; % NUM of Fitted Curves (num_fc)
        % record the num_fc at the beginning of each pixel
        rec_fc = num_fc;
        % initialize i_dense
        i_dense = 1;
        
        if recordSR
            % Initilization of for the clear observation recording
            i_start_segment = 1;
        end

        % while loop - process till the last clear observation - conse
        while i<= length(clrx)-conse
            % span of "i"
            i_span = i-i_start+1;
            % span of time (num of years)
            time_span = (clrx(i)-clrx(i_start))/num_yrs;
            % max time difference
            time_dif = max(clrx(i_start+1:i) - clrx(i_start:i-1));
            
            % basic requrirements: 1) enough observations; 2) enough time
            if i_span >= n_times*min_num_c && time_span >= mini_yrs
                % initializing model
                if BL_train == 0
                    % check max time difference
                    if time_dif > num_yrs
                        i = i+1;
                        i_start = i_start+1;
                        % i that is dense enough
                        i_dense = i_start;
                        continue
                    end
                    % Tmask: noise removal (good => 0 & noise => 1)
                    blIDs = autoTmask(clrx(i_start:i+conse),clry(i_start:i+conse,[num_B1,num_B2]),...
                        (clrx(i+conse)-clrx(i_start))/num_yrs,adj_rmse(num_B1),adj_rmse(num_B2),T_const);
                    
                    % IDs to be removed
                    IDs = i_start:i+conse;
                    rmIDs = IDs(blIDs(1:end-conse) == 1);
                    
                    % update i_span after noise removal
                    i_span = sum(~blIDs(1:end-conse));
                    
                    % check if there is enough observation
                    if i_span < n_times*min_num_c
                        % move forward to the i+1th clear observation
                        i = i+1;
                        % not enough clear observations
                        continue;
                        % check if there is enough time
                    else
                        % copy x & y
                        cpx = clrx;
                        cpy = clry;
                        cpwl = wl;
                        % remove noise pixels between i_start & i
                        cpx(rmIDs) = [];
                        cpy(rmIDs,:) = [];
                        cpwl(rmIDs) = [];
                        
                        % record i before noise removal
                        % This is very important as if model is not initialized
                        % the multitemporal masking shall be done again instead
                        % of removing outliers in every masking
                        i_rec = i;
                        
                        % update i afer noise removal (i_start stays the same)
                        i = i_start+i_span-1;
                        % update span of time (num of years)
                        time_span = (cpx(i)-cpx(i_start))/num_yrs;
                        
                        % check if there is enough time
                        if time_span < mini_yrs
                            % keep the original i
                            i = i_rec;
                            % move forward to the i+1th clear observation
                            i = i+1;
                            % not enough time
                            continue;
                            % Step 2: model fitting
                        else
                            % remove noise
                            clrx = cpx;
                            clry = cpy;
                            wl = cpwl;
                            % Step 2: model fitting
                            % initialize model testing variables
                            % defining computed variables
                            fit_cft = zeros(max_num_c,nbands+n_indices-2);
                            % rmse for each band
                            rmse = zeros(nbands+n_indices-2,1);
                            % value of differnce
                            v_dif = zeros(nbands+n_indices-2,1);
                            % record the diference in all bands
                            rec_v_dif = zeros(i-i_start+1,nbands+n_indices-2);
                            
                            % update number of coefficients
                            update_num_c = update_cft(i_span,n_times,min_num_c,mid_num_c,max_num_c,num_c);
                            
                            for i_B = 1:nbands+n_indices-2
                                % initial model fit
                                [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                                    autoTSFitWL(clrx(i_start:i),clry(i_start:i,i_B),update_num_c,wl(i_start:i));
                            end
                            
                            
                            % normalized to z-score
                            for i_B = B_detect
                                % minimum rmse
                                mini_rmse = max(adj_rmse(i_B),rmse(i_B));
                                
                                % compare the first clear obs
                                v_start = rec_v_dif(1,i_B)/mini_rmse;
                                % compare the last clear observation
                                v_end = rec_v_dif(end,i_B)/mini_rmse;
                                % anormalized slope values
                                v_slope = fit_cft(2,i_B)*(clrx(i)-clrx(i_start))/mini_rmse;
                                
                                % differece in model intialization
                                v_dif(i_B) = abs(v_slope) + max(abs(v_start),abs(v_end));
                            end
                            v_dif = norm(v_dif(B_detect))^2;
                            
                            % find stable start for each curve
                            if v_dif > T_cg
                                % start from next clear obs
                                i_start = i_start + 1;
                                % move forward to the i+1th clear observation
                                i = i + 1;
                                % keep all data and move to the next obs
                                continue
                            else
                                % model ready!
                                BL_train = 1;
                                % count difference of i for each iteration
                                i_count = 0;
                                % label 0 for computing the update frequency of model fit
                                num_skip = 0;
                                
                                % find the previous break point
                                if num_fc == rec_fc
                                    % first curve
                                    i_break = 1;
                                else
                                    % after the first curve
                                    i_break = find(clrx >= rec_cg(num_fc-1).t_break);
                                    i_break = i_break(1);
                                end
                                
                                if i_start > i_break
                                    % model fit at the beginning of the time series
                                    for i_ini = i_start-1:-1:i_break
                                        if i_start - i_break < conse
                                            ini_conse = i_start - i_break;
                                        else
                                            ini_conse = conse;
                                        end
                                        % value of difference for conse obs
                                        v_dif = zeros(ini_conse,nbands+n_indices-2);
                                        % record the magnitude of change
                                        v_dif_mag = v_dif;
                                        % record the date
                                        v_dif_clrx = zeros(ini_conse,1);
                                        % chagne vector magnitude
                                        vec_mag = zeros(ini_conse,1);
                                        
                                        for i_conse = 1:ini_conse
                                            for i_B = 1:nbands+n_indices-2
                                                % absolute difference
                                                v_dif_mag(i_conse,i_B) = clry(i_ini-i_conse+1,i_B)-autoTSPredWL(clrx(i_ini-i_conse+1),fit_cft(:,i_B),wl(i_ini-i_conse+1));
                                                % normalized to z-scores
                                                if sum(i_B == B_detect)
                                                    % minimum rmse
                                                    mini_rmse = max(adj_rmse(i_B),rmse(i_B));
                                                    
                                                    % z-scores
                                                    v_dif(i_conse,i_B) = v_dif_mag(i_conse,i_B)/mini_rmse;
                                                end
                                            end
                                            vec_mag(i_conse) = norm(v_dif(i_conse,B_detect))^2;
                                            v_dif_clrx(i_conse) = clrx(i_ini-i_conse+1);
                                        end
                                        
                                        % get the vec sign
                                        % vec_sign = abs(sum(sign(v_dif(:,B_detect)),1));
                                        max_angl = mean(angl(v_dif(:,B_detect)));
                                        
                                        % change detection
                                        if min(vec_mag) > T_cg && max_angl < nsign% change detected
                                            break
                                        elseif vec_mag(1) > Tmax_cg % false change
                                            % remove noise
                                            clrx(i_ini) = [];
                                            clry(i_ini,:) = [];
                                            wl(i_ini) = [];
                                            i=i-1; % stay & check again after noise removal
                                        end
                                        
                                        % update new_i_start if i_ini is not a confirmed break
                                        i_start = i_ini;
                                    end
                                end
                                % only fit first curve if 1) have more than
                                % conse obs 2) previous obs is less than a year
                                if num_fc == rec_fc && i_start - i_dense >= conse
                                    % defining computed variables
                                    fit_cft = zeros(max_num_c,nbands+n_indices-2);
                                    % rmse for each band
                                    rmse = zeros(nbands+n_indices-2,1);
                                    % start fit qa = 10
                                    qa = 10;
                                    
                                    update_num_c = update_cft(i_span,n_times,min_num_c,mid_num_c,max_num_c,num_c);
                                    
                                    for i_B=1:nbands+n_indices-2
                                        [fit_cft(:,i_B),rmse(i_B)] = ...
                                            autoTSFitWL(clrx(i_dense:i_start-1),clry(i_dense:i_start-1,i_B),update_num_c,wl(i_dense:i_start-1));
                                    end
                                    
                                    % record time of curve end
                                    rec_cg(num_fc).t_end = clrx(i_start-1);
                                    % record postion of the pixel
                                    rec_cg(num_fc).pos = (nrows-1)*ncols + i_ids;
                                    % record fitted coefficients
                                    rec_cg(num_fc).coefs = fit_cft;
                                    % record rmse of the pixel
                                    rec_cg(num_fc).rmse = rmse;
                                    % record break time
                                    rec_cg(num_fc).t_break = clrx(i_start);
                                    % record change probability
                                    rec_cg(num_fc).change_prob = 1;
                                    % record time of curve start
                                    rec_cg(num_fc).t_start = clrx(1);
                                    % record fit category
                                    rec_cg(num_fc).category = qa + min_num_c;
                                    % record number of observations
%                                     rec_cg(num_fc).num_obs = i_start - i_dense;
                                    rec_cg(num_fc).num_obs = i_start - 1;
                                    % record change magnitude
                                    rec_cg(num_fc).magnitude = - median(v_dif_mag,1);
                                    % record during-change
                                    rec_cg(num_fc).durchange = durchange(nbands, n_indices, v_dif_clrx, - v_dif_mag);

                                    if recordSR
                                        % record clear obs dates
%                                         rec_cg(num_fc).obs = clrx(1:i_start-1);
%                                         % record SR and indices
%                                         rec_cg(num_fc).bands = clry(1:i_start-1,:);
%                                         % record the tide information
%                                         rec_cg(num_fc).tide = wl(1:i_start-1);

                                        data = [clrx,clry,wl];
                                        rec_cg(num_fc).data = array2table(data(1:i_start-1,:),'VariableNames',{'sdate','blue','green','red','nir','swir1','swir2','tem','ndvi','mndvi','tcg','tcw','tcb','tide'});
                                    end
                                    % identified and move on for the next functional curve
                                    num_fc = num_fc + 1;


                                    if recordSR
                                        % Initilization of for the clear observation recording
                                        i_start_segment = i_start;
                                    end
                                end
                            end
                        end
                    end
                end
                
                % continuous monitoring started!!!
                if BL_train == 1
                    % all IDs
                    IDs = i_start:i;
                    % span of "i"
                    i_span = i-i_start+1;
                    % # of observations skipping to fit model
                    if num_skip == 0 % if is 0, not ready to skip for this observation, and thus to compute when we can start to skip
                        num_skip = fix(i_span*prct_update_model); % close to lower value
                    end
                    
                    % determine the time series model
                    update_num_c = update_cft(i_span,n_times,min_num_c,mid_num_c,max_num_c,num_c);
                    
                    % initial model fit when there are not many obs
                    if  i_count == 0 || i_span <= max_num_c*n_times
                        % update i_count at each interation
                        i_count = clrx(i)-clrx(i_start);
                        
                        % record previous end i
                        pre_end = i;
                        % label as 0, indicating for being not ready to skip next observations
                        num_skip = 0; 
                        
                        % defining computed variables
                        fit_cft = zeros(max_num_c,nbands+n_indices-2);
                        % rmse for each band
                        rmse = zeros(nbands+n_indices-2,1);
                        % record the diference in all bands
                        rec_v_dif = zeros(length(IDs),nbands+n_indices-2);
                        % normal fit qa = 0
                        qa = 0;
                        
                        for i_B=1:nbands+n_indices-2
                            [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                                autoTSFitWL(clrx(IDs),clry(IDs,i_B),update_num_c,wl(IDs));
                        end
                        
                        % updating information for the first iteration
                        % record time of curve start
                        rec_cg(num_fc).t_start = clrx(i_start);
                        % record time of curve end
                        rec_cg(num_fc).t_end = clrx(i);
                        % record break time
                        rec_cg(num_fc).t_break = 0; % no break at the moment
                        % record postion of the pixel
                        rec_cg(num_fc).pos = (nrows-1)*ncols+i_ids;
                        % record fitted coefficients
                        rec_cg(num_fc).coefs = fit_cft;
                        % record rmse of the pixel
                        rec_cg(num_fc).rmse = rmse;
                        % record change probability
                        rec_cg(num_fc).change_prob = 0;
                        % record number of observations
                        rec_cg(num_fc).num_obs = i-i_start+1;
                        % record fit category
                        rec_cg(num_fc).category = qa + update_num_c;
                        % record change magnitude
                        rec_cg(num_fc).magnitude = zeros(1,nbands+n_indices-2);
		                % record durchange
		                rec_cg(num_fc).durchange = zeros(2, nbands+n_indices-2);

                        if recordSR
                            % record clear obs dates
%                             rec_cg(num_fc).obs = clrx(i_start_segment:i);
%                             % record SR and indices
%                             rec_cg(num_fc).bands = clry(i_start_segment:i,:);
%                             % record the tide information
%                             rec_cg(num_fc).tide = wl(i_start_segment:i);

                            data = [clrx,clry,wl];
                            rec_cg(num_fc).data = array2table(data(i_start_segment:i,:),'VariableNames',{'sdate','blue','green','red','nir','swir1','swir2','tem','ndvi','mndvi','tcg','tcw','tcb','tide'});
                        end

                        % detect change
                        % value of difference for conse obs
                        v_dif = zeros(conse,nbands+n_indices-2);
                        % record the magnitude of change
                        v_dif_mag = v_dif;
                        v_dif_clrx = zeros(conse,1);
                        
                        vec_mag = zeros(conse,1);
                        
                        for i_conse = 1:conse
                            for i_B = 1:nbands+n_indices-2
                                % absolute difference
                                v_dif_mag(i_conse,i_B) = clry(i+i_conse,i_B)-autoTSPredWL(clrx(i+i_conse),fit_cft(:,i_B),wl(i+i_conse));

                                % normalized to z-scores
                                if sum(i_B == B_detect)
                                    % minimum rmse
                                    mini_rmse = max(adj_rmse(i_B),rmse(i_B));
                                    
                                    % z-scores
                                    v_dif(i_conse,i_B) = v_dif_mag(i_conse,i_B)/mini_rmse;
                                end
                            end
                            v_dif_clrx(i_conse) = clrx(i+i_conse);
                            vec_mag(i_conse) = norm(v_dif(i_conse,B_detect))^2;
                        end
                        % IDs that haven't updated
                        IDsOld = IDs;
                    else
                        if i - pre_end > num_skip % time to update the time series model
%                         if i - pre_end >= 3% % every 3 observations, update a model, version 13.04 recalled
                            % update i_count at each interation
                            i_count = clrx(i)-clrx(i_start);
                            
                            % record previous end i
                            pre_end = i;
                            % label as 0, indicating for being not ready to skip next observations
                            num_skip = 0; 
                        
                           % defining computed variables
                            fit_cft = zeros(max_num_c,nbands+n_indices-2);
                            % rmse for each band
                            rmse = zeros(nbands+n_indices-2,1);
                            % record the diference in all bands
                            rec_v_dif = zeros(length(IDs),nbands+n_indices-2);
                            % normal fit qa = 0
                            qa = 0;
                                                       
%                             % moving window with nyr = 5 start
%                             nyr = 5;
%                             id_nyr = find(clrx(IDs(end))-clrx(IDs) < nyr*num_yrs);
%                             if ~isempty(id_nyr) && length(id_nyr) >= max_num_c*n_times
%                                 IDs = IDs(id_nyr(1):end);
%                             end
%                             % end of moving window
                                                       
                            for i_B = 1:nbands+n_indices-2
                                [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                                    autoTSFitWL(clrx(IDs),clry(IDs,i_B),update_num_c,wl(IDs));
                            end
                            
                            % record fitted coefficients
                            rec_cg(num_fc).coefs = fit_cft;
                            % record rmse of the pixel
                            rec_cg(num_fc).rmse = rmse;
                            % record number of observations
                            rec_cg(num_fc).num_obs = i-i_start+1;
                            % record fit category
                            rec_cg(num_fc).category = qa + update_num_c;

                            % IDs that haven't updated
                            IDsOld = IDs;
                        end
                        
                        % record time of curve end
                        rec_cg(num_fc).t_end = clrx(i);
                        
                        if recordSR
%                             % record clear obs dates
%                             rec_cg(num_fc).obs = clrx(i_start_segment:i);
%                             % record SR and indices
%                             rec_cg(num_fc).bands = clry(i_start_segment:i,:);
%                             % record the tide information
%                             rec_cg(num_fc).tide = wl(i_start_segment:i);

                            data = [clrx,clry,wl];
                            rec_cg(num_fc).data = array2table(data(i_start_segment:i,:),'VariableNames',{'sdate','blue','green','red','nir','swir1','swir2','tem','ndvi','mndvi','tcg','tcw','tcb','tide'});
                        end

                        % use fixed number for RMSE computing
                        n_rmse = n_times*rec_cg(num_fc).category;
                        tmpcg_rmse = zeros(nbands+n_indices-2,1);
                        % better days counting for RMSE calculating
                        % relative days distance
                        d_rt = clrx(IDsOld) - clrx(i+conse);
                        d_yr = abs(round(d_rt/num_yrs)*num_yrs-d_rt);
                        
                        [~,sorted_indx] = sort(d_yr);
                        sorted_indx = sorted_indx(1:n_rmse);
                        
                        for i_B = B_detect
                            % temporally changing RMSE
                            tmpcg_rmse(i_B) = norm(rec_v_dif(IDsOld(sorted_indx)-IDsOld(1)+1,i_B))/...
                                sqrt(n_rmse-rec_cg(num_fc).category);
                        end
                        
                        % move the ith col to i-1th col
                        v_dif(1:conse-1,:) = v_dif(2:conse,:);
                        % only compute the difference of last consecutive obs
                        v_dif(conse,:) = 0;
                        % move the ith col to i-1th col
                        v_dif_mag(1:conse-1,:) = v_dif_mag(2:conse,:);
                        % record the magnitude of change of the last conse obs
                        v_dif_mag(conse,:) = 0;
                        % record the during change
                        v_dif_clrx(1:conse-1) = v_dif_clrx(2:conse);
                        % move the ith col to i-1th col
                        vec_mag(1:conse-1) = vec_mag(2:conse);
                        % change vector magnitude
                        vec_mag(conse) = 0;
                        
                        for i_B = 1:nbands+n_indices-2
                            % absolute difference
                            v_dif_mag(conse,i_B) = clry(i+conse,i_B)-autoTSPredWL(clrx(i+conse),fit_cft(:,i_B),wl(i+conse));
                            % normalized to z-scores
                            if sum(i_B == B_detect)
                                % minimum rmse
                                mini_rmse = max(adj_rmse(i_B),tmpcg_rmse(i_B));
                                
                                % z-scores
                                v_dif(conse,i_B) = v_dif_mag(conse,i_B)/mini_rmse;
                            end
                        end
                        v_dif_clrx(conse) = clrx(i+conse);
                        vec_mag(conse) = norm(v_dif(end,B_detect))^2;
                    end
                    
                    % sign of change vector
                    % vec_sign = abs(sum(sign(v_dif(:,B_detect)),1));
                    max_angl = mean(angl(v_dif(:,B_detect)));
                    
                    % change detection
                    if min(vec_mag) > T_cg && max_angl < nsign% change detected
                        % record break time
                        rec_cg(num_fc).t_break = clrx(i+1);
                        % record change probability
                        rec_cg(num_fc).change_prob = 1;
                        % record change magnitude
                        rec_cg(num_fc).magnitude = median(v_dif_mag,1);
                         
                        % update during change
                        rec_cg(num_fc).durchange = durchange(nbands,n_indices, v_dif_clrx, v_dif_mag);
                        

                        % identified and move on for the next functional curve
                        num_fc = num_fc + 1;
                        % start from i+1 for the next functional curve
                        i_start = i + 1;
                        % start training again
                        BL_train = 0;
                        
                        if recordSR
                            % Initilization of for the clear observation recording
                            i_start_segment = i_start;
                        end

                    elseif vec_mag(1) > Tmax_cg % false change
                        % remove noise
                        clrx(i+1) = [];
                        clry(i+1,:) = [];
                        wl(i+1) = [];
                        i=i-1; % stay & check again after noise removal
                    end
                end % end of continuous monitoring
            end % end of checking basic requrirements
            
            % move forward to the i+1th clear observation
            i=i+1;
        end % end of while iterative
        
        % Two ways for processing the end of the time series
        if BL_train == 1
            % 1) if no break find at the end of the time series
            % define probability of change based on conse
            for i_conse = conse:-1:1
                % sign of change vector
                % vec_sign = abs(sum(sign(v_dif(i_conse:conse,B_detect)),1));
                max_angl = mean(angl(v_dif(i_conse:conse,B_detect)));
                
                if vec_mag(i_conse) <= T_cg || max_angl >= nsign
                    % the last stable id
                    id_last = i_conse;
                    break;
                end
            end
            
            % update change probability
            rec_cg(num_fc).change_prob = (conse-id_last)/conse;
            % update end time of the curve
            rec_cg(num_fc).t_end=clrx(end-conse+id_last);
            
            if recordSR            
%                 % record clear obs dates
%                 rec_cg(num_fc).obs = clrx(i_start_segment:end);
%                 % record SR and indices
%                 rec_cg(num_fc).bands = clry(i_start_segment:end,:);
%                 % record the tide information
%                 rec_cg(num_fc).tide = wl(i_start_segment:end);

                data = [clrx,clry,wl];
                rec_cg(num_fc).data = array2table(data(i_start_segment:end,:),'VariableNames',{'sdate','blue','green','red','nir','swir1','swir2','tem','ndvi','mndvi','tcg','tcw','tcb','tide'});
            end
            

            if conse > id_last % > 1
                % update time of the probable change
                rec_cg(num_fc).t_break = clrx(end-conse+id_last+1);
                % update magnitude of change
                rec_cg(num_fc).magnitude = median(v_dif_mag(id_last+1:conse,:),1);
                % update during change
                rec_cg(num_fc).durchange = durchange(nbands, n_indices, v_dif_clrx, v_dif_mag);
            end
            

        elseif BL_train == 0
            % 2) if break find close to the end of the time series
            % Use [conse,min_num_c*n_times+conse) to fit curve
            
            if num_fc == rec_fc
                % first curve
                i_start = 1;
            else
                i_start = find(clrx >= rec_cg(num_fc-1).t_break);
                i_start = i_start(1);
            end
            
            % Tmask
            if length(clrx(i_start:end)) > conse
                blIDs = autoTmask(clrx(i_start:end),clry(i_start:end,[num_B1,num_B2]),...
                    (clrx(end)-clrx(i_start))/num_yrs,adj_rmse(num_B1),adj_rmse(num_B2),T_const);
                
                % update i_span after noise removal
                i_span = sum(~blIDs); %#ok<NASGU>
                
                IDs = i_start:length(clrx); % all IDs
                rmIDs = IDs(blIDs(1:end-conse) == 1); % IDs to be removed
                
                % remove noise pixels between i_start & i
                clrx(rmIDs) = [];
                clry(rmIDs,:) = [];
                wl(rmIDs) = [];
            end
            
            % enough data
            if length(clrx(i_start:end)) >= conse
                % defining computed variables
                fit_cft = zeros(max_num_c,nbands+n_indices-2);
                % rmse for each band
                rmse = zeros(nbands+n_indices-2,1);
                % end of fit qa = 20
                qa = 20;
                
                update_num_c = update_cft(i_span,n_times,min_num_c,mid_num_c,max_num_c,num_c);
                for i_B = 1:nbands+n_indices-2
                    [fit_cft(:,i_B),rmse(i_B)] = ...
                        autoTSFitWL(clrx(i_start:end),clry(i_start:end,i_B),update_num_c,wl(i_start:end));
                end
                
                % record time of curve start
                rec_cg(num_fc).t_start = clrx(i_start);
                % record time of curve end
                rec_cg(num_fc).t_end=clrx(end);
                % record break time
                rec_cg(num_fc).t_break = 0;
                % record postion of the pixel
                rec_cg(num_fc).pos = (nrows-1)*ncols+i_ids;
                % record fitted coefficients
                rec_cg(num_fc).coefs = fit_cft;
                % record rmse of the pixel
                rec_cg(num_fc).rmse = rmse;
                % record change probability
                rec_cg(num_fc).change_prob = 0;
                % record number of observations
                rec_cg(num_fc).num_obs = length(clrx(i_start:end));
                % record fit category
                rec_cg(num_fc).category = qa + min_num_c;
                % record change magnitude
                rec_cg(num_fc).magnitude = zeros(1,nbands+n_indices-2);
                % record durchange
                rec_cg(num_fc).durchange = zeros(2, nbands+n_indices-2);

                if recordSR
%                     % record clear obs dates
%                     rec_cg(num_fc).obs = clrx(i_start:end);
%                     % record SR and indices
%                     rec_cg(num_fc).bands = clry(i_start:end,:);
%                     % record the tide information
%                     rec_cg(num_fc).tide = wl(i_start:end);

                    data = [clrx,clry,wl];
                    rec_cg(num_fc).data = array2table(data(i_start:end,:),'VariableNames',{'sdate','blue','green','red','nir','swir1','swir2','tem','ndvi','mndvi','tcg','tcw','tcb','tide'});
                end
            end
        end
    end % end of if sum(idgood) statement    
   
end % end of for i_ids loop

end % end of function

% FUNCTION of updating the number of coeffs of the time series model
function update_num_c = update_cft(i_span,n_times,min_num_c,mid_num_c,max_num_c,num_c)

% determine the time series model
if i_span < mid_num_c*n_times
    % start with 5 coefficients model
    update_num_c = min(min_num_c,num_c);
elseif i_span < max_num_c*n_times
    % start with 7 coefficients model
    update_num_c = min(mid_num_c,num_c);
else
    % start with 9 coefficients model
    update_num_c =  min(max_num_c,num_c);
end

end

% FUNCTION of caculating included angle between ajacent pair of change vectors
function y = angl(v_dif)

[row,~] = size(v_dif);
y = zeros(row-1,1);

if row > 1
    for i = 1:row-1
        a = v_dif(i,:);
        b = v_dif(i+1,:);
        % y measures the opposite of cos(angle)
        y(i) = acos(a*b'/(norm(a)*norm(b)));
    end
else
    y = 0;
end

% convert angle from radiance to degree
y = y*180/pi;

end

% FUNCION of calculting dur-change trend with OLS fit
function v_dif_linear = durchange(nbands,n_indices, clrx, v_dif_mag)
    v_dif_linear = zeros(2,nbands+n_indices-2);
    for i_B = 1:nbands+n_indices-2
        p = polyfit(clrx, v_dif_mag(:,i_B),1);
        v_dif_linear(1, i_B) = p(1); % yfit =  p(1) * x + p(2)
        v_dif_pred = p(1) .* clrx + p(2);
        v_dif_linear(2, i_B) = sqrt(mean((v_dif_mag(:,i_B) - v_dif_pred).^2));  % Root Mean Squared Error
    end
end
